#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <Eigen/Cholesky>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/helper/AdvancedTimer.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 210600)
#include <sofa/helper/ScopedAdvancedTimer.h>
#endif
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/Mass.inl>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mass {

template<typename Element>
CaribouMass<Element>::CaribouMass()
: d_topology_container(initLink(
        "topology",
        "Topology container containing the elements on which this mass will be computed."))
, d_lumped(initData(
        &d_lumped,
        false,
        "lumped",
        "Whether or not the mass matrix should be lumped by scaling the diagonal entries such that the mass"
        "is constant per element. Note that the lumped matrix is always computed. But this parameter "
        "will determine if it (the lumped) matrix should be used to solve the acceleration (a = M^(-1).f)."))
, d_density(initData(
        &d_density,
        Real(1),
        "density",
        "Mass density of the material."))
{}

template<typename Element>
void CaribouMass<Element>::init() {
    using sofa::core::topology::BaseMeshTopology;
    using sofa::core::objectmodel::BaseContext;
    using CaribouTopology = SofaCaribou::topology::CaribouTopology<Element>;

    Inherit1::init();

    auto *context = this->getContext();

    if (!this->mstate) {
        msg_warning() << "No mechanical object found in the current context node. The data parameter "
                      << "'" << this->mstate.getName() << "' can be use to set the path to a mechanical "
                      << "object having a template of '" << DataTypes::Name() << "'";
        return;
    }

    // If not topology is specified, try to find one automatically in the current context
    if (not d_topology_container.get()) {
        // No topology specified. Try to find one suitable.
        auto caribou_containers = context->template getObjects<CaribouTopology>(BaseContext::Local);
        auto sofa_containers = context->template getObjects<BaseMeshTopology>(BaseContext::Local);
        std::vector<BaseMeshTopology *> sofa_compatible_containers;
        for (auto container : sofa_containers) {
            if (CaribouTopology::mesh_is_compatible(container)) {
                sofa_compatible_containers.push_back(container);
            }
        }
        if (caribou_containers.empty() and sofa_compatible_containers.empty()) {
            msg_warning() << "Could not find a topology container in the current context. "
                          << "Please add a compatible one in the current context or set the "
                          << "container's path using the '" << d_topology_container.getName()
                          << "' data parameter.";
        } else {
            if (caribou_containers.size() + sofa_compatible_containers.size() > 1) {
                msg_warning() << "Multiple topologies were found in the context node. "
                              << "Please specify which one contains the elements on "
                              << "which this force field will be applied "
                              << "by explicitly setting the container's path in the  '"
                              << d_topology_container.getName() << "' data parameter.";
            } else {
                // Prefer caribou's containers first
                if (not caribou_containers.empty()) {
                    d_topology_container.set(caribou_containers[0]);
                } else {
                    d_topology_container.set(sofa_compatible_containers[0]);
                }
            }

            msg_info() << "Automatically found the topology '" << d_topology_container.get()->getPathName()
                       << "'.";
        }
    }

    // Create a caribou internal Domain over the topology
    if (d_topology_container.get()) {
        auto sofa_topology = dynamic_cast<BaseMeshTopology *>(d_topology_container.get());
        auto caribou_topology = dynamic_cast<CaribouTopology *>(d_topology_container.get());
        if (sofa_topology) {
            // Initialize a new caribou topology from the SOFA topology
            p_topology = sofa::core::objectmodel::New<CaribouTopology>();
            p_topology->findData("indices")->setParent(CaribouTopology::get_indices_data_from(sofa_topology));
            p_topology->findData("position")->setParent(this->getMState()->findData("position"));
            p_topology->init();
        } else {
            // A Caribou topology already exists in the scene
            p_topology = caribou_topology;
        }

        if (number_of_elements() == 0) {
            msg_warning() << "No element found in the topology '" << d_topology_container.get()->getPathName() << "'";
        }
    }

    initialize_elements();
    assemble_mass_matrix();
}

template<typename Element>
template<typename Derived>
auto CaribouMass<Element>::canCreate(Derived *o, sofa::core::objectmodel::BaseContext *context,
                                     sofa::core::objectmodel::BaseObjectDescription *arg) -> bool {
    using namespace sofa::core::objectmodel;
    using CaribouTopology = SofaCaribou::topology::CaribouTopology<Element>;

    std::string requested_element_type = arg->getAttribute( "template", "");
    std::string this_element_type = Derived::templateName(o);

    // to lower
    std::string requested_element_type_lower = requested_element_type, this_element_type_lower = this_element_type;
    std::transform(requested_element_type.begin(), requested_element_type.end(), requested_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(this_element_type.begin(), this_element_type.end(), this_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });

    if (requested_element_type_lower == this_element_type_lower) {
        return Inherit1::canCreate(o, context, arg);
    }

    if (not requested_element_type.empty()) {
        arg->logError("Requested element type is not '"+this_element_type+"'.");
        return false;
    }

    std::string topology_path = arg->getAttribute("topology", "");
    if (not topology_path.empty()) {
        topology_path = topology_path.substr(1); // removes the "@"
        // Make sure the specified topology has elements of type Element
        auto topology = context->get<BaseObject>(topology_path);
        auto caribou_topology = dynamic_cast<CaribouTopology *>(topology);
        auto sofa_topology = dynamic_cast<sofa::core::topology::BaseMeshTopology *>(topology);
        if (not caribou_topology and (not sofa_topology or not CaribouTopology::mesh_is_compatible(sofa_topology))) {
            if (not sofa_topology) {
                arg->logError("The path '" + topology_path + "' doesn't point to either a SOFA mesh topology, or a Caribou topology.");
            } else {
                arg->logError("Cannot deduce the element type from the specified SOFA topology '" + topology->getPathName() + "'.");
            }
            return false;
        }
    } else {
        // Try to find a compatible topology in the current context
        BaseObject * topology = nullptr;
        auto objects = context->getObjects<BaseObject>(BaseContext::SearchDirection::Local);
        for (auto * object : objects) {
            auto caribou_topology = dynamic_cast<const CaribouTopology *>(object);
            auto sofa_topology = dynamic_cast<const sofa::core::topology::BaseMeshTopology *>(object);
            if (caribou_topology  or (sofa_topology and CaribouTopology::mesh_is_compatible(sofa_topology))) {
                topology = object;
                break;
            }
        }
        if (not topology) {
            arg->logError("Cannot find a topology in the current context from which the template '"+this_element_type+"' can be deduced.");
            return false;
        }

        if (Inherit1::canCreate(o, context, arg)) {
            arg->setAttribute("topology", "@" + topology->getPathName());
            return true;
        } else {
            return false;
        }
    }

    return Inherit1::canCreate(o, context, arg);
}

template<typename Element>
void CaribouMass<Element>::assemble_mass_matrix() {
    assemble_mass_matrix(*this->mstate->read (sofa::core::ConstVecCoordId::restPosition()));
}

template<typename Element>
void CaribouMass<Element>::assemble_mass_matrix(const sofa::Data<VecCoord> & x0) {
    using namespace sofa::core::objectmodel;

    const sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = x0;
    const auto nb_nodes = sofa_x0.size();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> X (sofa_x0.ref().data()->data(),  nb_nodes, Dimension);

    assemble_mass_matrix(X);
}

template<typename Element>
template<typename Derived>
void CaribouMass<Element>::assemble_mass_matrix(const Eigen::MatrixBase<Derived> & x0) {
    const auto density = d_density.getValue();
    if (density < std::numeric_limits<Real>::epsilon()) {
        return;
    }

    const auto nb_elements = this->number_of_elements();
    const auto nb_nodes = x0.rows();
    const auto nDofs = nb_nodes*Dimension;
    p_M.resize(nDofs, nDofs);
    p_Mdiag.setZero(nDofs);

    /// Triplets are used to store matrix entries before the call to 'compress'.
    /// Duplicates entries are summed up.
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(nDofs*3);

    sofa::helper::AdvancedTimer::stepBegin("CaribouMass::update_mass_matrix");
    for (int element_id = 0; element_id < static_cast<int>(nb_elements); ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = this->topology()->domain()->element_indices(element_id);

        // Assemble the element's mass matrix
        using Mass = Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodesPerElement*Dimension, NumberOfNodesPerElement*Dimension, Eigen::RowMajor>;
        Mass Me = Mass::Zero();
        Real element_mass = 0; ///< Used for the lumping algorithm

        for (const auto & gauss_node : gauss_nodes_of(element_id)) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Shape functions at gauss node
            const auto & N = gauss_node.N;

            // Element total mass (used for the lumping algorithm)
            element_mass += density*w*detJ;

            // Computation of the consistent mass sub-matrix M_ij
            for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
                const auto & Ni = N[i];

                // We now loop only on the upper triangular part of the
                // element stiffness matrix Ke since it is symmetric
                for (std::size_t j = i; j < NumberOfNodesPerElement; ++j) {
                    const auto & Nj = N[j];

                    // The 3x3 sub-matrix Mij is diagonal
                    Me.template block<3, 3>(i*Dimension, j*Dimension)
                      .template diagonal<0>()
                      += Vector<3>::Constant(density*Ni*Nj*w*detJ);
                }
            }
        }

        // Assembling the diagonal lumped mass matrix
        Real sum_of_nodal_masses = 0;
        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            sum_of_nodal_masses += Me(i*Dimension, i*Dimension);
        }
        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            // Node index of the ith node in the global stiffness matrix
            auto x = static_cast<int>(node_indices[i]*Dimension);
            // Factor that scales down the diagonal terms in such a way that the mass
            // is constant within the element (see Hinton et al. 1976 and Wriggers 2008)
            const auto scaling_factor = element_mass / sum_of_nodal_masses;
            const auto diagonal_mass = Me(i*Dimension, i*Dimension) * scaling_factor;
            for (int m = 0; m < Dimension; ++m) {
                p_Mdiag.diagonal()[x+m] += diagonal_mass;
            }
        }

        // Assembling the full consistent mass matrix
        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            // Node index of the ith node in the global stiffness matrix
            auto x = static_cast<int>(node_indices[i]*Dimension);
            for (std::size_t j = i; j < NumberOfNodesPerElement; ++j) {
                // Node index of the jth node in the global stiffness matrix
                auto y = static_cast<int>(node_indices[j]*Dimension);
                if (x > y) std::swap(x, y); // Fill-in the upper diagonal only
                for (int m = 0; m < Dimension; ++m) {
                    triplets.emplace_back(x+m, y+m, Me(i*Dimension+m,j*Dimension+m));
                }
            }
        }
    }
    p_M.setFromTriplets(triplets.begin(), triplets.end());
    sofa::helper::AdvancedTimer::stepEnd("CaribouMass::update_mass_matrix");
}

template <typename Element>
auto CaribouMass<Element>::get_gauss_nodes(const std::size_t & /*element_id*/, const Element & element) const -> GaussContainer {
    GaussContainer gauss_nodes {};
    if constexpr (NumberOfGaussNodesPerElement == caribou::Dynamic) {
        gauss_nodes.resize(element.number_of_gauss_nodes());
    }

    const auto nb_of_gauss_nodes = gauss_nodes.size();
    for (std::size_t gauss_node_id = 0; gauss_node_id < nb_of_gauss_nodes; ++gauss_node_id) {
        const auto & g = element.gauss_node(gauss_node_id);

        const auto J = element.jacobian(g.position);
        const auto detJ = std::abs(J.determinant());

        // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
        const Vector<NumberOfNodesPerElement> N = element.L(g.position);


        GaussNode & gauss_node = gauss_nodes[gauss_node_id];
        gauss_node.weight               = g.weight;
        gauss_node.jacobian_determinant = detJ;
        gauss_node.N                    = N;
    }

    return gauss_nodes;
}

template <typename Element>
void CaribouMass<Element>::initialize_elements()
{
    using namespace sofa::core::objectmodel;

    sofa::helper::AdvancedTimer::stepBegin("CaribouMass::initialize_elements");

    if (!this->mstate)
        return;

    // Resize the container of elements'quadrature nodes
    const auto nb_elements = this->number_of_elements();
    if (p_elements_quadrature_nodes.size() != nb_elements) {
        p_elements_quadrature_nodes.resize(nb_elements);
    }

    // Translate the Sofa's mechanical state vector to Eigen vector type
    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), sofa_x0.size(), Dimension);

    // Loop on each element and compute the shape functions and their derivatives for every of their integration points
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Get an Element instance from the Domain
        const auto initial_element = this->topology()->element(element_id);

        // Fill in the Gauss integration nodes for this element
        p_elements_quadrature_nodes[element_id] = get_gauss_nodes(element_id, initial_element);
    }

    // Compute the volume
    Real v = 0.;
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        for (const auto & gauss_node : gauss_nodes_of(element_id)) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            v += detJ*w;
        }
    }

    msg_info() << "Total mass of the geometry is " << v*d_density.getValue();

    sofa::helper::AdvancedTimer::stepEnd("CaribouMass::initialize_elements");
}

template<typename Element>
void CaribouMass<Element>::addForce(const sofa::core::MechanicalParams * /*mparams*/, CaribouMass::DataVecDeriv & d_f,
                                    const CaribouMass::DataVecCoord & d_x, const CaribouMass::DataVecDeriv & /*d_v*/) {
    const auto density = d_density.getValue();
    if (density < std::numeric_limits<Real>::epsilon()) {
        return;
    }

    const auto g_sofa = this->getContext()->getGravity();
    Eigen::Map<const Eigen::Matrix<Real, 3, 1>> g (g_sofa.data());

    // Map SOFA vectors to Eigen matrices
    const sofa::helper::ReadAccessor<DataVecCoord> sofa_x = d_x;
    sofa::helper::WriteAccessor<DataVecDeriv> sofa_f = d_f;
    const auto nb_nodes = sofa_x.size();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> x (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> f (&(sofa_f[0][0]),  nb_nodes, Dimension);

    // Start the timer
    sofa::helper::ScopedAdvancedTimer _t_ ("CaribouMass::addForce");

    // Assemble the force vector
    const auto nb_elements = this->number_of_elements();
    for (int element_id = 0; element_id < static_cast<int>(nb_elements); ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = this->topology()->domain()->element_indices(element_id);

        for (const auto & gauss_node : gauss_nodes_of(element_id)) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Shape functions at gauss node
            const auto & N = gauss_node.N;

            // Force at Gauss point
            const Eigen::Matrix<Real, 3, 1> force = g*density*w*detJ;

            // Force at the ith node
            for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
                const auto & Ni = N[i];
                f.row(node_indices[i]) += force*Ni;
            }
        }
    }
}

template<typename Element>
void CaribouMass<Element>::addMToMatrix(sofa::linearalgebra::BaseMatrix * matrix, SReal mFact, unsigned int & offset) {
    // Start the timer
    sofa::helper::ScopedAdvancedTimer _t_ ("CaribouMass::addMToMatrix");

    const bool lumped = d_lumped.getValue();

    if (lumped) {
        // Lumped diagonal mass matrix
        const auto n = p_Mdiag.rows();
        for (int i = 0; i < n; ++i) {
            const auto & v = p_Mdiag.diagonal()[i];
            matrix->add(offset + i, offset + i, v);
        }
    } else {
        // Sparse mass matrix
        for (int k = 0; k < p_M.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_M, k); it; ++it) {
                const auto i = it.row();
                const auto j = it.col();
                const auto v = it.value() * mFact;
                matrix->add(offset + i, offset + j, v);
                if (i != j)
                    matrix->add(offset + j, offset + i, v);
            }
        }
    }
}

template<typename Element>
void CaribouMass<Element>::addGravityToV(const sofa::core::MechanicalParams * mparams, CaribouMass::DataVecDeriv & d_v) {
    // Get an Eigen reference to v
    sofa::helper::WriteAccessor<DataVecDeriv> sofa_v = d_v;
    const auto nb_nodes = sofa_v.size();
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> v (&(sofa_v[0][0]),  nb_nodes, Dimension);

    // Get an Eigen reference to g
    const auto g_sofa = this->getContext()->getGravity();
    Eigen::Map<const Eigen::Matrix<Real, 3, 1>> g (g_sofa.data());

    // Get dt
    const auto dt = mparams->dt();

    // dt*g
    const Eigen::Matrix<Real, 1, 3> hg = (dt*g).transpose();
    sofa::helper::ScopedAdvancedTimer _t_ ("CaribouMass::addGravityToV");
    v.rowwise() += hg;
}

template<typename Element>
void CaribouMass<Element>::accFromF(const sofa::core::MechanicalParams * /*mparams*/, CaribouMass::DataVecDeriv & d_a,
                                    const CaribouMass::DataVecDeriv & d_f) {
    // Map SOFA vectors to Eigen matrices
    const sofa::helper::ReadAccessor<DataVecCoord> sofa_f = d_f;
    sofa::helper::WriteAccessor<DataVecDeriv> sofa_a = d_a;
    const auto nb_nodes = sofa_f.size();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> f (sofa_f.ref().data()->data(),  nb_nodes, Dimension);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> a (&(sofa_a[0][0]),  nb_nodes, Dimension);

    // a = M-1 f
    const bool lumped = d_lumped.getValue();
    if (lumped) {
        for (unsigned int i = 0; i < nb_nodes; ++i) {
            a.row(i) = f.row(i) / p_Mdiag.diagonal()[i*Dimension];
        }
    } else {
        Eigen::SimplicialCholesky<decltype(p_M), Eigen::Upper> solver;
        const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor> r = solver.compute(p_M).solve(f);
        a = r;
    }
}

template<typename Element>
void CaribouMass<Element>::addMDx(const sofa::core::MechanicalParams * /*mparams*/, CaribouMass::DataVecDeriv &d_f,
                                  const CaribouMass::DataVecDeriv &d_dx, SReal factor) {

    // Map SOFA vectors to Eigen matrices
    const sofa::helper::ReadAccessor<DataVecCoord> sofa_dx = d_dx;
    sofa::helper::WriteAccessor<DataVecDeriv> sofa_f = d_f;
    const auto nb_nodes = sofa_f.size();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> dx (sofa_dx.ref().data()->data(),  nb_nodes*Dimension);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> f (&(sofa_f[0][0]),  nb_nodes*Dimension);

    const bool lumped = d_lumped.getValue();

    if (lumped) {
        f = (p_Mdiag.diagonal().array() * dx.array()).matrix() * factor;
    } else {
        f = p_M * dx * factor;
    }
}

} // namespace SofaCaribou::mass
