#pragma once

#include "HyperelasticForcefield.h"
#include <sofa/helper/AdvancedTimer.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/visual/VisualParams.h>
#ifdef CARIBOU_WITH_OPENMP
#include <omp.h>
#endif

namespace SofaCaribou::forcefield {

template <typename Element>
HyperelasticForcefield<Element>::HyperelasticForcefield()
: d_topology_container(initLink(
    "topology",
    "Topology container containing the elements on which this forcefield will be applied."))
, d_material(initLink(
    "material",
    "Material used to compute the hyperelastic force field."))
, d_enable_multithreading(initData(&d_enable_multithreading,
    false,
    "enable_multithreading",
    "Enable the multithreading computation of the stiffness matrix. Only use this if you have a "
    "very large number of elements, otherwise performance might be worse than single threading."
    "When enabled, use the environment variable OMP_NUM_THREADS=N to use N threads."))
, d_drawScale(initData(&d_drawScale,
    0.85,
    "draw_scale",
    "Scaling factor for the drawing of elements (between 0 and 1). The factor allows to shrink "
    "the element relative to its center point when drawing it."))
{
}

template <typename Element>
void HyperelasticForcefield<Element>::init()
{
    Inherit::init();

    if (number_of_elements() == 0) {
        if (d_topology_container.get()) {
            msg_warning() << "No element found in the mesh topology '" << d_topology_container->getPathName() << "'.";
        } else if (not d_topology_container.get()) {
            // No topology specified. Try to find one suitable.
            auto containers = this->getContext()->template getObjects<sofa::core::topology::BaseMeshTopology>(
                BaseContext::Local);
            if (containers.empty()) {
                msg_warning() << "Could not find a topology container in the current context.";
            } else {
                std::vector<sofa::core::topology::BaseMeshTopology *> suitable_containers;
                for (const auto &container : containers) {
                    d_topology_container.set(container);
                    if (number_of_elements() > 0) {
                        suitable_containers.push_back(container);
                    }
                    d_topology_container.set(nullptr);
                }

                if (suitable_containers.empty()) {
                    msg_warning() << "Could not find a suitable topology container in the current context.";
                } else if (suitable_containers.size() > 1) {
                    msg_warning() <<
                                  "Multiple topology were found in the context node." <<
                                  " Please specify which one contains the elements on which this force field will be applied "
                                  <<
                                  "by explicitly setting the container's path in the  '"
                                  << d_topology_container.getName() << "'  parameter.";
                } else {
                    d_topology_container.set(suitable_containers[0]);
                    msg_info() << "Automatically found the topology '" << d_topology_container.get()->getPathName()
                               << "'.";
                }
            }
        }
    }

    // No material set, try to find one in the current context
    if (not d_material.get()) {
        auto materials = this->getContext()->template getObjects<material::HyperelasticMaterial<DataTypes>>(BaseContext::Local);
        if (materials.empty()) {
            msg_warning() << "Could not find an hyperelastic material in the current context.";
        } else if (materials.size() > 1) {
            msg_warning() << "Multiple materials were found in the context node. "   <<
                             "Please specify which one should be use by explicitly " <<
                             "setting the material's path in the '" << d_material.getName() << "' parameter.";
        } else {
            d_material.set(materials[0]);
            msg_info() << "Automatically found the material '" << d_material.get()->getPathName() << "'.";
        }
    }

    // Compute and store the shape functions and their derivatives for every integration points
    initialize_elements();

    // Update the stiffness matrix for every elements
    update_stiffness();
}

template <typename Element>
void HyperelasticForcefield<Element>::addForce(
    const MechanicalParams* mparams,
    Data<VecDeriv>& d_f,
    const Data<VecCoord>& d_x,
    const Data<VecDeriv>& d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    if (!this->mstate)
        return;

    const auto material = d_material.get();
    if (!material) {
        return;
    }

    // Update material parameters in case the user changed it
    material->before_update();

    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x = d_x;
    sofa::helper::WriteAccessor<Data<VecDeriv>> sofa_f = d_f;

    if (sofa_x.size() != sofa_f.size())
        return;
    const auto nb_nodes = sofa_x.size();
    const auto nb_elements = number_of_elements();

    if (nb_nodes == 0 || nb_elements == 0)
        return;

    if (p_elements_quadrature_nodes.size() != nb_elements)
        return;

    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> forces  (&(sofa_f[0][0]),  nb_nodes, Dimension);

    sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield::addForce");

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Fetch the node indices of the element
        const sofa::Index * node_indices = get_element_nodes_indices(element_id);

        // Fetch the initial and current positions of the element's nodes
        Matrix<NumberOfNodes, Dimension> current_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            current_nodes_position.row(i).noalias() = X.row(node_indices[i]);
        }

        // Compute the nodal forces
        Matrix<NumberOfNodes, Dimension> nodal_forces;
        nodal_forces.fill(0);

        for (GaussNode &gauss_node : p_elements_quadrature_nodes[element_id]) {

            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto & dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Deformation tensor at gauss node
            gauss_node.F.noalias() = current_nodes_position.transpose()*dN_dx;
            const auto & F = gauss_node.F;
            const auto J = F.determinant();

            // Right Cauchy-Green strain tensor at gauss node
            const Mat33 C = F.transpose() * F;

            // Second Piola-Kirchhoff stress tensor at gauss node
            const Mat33 S = material->PK2_stress(J, C);

            // Elastic forces w.r.t the gauss node applied on each nodes
            for (size_t i = 0; i < NumberOfNodes; ++i) {
                const auto dx = dN_dx.row(i).transpose();
                const Vector<Dimension> f_ = (detJ * w) * F*S*dx;
                for (size_t j = 0; j < Dimension; ++j) {
                    nodal_forces(i, j) += f_[j];
                }
            }
        }

        for (size_t i = 0; i < NumberOfNodes; ++i) {
            for (size_t j = 0; j < Dimension; ++j) {
                sofa_f[node_indices[i]][j] -= nodal_forces(i,j);
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield::addForce");

    K_is_up_to_date = false;
    eigenvalues_are_up_to_date = false;
}

template <typename Element>
void HyperelasticForcefield<Element>::addDForce(
    const MechanicalParams* mparams,
    Data<VecDeriv>& d_df,
    const Data<VecDeriv>& d_dx)
{
    if (not K_is_up_to_date) {
        update_stiffness();
    }

    auto kFactor = static_cast<Real> (mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue()));
    sofa::helper::ReadAccessor<Data<VecDeriv>> sofa_dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> sofa_df = d_df;

    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> DX   (&(sofa_dx[0][0]), sofa_dx.size()*3);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>>       DF   (&(sofa_df[0][0]), sofa_df.size()*3);

    sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield::addDForce");

    for (int k = 0; k < p_K.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_K, k); it; ++it) {
            const auto i = it.row();
            const auto j = it.col();
            const auto v = -1 * it.value() * kFactor;
            if (i != j) {
                DF[i] += v*DX[j];
                DF[j] += v*DX[i];
            } else {
                DF[i] += v*DX[i];
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield::addDForce");
}

template <typename Element>
void HyperelasticForcefield<Element>::addKToMatrix(
    sofa::defaulttype::BaseMatrix * matrix,
    SReal kFact, unsigned int & offset)
{
    if (not K_is_up_to_date) {
        update_stiffness();
    }

    sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield::addKToMatrix");

    // K is symmetric, so we only stored "one side" of the matrix.
    // But to accelerate the computation, coefficients were not
    // stored only in the upper or lower triangular part, but instead
    // in whatever triangular part (upper or lower) the first node
    // index of the element was. This means that a coefficient (i,j)
    // might be in the lower triangular part, while (k,l) is in the
    // upper triangular part. But no coefficient will be both in the
    // lower AND the upper part.

    for (int k = 0; k < p_K.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_K, k); it; ++it) {
            const auto i = it.row();
            const auto j = it.col();
            const auto v = -1 * it.value() * kFact;
            if (i != j) {
                matrix->add(offset+i, offset+j, v);
                matrix->add(offset+j, offset+i, v);
            } else {
                matrix->add(offset+i, offset+i, v);
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield::addKToMatrix");
}

template <typename Element>
SReal HyperelasticForcefield<Element>::getPotentialEnergy (
    const MechanicalParams* mparams,
    const Data<VecCoord>& d_x) const {
    SOFA_UNUSED(mparams);

    if (!this->mstate)
        return 0.;

    const auto material = d_material.get();
    if (!material) {
        return 0;
    }

    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();

    if (sofa_x.size() != sofa_x0.size() )
        return 0.;

    const auto nb_nodes = sofa_x.size();
    const auto nb_elements = number_of_elements();

    if (nb_nodes == 0 || nb_elements == 0)
        return 0;

    if (p_elements_quadrature_nodes.size() != nb_elements)
        return 0;

    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), nb_nodes, Dimension);

    SReal Psi = 0.;

    sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield::getPotentialEnergy");

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        // Fetch the node indices of the element
        const sofa::Index * node_indices = get_element_nodes_indices(element_id);

        // Fetch the initial and current positions of the element's nodes
        Matrix<NumberOfNodes, Dimension> initial_nodes_position;
        Matrix<NumberOfNodes, Dimension> current_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            initial_nodes_position.row(i).noalias() = X0.row(node_indices[i]);
            current_nodes_position.row(i).noalias() = X.row(node_indices[i]);
        }

        // Compute the nodal displacement
        Matrix<NumberOfNodes, Dimension> U {};
        for (size_t i = 0; i < NumberOfNodes; ++i) {
            const auto u = sofa_x[node_indices[i]] - sofa_x0[node_indices[i]];
            for (size_t j = 0; j < Dimension; ++j) {
                U(i, j) = u[j];
            }
        }

        // Compute the nodal forces

        for (const GaussNode & gauss_node : p_elements_quadrature_nodes[element_id]) {

            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto & dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Deformation tensor at gauss node
            const auto & F = caribou::mechanics::elasticity::strain::F(dN_dx, U);
            const auto J = F.determinant();

            // Strain tensor at gauss node
            const Mat33 C = F.transpose() * F;

            // Add the potential energy at gauss node
            Psi += (detJ * w) *  material->strain_energy_density(J, C);
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield::getPotentialEnergy");

    return Psi;
}

template <typename Element>
void HyperelasticForcefield<Element>::computeBBox(const sofa::core::ExecParams*, bool onlyVisible)
{
    if (!onlyVisible)  return;
    if (!this->mstate) return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = this->mstate->read(sofa::core::VecCoordId::position());

    static const Real max_real = std::numeric_limits<Real>::max();
    static const Real min_real = std::numeric_limits<Real>::lowest();
    Real maxBBox[3] = {min_real,min_real,min_real};
    Real minBBox[3] = {max_real,max_real,max_real};
    for (size_t i=0; i<x.size(); i++)
    {
        for (int c=0; c<3; c++)
        {
            if (x[i][c] > maxBBox[c]) maxBBox[c] = static_cast<Real>(x[i][c]);
            else if (x[i][c] < minBBox[c]) minBBox[c] = static_cast<Real>(x[i][c]);
        }
    }

    this->f_bbox.setValue(sofa::defaulttype::TBoundingBox<Real>(minBBox,maxBBox));
}

template <typename Element>
void HyperelasticForcefield<Element>::initialize_elements()
{
    sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield::initialize_elements");

    if (!this->mstate)
        return;

    // Resize the container of elements'quadrature nodes
    const auto nb_elements = number_of_elements();
    if (p_elements_quadrature_nodes.size() != nb_elements) {
        p_elements_quadrature_nodes.resize(nb_elements);
    }

    // Translate the Sofa's mechanical state vector to Eigen vector type
    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), sofa_x0.size(), Dimension);

    // Loop on each element and compute the shape functions and their derivatives for every of their integration points
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Fetch the node indices of the element
        const sofa::Index * node_indices = get_element_nodes_indices(element_id);
        Matrix<NumberOfNodes, Dimension> initial_nodes_position;

        // Fetch the initial positions of the element's nodes
        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            initial_nodes_position.row(i) = X0.row(node_indices[i]);
        }

        // Create an Element instance from the node positions
        const Element initial_element = Element(initial_nodes_position);

        // Fill in the Gauss integration nodes for this element
        p_elements_quadrature_nodes[element_id] = get_gauss_nodes(element_id, initial_element);
    }

    // Compute the volume
    Real v = 0.;
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        for (GaussNode &gauss_node : p_elements_quadrature_nodes[element_id]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            v += detJ*w;
        }
    }
    msg_info() << "Total volume of the geometry is " << v;

    sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield::initialize_elements");
}

template <typename Element>
void HyperelasticForcefield<Element>::update_stiffness()
{
    const auto material = d_material.get();
    const auto enable_multithreading = d_enable_multithreading.getValue();
    if (!material) {
        return;
    }

    // Update material parameters in case the user changed it
    material->before_update();

    static const auto Id = Mat33::Identity();
    const auto nb_elements = number_of_elements();

    const sofa::helper::ReadAccessor<Data<VecCoord>> X = this->mstate->readRestPositions();
    const auto nDofs = X.size() * 3;
    p_K.resize(nDofs, nDofs);

    ///< Triplets are used to store matrix entries before the call to 'compress'.
    /// Duplicates entries are summed up.
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(nDofs*24*2);

    sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield::update_stiffness");
    #pragma omp parallel for if (enable_multithreading)
    for (int element_id = 0; element_id < nb_elements; ++element_id) {
        // Fetch the node indices of the element
        const sofa::Index * node_indices = get_element_nodes_indices(element_id);

        using Stiffness = Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes*Dimension, NumberOfNodes*Dimension, Eigen::RowMajor>;
        Stiffness Ke = Stiffness::Zero();

        for (GaussNode &gauss_node : p_elements_quadrature_nodes[element_id]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            // Deformation tensor at gauss node
            const auto F = gauss_node.F;
            const auto J = F.determinant();

            // Right Cauchy-Green strain tensor at gauss node
            const Mat33 C = F.transpose() * F;

            // Second Piola-Kirchhoff stress tensor at gauss node
            const auto S = material->PK2_stress(J, C);

            // Jacobian of the Second Piola-Kirchhoff stress tensor at gauss node
            const auto D = material->PK2_stress_jacobian(J, C);

            // Computation of the tangent-stiffness matrix
            for (std::size_t i = 0; i < NumberOfNodes; ++i) {
                // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
                const Vec3 dxi = dN_dx.row(i).transpose();

                Matrix<6,3> Bi;
                Bi <<
                    F(0,0)*dxi[0],                 F(1,0)*dxi[0],                 F(2,0)*dxi[0],
                    F(0,1)*dxi[1],                 F(1,1)*dxi[1],                 F(2,1)*dxi[1],
                    F(0,2)*dxi[2],                 F(1,2)*dxi[2],                 F(2,2)*dxi[2],
                    F(0,0)*dxi[1] + F(0,1)*dxi[0], F(1,0)*dxi[1] + F(1,1)*dxi[0], F(2,0)*dxi[1] + F(2,1)*dxi[0],
                    F(0,1)*dxi[2] + F(0,2)*dxi[1], F(1,1)*dxi[2] + F(1,2)*dxi[1], F(2,1)*dxi[2] + F(2,2)*dxi[1],
                    F(0,0)*dxi[2] + F(0,2)*dxi[0], F(1,0)*dxi[2] + F(1,2)*dxi[0], F(2,0)*dxi[2] + F(2,2)*dxi[0];

                // The 3x3 sub-matrix Kii is symmetric, we only store its upper triangular part
                Mat33 Kii = (dxi.dot(S*dxi)*Id + Bi.transpose()*D*Bi) * detJ * w;
                Ke.template block<Dimension, Dimension>(i*Dimension, i*Dimension)
                    .template triangularView<Eigen::Upper>()
                    += Kii;

                // We now loop only on the upper triangular part of the
                // element stiffness matrix Ke since it is symmetric
                for (std::size_t j = i+1; j < NumberOfNodes; ++j) {
                    // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                    const Vec3 dxj = dN_dx.row(j).transpose();

                    Matrix<6,3> Bj;
                    Bj <<
                        F(0,0)*dxj[0],                 F(1,0)*dxj[0],                 F(2,0)*dxj[0],
                        F(0,1)*dxj[1],                 F(1,1)*dxj[1],                 F(2,1)*dxj[1],
                        F(0,2)*dxj[2],                 F(1,2)*dxj[2],                 F(2,2)*dxj[2],
                        F(0,0)*dxj[1] + F(0,1)*dxj[0], F(1,0)*dxj[1] + F(1,1)*dxj[0], F(2,0)*dxj[1] + F(2,1)*dxj[0],
                        F(0,1)*dxj[2] + F(0,2)*dxj[1], F(1,1)*dxj[2] + F(1,2)*dxj[1], F(2,1)*dxj[2] + F(2,2)*dxj[1],
                        F(0,0)*dxj[2] + F(0,2)*dxj[0], F(1,0)*dxj[2] + F(1,2)*dxj[0], F(2,0)*dxj[2] + F(2,2)*dxj[0];

                    // The 3x3 sub-matrix Kij is NOT symmetric, we store its full part
                    Mat33 Kij = (dxi.dot(S*dxj)*Id + Bi.transpose()*D*Bj) * detJ * w;
                    Ke.template block<Dimension, Dimension>(i*Dimension, j*Dimension)
                        .noalias() += Kij;
                }
            }
        }

        #pragma omp critical
        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            // Node index of the ith node in the global stiffness matrix
            const auto x = (node_indices[i]*Dimension);
            for (std::size_t m = 0; m < Dimension; ++m) {
                for (size_t n = m; n < Dimension; ++n) {
                    triplets.emplace_back(x+m, x+n, Ke(i*Dimension+m,i*Dimension+n));
                }
            }

            for (std::size_t j = i+1; j < NumberOfNodes; ++j) {
                // Node index of the jth node in the global stiffness matrix
                const auto y = (node_indices[j]*Dimension);
                for (std::size_t m = 0; m < Dimension; ++m) {
                    for (size_t n = 0; n < Dimension; ++n) {
                        triplets.emplace_back(x+m, y+n, Ke(i*Dimension+m,j*Dimension+n));
                    }
                }
            }
        }
    }
    p_K.setFromTriplets(triplets.begin(), triplets.end());
    sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield::update_stiffness");

    K_is_up_to_date = true;
    eigenvalues_are_up_to_date = false;
}

template <typename Element>
auto HyperelasticForcefield<Element>::get_gauss_nodes(const std::size_t & /*element_id*/, const Element & element) const -> GaussContainer {
    GaussContainer gauss_nodes {};
    if constexpr (NumberOfGaussNodes == caribou::Dynamic) {
        gauss_nodes.resize(element.number_of_gauss_nodes());
    }

    for (std::size_t gauss_node_id = 0; gauss_node_id < element.number_of_gauss_nodes(); ++gauss_node_id) {
        const auto & g = element.gauss_node(gauss_node_id);

        const auto J = element.jacobian(g.position);
        const Mat33 Jinv = J.inverse();
        const auto detJ = J.determinant();

        // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
        const Matrix<NumberOfNodes, Dimension> dN_dx =
            (Jinv.transpose() * element.dL(g.position).transpose()).transpose();


        GaussNode & gauss_node = gauss_nodes[gauss_node_id];
        gauss_node.weight               = g.weight;
        gauss_node.jacobian_determinant = detJ;
        gauss_node.dN_dx                = dN_dx;
        gauss_node.F = Mat33::Identity();
    }

    return gauss_nodes;
}

template <typename Element>
auto HyperelasticForcefield<Element>::eigenvalues() -> const Vector<Eigen::Dynamic> & {
    if (not eigenvalues_are_up_to_date) {
#ifdef EIGEN_USE_LAPACKE
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> k (K());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> eigensolver(k, Eigen::EigenvaluesOnly);
#else
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<Real>> eigensolver(K(), Eigen::EigenvaluesOnly);
#endif
        if (eigensolver.info() != Eigen::Success) {
            msg_error() << "Unable to find the eigen values of K.";
        }

        p_eigenvalues = eigensolver.eigenvalues();
        eigenvalues_are_up_to_date = true;
    }

    return p_eigenvalues;
}

template <typename Element>
auto HyperelasticForcefield<Element>::cond() -> Real {
    const auto & values = eigenvalues();
    const auto min = values.minCoeff();
    const auto max = values.maxCoeff();

    return min/max;
}

static const unsigned long long kelly_colors_hex[] =
{
    0xFFFFB300, // Vivid Yellow
    0xFF803E75, // Strong Purple
    0xFFFF6800, // Vivid Orange
    0xFFA6BDD7, // Very Light Blue
    0xFFC10020, // Vivid Red
    0xFFCEA262, // Grayish Yellow
    0xFF817066, // Medium Gray

    // The following don't work well for people with defective color vision
    0xFF007D34, // Vivid Green
    0xFFF6768E, // Strong Purplish Pink
    0xFF00538A, // Strong Blue
    0xFFFF7A5C, // Strong Yellowish Pink
    0xFF53377A, // Strong Violet
    0xFFFF8E00, // Vivid Orange Yellow
    0xFFB32851, // Strong Purplish Red
    0xFFF4C800, // Vivid Greenish Yellow
    0xFF7F180D, // Strong Reddish Brown
    0xFF93AA00, // Vivid Yellowish Green
    0xFF593315, // Deep Yellowish Brown
    0xFFF13A13, // Vivid Reddish Orange
    0xFF232C16, // Dark Olive Green
};

template <typename Element>
void HyperelasticForcefield<Element>::draw(const sofa::core::visual::VisualParams *vparams) {
    using Color = sofa::helper::types::RGBAColor;
    using Face = typename caribou::geometry::traits<Element>::BoundaryElementType;
    constexpr static auto NumberOfFaces = caribou::geometry::traits<Element>::NumberOfBoundaryElementsAtCompileTime;
    constexpr static auto NumberOfNodesPerFaces = caribou::geometry::traits<Face>::NumberOfNodesAtCompileTime;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    const auto nb_elements = number_of_elements();

    if (nb_elements == 0)
        return;

    vparams->drawTool()->saveLastState();
    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);
    vparams->drawTool()->disableLighting();

    const double & scale = d_drawScale.getValue();
    const VecCoord& sofa_x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.data()->data(),  sofa_x.size(), Dimension);

    std::vector< sofa::defaulttype::Vec<Dimension, Real> > faces_points[NumberOfFaces];
    for (std::size_t face_id = 0; face_id < NumberOfFaces; ++face_id) {
        faces_points[face_id].reserve(nb_elements*NumberOfNodesPerFaces);
    }

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        // Fetch the node indices of the element
        const sofa::Index * node_indices = get_element_nodes_indices(element_id);
        Matrix<NumberOfNodes, Dimension> element_nodes_position;

        // Fetch the initial positions of the element's nodes
        for (std::size_t node_id = 0; node_id < NumberOfNodes; ++node_id) {
            element_nodes_position.row(node_id) = X.row(node_indices[node_id]);
        }

        // Create an Element instance from the node positions
        const Element e = Element(element_nodes_position);

        // Scale down the element around its center point
        const auto c = e.center();
        for (std::size_t node_id = 0; node_id < NumberOfNodes; ++node_id) {
            const auto & p = element_nodes_position.row(node_id).transpose();
            element_nodes_position.row(node_id) = (c + (p - c)*scale).transpose();
        }

        // Push the faces scaled-down nodes
        const auto face_node_indices = e.boundary_elements_node_indices();
        for (std::size_t face_id = 0; face_id < NumberOfFaces; ++face_id) {
            for (std::size_t face_node_id = 0; face_node_id < NumberOfNodesPerFaces; ++face_node_id) {
                const auto & p = element_nodes_position.row(face_node_indices[face_id][face_node_id]);
                if constexpr (Dimension == 2) {
                    faces_points[face_id].emplace_back(p[0], p[1]);
                } else if constexpr (Dimension == 3) {
                    faces_points[face_id].emplace_back(p[0], p[1], p[2]);
                }
            }
        }
    }

    for (std::size_t face_id = 0; face_id < NumberOfFaces; ++face_id) {
        const auto &hex_color = kelly_colors_hex[face_id % 20];
        const Color face_color(
            static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(16))) / 255.),
            static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(8) )) / 255.),
            static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(0) )) / 255.),
            static_cast<float> (1));

        if (NumberOfNodesPerFaces == 3) {
            vparams->drawTool()->drawTriangles(faces_points[face_id], face_color);
        } else if (NumberOfNodesPerFaces == 4) {
            vparams->drawTool()->drawQuads(faces_points[face_id], face_color);
        } else {
            throw std::runtime_error("Drawing an unsupported face type");
        }
    }

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

template <typename Element>
auto HyperelasticForcefield<Element>::canCreate(HyperelasticForcefield<Element>* o, BaseContext* context, BaseObjectDescription* arg) -> bool {
    std::string requested_element_type = arg->getAttribute( "template", "");
    std::string this_element_type = templateName(o);

    // to lower
    std::string requested_element_type_lower = requested_element_type, this_element_type_lower = this_element_type;
    std::transform(requested_element_type.begin(), requested_element_type.end(), requested_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(this_element_type.begin(), this_element_type.end(), this_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });

    if (requested_element_type_lower == this_element_type_lower) {
        return Inherit::canCreate(o, context, arg);
    }

    if (not requested_element_type.empty()) {
        arg->logError("Requested element type is not '"+templateName(o)+"'.");
        return false;
    }

    std::string topology_path = arg->getAttribute("topology", "");
    if (not topology_path.empty()) {
        topology_path = topology_path.substr(1); // removes the "@"
        // Make sure the specified topology has elements of type Element
        auto topology = context->get<sofa::core::topology::BaseMeshTopology>(topology_path);
        if (not topology or not mesh_is_compatible(topology)) {
            arg->logError("Cannot deduce the element type from the specified mesh topology '" + topology_path + "'. Add template=\""+this_element_type+"\" to use it.");
            return false;
        }
    } else {
        // Try to find a compatible topology in the current context
        sofa::core::topology::BaseMeshTopology * topology = nullptr;
        const auto topologies = context->getObjects<sofa::core::topology::BaseMeshTopology>(
            BaseContext::SearchDirection::Local);
        for (const auto t : topologies) {
            if (mesh_is_compatible(t)) {
                topology = t;
                break;
            }
        }
        if (not topology) {
            arg->logError("Cannot find a topology in the current context from which the template '"+this_element_type+"' can be deduced.");
            return false;
        }

        if (Inherit::canCreate(o, context, arg)) {
            arg->setAttribute("topology", "@" + topology->getPathName());
            return true;
        } else {
            return false;
        }
    }

    return Inherit::canCreate(o, context, arg);
}

} // namespace SofaCaribou::forcefield
