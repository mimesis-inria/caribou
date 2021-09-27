#include <SofaCaribou/config.h>
#include <SofaCaribou/Forcefield/TractionForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION <= 201299)
namespace sofa::helper::visual { using DrawTool = sofa::core::visual::DrawTool; }
#endif

namespace SofaCaribou::forcefield {

template <typename Element>
TractionForcefield<Element>::TractionForcefield()
    // Inputs
    : d_traction(initData(&d_traction,
            "traction",
            "Tractive force per unit area (if an incremental load is set by the slope parameter, this is the final load "
            "reached after all increments)."))
    , d_slope(initData(&d_slope,
            (Real) 0,
            "slope",
            "Slope of load increment, the resulting tractive force will be p^t = p^{t-1} + p*slope where p is the "
            "traction force passed as a data and p^t is the traction force applied at time step t. "
            "If slope = 0, the traction will be constant."))
    , d_number_of_steps_before_increment(initData(&d_number_of_steps_before_increment,
            (unsigned int) 1,
            "number_of_steps_before_increment",
            "Number of time steps to wait before adding an increment. "
            "This can be used to simulate Newton-Raphson solving process where the time steps are the Newton iterations."))

    // Outputs
    , d_nodal_forces(initData(&d_nodal_forces,
            "nodal_forces",
            "Current nodal forces from the applied traction", true, true))
    , d_total_load(initData(&d_total_load,
            (Real) 0,
            "total_load",
            "Accumulated load applied on all the surface area."))
{
    this->f_listening.setValue(true);
}

template <typename Element>
void TractionForcefield<Element>::init()
{
    Inherit::init();

    // Initialization from the data inputs
    // @todo(jnbrunet2000@gmail.com) : this part should be moved to a real data update function
    const auto rest_positions = this->getMState()->readRestPositions();
    sofa::helper::WriteOnlyAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;

    if (d_slope.getValue() > -EPSILON && d_slope.getValue() < EPSILON)
        m_traction_is_constant = true;

    nodal_forces.resize(rest_positions.size());

    // Do an increment on the first step of the simulation
    m_number_of_steps_since_last_increment = d_number_of_steps_before_increment.getValue();

    // Compute and store the shape functions and their derivatives for every integration points
    initialize_elements();

    // If all the load is applied at once, increment now the load since it will be used for the computation of the first time step
    if (m_traction_is_constant)
        increment_load(d_traction.getValue());
}

template <typename Element>
void TractionForcefield<Element>::reset()
{
    // Do an increment on the first step of the simulation
    m_number_of_steps_since_last_increment = d_number_of_steps_before_increment.getValue();

    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;

    nodal_forces.clear();
    d_total_load.setValue(0.);
    m_current_traction = Deriv();
}

template <typename Element>
void TractionForcefield<Element>::handleEvent(sofa::core::objectmodel::Event* event)
{
    if (!sofa::simulation::AnimateBeginEvent::checkEventType(event))
        return;

    if (m_traction_is_constant)
        return;

    const Deriv & maximum_traction_to_apply = d_traction.getValue();
    const Real &  slope = d_slope.getValue();
    const unsigned int number_of_steps_before_increment = d_number_of_steps_before_increment.getValue();
    unsigned int & number_of_steps_since_last_increment = m_number_of_steps_since_last_increment;

    // If we've reach the total traction threshold, do not increment the load again
    if (m_current_traction.norm() >= maximum_traction_to_apply.norm())
        return;

    // Update the counter of steps since the last increment was done
    number_of_steps_since_last_increment++;

    // We did not reach the amount of steps required, skip the increment
    if (number_of_steps_since_last_increment < number_of_steps_before_increment)
        return;

    // We reached the amount of steps required, reset the counter as we will do an increment right now
    number_of_steps_since_last_increment = 0;

    // Compute the increment value from the given slope
    Deriv increment = maximum_traction_to_apply * slope;

    // Add the increment to the current applied load to get an idea of the traction we would apply
    Deriv traction_to_apply = m_current_traction + increment;

    // If this traction is greater than the goal traction, then only increment with the remaining load to get exactly the
    // prescribed goal load
    if (traction_to_apply.norm() > maximum_traction_to_apply.norm())
        increment = maximum_traction_to_apply - m_current_traction;

    m_current_traction += increment;

    // Finally, apply the increment
    increment_load(increment);
}

template <typename Element>
void TractionForcefield<Element>::increment_load(Deriv traction_increment_per_unit_area)
{
    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;
    sofa::helper::WriteAccessor<Data<Real>> current_load = d_total_load;

    const auto rest_positions = this->getMState()->readRestPositions();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0 (rest_positions.ref().data()->data(),  rest_positions.size(), Dimension);

    msg_info() << "Incrementing the load by " << traction_increment_per_unit_area << " (tractive force per unit area).";

    Deriv load;
    const auto nb_elements = this->number_of_elements();
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Fetch the node indices of the element
        auto node_indices = this->topology()->domain()->element_indices(element_id);

        // Fetch the initial positions of the element's nodes
        Matrix<NumberOfNodesPerElement, Dimension> initial_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            initial_nodes_position.row(i).noalias() = X0.row(node_indices[i]);
        }

        // Integration of the traction increment over the element.
        for (GaussNode &gauss_node : p_elements_quadrature_nodes[element_id]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Shape values at each nodes evaluated at the gauss point position
            const auto N = gauss_node.N;

            // Traction evaluated at the gauss point position
            const auto F = traction_increment_per_unit_area * w * detJ;

            // Tractive forces w.r.t the gauss node applied on each nodes
            for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
                nodal_forces[node_indices[i]] += F*N[i];
                load += F*N[i];
            }
        }
    }

    current_load += load.norm();
}

template <typename Element>
void TractionForcefield<Element>::addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& /*d_v*/)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_x);
    sofa::helper::AdvancedTimer::stepBegin("TractionForce::addForce");
    sofa::helper::ReadAccessor<Data<Real>> current_load = d_total_load;
    sofa::helper::ReadAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    for (size_t i = 0; i < f.size(); ++i)
        f[i] += nodal_forces[i];

    sofa::helper::AdvancedTimer::valSet("load", current_load);
    sofa::helper::AdvancedTimer::stepEnd("TractionForce::addForce");
}

template<typename Element>
void TractionForcefield<Element>::initialize_elements() {
    using namespace sofa::core::objectmodel;

    sofa::helper::ScopedAdvancedTimer _t_ ("TractionForcefield::initialize_elements");

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

    // Compute the surface area
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
    msg_info() << "Total surface of the geometry is " << v;
}

template<typename Element>
auto
TractionForcefield<Element>::get_gauss_nodes(const size_t & /*element_id*/, const Element &element) const -> GaussContainer {
    GaussContainer gauss_nodes {};
    if constexpr (NumberOfGaussNodesPerElement == caribou::Dynamic) {
        gauss_nodes.resize(element.number_of_gauss_nodes());
    }

    const auto nb_of_gauss_nodes = gauss_nodes.size();
    for (std::size_t gauss_node_id = 0; gauss_node_id < nb_of_gauss_nodes; ++gauss_node_id) {
        const auto & g = element.gauss_node(gauss_node_id);

        const auto J = element.jacobian(g.position);
        Real detJ;
        if constexpr (Element::Dimension == 3 && Element::CanonicalDimension == 2) {
            detJ = std::abs(J.col(0).cross(J.col(1)).norm());
        } else {
            detJ = std::abs(J.determinant());
        }

        // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
        const Vector<NumberOfNodesPerElement> N = element.L(g.position);


        GaussNode & gauss_node = gauss_nodes[gauss_node_id];
        gauss_node.weight               = g.weight;
        gauss_node.jacobian_determinant = detJ;
        gauss_node.N                    = N;
    }

    return gauss_nodes;
}

} // namespace SofaCaribou::forcefield
