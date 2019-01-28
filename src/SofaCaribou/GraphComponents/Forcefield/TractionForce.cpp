#include "TractionForce.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/AdvancedTimer.h>

#include <Caribou/config.h>
#include <Caribou/Geometry/Triangle.h>
#include <SofaCaribou/Traits.h>

namespace SofaCaribou {
namespace GraphComponents {
namespace forcefield {

template<class DataTypes>
TractionForce<DataTypes>::TractionForce()
    // Inputs
    : d_traction(initData(&d_traction,
            "traction",
            "Tractive force per unit area (if an incremental load is set by the slope parameter, this is the final load "
            "reached after all increments)."))
    , d_triangles(initData(&d_triangles,
            "triangles",
            "List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...])."))
    , d_slope(initData(&d_slope,
            (Real) 0,
            "slope",
            "Slope of load increment, the resulting tractive force will be p^t = p^{t-1} + p*slope. "
            "If slope = 0, the traction will be constant."))
    , d_triangleContainer(initLink(
            "triangle_container",
            "Triangle set topology container that contains the triangle indices."))
    , d_mechanicalState(initLink(
            "state",
            "Mechanical state that contains the positions of the surface elements."))
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

template<class DataTypes>
void TractionForce<DataTypes>::init()
{
    // If no triangle container specified, but a link is set for the triangles, use the linked container
    if ( !  d_triangleContainer.get()
         && d_triangles.getParent()
         && dynamic_cast<sofa::component::topology::TriangleSetTopologyContainer*> (d_triangles.getParent()->getOwner())) {

        d_triangleContainer.set(dynamic_cast<sofa::component::topology::TriangleSetTopologyContainer*> (d_triangles.getParent()->getOwner()));
    }

    // If no triangle indices set, get the ones from the container (if one was provided, or automatically found in the context)
    if (d_triangles.getValue().empty()) {
        if (!d_triangleContainer.get()) {
            auto container = this->getContext()->template get<sofa::component::topology::TriangleSetTopologyContainer>();
            if (container) {
                d_triangleContainer.set(container);
                d_triangles.setParent(&d_triangleContainer.get()->getTriangleDataArray());
                if (d_triangles.getValue().empty()) {
                    msg_error() << "A triangle topology container was found, but contained an empty set of triangles.";
                }
            } else
                msg_error() << "A set of triangles or a triangle topology containing a set of triangles must be provided.";
        }
    }

    // If no mechanical state specified, try to find one in the context
    if (! d_mechanicalState.get()) {
        auto state = this->getContext()->template get<sofa::core::behavior::MechanicalState<DataTypes>>();
        if (state)
            d_mechanicalState.set(state);
        else
            msg_error() << "No mechanical state provided or found in the current context.";
    }


    // Initialization from the data inputs
    // @todo(jnbrunet2000@gmail.com) : this part should be moved to a real data update function
    const auto rest_positions = d_mechanicalState.get()->readRestPositions();
    sofa::helper::WriteOnlyAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;

    if (d_slope.getValue() > -EPSILON && d_slope.getValue() < EPSILON)
        m_traction_is_constant = true;

    nodal_forces.resize(rest_positions.size());

    // Do an increment on the first step of the simulation
    m_number_of_steps_since_last_increment = d_number_of_steps_before_increment.getValue();

    // If all the load is applied at once, increment now the load since it will be used for the computation of the first time step
    if (m_traction_is_constant)
        increment_load(d_traction.getValue());
}

template<class DataTypes>
void TractionForce<DataTypes>::reset()
{
    // Do an increment on the first step of the simulation
    m_number_of_steps_since_last_increment = d_number_of_steps_before_increment.getValue();

    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;

    nodal_forces.clear();
    d_total_load.setValue(0.);
    m_current_traction = Deriv();
}

template<class DataTypes>
void TractionForce<DataTypes>::handleEvent(sofa::core::objectmodel::Event* event)
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

template<class DataTypes>
void TractionForce<DataTypes>::increment_load(Deriv traction_increment_per_unit_area)
{
    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;
    sofa::helper::WriteAccessor<Data<Real>> current_load = d_total_load;

    const auto & triangles = d_triangles.getValue();
    const auto rest_positions = d_mechanicalState.get()->readRestPositions();

    msg_info() << "Incrementing the load by " << traction_increment_per_unit_area << " (tractive force per unit area).";

    Deriv load;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto & triangle_node_indices = triangles[i];
        const auto & p1 = rest_positions[triangle_node_indices[0]];
        const auto & p2 = rest_positions[triangle_node_indices[1]];
        const auto & p3 = rest_positions[triangle_node_indices[2]];

        const auto triangle = caribou::geometry::Triangle<3>(p1, p2, p3);

        // Integration of the traction increment over the element.
        const auto area = triangle.area();
        const Deriv integrated_force_increment = area * traction_increment_per_unit_area;

        nodal_forces[triangle_node_indices[0]] += integrated_force_increment / (Real) 3;
        nodal_forces[triangle_node_indices[1]] += integrated_force_increment / (Real) 3;
        nodal_forces[triangle_node_indices[2]] += integrated_force_increment / (Real) 3;

        load += integrated_force_increment;
    }

    current_load += load.norm();
}

template<class DataTypes>
void TractionForce<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& /*d_v*/)
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

template<class DataTypes>
void TractionForce<DataTypes>::draw(const sofa::core::visual::VisualParams* vparams)
{
    using Color = sofa::core::visual::DrawTool::RGBAColor;
    using Vector3 = sofa::core::visual::DrawTool::Vector3;

    if (! vparams->displayFlags().getShowForceFields())
        return;

    const auto & triangles = d_triangles.getValue();
    const auto positions = d_mechanicalState.get()->readPositions();

    sofa::helper::vector<Vector3> points (triangles.size() * 2);

    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto & triangle_node_indices = triangles[i];
        const auto & p1 = positions[triangle_node_indices[0]];
        const auto & p2 = positions[triangle_node_indices[1]];
        const auto & p3 = positions[triangle_node_indices[2]];

        const auto triangle = caribou::geometry::Triangle<3>(p1, p2, p3);

        const auto c = triangle.center();
        const Vector3 center(c[0], c[1], c[2]);

        const auto n = triangle.normal();
        const Vector3 normal(n[0], n[1], n[2]);

        points[2*i] = center;
        points[2*i + 1] = center + normal;
    }

    vparams->drawTool()->drawLines(points, 1.f, Color(1, 0, 0, 1));
}


SOFA_DECL_CLASS(TractionForce)
static int TractionForceClass = sofa::core::RegisterObject("Traction forcefield.")
                                          .add< TractionForce<sofa::defaulttype::Vec3dTypes> >(true)
;

} // namespace forcefield
} // namespace GraphComponents
} // namespace SofaCaribou