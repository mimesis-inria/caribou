#include "PressureForcefield.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/AdvancedTimer.h>

#include "../Helper/Triangle.h"

#define EPSILON 0.000000001

namespace sofa
{

namespace caribou
{

namespace forcefield
{

template<class DataTypes>
PressureForcefield<DataTypes>::PressureForcefield()
    // Inputs
    : d_pressure(initData(&d_pressure, "pressure", "Pressure force per unit area"))
    , d_triangles(initData(&d_triangles, "triangles", "List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...]) "))
    , d_slope(initData(&d_slope, (Real) 0, "slope",
                       "Slope of pressure increment, the resulting pressure will be p^t = p^{t-1} + p*slope. "
                       "If slope = 0, the pressure will be constant"))
    , d_triangleContainer(
        initLink("triangle_container", "Triangle set topology container that contains the triangle indices"))
    , d_mechanicalState(
        initLink("state", "Mechanical state that contains the triangle positions"))
    , d_number_of_steps_before_increment(
        initData(&d_number_of_steps_before_increment, (unsigned int) 1, "number_of_steps_before_increment",
                 "Number of steps to wait before adding an increment. This can be used to simulate a Newton-Raphson solver."))

    // Outputs
    , d_nodal_forces(
        initData(&d_nodal_forces, "nodal_forces", "Current nodal forces from the applied pressure", true, true))
    , d_current_load(initData(&d_current_load, (Real) 0, "current_load", "Current total load applied"))
{
    this->f_listening.setValue(true);
}

template<class DataTypes>
void PressureForcefield<DataTypes>::init()
{
    core::behavior::ForceField<DataTypes>::init();

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
        m_pressure_is_constant = true;

    nodal_forces.resize(rest_positions.size());

    // Force an increment on the first step of the simulation
    m_number_of_steps_since_last_increment = d_number_of_steps_before_increment.getValue();

    if (m_pressure_is_constant)
        addNodalPressures(d_pressure.getValue());
}

template<class DataTypes>
void PressureForcefield<DataTypes>::reset()
{
    core::behavior::ForceField<DataTypes>::reset();
    // Force an increment on the first step of the simulation
    m_number_of_steps_since_last_increment = d_number_of_steps_before_increment.getValue();

    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;

    nodal_forces.clear();
    d_current_load.setValue(0.);
    m_current_pressure = Deriv();
}

template<class DataTypes>
void PressureForcefield<DataTypes>::handleEvent(sofa::core::objectmodel::Event* event)
{
    if (!sofa::simulation::AnimateBeginEvent::checkEventType(event))
        return;

    if (m_pressure_is_constant)
        return;

    const Deriv & pressure = d_pressure.getValue();
    const Real &  slope = d_slope.getValue();
    const unsigned int number_of_steps_before_increment = d_number_of_steps_before_increment.getValue();
    unsigned int & number_of_steps_since_last_increment = m_number_of_steps_since_last_increment;

    // If we've reach the total threshold, stop adding pressure
    if (m_current_pressure.norm() >= pressure.norm())
        return;

    // Update the counter of steps since the last increment was done
    number_of_steps_since_last_increment++;

    if (number_of_steps_since_last_increment < number_of_steps_before_increment)
        return;

    // Reset the counter as we are doing an increment right now
    number_of_steps_since_last_increment = 0;

    Deriv increment = pressure * slope;
    Deriv pressure_to_apply = m_current_pressure + increment;

    if (pressure_to_apply.norm() > pressure.norm())
        increment = pressure - m_current_pressure;

    m_current_pressure += increment;

    addNodalPressures(increment);
}

template<class DataTypes>
void PressureForcefield<DataTypes>::addNodalPressures(Deriv pressure)
{
    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;
    sofa::helper::WriteAccessor<Data<Real>> current_load = d_current_load;

    const auto & triangles = d_triangles.getValue();
    const auto rest_positions = d_mechanicalState.get()->readRestPositions();

    msg_info() << "Adding a pressure increment of " << pressure;

    Deriv load;
    Real total_area = 0;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto & triangle = triangles[i];
        const auto & p1 = rest_positions[triangle[0]];
        const auto & p2 = rest_positions[triangle[1]];
        const auto & p3 = rest_positions[triangle[2]];

        const Real area = caribou::helper::triangle::area(p1, p2, p3);
        const Deriv force = area * pressure;

        total_area += area;

        nodal_forces[triangle[0]] += force / (Real) 3;
        nodal_forces[triangle[1]] += force / (Real) 3;
        nodal_forces[triangle[2]] += force / (Real) 3;

        load += force;
    }

    current_load += load.norm();
}

template<class DataTypes>
void PressureForcefield<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& /*d_v*/)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_x);
    sofa::helper::AdvancedTimer::stepBegin("PressureForcefield::addForce");
    sofa::helper::ReadAccessor<Data<Real>> current_load = d_current_load;
    sofa::helper::ReadAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    for (size_t i = 0; i < f.size(); ++i)
        f[i] += nodal_forces[i];

    sofa::helper::AdvancedTimer::valSet("load", current_load);
    sofa::helper::AdvancedTimer::stepEnd("PressureForcefield::addForce");
}

template<class DataTypes>
void PressureForcefield<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    using Color = core::visual::DrawTool::RGBAColor;
    using Vector3 = core::visual::DrawTool::Vector3;

    if (! vparams->displayFlags().getShowForceFields())
        return;

    const auto & triangles = d_triangles.getValue();
    const auto positions = d_mechanicalState.get()->readPositions();

    sofa::helper::vector<Vector3> points (triangles.size() * 2);

    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto & triangle = triangles[i];
        const auto & p1 = positions[triangle[0]];
        const auto & p2 = positions[triangle[1]];
        const auto & p3 = positions[triangle[2]];
        const auto   center = caribou::helper::triangle::center(p1, p2, p3);
        const auto   normal = caribou::helper::triangle::normal(p1, p2, p3);

        points[2*i] = center;
        points[2*i + 1] = center + normal;
    }

    vparams->drawTool()->drawLines(points, 1.f, Color(1, 0, 0, 1));
}


SOFA_DECL_CLASS(PressureForcefield)
static int PressureForceFieldClass = core::RegisterObject("TrianglePressure")
                                          .add< PressureForcefield<sofa::defaulttype::Vec3dTypes> >(true)
;


} // namespace forcefield

} // namespace caribou

} // namespace sofa