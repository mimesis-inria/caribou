#include "IBMForcefield.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa
{

using namespace core::objectmodel;
using namespace core::behavior;
using namespace component::topology;
using namespace defaulttype;

namespace caribou
{

namespace forcefield
{

template<class DataTypes>
IBMForcefield<DataTypes>::IBMForcefield()
        : d_youngModulus(initData(&d_youngModulus, Real(1000), "youngModulus", "[IN] Young's modulus of the material"))
        , d_poissonRatio(initData(&d_poissonRatio, Real(0.3),  "poissonRatio", "[IN] Poisson's ratio of the material"))
        , d_points_flags(initData(&d_points_flags, "points_flags", "[IN] Flags (Outside -1, Boundary 0 or Inside 1) of the points"))
        , d_hexahedrons_flags(initData(&d_hexahedrons_flags, "hexahedrons_flags", "[IN] Flags (Outside -1, Boundary 0 or Inside 1) of the hexahedrons"))
        , d_container(initLink("container", "Topology container that contains the hexahedrons"))
{
}

template<class DataTypes>
void IBMForcefield<DataTypes>::init()
{
    Inherit::init();

    if (!this->getMState()) {
        msg_error() << "No mechanical state provided.";
        return;
    }

    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Flag>>> points_flags = d_points_flags;
    sofa::helper::ReadAccessor<Data<VecCoord>> positions = this->getMState()->readPositions();

    if (points_flags.size() != positions.size()) {
        msg_error() << "The points flags must match exactly the number of nodes in the mechanical object.";
        return;
    }
}

template<class DataTypes>
void IBMForcefield<DataTypes>::reset()
{
    Inherit::reset();
}

template<class DataTypes>
void IBMForcefield<DataTypes>::addForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_f);
    SOFA_UNUSED(d_x);
    SOFA_UNUSED(d_v);
}

template<class DataTypes>
void IBMForcefield<DataTypes>::addDForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_df, const Data<VecDeriv>& d_dx)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_df);
    SOFA_UNUSED(d_dx);
}

template<class DataTypes>
void IBMForcefield<DataTypes>::addKToMatrix(BaseMatrix * matrix, SReal kFact, unsigned int &offset)
{
    SOFA_UNUSED(offset);
    SOFA_UNUSED(matrix);
    SOFA_UNUSED(kFact);
}

template<class DataTypes>
void IBMForcefield<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    SOFA_UNUSED(vparams);
}

SOFA_DECL_CLASS(IBMForcefield)
static int IBMForcefieldClass = core::RegisterObject("Caribou IBM Forcefield")
        .add< IBMForcefield<sofa::defaulttype::Vec3dTypes> >(true)
;

} // namespace forcefield

} // namespace caribou

} // namespace sofa