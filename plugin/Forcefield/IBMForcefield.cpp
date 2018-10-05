#include "IBMForcefield.h"
#include "../Helper/MirtichIntegration.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace sofa
{

using namespace core::objectmodel;
using namespace core::behavior;
using namespace core::topology;
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
        , d_hexahedrons(initData(&d_hexahedrons, "hexahedrons",
                "[IN] List of hexahedrons [n1H1, n2H1, n3H1, ..., n8Hm] where niHj is the node id of vertex i (0 to 7) of the hexahedron j (0 to m)."))
        , d_hexahedrons_flags(initData(&d_hexahedrons_flags, "hexahedrons_flags", "[IN] Flags (Outside -1, Boundary 0 or Inside 1) of the hexahedrons"))
        , d_triangle_positions(initData(&d_triangle_positions, "triangle_positions", "[IN] List of triangles vertices positions."))
        , d_triangles(initData(&d_triangles, "triangles",
                "[IN] List of triangles per hexahedron (they should represent both the cut surface and the hexa's subsurfaces under the cut).", false))
{
            f_listening.setValue(true);
}

template<class DataTypes>
void IBMForcefield<DataTypes>::init()
{
    Inherit::init();

    if (!this->getMState()) {
        msg_error() << "No mechanical state provided.";
        return;
    }

    // If no link to an hexahedron list is set and no hexahedron were manually added to the input, try to find a
    // list of hexahedron in the current context.
    if (!d_hexahedrons.getParent() && d_hexahedrons.getValue().empty()) {
        auto topology = this->getContext()->template get<BaseMeshTopology>();
        if (topology && topology->findData("hexahedra")) {
            d_hexahedrons.setParent(topology->findData("hexahedra"));
            msg_info() << "Hexahedron list automatically linked to the mesh topology '"
                << topology->getPathName() << "'";
            if (topology->getNbHexahedra() == 0) {
                msg_warning()
                    << "Hexahedron list was automatically linked to the mesh topology '"
                    << topology->getPathName() << "' but it contain an empty set of hexahedrons.";
            }
        } else {
            msg_error() << "No hexahedrons list was provided and no mesh topology was found in the current context.";
        }
    }

    // If no hexahedron flags set, try to find  CutGridEngine in the current context
    if (!d_hexahedrons_flags.getParent() && d_hexahedrons_flags.getValue().empty()) {
        auto engine = this->getContext()->template get<engine::CutGridEngine>();
        if (engine && engine->findData("hexahedrons_flags")) {
            d_hexahedrons_flags.setParent(engine->findData("hexahedrons_flags"));
            msg_info() << "Hexahedron flag list automatically linked to the cut grid engine '"
                << engine->getPathName() << "'";
        }
    }

    // If no triangles set, try to find a CutGridEngine in the current context
    if (!d_triangle_positions.getParent() || !d_triangles.getParent()) {
        auto engine = this->getContext()->template get<engine::CutGridEngine>();
        if (engine && engine->findData("triangle_positions") && engine->findData("triangles")) {
            d_triangle_positions.setParent(engine->findData("triangle_positions"));
            d_triangles.setParent(engine->findData("triangles"));
            msg_info() << "Hexahedron triangles list automatically linked to the cut grid engine '"
                       << engine->getPathName() << "'";
        }
    }

    // Compute the total volume
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Hexahedron>>> hexahedrons = d_hexahedrons;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Flag>>> hexahedrons_flags = d_hexahedrons_flags;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Coord>>> triangle_positions = d_triangle_positions;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<sofa::helper::vector<Triangle>>>> triangles = d_triangles;
    sofa::helper::ReadAccessor<Data<VecCoord>> x = this->getMState()->readRestPositions();

    Real mass = 0.;
    for (HexahedronID i = 0; i < hexahedrons.size(); ++i) {
        if (hexahedrons_flags[i] == Flag::Outside)
            continue;

        if (triangles[i].empty()) {
            std::array<Coord, 8> nodes = {{
                x[hexahedrons[i][0]], x[hexahedrons[i][1]], x[hexahedrons[i][2]], x[hexahedrons[i][3]],
                x[hexahedrons[i][4]], x[hexahedrons[i][5]], x[hexahedrons[i][6]], x[hexahedrons[i][7]]
            }};

            mass += (nodes[1]-nodes[0]).norm() * (nodes[3]-nodes[0]).norm() * (nodes[4]-nodes[0]).norm();
        } else {

        }

    }

    std::cout << "MASS IS " << std::to_string(mass) << std::endl;

}

template<class DataTypes>
void IBMForcefield<DataTypes>::reinit()
{
    Inherit::reinit();
    update();
}

template<class DataTypes>
void IBMForcefield<DataTypes>::update()
{

    if (   !d_hexahedrons_flags.isDirty()
        && !d_hexahedrons.isDirty()
        && !d_triangle_positions.isDirty()
        && !d_triangles.isDirty()) {
        cleanDirty();
        return;
    }

    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Hexahedron>>> hexahedrons = d_hexahedrons;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Flag>>> hexahedrons_flags = d_hexahedrons_flags;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Coord>>> triangle_positions = d_triangle_positions;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<sofa::helper::vector<Triangle>>>> triangles = d_triangles;

    cleanDirty();

    sofa::helper::WriteOnlyAccessor<Data<sofa::helper::vector<std::array<Real, 84>>>> integrated_monomials = d_integrated_monomials;
    integrated_monomials.resize(hexahedrons.size());




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