#ifndef CARIBOU_ENGINE_CUTGRIDENGINE_H
#define CARIBOU_ENGINE_CUTGRIDENGINE_H

#include <SofaBaseTopology/MeshTopology.h>
#include <SofaBaseTopology/GridTopology.h>
#include <sofa/core/DataEngine.h>
#include <sofa/helper/OptionsGroup.h>

namespace sofa
{

using namespace core;
using namespace core::objectmodel;
using namespace core::behavior;
using namespace component::topology;

namespace caribou
{

namespace engine
{

class CutGridEngine : public DataEngine
{
public:
    SOFA_CLASS(CutGridEngine, DataEngine);

    //### TYPE DEFINITIONS ### <editor-fold desc="Types definitions">
    using Inherit = DataEngine;
    using DataTypes =  sofa::defaulttype::Vec3dTypes;
    using Real = typename DataTypes::Real;
    using Coord = typename DataTypes::Coord;
    using VecCoord = typename DataTypes::VecCoord;
    using VecReal = typename DataTypes::VecReal;
    using HexaID = size_t;
    using Hexa = sofa::core::topology::Topology::Hexa;
    using Triangle = sofa::core::topology::Topology::Triangle;
    using PointID = sofa::core::topology::Topology::PointID;
    using Color = sofa::defaulttype::Vec4f;

    template <typename T>
    using Link = SingleLink<CutGridEngine, T, BaseLink::FLAG_STRONGLINK>;
    //######################## </editor-fold>

    enum class Flag {
        Outside = -1,
        Boundary = 0,
        Inside = 1
    };

    void init() override;
    void update() override;
    void draw(const core::visual::VisualParams* vparams) override;

    /// Input stream
    inline friend std::istream& operator>> ( std::istream& in, Flag& flag )
    {
        std::underlying_type<Flag>::type val;
        in >> val;
        flag = static_cast<Flag>(val);
        return in;
    }

    /// Output stream
    inline friend std::ostream& operator << (std::ostream& os, const Flag& flag)
    {
        os << static_cast<std::underlying_type<Flag>::type>(flag);
        return os;
    }

protected:
    CutGridEngine();

    // Inputs
    Link<GridTopology> d_grid_topology;
    Link<MeshTopology> d_surface_topology;
    Data<sofa::helper::OptionsGroup> d_showHexahedrons;
    Data<bool> d_showTriangles;
    Data<Color> d_showInsideColor;
    Data<Color> d_showOutsideColor;
    Data<Color> d_showBoundaryColor;
    Data<Color> d_showTrianglesColor;

    // Outputs
    Data<sofa::helper::vector<Flag>> d_points_flags;
    Data<sofa::helper::vector<Flag>> d_hexahedrons_flags;
    Data<sofa::helper::vector<Coord>> d_triangle_positions;
    Data<sofa::helper::vector<sofa::helper::vector<Triangle>>> d_triangles;

private:

//    Coord isovalue(const Coord & p) const;
    Flag getFlag(const Coord & p) const;
    Flag getFlag(const Hexa & hexa) const;

};

} // namespace engine

} // namespace caribou

} // namespace sofa
#endif //CARIBOU_ENGINE_CUTGRIDENGINE_H
