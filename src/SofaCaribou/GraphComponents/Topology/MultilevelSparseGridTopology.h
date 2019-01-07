#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_H
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/DDGNode.h>

#include <Caribou/config.h>
#include <Caribou/Topology/Engine/Grid/Grid.h>
#include <Caribou/Topology/Engine/Grid/Cell.h>

#include <array>
#include <memory>

namespace SofaCaribou {
namespace GraphComponents {
namespace topology {

using namespace sofa::core::objectmodel;

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

template <typename VecType>
class MultilevelSparseGridTopology : public virtual BaseObject, public virtual DDGNode
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(MultilevelSparseGridTopology, VecType), BaseObject);

    static constexpr char Dimension = VecType::spatial_dimensions;
    using CellType = caribou::topology::engine::Cell<Dimension>;
    using GridType = caribou::topology::engine::Grid<3, CellType>;
    using Index = size_t;
    using VecFloat = caribou::algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using VecInt   = caribou::algebra::Vector<Dimension, Index>;
    using Int   = typename VecInt::ValueType;
    using Float = typename VecFloat::ValueType;

    using SofaVecInt = sofa::defaulttype::Vec<VecType::spatial_dimensions>;

    using BaseObject::className;
    using BaseObject::templateName;
    using BaseObject::namespaceName;
    using BaseObject::shortName;
    using BaseObject::dynamicCast;

    MultilevelSparseGridTopology();
    void init() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;



    void update() override {};
    void setDirtyValue(const sofa::core::ExecParams* params = 0) override {
        onUpdate();

        // the input needs to be inform their output
        // are not dirty, to be sure they will call setDirtyValue when they are modified
        cleanDirtyOutputsOfInputs(params);
    };

    /// This method is needed by DDGNode
    Base* getOwner() const override { return nullptr; };

    /// This method is needed by DDGNode
    BaseData* getData() const override { return nullptr; };

    /// This method is needed to avoid ambiguity between DDGNode and BaseObject
    const std::string & getName() const override {
        return BaseObject::getName();
    }

private:
    void onUpdate();

private:
    Data<SofaVecInt> d_n;
    Data<unsigned char> d_number_of_subdivision;

    std::unique_ptr<GridType> p_grid;

};

extern template class MultilevelSparseGridTopology<Vec2Types>;
extern template class MultilevelSparseGridTopology<Vec3Types>;

} // namespace topology
} // namespace GraphComponents
} // namespace SofaCaribou

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_H
