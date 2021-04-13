#pragma once

#include <sofa/version.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace SofaCaribou::topology {

template <typename DataTypes>
class CaribouMesh : public sofa::core::objectmodel::BaseObject {
public:
    SOFA_CLASS(SOFA_TEMPLATE(CaribouMesh, DataTypes), BaseObject);
};

} // namespace SofaCaribou::topology