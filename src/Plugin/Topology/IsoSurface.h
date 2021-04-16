#pragma once

#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/objectmodel/BaseObject.h>
DISABLE_ALL_WARNINGS_END

#include <Eigen/Dense>

namespace SofaCaribou::topology {

template<class DataTypes>
class IsoSurface : public sofa::core::objectmodel::BaseObject {
public:
    SOFA_CLASS(SOFA_TEMPLATE(IsoSurface, DataTypes), sofa::core::objectmodel::BaseObject);
    static constexpr unsigned char Dimension = DataTypes::spatial_dimensions;
    using Real = typename DataTypes::Real;
    using Coord = typename DataTypes::Coord;

    /*!
     * Get the iso value at the given world coordinates
     */
    inline Real iso_value(const Coord & x) const {
        Eigen::Map<Eigen::Matrix<Real, Dimension, 1>> mapped_x(&x[0]);
        iso_value(mapped_x);
    };

    /*!
     * Get the iso value at the given world coordinates
     */
    virtual Real iso_value(const Eigen::Matrix<Real, Dimension, 1> & x) const = 0;
};

}
