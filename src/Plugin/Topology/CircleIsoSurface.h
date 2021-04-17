#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/IsoSurface.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::topology {

class CircleIsoSurface : public IsoSurface<sofa::defaulttype::Vec2Types> {
    template < class T = void* >
    using Data = sofa::core::objectmodel::Data<T>;
    using Base = IsoSurface<sofa::defaulttype::Vec2Types>;
    using Base::Coord;
    using Base::Real;
public:

    CircleIsoSurface()
    : p_radius(initData(&p_radius, 1., "radius", "Radius of the circle."))
    , p_center(initData(&p_center, Coord(0, 0), "center", "Coordinates at the center of the circle."))
    {}

    inline Real iso_value(const Eigen::Matrix<Real, Dimension, 1> & x) const final {
        const auto & r = p_radius.getValue();
        Eigen::Map<const Eigen::Matrix<Real, Dimension, 1>> c (p_center.getValue().data());

        const auto d = (x-c).dot(x-c);
        return d - r*r;
    }

private:
    Data<Real>  p_radius;
    Data<Coord> p_center;
};

}
