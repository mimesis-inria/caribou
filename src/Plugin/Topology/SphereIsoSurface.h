#pragma once

#include <sofa/defaulttype/VecTypes.h>

#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/IsoSurface.h>

namespace SofaCaribou::topology {

class CARIBOU_API SphereIsoSurface : public IsoSurface<sofa::defaulttype::Vec3Types> {
    template < class T = void* >
    using Data = sofa::core::objectmodel::Data<T>;
    using Base = IsoSurface<sofa::defaulttype::Vec3Types>;
    using Base::Coord;
    using Base::Real;
public:

    SphereIsoSurface()
    : p_radius(initData(&p_radius, 1., "radius", "Radius of the sphere."))
    , p_center(initData(&p_center, Coord(0, 0, 0), "center", "Coordinates at the center of the sphere."))
    {}

    inline Real iso_value(const Eigen::Matrix<Real, 3, 1> & x) const final {
        const auto & r = p_radius.getValue();
        Eigen::Map<const Eigen::Matrix<Real, 3, 1>> c (p_center.getValue().data());

        const auto d = (x-c).dot(x-c);
        return d - r*r;
    }

private:
    Data<Real>  p_radius;
    Data<Coord> p_center;
};

}
