#pragma once

#include <sofa/defaulttype/VecTypes.h>

#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/IsoSurface.h>

namespace SofaCaribou::topology {

class CARIBOU_API CylinderIsoSurface : public IsoSurface<sofa::defaulttype::Vec3Types> {
    template < class T = void* >
    using Data = sofa::core::objectmodel::Data<T>;
    using Base = IsoSurface<sofa::defaulttype::Vec3Types>;
    using Base::Coord;
    using Base::Real;
public:

    CylinderIsoSurface()
    : p_radius(initData(&p_radius, 1., "radius", "Radius of the cylinder."))
    , p_length(initData(&p_length, 5., "length", "Length of the cylinder."))
    , p_center(initData(&p_center, Coord(0, 0, 0), "center", "Coordinates at the center of the cylinder."))
    {}

    inline Real iso_value(const Eigen::Matrix<Real, 3, 1> & x) const final {
        const auto & r = p_radius.getValue();
        Eigen::Map<const Eigen::Matrix<Real, 3, 1>> c (p_center.getValue().data());

        const auto d = x.block<2,1>(0,0) - c.block<2,1>(0,0);
        return d.squaredNorm() - r*r;
    }

private:
    Data<Real>  p_radius;
    Data<Real>  p_length;
    Data<Coord> p_center;
};

}
