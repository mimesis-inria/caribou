#ifndef CARIBOU_HELPER_TRIANGLE_H
#define CARIBOU_HELPER_TRIANGLE_H

namespace sofa {

namespace caribou {

namespace helper {

namespace triangle {

template<class Coord, class Real = double>
inline Coord center(const Coord & p1, const Coord & p2, const Coord & p3) noexcept {
    return (p1 + p2 + p3)/ (Real) 3.0;
};

template<class Coord, class Real = double>
inline Coord normal(const Coord & p1, const Coord & p2, const Coord & p3) noexcept {
    Coord d = (p2-p1).cross(p3-p1);
    Real  l = d.norm();
    if (l > 0)
        return d/l;
    else
        return Coord();
};

template<class Coord, class Real = double>
inline Real area(const Coord & p1, const Coord & p2, const Coord & p3) noexcept {
    return ((p2-p1).cross(p3-p1)).norm() * (Real) 0.5;
};

} // namespace triangle

} // namespace helper

} // namespace caribou

} // namespace sofa

#endif //CARIBOU_TRIANGLE_H
