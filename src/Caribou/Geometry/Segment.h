#ifndef CARIBOU_GEOMETRY_SEGMENT_H
#define CARIBOU_GEOMETRY_SEGMENT_H


#include <Caribou/Geometry/Entity.h>
#include <Caribou/Geometry/Point.h>

namespace caribou
{
namespace geometry
{

/**
 * A segment (between 2 nodes) in space (independent of the space dimension).
 * @tparam Dim Dimension of the current space (default to 3D).
 * @tparam Real Type of the floating values.
 */
template<int Dim=3, typename Data=BaseData, typename Real=float, typename TPoint=Point<Dim, BaseData, Real>>
class Segment : public Entity<Data>
{
public:
    Segment() = delete;
    Segment(const TPoint & p1, const TPoint & p2) : nodes({p1, p2}) {};
    Segment(std::array<TPoint, 2> l) : nodes(l) {}
    Segment(const Segment<Dim, Data, Real, TPoint> & s) {
        std::copy(std::begin(s.nodes), std::end(s.nodes), std::begin(nodes));
    }

    bool operator==(const Segment<Dim, BaseData, Real, TPoint> & s) const {
        return (
                this->data == s.data &&
                std::equal(std::begin(nodes), std::end(nodes), std::begin(s.nodes))
        );
    }

    bool operator!=(const Segment<Dim, Data, Real> & s) const {
        return not (*this == s);
    }

    TPoint & operator[] (std::size_t x) {
        return nodes[x];
    }

protected:
    std::array<TPoint, 2> nodes;
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_SEGMENT_H
