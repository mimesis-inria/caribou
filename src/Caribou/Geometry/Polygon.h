#ifndef CARIBOU_GEOMETRY_POLYGON_H
#define CARIBOU_GEOMETRY_POLYGON_H
#include <cassert>

#include <Caribou/Geometry/Entity.h>
#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>

namespace caribou
{
namespace geometry
{

/**
 * A polygon (of n nodes) in space (independent of the space dimension).
 */
template<size_t NNodes, typename TVector>
class Polygon : public Entity
{
public:
    static constexpr size_t Dimension = TVector::Dimension;
    static constexpr size_t NumberOfSegments = NNodes;
    static constexpr size_t NumberOfNodes = NNodes;

    using VectorType = TVector;
    using PointType = Point<Dimension, VectorType>;

    static_assert(Dimension >= 2, "A polygon must be defined on space of dimension 2 or greater.");
    static_assert(NNodes >= 3, "A polygon must have at least three nodes.");

    Polygon() = default;

    Polygon(PointType const (&n)[NumberOfNodes]) : Entity(), nodes(n) {}
    Polygon(const std::array<PointType, NumberOfNodes> & n) : Entity(), nodes(n) {}
    Polygon(const std::initializer_list<PointType > & il) : Entity() {
        std::copy(std::begin(il), std::end(il), std::begin(nodes));
    }

    template<typename TOtherVector>
    inline bool operator==(const Polygon<NNodes, TOtherVector> & other) const {
        using OtherPointType = typename Polygon<NNodes, TOtherVector>::PointType;
        return (
                std::equal(std::begin(nodes), std::end(nodes), std::begin(other.nodes), [](const PointType& a, const OtherPointType &b) -> bool {
                    return a == b;
                })
        );
    }

    template<typename TOtherPoint>
    inline bool operator!=(const Polygon<NNodes, TOtherPoint> & other) const {
        return not (*this == other);
    }

    inline PointType & operator[] (std::size_t x) {
        return nodes[x];
    }

    inline const PointType & operator[] (std::size_t x) const {
        return nodes[x];
    }

    std::array<PointType, NumberOfNodes> nodes;
};

//
// Tool functions to create a polygon with a list of segments
//
template<size_t NSegments, typename VectorType>
Polygon<NSegments, VectorType>
make_polygon_from_segments(Segment<VectorType> const (&segments)[NSegments])
{
    assert(
            segments[0][0] == segments[NSegments-1][1] &&
            "The segments of a polygon must form a closed loop (the first node of the first segment must be "
            "equal to the last node of the last segment)."
    );

    constexpr size_t Dimension = VectorType::Dimension;
    using PointType = Point<Dimension, VectorType>;

    std::array<PointType, NSegments> points;
    for (size_t i = 0; i < NSegments; ++i)
        points[i] = segments[i][0];

    return Polygon<NSegments, VectorType>(points);
}

template<class VectorType, typename... Args>
auto make_polygon (Segment<VectorType> & arg, Args&&... args) {
    return make_polygon_from_segments<sizeof...(args)+1>({ arg, std::forward<Args>(args)... });
}

//
// Tool functions to create a polygon with a list of position vectors
//
template<size_t Dimension, typename TReal, typename... Args>
auto make_polygon (const Point<Dimension, algebra::Vector<Dimension, TReal>> & arg, Args&&... args)
{
    return Polygon<sizeof...(args)+1, algebra::Vector<Dimension, TReal>>({ arg, std::forward<Args>(args)... });
}

template<size_t Dimension, size_t NumberOfNodes, typename ValueType>
auto make_polygon(ValueType const (&arg)[NumberOfNodes][Dimension])
{
    Polygon<NumberOfNodes, caribou::algebra::Vector<Dimension,ValueType>> p;
    for (size_t i = 0; i < NumberOfNodes; ++i) {
        p[i] = arg[i];
    }
    return p;
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_POLYGON_H
