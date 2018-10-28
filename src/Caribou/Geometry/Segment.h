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
 * @tparam TPoint The type of Point we want to use for our segment
 */
template<class TPoint>
class Segment : public Entity
{
public:
    typedef TPoint PointType;
    static constexpr int Dimension = PointType::Dimension;

    Segment() = default;
    Segment(const TPoint & p1, const TPoint & p2) : Entity() , nodes({p1, p2}){};
    Segment(const Segment & s) : Entity() {
        std::copy(std::begin(s.nodes), std::end(s.nodes), std::begin(nodes));
    }

    /** Const accessor to the first node of the segment */
    inline const TPoint &
    first_node () const
    { return nodes[0]; }

    /** Accessor to the first node of the segment */
    inline TPoint &
    first_node ()
    { return nodes[0]; }

    /** Const accessor to the second node of the segment */
    inline const TPoint &
    second_node () const
    { return nodes[1]; }

    /** Accessor to the second node of the segment */
    inline TPoint &
    second_node ()
    { return nodes[1]; }

    /**
     * Assignment operator from an other segment (this = other;)
     * @tparam TOtherPoint The type of Point the other segment
     * @param other The other segment from which data we want to copy inside this one
     * @return This segment
     */
    template<class TOtherPoint>
    inline Segment<TPoint>
    &operator=(const Segment<TOtherPoint> & other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        std::copy(std::begin(other.nodes), std::end(other.nodes), std::begin(nodes));

        return *this;
    }

    /**
     * Equality operator (this == other)
     * @tparam TOtherPoint The type of Point the other segment
     * @param other The other segment against who we want to compare ourself
     * @return True if other is equal to this segment, false otherwise
     */
    template<class TOtherPoint>
    inline bool
    operator==(const Segment<TOtherPoint> & other) const
    {
        return (
                std::equal(std::begin(nodes), std::end(nodes), std::begin(other.nodes))
        );
    }

    /**
     * see operator==
     */
    template<class TOtherPoint>
    inline bool
    operator!=(const Segment<TOtherPoint>& other) const
    {
        return not (*this == other);
    }


    /**
     * Subscript operator (this[index])
     * @param index The index (0 or 1) of the point we want to get
     * @return A reference to the node at position index
     */
    inline TPoint &
    operator[] (std::size_t index)
    {
        return nodes[index];
    }

    /**
     * Subscript operator (this[index])
     * @param index The index (0 or 1) of the point we want to get
     * @return A const reference to the node at position index
     */
    inline const TPoint &
    operator[] (std::size_t x) const
    {
        return nodes[x];
    }

    std::array<TPoint, 2> nodes;
};

/**
 * Tool function to create a new segment between two points.
 *
 * Example:
 *
 * \code{.cpp}
 * // Create a segment between 2 points
 * auto p1 = make_point(0,0,0)
 * auto p2 = make_point(1,1,1)
 *
 * auto s1 = make_segment(p1, p2)
 * \endcode
 *
 * @tparam TPoint The type of Point we want to add to our segment. It is usually automatically deducted from the parameters.
 * @param p1 First node of the segment
 * @param p2 Second node of the segment
 *
 * @return A newly created segment between p1 and p2
 */
template<typename TPoint>
Segment<TPoint> make_segment(const TPoint & p1, const TPoint & p2) {
    return Segment<TPoint>(p1, p2);
}

template<size_t N, typename ValueType>
Segment<Point<N, caribou::algebra::Vector<N, ValueType>>>
make_segment(ValueType const (&a1)[N], ValueType const (&a2)[N])
{
    Point <N, caribou::algebra::Vector<N, ValueType>> p1(a1);
    Point <N, caribou::algebra::Vector<N, ValueType>> p2(a2);

    return make_segment(p1, p2);
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_SEGMENT_H
