#ifndef CARIBOU_GEOMETRY_SEGMENT_H
#define CARIBOU_GEOMETRY_SEGMENT_H

#include <Caribou/Geometry/Point.h>

#ifdef CARIBOU_USE_DOUBLE
#define FLOATING_POINT_TYPE double
#else
#define FLOATING_POINT_TYPE float
#endif

namespace caribou
{
namespace geometry
{

/**
 * A segment (between 2 nodes) in space (independent of the space dimension).
 * @tparam TPoint The type of Point we want to use for our segment
 */
template <size_t Dim>
class Segment
{
public:
    static constexpr size_t Dimension = Dim;

    using PointType = Point<Dimension>;
    using VectorType = typename PointType::VectorType ;

    Segment() = default;
    Segment(const PointType & p1, const PointType & p2) : nodes({p1, p2}){};
    Segment(const Segment & s) {
        std::copy(std::begin(s.nodes), std::end(s.nodes), std::begin(nodes));
    }

    /** Const accessor to the first node of the segment */
    inline const PointType &
    first_node () const
    { return nodes[0]; }

    /** Accessor to the first node of the segment */
    inline PointType &
    first_node ()
    { return nodes[0]; }

    /** Const accessor to the second node of the segment */
    inline const PointType &
    second_node () const
    { return nodes[1]; }

    /** Accessor to the second node of the segment */
    inline PointType &
    second_node ()
    { return nodes[1]; }

    /** Get the direction vector (from p1 to p2) **/
    inline VectorType
    direction() const
    {
        return second_node() - first_node();
    }

    /** Get the unit vector (aka normalized direction vector)  **/
    inline VectorType
    unit() const
    {
        VectorType dir = direction();
        return dir / dir.length();
    }

    /** Get the reversed segment (from p2 to p1) **/
    inline Segment<Dimension>
    reversed() const
    {
        return Segment<Dimension>(second_node(), first_node());
    }

    /** Get the length of the segment **/
    inline typename VectorType::ValueType
    length() const
    {
        return direction().length();
    }

    /** Test if the segment contains the position x **/
    inline bool
    contains(const VectorType & x, FLOATING_POINT_TYPE tolerance = 0.00001 ) const
    {
        const auto & p0 = first_node();
        const auto & p1 = second_node();

        if (((x-p0)^(x-p1)).length_squared() > tolerance*tolerance) {
            // The point isn't on the line made by the segment
            return false;
        }

        const auto c = (p1 - p0) * (x - p0);
        if (c < 0 || c > (p1 - p0) * (p1 - p0)) {
            // The point is on the line, but isn't between the two nodes of the segment
            return false;
        }

        return true;
    }

    /**
     * Assignment operator from an other segment (this = other;)
     * @param other The other segment from which data we want to copy inside this one
     * @return This segment
     */
    inline Segment &
    operator=(const Segment<Dimension> & other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        std::copy(std::begin(other.nodes), std::end(other.nodes), std::begin(nodes));

        return *this;
    }

    /**
     * Equality operator (this == other)
     * @param other The other segment against who we want to compare ourself
     * @return True if other is equal to this segment, false otherwise
     */
    inline bool
    operator==(const Segment & other) const
    {
        return (
                std::equal(std::begin(nodes), std::end(nodes), std::begin(other.nodes))
        );
    }

    /**
     * see operator==
     */
    inline bool
    operator!=(const Segment & other) const
    {
        return not (*this == other);
    }


    /**
     * Subscript operator (this[index])
     * @param index The index (0 or 1) of the point we want to get
     * @return A reference to the node at position index
     */
    inline PointType &
    operator[] (std::size_t index)
    {
        return nodes[index];
    }

    /**
     * Subscript operator (this[index])
     * @param index The index (0 or 1) of the point we want to get
     * @return A const reference to the node at position index
     */
    inline const PointType &
    operator[] (std::size_t x) const
    {
        return nodes[x];
    }

    std::array<PointType, 2> nodes;
};

/**
 * Create a new segment between two points.
 *
 * Example:
 *
 * \code{.cpp}
 * auto p1 = make_point(0,0,0)
 * auto p2 = make_point(1,1,1)
 *
 * auto s1 = make_segment(p1, p2)
 * \endcode
 *
 * @tparam Dim The dimension of the point
 * @param p1 First node of the segment
 * @param p2 Second node of the segment
 *
 * @return A newly created segment between p1 and p2
 */
template <size_t Dim>
Segment<Dim> make_segment(const Point<Dim> & p1, const Point<Dim> & p2) {
    return Segment<Dim>(p1, p2);
}

template<size_t N, typename ValueType>
Segment<N>
make_segment(ValueType const (&a1)[N], ValueType const (&a2)[N])
{
    Point <N> p1(a1);
    Point <N> p2(a2);

    return make_segment(p1, p2);
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_SEGMENT_H
