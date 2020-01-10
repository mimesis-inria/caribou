#ifndef CARIBOU_GEOMETRY_TRAITS_H
#define CARIBOU_GEOMETRY_TRAITS_H

#include <Caribou/Traits.h>


namespace caribou {
template <template<size_t, typename> typename ElementType, size_t Dim, typename CanonicalElementType>
struct traits<ElementType<Dim, CanonicalElementType>>
{
    using LocalCoordinates = typename CanonicalElementType::LocalCoordinates;
    using WorldCoordinates = typename ElementType<Dim, CanonicalElementType>::WorldCoordinates;
    enum {
        Dimension = Dim,
        NumberOfNodes = ElementType<Dim, CanonicalElementType>::NumberOfNodes,
        NumberOfGaussNodes = ElementType<Dim, CanonicalElementType>::number_of_gauss_nodes,
    };
};

template <template<typename> typename ElementType, typename CanonicalElementType>
struct traits<ElementType<CanonicalElementType>> {
    using LocalCoordinates = typename CanonicalElementType::LocalCoordinates;
    using WorldCoordinates = typename ElementType<CanonicalElementType>::WorldCoordinates;
    enum {
        Dimension = 3,
        NumberOfNodes = ElementType<CanonicalElementType>::NumberOfNodes,
        NumberOfGaussNodes = ElementType<CanonicalElementType>::number_of_gauss_nodes,
        NumberOfFaces = ElementType<CanonicalElementType>::number_of_faces,
    };
};
} // namespace caribou::geometry

#endif //CARIBOU_GEOMETRY_TRAITS_H
