#ifndef CARIBOU_TOPOLOGY_DOMAIN_H
#define CARIBOU_TOPOLOGY_DOMAIN_H

#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Caribou/Constants.h>
#include <Caribou/Topology/BaseDomain.h>
#include <Caribou/Geometry/Element.h>

#include <vector>
#include <array>

namespace caribou::topology {

    // Forward declaration
    template <UNSIGNED_INTEGER_TYPE Dimension>
    class Mesh;

    /*!
     * The domain storage is used as a way to personalize the storage type of the domain based on the Element type.
     * It is empty by default, but can be changed by template specialization of a given Element type.
     * @tparam Element
     */
    template <typename Element>
    class DomainStorage {};

    /*!
     * A domain implements the BaseDomain interface with a given element type.
     *
     * @tparam Element See caribou::geometry::Element
     * @tparam NodeIndex The type of integer used for a node index
     */
    template <typename Element, typename NodeIndex = UNSIGNED_INTEGER_TYPE>
    class Domain final : public BaseDomain, private DomainStorage<Element> {
        friend Mesh<geometry::traits<Element>::Dimension>;
    public:
        static constexpr INTEGER_TYPE Dimension = geometry::traits<Element>::Dimension;

        /*!
         * The type of container that stores the node indices of all the elements of the domain.
         *
         * The indices should be stored in a dense matrix where each row represent an element, and the columns
         * are the node indices of the element row.
         */
        using ElementsIndices = Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime, Eigen::RowMajor>;

        /*!
         * The type of container that stores the node indices of a single element.
         */
        using ElementIndices = Eigen::Matrix<NodeIndex, geometry::traits<Element>::NumberOfNodesAtCompileTime, 1>;

        /*! Empty constructor is prohibited */
        Domain() = delete;

        /*! Copy constructor */
        Domain(const Domain<Element> & other) noexcept
        : DomainStorage<Element>(other), p_mesh(other.p_mesh), p_buffer(other.p_buffer), p_elements(other.p_elements) {}

        /*! Move constructor */
        Domain(Domain && other) noexcept {
            swap(*this, other);
        }

        /*! copy-and-swap assigment (valid for both copy and move assigment) */
        auto operator=(Domain other) noexcept -> Domain & {
            swap(*this, other);
            return *this;
        }

        /*! Destructor */
        ~Domain() final = default;

        /*!
         * \copydoc caribou::topology::BaseDomain::canonical_dimension
         */
        [[nodiscard]]
        auto canonical_dimension() const -> UNSIGNED_INTEGER_TYPE final {
            return geometry::traits<Element>::CanonicalDimension;
        }

        /*!
         * \copydoc caribou::topology::BaseDomain::number_of_nodes_per_elements
         */
        [[nodiscard]]
        auto number_of_nodes_per_elements() const -> UNSIGNED_INTEGER_TYPE final;

        /*!
         * \copydoc caribou::topology::BaseDomain::number_of_elements
         */
        [[nodiscard]]
        auto number_of_elements() const -> UNSIGNED_INTEGER_TYPE final;

        /*!
         * Construct an element of the domain
         * @param element_id The id of the element in this domain.
         * @param positions An eigen matrix containing all the positions of the mesh. The matrix should have N rows
         *                  and D columns for a mesh of dimension D containing N nodes.
         * @return A new Element instance from the domain
         */
         template<typename EigenMatrix>
        inline auto element(const UNSIGNED_INTEGER_TYPE & element_id, const Eigen::DenseBase<EigenMatrix> & positions) const -> Element;

        /*!
        * Construct an element of the domain using the positions vector of the associated mesh.
        * @param element_id The id of the element in this domain.
        * @return A new Element instance from the domain
        */
        inline auto element(const UNSIGNED_INTEGER_TYPE & element_id) const -> Element;

        /*!
         * Get the indices of an element in the domain.
         * @param index The index of the element.
         * @return An Eigen vector containing the element indices.
         *
         */
        inline auto element_indices(const UNSIGNED_INTEGER_TYPE & index) const;

        /*!
         * Get the mesh associated to this domain.
         */
        inline auto mesh() const -> const Mesh<Dimension> & {
            return *p_mesh;
        }

        friend void swap(Domain<Element> & first, Domain<Element>& second) noexcept
        {
            // enable ADL
            using std::swap;
            swap(first.p_buffer, second.p_buffer);
            swap(first.p_elements, second.p_elements);
        }

    protected:

        /*!
         * Construct the domain from an array of indices.
         *
         * \note The indices are copied.
         */
        Domain(const Mesh<Dimension> * mesh, const ElementsIndices & elements)
            : p_mesh(mesh), p_buffer(elements), p_elements(p_buffer.data(), p_buffer.rows(), p_buffer.cols(), {p_buffer.cols(), 1}) {}

        /*!
         * Construct the domain from an array of indices.
         *
         * \note The indices are copied.
         */
        Domain(const Mesh<Dimension> * mesh, ElementsIndices & elements)
            : p_mesh(mesh), p_buffer(elements), p_elements(p_buffer.data(), p_buffer.rows(), p_buffer.cols(), {p_buffer.cols(), 1}) {}

        /*!
         * Construct the domain from a reference to an external array of indices.
         *
         * \note The indices are NOT copied. If the external array of indices is deleted, or changes, the behavior
         *       of the domain is undefined.
         */
        Domain(const Mesh<Dimension> * mesh, const ElementsIndices * elements)
            : p_mesh(mesh), p_buffer(), p_elements(elements->data(), elements->rows(), elements->cols(), {elements->cols(), 1}) {}

        /*!
         * Construct the domain from a reference to an external array of indices.
         *
         * \note The indices are NOT copied. If the external array of indices is deleted, or changes, the behavior
         *       of the domain is undefined.
         */
        explicit Domain(const Mesh<Dimension> * mesh, const NodeIndex * data, const Eigen::Index & number_of_elements, const Eigen::Index & number_of_nodes_per_elements)
            : p_mesh(mesh), p_buffer(), p_elements(data, number_of_elements, number_of_nodes_per_elements, {number_of_nodes_per_elements, 1}) {}

        /*!
         * Construct the domain from a reference to an array of indices.
         *
         * \note The indices are NOT copied. If the external array of indices is deleted, or changes, the behavior
         *       of the domain is undefined.
         */
        explicit Domain(const Mesh<Dimension> * mesh, const NodeIndex * data, const Eigen::Index & number_of_elements, const Eigen::Index & number_of_nodes_per_elements, Eigen::Index outer_stride, Eigen::Index inner_stride)
            : p_mesh(mesh), p_buffer(), p_elements(data, number_of_elements, number_of_nodes_per_elements, {outer_stride, inner_stride}) {}

        /// The mesh associated with this domain.
        const Mesh<Dimension> * p_mesh;

        /// Buffer containing the element indices. This buffer is used when the domain is constructed by copying
        /// an array of element indices. If the domain is constructed from a mapped buffer (ie, the indices are
        /// stored outside of the domain instance), then this buffer is empty.
        ElementsIndices p_buffer;

        /// Actual pointer to the element indices. When the domain is constructed by copying an array of element
        /// indices, this map points to p_buffer. Else, it points to an external buffer.
        Eigen::Map<const ElementsIndices, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> p_elements;
    };

    // -------------------------------------------------------------------------------------
    //                                    IMPLEMENTATION
    // -------------------------------------------------------------------------------------
    // Implementation of methods that can be specialized later for an explicit Element type
    // -------------------------------------------------------------------------------------

    template <typename Element, typename NodeIndex>
    auto Domain<Element, NodeIndex>::number_of_nodes_per_elements() const -> UNSIGNED_INTEGER_TYPE {
        return p_elements.cols();
    }

    template <typename Element, typename NodeIndex>
    auto Domain<Element, NodeIndex>::number_of_elements() const -> UNSIGNED_INTEGER_TYPE {
        return p_elements.rows();
    }

    template <typename Element, typename NodeIndex>
    inline auto Domain<Element, NodeIndex>::element_indices(const UNSIGNED_INTEGER_TYPE & index) const {
        caribou_assert(index < number_of_elements() and "Trying to get an element that does not exists.");

        return Eigen::Map<const Eigen::Matrix<NodeIndex, geometry::traits<Element>::NumberOfNodesAtCompileTime, 1>, Eigen::Unaligned, Eigen::Stride<1, Eigen::Dynamic>>(
            p_elements.row(index).data(), {1, p_elements.innerStride()}
        );
    }

    template <typename Element, typename NodeIndex>
    template<typename EigenMatrix>
    inline auto Domain<Element, NodeIndex>::element(const UNSIGNED_INTEGER_TYPE & element_id, const Eigen::DenseBase<EigenMatrix> & positions) const -> Element {
        caribou_assert(element_id < number_of_elements() &&
                       "Trying to get the element #"+std::to_string(element_id) + ", but the domain only has " +
                       std::to_string(number_of_elements()) + " elements."
        );

        using NodeMatrix = typename geometry::Element<Element>::template Matrix<geometry::traits<Element>::NumberOfNodesAtCompileTime, Dimension>;
        NodeMatrix node_positions;
        if constexpr (geometry::traits<Element>::NumberOfNodesAtCompileTime == caribou::Dynamic) {
            node_positions.resize(number_of_nodes_per_elements(), Dimension);
        }

        const auto node_indices = element_indices(element_id);
        for (std::size_t i = 0; i < node_indices.size(); ++i) {
            node_positions.row(i) = positions.row(node_indices[i]);
        }

        return Element(node_positions);
    }

template <typename Element, typename NodeIndex>
inline auto Domain<Element, NodeIndex>::element(const UNSIGNED_INTEGER_TYPE & element_id) const -> Element {
    caribou_assert(element_id < number_of_elements() &&
                   ("Trying to get the element #"+std::to_string(element_id) + ", but the domain only has " +
                   std::to_string(number_of_elements()) + " elements.").c_str()
    );

    return Element(mesh().positions(element_indices(element_id)));
}
}

#endif //CARIBOU_TOPOLOGY_DOMAIN_H
