#pragma once

#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Topology/BaseDomain.h>
#include <Caribou/Topology/BarycentricContainer.h>
#include <Caribou/Geometry/Element.h>

#include <Eigen/Core>
#include <vector>
#include <array>

namespace caribou::topology {

    /*!
     * The domain storage is used as a way to personalize the storage type of the domain based on the Element type.
     * It is empty by default, but can be changed by template specialization of a given Element type.
     * @tparam Element
     */
    template <typename Element>
    class DomainStorage {};

    /*!
     * A Domain is a subspace of a Mesh containing a set of points and the topological relation between them. It does not
     * contain any world positions of the points, but only their connectivity.
     *
     * The Domain class supports either internal storing of the node connectivity, or external storing (
     * for example when the vector of node indices for every elements are stored externally).
     *
     * In a Domain, all the elements are of the same type. For example, a Domain can not contain both hexahedrons and
     * tetrahedrons.
     *
     * A Domain can only reside inside one and only one Mesh. In fact, only a Mesh can create a Domain instance. The Mesh
     * will typically contain one or more domains.
     *
     * Example of a domain that stores internally its connectivity:
     * \code{.cpp}
     * // We supposed the Mesh containing the position of the nodes have been created before.
     * Mesh<_3D> * mesh = get_mesh();
     *
     * // Set the node connectivity of 4 triangles (each having 3 nodes).
     * Eigen::Matrix<unsigned int, 4, 3> indices;
     * indices << 0, 1, 3, // Triangle 1
     *            1, 4, 5, // Triangle 2
     *            8, 3, 1, // Triangle 3
     *            9, 5, 1; // Triangle 4
     *
     * // Here the indices array will be copied into the domain. It can therefore
     * // be safely deleted once the domain has been added to the mesh
     * mesh->add_domain<Triangle<_3D, Linear>>(indices);
     * \endcode
     *
     * Example of a domain that uses the connectivity stored externally:
     * \code{.cpp}
     * // We supposed the mesh containing the position of the nodes have been created before.
     * Mesh<_3D> * mesh = get_mesh();
     *
     * // Set the node connectivity of 4 triangles (each having 3 nodes).
     * unsigned int indices[8] = {0, 1, 3,  // Triangle 1
     *                            1, 4, 5,  // Triangle 2
     *                            8, 3, 1,  // Triangle 3
     *                            9, 5, 1}; // Triangle 4
     *
     * // Here the indices array will NOT be copied into the domain.
     * // Hence, it must remain valid for the entire lifetime of the domain.
     * mesh->add_domain<Triangle<_3D, Linear>>(indices, 4, 3);
     * \endcode
     *
     * \note More examples can be found in the file src/Caribou/Topology/test/test_domain.cpp
     *
     * @tparam Element See caribou::geometry::Element
     * @tparam NodeIndex The type of integer used for a node index
     */
    template <typename Element, typename NodeIndex = UNSIGNED_INTEGER_TYPE>
    class Domain : public BaseDomain, private DomainStorage<Element> {
    public:
        static constexpr INTEGER_TYPE Dimension = geometry::traits<Element>::Dimension;
        using ElementType = Element;
        using NodeIndexType = NodeIndex;

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
        Domain(const Domain & other) noexcept
        : DomainStorage<Element>(other), p_buffer(other.p_buffer), p_elements(other.p_elements) {}

        /*! Move constructor */
        Domain(Domain && other) noexcept : p_elements {nullptr, 0, 0} {
            swap(*this, other);
        }

        /*! Destructor */
        ~Domain() override = default;

        /*!
         * \copydoc caribou::topology::BaseDomain::clone
         */
        BaseDomain * clone() const override {
            return new Domain(*this); // Will call the copy constructor
        }

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
         * \copydoc caribou::topology::BaseDomain::mesh
         */
        inline auto mesh() const -> const BaseMesh * override { return p_mesh; }

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
        virtual inline auto element(const UNSIGNED_INTEGER_TYPE & element_id) const -> Element;

        /*!
         * Get the indices of an element in the domain.
         * @param index The index of the element.
         * @return An Eigen vector containing the element indices.
         *
         */
        inline auto element_indices(const UNSIGNED_INTEGER_TYPE & index) const;

        /**
         * Embed a set of nodes (in world coordinates) inside this domain. This will return a BarycentricContainer that
         * can be used to interpolate field values on these embedded nodes.
         * @tparam Derived NXD Eigen matrix representing the D dimensional coordinates of the N embedded points.
         * @param points The positions (in world coordinates) of the nodes embedded in this domain.
         * @return A BarycentricContainer instance.
         */
        template <typename Derived>
        inline auto embed(const Eigen::MatrixBase<Derived> & points) const -> BarycentricContainer<Domain> {
            return {this, points};
        }

        /*!
         * Construct the domain from an array of indices.
         *
         * \note The indices are copied.
         */
        Domain(const BaseMesh * mesh, const ElementsIndices & elements)
            : p_mesh(mesh), p_buffer(elements), p_elements(p_buffer.data(), p_buffer.rows(), p_buffer.cols(), {p_buffer.cols(), 1}) {}

        /*!
         * Construct the domain from an array of indices.
         *
         * \note The indices are copied.
         */
        Domain(const BaseMesh * mesh, ElementsIndices & elements)
            : p_mesh(mesh), p_buffer(elements), p_elements(p_buffer.data(), p_buffer.rows(), p_buffer.cols(), {p_buffer.cols(), 1}) {}

        /*!
         * Construct the domain from a reference to an external array of indices.
         *
         * \note The indices are NOT copied. If the external array of indices is deleted, or changes, the behavior
         *       of the domain is undefined.
         */
        Domain(const BaseMesh * mesh, const ElementsIndices * elements)
            : p_mesh(mesh), p_buffer(), p_elements(elements->data(), elements->rows(), elements->cols(), {elements->cols(), 1}) {}

        /*!
         * Construct the domain from a reference to an external array of indices.
         *
         * \note The indices are NOT copied. If the external array of indices is deleted, or changes, the behavior
         *       of the domain is undefined.
         */
        Domain(const BaseMesh * mesh, const NodeIndex * data, const Eigen::Index & number_of_elements, const Eigen::Index & number_of_nodes_per_elements)
            : p_mesh(mesh), p_buffer(), p_elements(data, number_of_elements, number_of_nodes_per_elements, {number_of_nodes_per_elements, 1}) {}

        /*!
         * Construct the domain from a reference to an array of indices.
         *
         * \note The indices are NOT copied. If the external array of indices is deleted, or changes, the behavior
         *       of the domain is undefined.
         */
        Domain(const BaseMesh * mesh, const NodeIndex * data, const Eigen::Index & number_of_elements, const Eigen::Index & number_of_nodes_per_elements, Eigen::Index outer_stride, Eigen::Index inner_stride)
            : p_mesh(mesh), p_buffer(), p_elements(data, number_of_elements, number_of_nodes_per_elements, {outer_stride, inner_stride}) {}

        /*!
         * Attach this Domain to the given Mesh.
         */
        inline void attach_to(const BaseMesh * mesh) override {
            // Trying to attach the same mesh as the one which is already attached to this domain
            if (mesh == p_mesh) {
                return;
            }

            // There is already a Mesh instance attached to this Domain, and this domain can still be found in the
            // list of domains of the mesh, let's throw an exception
            if (p_mesh != nullptr)  {
                throw std::logic_error("Trying to add a domain that is already attached to another Mesh instance.");
            }

            p_mesh = mesh;
        }

    protected:

        friend void swap(Domain & first, Domain& second) noexcept
        {
            // enable ADL
            using std::swap;
            swap(first.p_buffer, second.p_buffer);

            // Map isn't swappable, we need to do it manually using placement new operators
            using MapType = decltype(p_elements);
            MapType temp (first.p_elements.data(), first.p_elements.rows(), first.p_elements.cols(), {first.p_elements.outerStride(), first.p_elements.innerStride()});
            new (&first.p_elements) MapType (second.p_elements.data(), second.p_elements.rows(), second.p_elements.cols(), {second.p_elements.outerStride(), second.p_elements.innerStride()});
            new (&second.p_elements) MapType (temp.data(), temp.rows(), temp.cols(), {temp.outerStride(), temp.innerStride()});
        }

        /// Pointer to the associated Mesh. This Domain instance should be part of one Mesh, where the latter could
        /// contains many Domain instances of different or same element types.
        const BaseMesh * p_mesh;

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
        caribou_assert(element_id < number_of_elements());

        using NodeMatrix = typename geometry::Element<Element>::template Matrix<geometry::traits<Element>::NumberOfNodesAtCompileTime, Dimension>;
        NodeMatrix node_positions;
        if constexpr (geometry::traits<Element>::NumberOfNodesAtCompileTime == caribou::Dynamic) {
            node_positions.resize(number_of_nodes_per_elements(), Dimension);
        }

        const auto node_indices = element_indices(element_id);
        for (std::size_t i = 0; i < static_cast<std::size_t>(node_indices.size()); ++i) {
            auto p1 = node_positions.row(i);
            auto p2 = this->mesh()->node(node_indices[i]);
            p1[0] = p2[0];
            if constexpr (Dimension > 1) {
                p1[1] = p2[1];
            }
            if constexpr(Dimension > 2) {
                p1[2] = p2[2];
            }
        }

        return Element(node_positions);
    }
}
