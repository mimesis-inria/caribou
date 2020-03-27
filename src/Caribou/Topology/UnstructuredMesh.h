#ifndef CARIBOU_TOPOLOGY_UNSTRUCTUREDMESH_H
#define CARIBOU_TOPOLOGY_UNSTRUCTUREDMESH_H

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Caribou/Topology/Mesh.h>

#include <vector>
#include <string>
#include <algorithm>

namespace caribou::topology {


    template <UNSIGNED_INTEGER_TYPE _Dimension>
    class UnstructuredMesh : public Mesh<UnstructuredMesh<_Dimension>> {
    public:
        friend class Mesh<UnstructuredMesh<_Dimension>>;

        using Parent = Mesh<UnstructuredMesh<_Dimension>>;

        template<typename Element>
        using Domain = Domain<Element>;

        static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
        using WorldCoordinates = std::array<FLOATING_POINT_TYPE, Dimension>;

        struct DomainRange;

        /*!
         * Default constructor.
         * Initializes an empty unstructured mesh.
         */
        UnstructuredMesh() = default;

        /**
         * Construct the unstructured mesh with a set of point positions from a position vector
         * @param positions
         */
        explicit UnstructuredMesh(const std::vector<WorldCoordinates> & positions) : p_positions(positions) {}

        /*! copy-and-swap move constructor */
        UnstructuredMesh(UnstructuredMesh<Dimension> && other) noexcept {
            swap(*this, other);
        }

        /*! copy-and-swap assigment */
        auto operator=(UnstructuredMesh<Dimension> other) noexcept -> UnstructuredMesh &
        {
            swap(*this, other);
            return *this;
        }

        /*!
         * Get the number of nodes of the mesh.
         */
        [[nodiscard]] inline auto number_of_nodes() const -> UNSIGNED_INTEGER_TYPE final {return p_positions.size();};

        /*!
         * Adds a node to the mesh.
         * @param position The position coordinates of the node.
         */
        inline void add_node(const WorldCoordinates position) {
            p_positions.push_back(position);
        }

        /*!
         * Adds nodes to the mesh.
         * @param position The position coordinates of the node.
         */
        inline void add_nodes(const std::vector<WorldCoordinates> & positions) {
            p_positions.insert(p_positions.end(), positions.begin(), positions.end());
        }

        friend void swap(UnstructuredMesh<Dimension> & first, UnstructuredMesh<Dimension>& second) noexcept
        {
            // enable ADL
            using std::swap;
            swap(first.p_positions, second.p_positions);
        }

    private:
        /*!
         * \copydoc caribou::topology::Mesh::position
         */
        [[nodiscard]]
        inline auto _position(UNSIGNED_INTEGER_TYPE index) const -> WorldCoordinates {
            return p_positions[index];
        }

        /*!
         * \copydoc caribou::topology::Mesh::_positions
         */
        template <typename IntegerType, std::size_t N>
        inline auto _positions(const IntegerType(&indices)[N]) const {
            std::array<WorldCoordinates, N> positions {};
            for (std::size_t i = 0; i < N; ++i) {
                positions[i] = this->p_positions[indices[i]];
            }
            return positions;
        }

        /*!
         * \copydoc caribou::topology::Mesh::_positions
         */
        template <typename IntegerType, std::size_t N>
        inline auto _positions(const std::array<IntegerType, N> & indices) const {
            std::array<WorldCoordinates, N> positions {};
            for (std::size_t i = 0; i < N; ++i) {
                positions[i] = this->p_positions[indices[i]];
            }
            return positions;
        }

        /*!
         * \copydoc caribou::topology::Mesh::_positions
         */
        template <typename IntegerType>
        inline auto _positions(const std::vector<IntegerType> & indices) const {
            std::vector<WorldCoordinates> positions {};
            positions.reserve(indices.size());
            for (const auto & index : indices) {
                positions.emplace_back(this->p_positions[index]);
            }
            return positions;
        }

    private:
        std::vector<WorldCoordinates> p_positions;
    };

    extern template class UnstructuredMesh<1>;
    extern template class UnstructuredMesh<2>;
    extern template class UnstructuredMesh<3>;
}

#endif //CARIBOU_TOPOLOGY_UNSTRUCTUREDMESH_H
