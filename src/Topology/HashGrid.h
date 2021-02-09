#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Element.h>
#include <unordered_map>
#include <vector>
#include <set>
#include <Eigen/Core>
#include <bitset>

namespace caribou::topology {

struct CwiseRound {
    auto operator()(const FLOATING_POINT_TYPE& x) const -> FLOATING_POINT_TYPE { return std::round(x); }
};
struct CwiseFloor {
    auto operator()(const FLOATING_POINT_TYPE& x) const -> FLOATING_POINT_TYPE { return std::floor(x); }
};

template <typename Element>
class HashGrid {
public:

    static constexpr UNSIGNED_INTEGER_TYPE Dimension = caribou::geometry::traits<Element>::Dimension;

    using Index = UNSIGNED_INTEGER_TYPE;
    using GridCoordinates = Eigen::Matrix<INTEGER_TYPE, Dimension, 1>;
    using WorldCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;
    using VecFloat = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    HashGrid (const FLOATING_POINT_TYPE & cell_size) : p_cell_size(cell_size) {}

    HashGrid (const FLOATING_POINT_TYPE & cell_size, const UNSIGNED_INTEGER_TYPE & number_of_elements) : p_cell_size(cell_size) {
        p_hash_table.reserve(2*number_of_elements);
    }

    /**
     * Add an element's data to the hash grid
     * @param e The element defined by its nodes.
     * @param id The element's identifier (usually its index within a topology).
     */
    inline
    void add(const Element & e, const Index & id) {
        const auto min = (e.nodes().colwise().minCoeff() * 1/p_cell_size).unaryExpr(CwiseFloor()).eval();
        const auto max = (e.nodes().colwise().maxCoeff() * 1/p_cell_size).unaryExpr(CwiseFloor()).eval();

        for (INTEGER_TYPE i = min[0]; i <= max[0]; ++i) {
            if constexpr (Dimension == 1) {
                p_hash_table[GridCoordinates {i}].emplace_back(id);
            } else {
                for (INTEGER_TYPE j = min[1]; j <= max[1]; ++j) {
                    if constexpr (Dimension == 2) {
                        p_hash_table[GridCoordinates {i, j}].emplace_back(id);
                    } else {
                        for (INTEGER_TYPE k = min[2]; k <= max[2]; ++k) {
                            p_hash_table[GridCoordinates {i, j, k}].emplace_back(id);
                        }
                    }
                }
            }
        }
    }

    /**
     * Get all the data of all elements that are very close to the point p. Note that the returned elements do not
     * ensure that the point p resides inside of them. One has to further check each ones of them with an
     * intersection test.
     */
    inline
    auto get(const WorldCoordinates & p) const -> std::set<Index> {
        const VecFloat absolute = p / p_cell_size;
        const VecFloat rounded = absolute.unaryExpr(CwiseRound());
        const VecFloat distance = absolute - rounded;

        auto close_to_axis = std::bitset<Dimension>();

        std::vector<GridCoordinates> cells;
        cells.reserve((unsigned) 1<<Dimension);

        // Let's find the axis for which our coordinate is very close to
        for (UNSIGNED_INTEGER_TYPE axis = 0; axis < Dimension; ++axis) {
            if (distance[axis]*distance[axis] < EPSILON*EPSILON )
                close_to_axis[axis] = true;
        }

        if (close_to_axis.none()) {
            // We are not near any axis, which means we are well inside a cell's boundaries
            cells.emplace_back(absolute.unaryExpr(CwiseFloor()). template cast<INTEGER_TYPE>());
        } else {
            const VecFloat d = VecFloat::Constant(1 / 2.);

            std::array<std::vector<INTEGER_TYPE>, Dimension> axis_indices;
            for (auto & indices  : axis_indices) {
                indices.reserve(Dimension);
            }

            for (UNSIGNED_INTEGER_TYPE axis = 0; axis < Dimension; ++axis) {
                if (close_to_axis[axis]) {
                    axis_indices[axis].emplace_back(floor(rounded[axis]-d[axis]));
                    axis_indices[axis].emplace_back(floor(rounded[axis]+d[axis]));
                } else {
                    axis_indices[axis].emplace_back(floor(absolute[axis]));
                }
            }

            for (const INTEGER_TYPE & i : axis_indices[0]) {
                if constexpr (Dimension == 1) {
                    cells.emplace_back(GridCoordinates{i});
                } else {
                    for (const INTEGER_TYPE & j : axis_indices[1]) {
                        if constexpr (Dimension == 2) {
                            cells.emplace_back(GridCoordinates{i, j});
                        } else {
                            for (const INTEGER_TYPE & k : axis_indices[2]) {
                                cells.emplace_back(GridCoordinates{i, j, k});
                            }
                        }
                    }
                }
            }
        }

        std::set<Index> elements;
        for (const auto & cell_coordinates : cells) {
            const auto & iter = p_hash_table.find(cell_coordinates);
            if (iter != p_hash_table.end()) {
                for (const auto &element_data : iter->second) {
                    elements.emplace(element_data);
                }
            }
        }

        return elements;
    }

private:

    struct HashFunction
    {
        auto operator()(const GridCoordinates & coordinates) const -> std::size_t
        {
            // We use the large prime numbers proposed in paper:
            // M.Teschner et al "Optimized Spatial Hashing for Collision Detection of Deformable Objects" (2003)
            INTEGER_TYPE h = 73856093 * coordinates[0] ;
            if constexpr (Dimension > 1)
                h ^= 19349663*coordinates[1];
            if constexpr (Dimension > 2)
                h ^= 83492791*coordinates[2];

            return h;
        }
    };

    struct HashEqual
    {
        auto operator()(const GridCoordinates &a, const GridCoordinates &b) const -> bool
        {
            return ((a - b).norm() == 0);
        }
    };


    FLOATING_POINT_TYPE p_cell_size;
    std::unordered_map<GridCoordinates, std::vector<Index>, HashFunction, HashEqual> p_hash_table;
};

} // namespace caribou::topology