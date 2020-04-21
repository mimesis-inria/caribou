#ifndef CARIBOU_TOPOLOGY_MESH_H
#define CARIBOU_TOPOLOGY_MESH_H

#include <Caribou/config.h>
#include <Caribou/Topology/BaseMesh.h>
#include <Caribou/Topology/Domain.h>

#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace caribou::topology {

    template <UNSIGNED_INTEGER_TYPE _Dimension>
    class Mesh : public BaseMesh {
    public:

        static constexpr INTEGER_TYPE Dimension = _Dimension;
        using WorldCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

        /*!
         * Default constructor.
         * Initializes an empty unstructured mesh.
         */
        Mesh() = default;

        /**
         * Construct the unstructured mesh with a set of point positions from a position vector.
         * @param positions
         * \note A copy of the full position vector is made.
         */
        explicit Mesh(const std::vector<WorldCoordinates> & positions) : p_positions(positions) {}

        /*! Copy constructor */
        Mesh(const Mesh<_Dimension> & other)
            : p_positions(other.p_positions) {}

        /*! Move constructor */
        Mesh(Mesh && other) noexcept {
            swap(*this, other);
        }

        /*! copy-and-swap assigment (valid for both copy and move assigment) */
        auto operator=(Mesh other) noexcept -> Mesh &
        {
            swap(*this, other);
            return *this;
        }

        /*!
         * \copydoc caribou::topology::BaseMesh::dimension
         */
        [[nodiscard]]
        inline auto dimension() const -> UNSIGNED_INTEGER_TYPE final {
            return Dimension;
        };

        /*!
         * \copydoc caribou::topology::BaseMesh::number_of_domains
         */
        [[nodiscard]]
        inline auto number_of_domains() const -> UNSIGNED_INTEGER_TYPE final {
            return p_domains.size();
        }

        /*!
         * Get the number of nodes of the mesh.
         */
        [[nodiscard]]
        inline auto number_of_nodes() const -> UNSIGNED_INTEGER_TYPE final {return p_positions.size();};

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

        /*!
         * Get the list of domains of the mesh.
         */
        [[nodiscard]]
        inline auto domains() const -> std::vector<std::pair<std::string, const BaseDomain *>> {
            std::vector<std::pair<std::string, const BaseDomain *>> domains;
            domains.reserve(p_domains.size());
            for (const auto & d : p_domains) {
                domains.emplace_back(d.first, d.second.get());
            }
            return domains;
        }

        /*!
         * Get the ith domain of the mesh.
         * /note Undefined behavior is found if the given index is outside the size of domains in the mesh
         * */
        [[nodiscard]]
        inline auto domain(const UNSIGNED_INTEGER_TYPE & i) const -> const BaseDomain * {
            caribou_assert((i < p_domains.size()) &&
                           ("Trying to get the " + std::to_string(i+1) + "th domain while the mesh has only " +
                std::to_string(p_domains.size()) + " domains.").c_str()
            );

            return p_domains[i].second.get();
        }

        /*! Get a domain of the mesh from its name. Returns null if no domain has the given name. */
        [[nodiscard]]
        inline auto domain(const std::string & name) const -> const BaseDomain * {
            for (const auto & d : p_domains) {
                if (d.first == name) {
                    return d.second.get();
                }
            }

            return nullptr;
        }

        /*!
         * Adds a domain with the Element type supplied as a template parameter.
         * @tparam Element The type of Element
         * @param name
         * @return A pointer to the newly created domain, or null if the creation failed
         *
         * \note The domain is managed by this mesh instance, and will be freed upon the deletion of the mesh.
         */
        template<typename Element, typename... Args>
        inline auto add_domain(const std::string & name, Args ...args) -> Domain<Element> * {
            static_assert(geometry::traits<Element>::Dimension == Dimension, "The dimension of the mesh doesn't match the dimension of the element type.");
            for (const auto & d : p_domains) {
                if (d.first == name) {
                    throw std::domain_error(std::string("A domain with the name ") + name + " already exists.");
                }
            }

            auto domain_ptr = new Domain<Element>(this, std::forward<Args>(args)...);
            auto res = p_domains.emplace(p_domains.end(), name, std::unique_ptr<BaseDomain>(static_cast<BaseDomain *>(domain_ptr)));
            if (res->second) {
                return domain_ptr;
            } else {
                delete domain_ptr;
                return static_cast<Domain<Element>*> (nullptr);
            }
        }

        /*!
         * Adds a domain with the Element type supplied as a template parameter.
         * @tparam Element The type of Element
         * @param name
         * @return A pointer to the newly created domain, or null if the creation failed
         *
         * \note The domain is managed by this mesh instance, and will be freed upon the deletion of the mesh.
         */
        template<typename Element, typename... Args>
        inline auto add_domain(const char * name, Args ...args) -> Domain<Element> * {
            return add_domain<Element, Args...>(std::string(name), std::forward<Args>(args)...);
        }

        /*!
         * Adds a domain with the Element type supplied as a template parameter.
         * @tparam Element The type of Element
         * @return A pointer to the newly created domain, or null if the creation failed
         *
         * \note The domain is managed by this mesh instance, and will be freed upon the deletion of the mesh.
         */
        template<typename Element, typename... Args>
        inline auto add_domain(Args ...args) -> Domain<Element> * {
            static_assert(geometry::traits<Element>::Dimension == Dimension, "The dimension of the mesh doesn't match the dimension of the element type.");

            auto domain_ptr = new Domain<Element>(this, std::forward<Args>(args)...);

            std::ostringstream ss;
            ss << "domain_" << domain_ptr;

            auto res = p_domains.emplace(p_domains.end(), ss.str(), std::unique_ptr<BaseDomain>(static_cast<BaseDomain *>(domain_ptr)));
            if (res->second) {
                return domain_ptr;
            } else {
                delete domain_ptr;
                return static_cast<Domain<Element>*> (nullptr);
            }
        }

        /*!
         * Remove and destruct a domain of the mesh.
         * @param name Name of the domain
         * @return bool True on success, false otherwise
         */
        inline auto remove_domain(const std::string & name) -> bool {
            auto it = p_domains.begin();
            while (it != p_domains.end() and it->first != name);
            if (it != p_domains.end()) {
                p_domains.erase(it);
                return true;
            }
            return false;
        }

        /*!
         * Remove and destruct a domain of the mesh.
         * @param domain Pointer to the domain to remove.
         * @return bool True on success, false otherwise
         */
        inline auto remove(const BaseDomain * domain) -> bool {
            auto it = p_domains.begin();
            while (it != p_domains.end() and it->second.get() != domain);
            if (it != p_domains.end()) {
                p_domains.erase(it);
                return true;
            }
            return false;
        }

        /*!
        * Get the position coordinates of a node from its index.
        */
        [[nodiscard]]
        inline auto position(UNSIGNED_INTEGER_TYPE index) const -> const WorldCoordinates & {
            caribou_assert(index < p_positions.size() && "Trying to access the position vector at a node index greater "
                                                         "than the number of nodes in the mesh.");
            return p_positions[index];
        }

        /*!
         * Get an array of positions from the given indices.
         * @tparam IntegerType Integer type of the indices passed
         * @param indices [IN] The indices from which the positions array will be filled
         * @return The position coordinates of every nodes queried as an Eigen dense matrix.
         *
         * Example:
         * \code{.cpp}
         * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, 2, Dimension, Eigen::RowMajor>;
         * Nodes positions = mesh.positions({55, 62});
         * auto p55 = positions.row(0); // Position vector of the node id 55
         * \endcode
         *
         */
        template <typename IntegerType, std::size_t N>
        inline auto positions(const IntegerType(&indices)[N]) const {
            Eigen::Matrix<FLOATING_POINT_TYPE, N, Dimension, Eigen::RowMajor> positions;
            for (std::size_t i = 0; i < N; ++i) {
                positions.row(i) = p_positions[indices[i]];
            }
            return positions;
        }

        /*!
         * Get an array of positions from the given indices.
         * @tparam IntegerType Integer type of the indices passed
         * @param indices [IN] The indices from which the positions array will be filled
         * @return The position coordinates of every nodes queried as an Eigen dense matrix.
         *
         * Example:
         * \code{.cpp}
         * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, 2, Dimension, Eigen::RowMajor>;
         * std::array<int, 2> indices = {55, 62};
         * Nodes positions = mesh.positions(indices);
         * auto p55 = positions.row(0); // Position vector of the node id 55
         * \endcode
         *
         */
        template <typename IntegerType, std::size_t N>
        inline auto positions(const std::array<IntegerType, N> & indices) const {
            Eigen::Matrix<FLOATING_POINT_TYPE, N, Dimension, Eigen::RowMajor> positions;
            for (std::size_t i = 0; i < N; ++i) {
                positions.row(i) = p_positions[indices[i]];
            }
            return positions;
        }

        /*!
         * Get an array of positions from the given indices.
         * @tparam IntegerType Integer type of the indices passed
         * @param indices [IN] he indices from which the positions array will be filled
         * @return The position coordinates of every nodes queried as an Eigen dense matrix.
         *
         * Example:
         * \code{.cpp}
         * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Dimension, Eigen::RowMajor>;
         * std::vector<int> indices = {55, 62};
         * Nodes positions = mesh.positions(indices);
         * auto p55 = positions.row(0); // Position vector of the node id 55
         * \endcode
         *
         */
        template <typename IntegerType>
        inline auto positions(const std::vector<IntegerType> & indices) const {
            const auto number_of_indices = indices.size();
            Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Dimension, Eigen::RowMajor> positions;
            positions.resize(number_of_indices, Dimension);
            for (std::size_t i = 0; i < number_of_indices; ++i) {
                positions.row(i) = p_positions[indices[i]];
            }
            return positions;
        }

        /*!
         * Get an array of positions from the given indices.
         * @tparam EigenDerived The type of Eigen matrix/array containing the indices.
         * @param indices [IN] he indices from which the positions array will be filled
         * @return The position coordinates of every nodes queried as an Eigen dense matrix.
         *
         * Example:
         * \code{.cpp}
         * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Dimension, Eigen::RowMajor>;
         * Eigen::VectorXi indices(2);
         * indices << 55, 62;
         * Nodes positions = mesh.positions(indices);
         * auto p55 = positions.row(0); // Position vector of the node id 55
         * \endcode
         *
         */
        template <typename EigenDerived>
        inline auto positions(const Eigen::EigenBase<EigenDerived> & indices) const {
            static_assert(EigenDerived::RowsAtCompileTime == 1 or EigenDerived::ColsAtCompileTime == 1,
                "Only vector type matrix can be used as the indices container.");
            const auto number_of_indices = indices.size();
            Eigen::Matrix<FLOATING_POINT_TYPE,
                          EigenDerived::RowsAtCompileTime == 1 ? EigenDerived::ColsAtCompileTime : EigenDerived::RowsAtCompileTime,
                          Dimension,
                          Eigen::RowMajor> positions;
            positions.resize(number_of_indices, Dimension);
            for (Eigen::Index i = 0; i < number_of_indices; ++i) {
                positions.row(i) = p_positions[indices.derived()[i]];
            }
            return positions;
        }

        /*! Swap the data of two meshes */
        friend void swap(Mesh<Dimension> & first, Mesh<Dimension>& second) noexcept
        {
            // enable ADL
            using std::swap;
            swap(first.p_positions, second.p_positions);
        }

    private:
        std::vector<WorldCoordinates> p_positions;
        std::vector<std::pair<std::string, std::unique_ptr<BaseDomain>>> p_domains;
    };
}

#endif //CARIBOU_TOPOLOGY_MESH_H
