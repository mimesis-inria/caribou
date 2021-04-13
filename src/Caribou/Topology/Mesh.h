#pragma once

#include <Caribou/config.h>
#include <Caribou/Topology/BaseMesh.h>
#include <Caribou/Topology/Domain.h>

#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <algorithm>

namespace caribou::topology {

template<typename Derived>
struct BaseEigenNodesHolder {
    using Scalar = typename Eigen::MatrixBase<Derived>::Scalar ;
    using MatrixType = Derived;

    /*! Copy constructor */
    BaseEigenNodesHolder(const BaseEigenNodesHolder & other)
            : p_nodes(other.p_nodes)
    {}

    /*! Move constructor */
    BaseEigenNodesHolder(BaseEigenNodesHolder && other) noexcept
    : BaseEigenNodesHolder() {
        this->p_nodes.swap(other.p_nodes);
    }

    /*! copy-and-swap assigment (valid for both copy and move assigment) */
    auto operator=(BaseEigenNodesHolder other) noexcept -> BaseEigenNodesHolder & {
        this->p_nodes.swap(other.p_nodes);
        return *this;
    }

    template<typename Index>
    auto node(Index && index) const -> auto {return this->p_nodes.row(index);}

    template<typename Index>
    auto node(Index && index) -> auto {return this->p_nodes.row(index);}

    template<typename Size1>
    auto resize(Size1 && n) -> auto {return this->p_nodes.resize(std::forward<Size1>(n), this->p_nodes.cols());}

    auto size() const -> auto {return this->p_nodes.rows();}

    BaseEigenNodesHolder()
            : p_nodes()
    {}

protected:
    MatrixType p_nodes;
};

template <typename T>
struct EigenNodesHolder {};

template<typename Scalar_t, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct EigenNodesHolder<Eigen::Matrix<Scalar_t, Rows, Cols, Options, MaxRows, MaxCols>>
: public BaseEigenNodesHolder<Eigen::Matrix<Scalar_t, Rows, Cols, Options, MaxRows, MaxCols>>
{
    using Base = BaseEigenNodesHolder<Eigen::Matrix<Scalar_t, Rows, Cols, Options, MaxRows, MaxCols>>;
    using Scalar = typename Base::Scalar;
    using MatrixType = typename Base::MatrixType;

    using Base::Base;
};

// todo(jnbrunet2000@gmail.com): Implement a EigenNodesHolder<Eigen::Map<...>> container type to enable
//                               references to external buffers.

    /**
     * The Mesh class represents a collection of polygonal domains (see caribou::topology::Domain) and
     * holds the position of their vertices. Hence, the indices of the nodes of each domains are relative
     * to the vector of positions contained in this Mesh. The mesh is therefore responsible to manage either
     * an internal buffer containing all the nodes, or can be link to an external buffer allocated somewhere.
     * It does this by allowing one to specify the node container type to be used for this buffer. Hence, to get
     * a copy of the position vector, the template parameter ContainerType can be set to a standard dense
     * Eigen::Matrix (the default when no template argument is specified). If, instead, one has already an
     * allocated buffer somewhere, this template parameter can be set to a type that can handle external buffers.
     *
     * Example of a the construction of a 3D mesh where the position of the nodes are copied from an std::vector
     * into the internal buffer of the mesh (i.e. a copy of all the positions is made).
     * \code{.cpp}
     * using Coordinates = Mesh<3>::WorldCoordinates;
     * std::vector<Coordinates> positions = {{0,0,0}, {1,1,1}};
     * Mesh<3> mesh (initial_positions);
     * std::cout << "The mesh has " << mesh.number_of_nodes() << " nodes\n"; // The mesh has 2 nodes
     * std::cout << "First node is " << mesh.position(0) << "\n"; // First node is [0,0,0]
     * \endcode
     *
     * @tparam WorldDimension The dimension (1D, 2D or 3D) of the node positions.
     * @tparam NodeContainerType  Holder type that contains the nodes of the mesh. Let A be a NodeContainerType,
     * it must follow these requirements:
     * - Must be default constructable
     * <table>
     * <caption id="mesh_node_container_inner_types">Inner types</caption>
     * <tr><th>Type-id</th><th>Description</th><th>Requirements</th></tr>
     * <tr><td>A::Scalar</td><td>Scalar type of a position coordinate</td><td></td></tr>
     * </table>
     *
     * <table>
     * <caption id="mesh_node_container_storage">Storage operations</caption>
     * <tr><th>Expression</th><th>Return type</th><th>Description</th><th>Requirements</th></tr>
     * <tr><td>a.resize(n)</td><td>(not used)</td><td>Resize the container 'a' to 'n' nodes.</td><td>If the container
     * already had storage allocated for a number of nodes greater than 'n', the 'n'th nodes must still be accessible to their previous index before this call.</td></tr>
     * <tr><td>a.node(i)</td><td>[const] NodeType [&]</td><td>Obtain a reference or a copy to the ith node</td><td>NodeType
     * must implement the [](IntegerType j) operator which should return a reference or a copy to the jth scalar coordinates
     * of the node</td></tr>
     * <tr><td>a.size()</td><td>unsigned integer type</td><td>Get the number of nodes actually contained inside this container</td>
     * <td></td></tr>
     * </table>
     */

    template <
        unsigned int WorldDimension,
        typename NodeContainerType = EigenNodesHolder<Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, WorldDimension, (WorldDimension>1?Eigen::RowMajor:Eigen::ColMajor)>>
    >
    class Mesh : public BaseMesh {
    public:
        /**
         * Holder type that contains the nodes of the mesh
         */
        using NodeContainer_t = std::decay_t<NodeContainerType>;
        static_assert(WorldDimension == 1 or WorldDimension == 2 or WorldDimension == 3, "The world dimension must be 1, 2 or 3.");

        using Self = Mesh<WorldDimension, NodeContainer_t>;
        static constexpr INTEGER_TYPE Dimension = WorldDimension;
        using Real = typename NodeContainer_t::Scalar;
        using WorldCoordinates = Eigen::Matrix<Real, 1, Dimension>;

        template <typename Element, typename NodeIndex = UNSIGNED_INTEGER_TYPE>
        using Domain = Domain<Element, NodeIndex>;

        /*!
         * Default constructor for an empty mesh.
         */
        Mesh() : p_nodes {}, p_domains {} {}

        /*!
         * Virtual default destructor
         */
        virtual ~Mesh() = default;

        /**
         * Construct the unstructured mesh with a set of point positions from a position vector.
         * @param positions A reference to a std::vector containing the position of the nodes
         *
         * @note A copy is made of all nodes position vectors is made.
         */
        explicit Mesh(const std::vector<WorldCoordinates> & positions)
        {
            const auto n = positions.size();
            p_nodes.resize(n);
            for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
                auto node = this->p_nodes.node(i);
                for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                    node[j] = positions[i][j];
                }
            }
        }

        /**
         * Construct the unstructured mesh with an Eigen matrix containing the position vector
         * (NxD with N nodes of D world dimension).
         *
         * @param positions A reference to a NxD matrix containing the position vector
         */
        template <typename Derived>
        explicit Mesh(const Eigen::MatrixBase<Derived> & positions)
        {
            static_assert(
                Eigen::MatrixBase<Derived>::ColsAtCompileTime == Dimension or Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic,
                "The number of columns at compile time should match the Dimension of the mesh, or by dynamic (known at compile time)."
            );

            // The number of columns must equal the World dimension
            caribou_assert(positions.cols() == Dimension);

            // Do the copy
            const auto n = positions.rows();
            p_nodes.resize(n);
            for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
                auto node = this->p_nodes.node(i);
                for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                    node[j] = positions(i,j);
                }
            }
        }

        /**
         * Construct the unstructured mesh with an Eigen matrix containing the position vector
         * (NxD with N nodes of D world dimension).
         *
         * @param positions A reference to a NxD matrix containing the position vector
         * @param copy If true, a copy of the position is made into an internal buffer. If false,
         *             no copy is made and the user must make sure that the external buffer specified
         *             by the positions parameter will exists as long as this Mesh will.
         *             Default to true, i.e. a copy will be made.
         */
//        explicit Mesh(const MatrixType & positions, bool copy = true)
//        : p_positions(nullptr, positions.rows(), positions.cols(), StrideType(positions.outerStride(), positions.innerStride()))
//        {
//            const auto number_of_nodes = positions.rows();
//            if (copy) {
//                p_buffer = positions;
//                new (&p_positions) PositionsReference (
//                    p_buffer.data(), number_of_nodes, Dimension,
//                    StrideType(positions.outerStride(), positions.innerStride())
//                );
//            } else {
//                new (&p_positions) PositionsReference (
//                    positions.data(), number_of_nodes, Dimension,
//                    StrideType(positions.outerStride(), positions.innerStride())
//                );
//            }
//        }

        /*! Copy constructor */
        Mesh(const Mesh & other)
        : p_nodes (other.p_nodes), p_domains {}
        {
            // Copy domains
            p_domains.reserve(other.p_domains.size());
            for (const auto & p : other.p_domains) {
                const std::string & domain_name = p.first;
                BaseDomain * cloned_domain_ptr  = p.second->clone();
                p_domains.emplace_back(domain_name, cloned_domain_ptr);
            }
        }

        /*! Move constructor */
        Mesh(Mesh && other) noexcept
        : Mesh() {
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
        inline auto number_of_nodes() const -> UNSIGNED_INTEGER_TYPE final {return p_nodes.size();};

        /*!
         * \copydoc caribou::topology::BaseMesh::node
         */
        [[nodiscard]]
        inline auto node(const UNSIGNED_INTEGER_TYPE & node_id) const -> Eigen::Vector3d {
            const auto p = position(node_id);
            Eigen::Vector3d n;
            n[0] = p[0];
            if constexpr (Dimension > 1) {
                n[1] = p[1];
            }
            if constexpr (Dimension > 2) {
                n[2] = p[2];
            }
            return n;
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
            caribou_assert(
                (i < p_domains.size()) &&
                "Trying to get a domain outside from a domain id larger than the number of domains."
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
        template<typename Element, typename NodeIndex, typename... Args>
        inline
        typename std::enable_if_t<std::is_integral_v<NodeIndex>, Domain<Element, NodeIndex> *>
        add_domain(const std::string & name, Args ...args) {
            static_assert(geometry::traits<Element>::Dimension == Dimension, "The dimension of the mesh doesn't match the dimension of the element type.");
            for (const auto & d : p_domains) {
                if (d.first == name) {
                    throw std::domain_error(std::string("A domain with the name ") + name + " already exists.");
                }
            }

            auto domain_ptr = new Domain<Element, NodeIndex>(this, std::forward<Args>(args)...);
            add_domain(domain_ptr, name);

            return domain_ptr;
        }

        /*!
         * Adds a domain with the Element type supplied as a template parameter.
         * @tparam Element The type of Element
         * @tparam NodeIndex The integer type for stored node indices
         * @return A pointer to the newly created domain, or null if the creation failed
         *
         * \note The domain is managed by this mesh instance, and will be freed upon the deletion of the mesh.
         */
        template<typename Element, typename NodeIndex, typename... Args>
        inline
        typename std::enable_if_t<std::is_integral_v<NodeIndex>, Domain<Element, NodeIndex> *>
        add_domain(Args ...args) {
            static_assert(geometry::traits<Element>::Dimension == Dimension, "The dimension of the mesh doesn't match the dimension of the element type.");

            auto domain_ptr = new Domain<Element, NodeIndex>(this, std::forward<Args>(args)...);
            add_domain(domain_ptr);

            return domain_ptr;
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
            return add_domain<Element, UNSIGNED_INTEGER_TYPE, Args...>(name, std::forward<Args>(args)...);
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

            add_domain(domain_ptr);
            return domain_ptr;
        }

        /*!
         * \copydoc caribou::topology::BaseMesh::add_domain
         */
        inline auto add_domain(BaseDomain * domain, const std::string & name) -> BaseDomain * override {
            if (domain->mesh() == this) {
                // Trying to attach the same mesh as the one which is already attached to this domain
                using namespace std;
                auto it = find_if(begin(p_domains), end(p_domains), [&domain](const auto & d) {
                    return (d.second.get() == domain);
                });
                if (it != end(p_domains)) {
                    if ((*it).first != name) {
                        // The domain is already linked to this mesh, but with a different name
                        throw std::logic_error("Trying to add a domain that is already attached to this mesh but with a different name.");
                    }

                    // The domain is already linked to this mesh, with the same name, let's simply return it
                    return domain;
                }
            } else if (domain->mesh() != nullptr) {
                // The domain is already attached to another Mesh
                throw std::logic_error("Trying to add a domain that is already attached to another Mesh instance.");
            } else {
                domain->attach_to(this);
            }

            p_domains.template emplace_back(name, domain);
            return domain;
        }

        /*!
         * \copydoc caribou::topology::BaseMesh::add_domain
         */
        inline auto add_domain(BaseDomain * domain) -> BaseDomain * override {
            std::ostringstream ss;
            ss << "domain_" << domain;
            return add_domain(domain, ss.str());
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
        inline auto position(UNSIGNED_INTEGER_TYPE index) const {
            caribou_assert(static_cast<Eigen::Index>(index) < p_nodes.size()
                           && "Trying to access the position vector at a node index greater "
                              "than the number of nodes in the mesh.");
            return p_nodes.node(index);
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
            Eigen::Matrix<Real, N, Dimension, (Dimension>1?Eigen::RowMajor:Eigen::ColMajor)> positions;
            for (std::size_t i = 0; i < N; ++i) {
                caribou_assert(indices[i] < static_cast<IntegerType>(number_of_nodes()));
                positions.row(i) = p_nodes.node(indices[i]);
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
            return positions<IntegerType, N>(indices.data());
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
            Eigen::Matrix<Real, Eigen::Dynamic, Dimension, (Dimension == 1) ? Eigen::ColMajor : Eigen::RowMajor> positions;
            positions.resize(number_of_indices, Dimension);
            for (std::size_t i = 0; i < number_of_indices; ++i) {
                positions.row(i) = p_nodes.node(indices[i]);
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
            Eigen::Matrix<Real,
                          EigenDerived::RowsAtCompileTime == 1 ? EigenDerived::ColsAtCompileTime : EigenDerived::RowsAtCompileTime,
                          Dimension,
                          (Dimension>1 ? Eigen::RowMajor : Eigen::ColMajor)> positions;
            positions.resize(number_of_indices, Dimension);
            for (Eigen::Index i = 0; i < number_of_indices; ++i) {
                positions.row(i) = p_nodes.node(indices.derived()[i]);
            }
            return positions;
        }

        /*! Swap the data of two meshes */
        friend void swap(Mesh & first, Mesh& second) noexcept
        {
            // enable ADL
            using std::swap;

            swap(first.p_nodes, second.p_nodes);
            swap(first.p_domains, second.p_domains);
        }

    private:

        /**
         * Container of the mesh node positions
         */
        NodeContainer_t p_nodes;

        /**
         * The list of domains of the mesh.
         */
        std::vector<std::pair<std::string, std::unique_ptr<BaseDomain>>> p_domains;
    };
} // namespace caribou::topology


namespace caribou::topology {
    // -------------------------------------------------------------------------------------
    //                                    IMPLEMENTATION
    // -------------------------------------------------------------------------------------
    // Implementation of methods that can be specialized later for an explicit container type
    // -------------------------------------------------------------------------------------


    // -------------------------------------------------------------------------------------
    //                               TEMPLATE DEDUCTION GUIDES
    // -------------------------------------------------------------------------------------
    // Template deduction guides that can help the compiler to automatically determine the
    // template parameters of the Mesh from the constructor used.
    // -------------------------------------------------------------------------------------
    template <int WorldDimension, typename Real>
    Mesh(const std::vector<Eigen::Matrix<Real, 1, WorldDimension>> &) ->
    Mesh <
        WorldDimension,
        EigenNodesHolder<Eigen::Matrix<Real, Eigen::Dynamic, WorldDimension, (WorldDimension>1?Eigen::RowMajor:Eigen::ColMajor)>>
    >;

    template <typename Derived>
    Mesh(const Eigen::MatrixBase<Derived> &) ->
    Mesh <
        Eigen::MatrixBase<Derived>::ColsAtCompileTime,
        EigenNodesHolder<
            Eigen::Matrix<
                typename Eigen::MatrixBase<Derived>::Scalar,
                Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                (Eigen::MatrixBase<Derived>::ColsAtCompileTime>1?Eigen::RowMajor:Eigen::ColMajor)
            >
        >
    >;
} // namespace caribou::topology
