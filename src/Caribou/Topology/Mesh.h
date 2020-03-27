#ifndef CARIBOU_TOPOLOGY_MESH_H
#define CARIBOU_TOPOLOGY_MESH_H

#include <Caribou/config.h>
#include <Caribou/Topology/BaseMesh.h>
#include <Caribou/Topology/Domain.h>

#include <array>
#include <memory>
#include <vector>

namespace caribou::topology {

    template <typename Derived>
    class Mesh : public BaseMesh {
    public:

        /*!
         * \copydoc caribou::topology::BaseMesh::dimension
         */
        [[nodiscard]] inline auto dimension() const -> UNSIGNED_INTEGER_TYPE final {
            return Derived::Dimension;
        };

        /*!
         * \copydoc caribou::topology::BaseMesh::number_of_domains
         */
        [[nodiscard]]
        inline auto number_of_domains() const -> UNSIGNED_INTEGER_TYPE final {
            return p_domains.size();
        }

        /*!
         * Get the list of domains of the mesh.
         */
        inline std::vector<std::pair<std::string, const BaseDomain *>> domains() const {
            std::vector<std::pair<std::string, const BaseDomain *>> domains;
            domains.reserve(p_domains.size());
            for (const auto & d : p_domains) {
                domains.emplace_back(d.first, d.second.get());
            }
            return domains;
        }

        /*!
         * Adds a domain with the Element type supplied as a template parameter.
         * @tparam Element The type of Element
         * @param name
         * @return A pointer to the newly created domain, or null if the creation failed
         *
         * \note The domain is managed by this mesh instance, and will be freed upon the deletion of the mesh.
         */
        template<typename Element>
        inline auto add_domain(const std::string & name) -> Domain<Element> * {
            for (const auto & d : p_domains) {
                if (d.first == name) {
                    throw std::domain_error(std::string("A domain with the name ") + name + " already exists.");
                }
            }

            auto domain_ptr = new Domain<Element>();
            auto res = p_domains.emplace(p_domains.end(), name, std::unique_ptr<BaseDomain>(static_cast<BaseDomain *>(domain_ptr)));
            if (res->second) {
                return domain_ptr;
            } else {
                delete domain_ptr;
                return static_cast<Domain<Element>*> (nullptr);
            }
        }

        /*!
        * Get the position coordinates of a node from its index.
        */
        [[nodiscard]] inline auto position(UNSIGNED_INTEGER_TYPE index) const {
            return self()._position(index);
        }

        /*!
         * Get a position vector from the given indices.
         * @tparam IntegerType Integer type of the indices passed
         * @param indices [IN] The indices from which the positions vector will be filled
         * @return The position coordinates of every nodes queried.
         *
         * Example:
         * \code{.cpp}
         * auto positions = mesh.positions({1, 2});
         * \endcode
         *
         */
        template <typename IntegerType, std::size_t N>
        inline auto positions(const IntegerType(&indices)[N]) const {
            return self()._positions(indices);
        }

        /*!
         * Get a position vector or array from the given indices as a vector or array, respectively.
         * @param indices [IN] The indices from which the positions vector or array will be filled
         * @return The position coordinates of every nodes queried.
         *
         * Example from an array of indices:
         * \code{.cpp}
         * std::array<int, 2> indices = {1, 2};
         * auto positions = mesh.positions(indices);
         * \endcode
         *
         * Example from a vector of indices:
         * \code{.cpp}
         * std::vector<int> indices = {1, 2};
         * auto positions = mesh.positions(indices);
         * \endcode
         *
         */
        template <typename T>
        inline auto positions(T && indices) const {
            return self()._positions(std::forward<T>(indices));
        }

    private:
        auto self () const -> const Derived &
        {
            return static_cast<const Derived &>(*this);
        }


        std::vector<std::pair<std::string, std::unique_ptr<BaseDomain>>> p_domains;
    };
}

#endif //CARIBOU_TOPOLOGY_MESH_H
