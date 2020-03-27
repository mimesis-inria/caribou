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
    namespace internal {
        template <typename Element, class Enable = void>
        class ElementStorage;
    }

    /*!
     * A domain implements the BaseDomain interface with a given element type.
     *
     * @tparam Element See caribou::geometry::is_an_element
     */
    template <typename Element>
    class Domain : public BaseDomain, private internal::ElementStorage<Element> {
    public:
        using Storage = internal::ElementStorage<Element>;
        /*!
         * If the Element type as a static number of nodes at compile time, than Indices is
         *  std::vector<std::array<unsigned int, Element::NumberOfNodesAtCompileTime>>
         * else, if the Element type as a dynamic number of nodes (known at runtime), than Indices is
         *  std::vector<std::vector<unsigned int>>
         */
        using Indices = typename Storage::Indices;

        /*! Destructor */
        ~Domain() final = default;

        /*! Empty constructor */
        explicit Domain() {}

        explicit Domain(const std::vector<Indices> & elements)
        : Storage::p_elements(elements) {}

        /*! copy-and-swap move constructor */
        Domain(Domain<Element> && other) noexcept {
            swap(*this, other);
        }

        /*! copy-and-swap assigment */
        auto operator=(Domain<Element> other) noexcept -> Domain & {
            swap(*this, other);
            return *this;
        }

        /*!
         * \copydoc caribou::topology::BaseDomain::canonical_dimension
         */
        [[nodiscard]] auto canonical_dimension() const -> UNSIGNED_INTEGER_TYPE final {
            return geometry::traits<Element>::CanonicalDimension;
        }

        /*!
         * \copydoc caribou::topology::BaseDomain::number_of_nodes_per_elements
         */
        [[nodiscard]] auto number_of_nodes_per_elements() const -> INTEGER_TYPE final {
            return geometry::traits<Element>::NumberOfNodesAtCompileTime;
        }

        /*!
         * \copydoc caribou::topology::BaseDomain::number_of_nodes_of_element
         */
        [[nodiscard]] auto number_of_nodes_of_element(const UNSIGNED_INTEGER_TYPE & element_id) const -> UNSIGNED_INTEGER_TYPE final {
            caribou_assert(element_id >= 0 and element_id < number_of_elements());
            if constexpr (geometry::traits<Element>::NumberOfNodesAtCompileTime == caribou::Dynamic) {
                return this->p_elements[element_id].size();
            } else {
                return static_cast<UNSIGNED_INTEGER_TYPE>(geometry::traits<Element>::NumberOfNodesAtCompileTime);
            }
        }

        /*!
         * \copydoc caribou::topology::BaseDomain::number_of_elements
         */
        [[nodiscard]] auto number_of_elements() const -> UNSIGNED_INTEGER_TYPE final {
            return this->p_elements.size();
        };

        /*!
         * \copydoc caribou::topology::BaseDomain::add_element
         */
        void add_element(const UNSIGNED_INTEGER_TYPE * node_indices, UNSIGNED_INTEGER_TYPE number_of_nodes) final {
            if (geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic and
                geometry::traits<Element>::NumberOfNodesAtCompileTime != number_of_nodes) {
                throw std::runtime_error(
                    "Trying to add an element of " + std::to_string(geometry::traits<Element>::NumberOfNodesAtCompileTime) +
                    " nodes with an array containing " + std::to_string(number_of_nodes) + " node indices.");
            } else if (number_of_nodes == 0) {
                throw std::runtime_error("Trying to add an element with an empty set of node indices.");
            }

            auto last_inserted = this->p_elements.insert(this->p_elements.end(), Indices());
            for (UNSIGNED_INTEGER_TYPE i = 0; i < number_of_nodes; ++i) {
                auto & v = *last_inserted;
                v[i] = node_indices[i];
            }
        }

        /*!
         * Add an element to the domain.
         * @param indices The indices of its nodes in the mesh's position vector.
         */
        inline void add(const Indices & element) {
            this->p_elements.push_back(element);
        }

        /*!
         * Add elements to the domain.
         * @param indices The indices of the nodes of the added elements.
         */
        inline void add(const std::vector<Indices> & elements) {
            this->p_elements.insert(this->p_elements.begin(), elements.begin(), elements.end());
        }

        /*!
         * Get the indices of an element in the domain.
         * @param index The index of the element.
         * @return The indices of its nodes in the mesh's position vector as either a vector or an array:
         *   1. If the element type as a static number of nodes, then a const std::array<uint, number_of_nodes> & is returned
         *   2. If the element type as a dynamic number of nodes, then a const std::vector<uint> & is returned
         *
         */
        inline auto element_indices(const UNSIGNED_INTEGER_TYPE & index) const -> const Indices & {
            caribou_assert(index < this->p_elements.size() and "Trying to get an element that does not exists.");

            return this->p_elements[index];
        }
    };

    namespace internal {
        template<typename Element, class Enable>
        class ElementStorage {
        protected:
            using Indices = std::vector<UNSIGNED_INTEGER_TYPE>;
            std::vector<Indices> p_elements;
        };

        template<typename Element>
        class ElementStorage<Element, CLASS_REQUIRES(geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic)> {
        protected:
            using Indices = std::array<UNSIGNED_INTEGER_TYPE, geometry::traits<Element>::NumberOfNodesAtCompileTime>;
            std::vector<Indices> p_elements;
        };
    }
}

#endif //CARIBOU_TOPOLOGY_DOMAIN_H
