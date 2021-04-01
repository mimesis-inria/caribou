#pragma once

#include <Caribou/Topology/Domain.h>

namespace caribou::topology {
/**
 * Implementation of a Domain<Element, NodeIndex> that can construct a geometrical element
 * using the position of the nodes of this Mesh.
 *
 * @sa class caribou::topology::Domain
 */
template <unsigned int WorldDimension, typename NodeContainerType>
template <typename Element, typename NodeIndex>
class Mesh<WorldDimension, NodeContainerType>::MeshDomain : public ::caribou::topology::Domain <Element, NodeIndex> {
friend Mesh;
public:

    // Import base constructors and aliases
    using Base = ::caribou::topology::Domain <Element, NodeIndex>;
    using ElementsIndices = typename Base::ElementsIndices;

    /*! copy-and-swap assigment (valid for both copy and move assigment) */
    auto operator=(MeshDomain other) noexcept -> MeshDomain & {
        this->swap(*this, other);
        return *this;
    }

    inline auto element(const UNSIGNED_INTEGER_TYPE &element_id) const -> Element override {
        caribou_assert(element_id < this->number_of_elements() &&
                       ("Trying to get the element #" + std::to_string(element_id) + ", but the domain only has " +
                        std::to_string(this->number_of_elements()) + " elements.").c_str()
        );

        return Element(p_mesh->positions(this->element_indices(element_id)));
    }

    BaseDomain * clone() const override {
        return new MeshDomain(*this); // Will call the copy constructor
    }

    inline auto mesh() const -> const BaseMesh * override { return p_mesh; }

private:
    MeshDomain(Mesh * mesh) : p_mesh(mesh) {}

    MeshDomain(Mesh * mesh, const ElementsIndices & elements)
        : Base (elements), p_mesh(mesh) {}

    MeshDomain(Mesh * mesh, ElementsIndices & elements)
        : Base (elements), p_mesh(mesh) {}

    MeshDomain(Mesh * mesh, const ElementsIndices * elements)
        : Base (elements), p_mesh(mesh) {}

    MeshDomain(Mesh * mesh, const NodeIndex * data, const Eigen::Index & number_of_elements, const Eigen::Index & number_of_nodes_per_elements)
        : Base (data, number_of_elements, number_of_nodes_per_elements), p_mesh(mesh) {}

    MeshDomain(Mesh * mesh, const NodeIndex * data, const Eigen::Index & number_of_elements, const Eigen::Index & number_of_nodes_per_elements, Eigen::Index outer_stride, Eigen::Index inner_stride)
        : Base (data, number_of_elements, number_of_nodes_per_elements, outer_stride, inner_stride), p_mesh(mesh) {}


    /// The mesh associated with this domain.
    const Mesh * p_mesh;
};
} // namespace caribou::topology

