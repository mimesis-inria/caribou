#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Topology/FictitiousGrid.h>

#include <sofa/helper/OptionsGroup.h>


namespace caribou::geometry {
// --------------------------------------------------------------------------------
// Subdivided Gauss Hexahedron Element
//---------------------------------------------------------------------------------
// Regular hexahedron that is cut by a boundary, hence it is only partially filled.
// The set of Gauss integration nodes are obtained by recursively subdividing each
// cells into 8 sub-cells. The Gauss nodes of the leaf sub-cells are used.
//---------------------------------------------------------------------------------
struct SubdividedGaussHexahedron;
template<>
struct traits<SubdividedGaussHexahedron> : traits<caribou::geometry::Hexahedron<caribou::Linear>> {
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = caribou::Dynamic;
};

struct SubdividedGaussHexahedron : public caribou::geometry::Hexahedron<caribou::Linear> {
    using Base = caribou::geometry::Hexahedron<caribou::Linear>;
    using Base::Base;
};

// --------------------------------------------------------------------------------
// Subdivided Volume Hexahedron Element
//---------------------------------------------------------------------------------
// Regular hexahedron that is cut by a boundary, hence it is only partially filled.
// The set of Gauss integration nodes are the same as a regular hexahedral elements.
// The weight of each Gauss integration nodes however are computed by recursively
// subdividing each cells into 8 sub-cells.
//---------------------------------------------------------------------------------
struct SubdividedVolumeHexahedron;
template<>
struct traits<SubdividedVolumeHexahedron> : traits<caribou::geometry::Hexahedron<caribou::Linear>> {};

struct SubdividedVolumeHexahedron : public caribou::geometry::Hexahedron<caribou::Linear> {
    using Base = caribou::geometry::Hexahedron<caribou::Linear>;
    using Base::Base;
};
} // namespace caribou::geometry

namespace SofaCaribou::forcefield {

/**
 * Hyperelastic forcefield for cut hexahedral elements.
 * @tparam Element Can be either SubdividedGaussHexahedron or SubdividedVolumeHexahedron
 */
template <typename Element>
class FictitiousGridHyperelasticForcefield : public HyperelasticForcefield<Element> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(FictitiousGridHyperelasticForcefield, Element), SOFA_TEMPLATE(HyperelasticForcefield, Element));

    using Base = HyperelasticForcefield<Element>;
    using FictitiousGrid = SofaCaribou::topology::FictitiousGrid<sofa::defaulttype::Vec3Types>;

    using GaussContainer = typename Base::GaussContainer;
    using Index = typename Base::Index;

    template <typename ObjectType>
    using Link = SingleLink<FictitiousGridHyperelasticForcefield<Element>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Constructor
    FictitiousGridHyperelasticForcefield()
    : Base::Base()
    , d_grid(initLink("fictitious_grid", "Fictitious grid containing the elements on which this forcefield will be applied."))
    , d_integration_method(initData(&d_integration_method,
        "integration_method",
        R"(
                Integration method used to integrate the stiffness matrix.

                Methods are:
                  SubdividedVolume: Hexas are recursively subdivided into cubic subcells and these subcells are used to
                                    compute the inside volume of the regular hexa's gauss points.
                  SubdividedGauss:  Hexas are recursively subdivided into cubic subcells and these subcells are used to
                                    add new gauss points. Gauss points outside of the boundary are ignored.
        )",
        true /*displayed_in_GUI*/, true /*read_only_in_GUI*/))
    {
        d_integration_method.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
           "SubdividedVolume", "SubdividedGauss"
        }));
        sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> integration_method = d_integration_method;
        integration_method->setSelectedItem((unsigned int) 0);
    }

    void init() override {
        if (d_grid.get() and number_of_elements() == 0) {
            msg_warning() << "No element found in the grid '" << d_grid->getPathName() << "'.";
        } else if (not d_grid.get()) {
            // No topology specified. Try to find one suitable.
            auto containers = this->getContext()->template getObjects<FictitiousGrid>(BaseContext::Local);
            if (containers.empty()) {
                msg_warning() << "Could not find a fictitious grid in the current context.";
            } else {
                std::vector<FictitiousGrid *> suitable_containers;
                for (const auto & container : containers) {
                    d_grid.set(container);
                    if (number_of_elements() > 0) {
                        suitable_containers.push_back(container);
                    }
                    d_grid.set(nullptr);
                }

                if (suitable_containers.empty()) {
                    msg_warning() << "Could not find a suitable fictitious grid in the current context, they are all empty.";
                } else if (suitable_containers.size() > 1) {
                    d_grid.set(suitable_containers[0]);
                    msg_warning() <<
                                  "Multiple fictitious grid were found in the context node. The first one will be taken ('" <<  d_grid.get()->getPathName() << "'). " <<
                                  "If it isn't the correct one, please specify which one contains the elements on which this force field will be applied " <<
                                  "by explicitly setting the container's path in the  '" << d_grid.getName() << "'  parameter.";
                    msg_info() << "Automatically found the topology '" << d_grid.get()->getPathName() << "'.";
                } else {
                    d_grid.set(suitable_containers[0]);
                    msg_info() << "Automatically found the topology '" << d_grid.get()->getPathName() << "'.";
                }
            }
        }

        Base::init();
    }

    static auto templateName(const FictitiousGridHyperelasticForcefield<Element>* = nullptr) -> std::string;

    static auto canCreate(FictitiousGridHyperelasticForcefield<Element>*, BaseContext*, BaseObjectDescription* arg) -> bool;

    /** Get the number of elements contained in this field **/
    [[nodiscard]]
    inline auto number_of_elements() const -> std::size_t override {
        if (d_grid.get()) {
            return d_grid->number_of_cells();
        }
        return 0;
    }

private:
    // These private methods are implemented but can be overridden

    /** Get the element nodes indices relative to the state vector */
    auto get_element_nodes_indices(const std::size_t & element_id) const -> const Index * override {
        if (not d_grid.get() or d_grid->number_of_cells() == 0) {
            return nullptr;
        }

        return &(d_grid->get_node_indices_of(element_id)[0]);
    }

    /**
     * Return true if the mesh topology is compatible with the type Element.
     *
     * Always false since this forcefield needs a Fictitious Grid and not a Sofa mesh topology.
     */
    inline
    static auto mesh_is_compatible(const sofa::core::topology::BaseMeshTopology *) -> bool {
        return false;
    }

    /** Get the set of Gauss integration nodes of the given element */
    auto get_gauss_nodes(const std::size_t & /*element_id*/, const Element & /*element*/) const -> GaussContainer override {
        return {};
    }

    // Data members
    Link<FictitiousGrid> d_grid;
    Data< sofa::helper::OptionsGroup > d_integration_method;
};


} // namespace SofaCaribou::forcefield