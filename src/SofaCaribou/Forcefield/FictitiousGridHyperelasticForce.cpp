#include <SofaCaribou/Forcefield/FictitiousGridHyperelasticForce.h>
#include <SofaCaribou/Forcefield/HyperelasticForcefield.inl>

#include <sofa/core/ObjectFactory.h>

namespace SofaCaribou::forcefield {
using namespace caribou::geometry;

// ---------------------------
// Subdivided Gauss Hexahedron
// ---------------------------
template<>
 auto FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>::templateName(const FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>*) -> std::string {
    return "SubdividedGaussHexahedron";
}

template<>
auto FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>::canCreate(FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>*, BaseContext*, BaseObjectDescription* arg) -> bool {
    std::string method = arg->getAttribute("integration_method", "SubdividedVolume");
    return (method == "SubdividedGauss");
}

template <>
auto FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>::get_gauss_nodes(const std::size_t & element_id, const SubdividedGaussHexahedron & e) const -> GaussContainer {
    // GaussContainer is a std::vector<GaussNode>
    GaussContainer gauss_nodes {};
    const auto g_nodes = d_grid->get_gauss_nodes_of_cell(element_id);
    gauss_nodes.resize(gauss_nodes.size());
    for (std::size_t gauss_node_id = 0; gauss_node_id < gauss_nodes.size(); ++gauss_node_id) {
        const auto & gauss_node = g_nodes[gauss_node_id].first;
        const auto & gauss_weight = g_nodes[gauss_node_id].second;
        const auto J = e.jacobian(gauss_node);
        const Mat33 Jinv = J.inverse();

        gauss_nodes[gauss_node_id].weight = gauss_weight;
        gauss_nodes[gauss_node_id].jacobian_determinant = 1;
        gauss_nodes[gauss_node_id].dN_dx  = (Jinv.transpose() * e.dL(gauss_node).transpose()).transpose();
        gauss_nodes[gauss_node_id].F      = Mat33::Identity();
    }
    return gauss_nodes;
}


// ----------------------------
// Subdivided Volume Hexahedron
// ----------------------------
template<>
auto FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>::templateName(const FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>*) -> std::string {
    return "SubdividedVolumeHexahedron";
}

template<>
auto FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>::canCreate(FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>*, BaseContext*, BaseObjectDescription* arg) -> bool {
    std::string method = arg->getAttribute("integration_method", "SubdividedVolume");
    return (method == "SubdividedVolume");
}

template <>
auto FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>::get_gauss_nodes(const std::size_t & element_id, const SubdividedVolumeHexahedron & e) const -> GaussContainer {
    // GaussContainer is a std::array<GaussNode, 8>
    GaussContainer gauss_nodes {};
    const auto g_nodes = d_grid->get_gauss_nodes_of_cell(element_id, 0);
    for (std::size_t gauss_node_id = 0; gauss_node_id < gauss_nodes.size(); ++gauss_node_id) {
        const auto & gauss_node = g_nodes[gauss_node_id].first;
        const auto & gauss_weight = g_nodes[gauss_node_id].second;
        const auto J = e.jacobian(gauss_node);
        const Mat33 Jinv = J.inverse();

        gauss_nodes[gauss_node_id].weight = gauss_weight;
        gauss_nodes[gauss_node_id].jacobian_determinant = 1;
        gauss_nodes[gauss_node_id].dN_dx  = (Jinv.transpose() * e.dL(gauss_node).transpose()).transpose();
    }
    return gauss_nodes;
}

template class FictitiousGridHyperelasticForcefield<caribou::geometry::SubdividedVolumeHexahedron>;
template class FictitiousGridHyperelasticForcefield<caribou::geometry::SubdividedGaussHexahedron>;

static int FictitiousGridHyperelasticForceClass = RegisterObject("Caribou Fictitious grid hyperelastic FEM Forcefield")
    .add< FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron> >(true)
    .add< FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>  >()
;

} // namespace SofaCaribou::forcefield