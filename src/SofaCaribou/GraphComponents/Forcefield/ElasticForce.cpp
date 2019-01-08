#include <numeric>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/Node.h>
#include <SofaBaseTopology/EdgeSetTopologyContainer.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <SofaBaseTopology/QuadSetTopologyContainer.h>
#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>

#include "ElasticForce.h"


namespace SofaCaribou::GraphComponents::forcefield {

using namespace sofa::component::topology;

template<class DataTypes>
ElasticForce<DataTypes>::ElasticForce()
: d_youngModulus(initData(&d_youngModulus,
        Real(1000), "youngModulus", "Young's modulus of the material", true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_poissonRatio(initData(&d_poissonRatio,
        Real(0.3),  "poissonRatio", "Poisson's ratio of the material", true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_element_type(initData(&d_element_type,
        "element_type","The type of elements on which this force field will be computed. If it is not specified, the "
                       "type will be deduced from the topological container. Element types are "
                       "\"Edge\", \"Triangle\", \"Quad\", \"Tetrahedron\" or \"Hexahedron\"."
                       , true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_topology_container(initLink(
        "topology_container", "Topology that contains the elements on which this force will be computed."))
{
    sofa::helper::WriteAccessor<Data<sofa::helper::OptionsGroup>> elementTypeOptions = d_element_type;
    elementTypeOptions->setNames(5,"Automatic", "Edge", "Triangle", "Quad", "Tetrahedron", "Hexahedron");
    elementTypeOptions->setSelectedItem(0);
}

template<class DataTypes>
void ElasticForce<DataTypes>::init()
{
    Inherit::init();
    if (not d_topology_container.get()) {
        auto containers = this->getContext()->template getObjects<TopologyContainer>(BaseContext::Local);
        auto node = static_cast<const sofa::simulation::Node *> (this->getContext());
        if (containers.empty()) {
            msg_error() << "No topology container were found in the context node '" << node->getPathName() << "'.";
        } else if (containers.size() > 1) {
            msg_error() <<
            "Multiple topology containers were found in the node '" << node->getPathName() << "'." <<
            " Please specify which one contains the elements on which this force field will be computed " <<
            "by explicitly setting the container's path in the 'topology_container' parameter.";
        } else {
            d_topology_container.set(containers[0]);
            msg_info() << "Automatically found the topology container '" << d_topology_container.get()->getPathName() << "'.";
        }
    }

    if (!this->mstate.get())
        msg_error() << "No mechanical state set. Please add a mechanical object in the current node or specify a path "
                       "to one in the 'mstate' parameter.";


    reinit();
}

template<class DataTypes>
void ElasticForce<DataTypes>::reinit()
{
    sofa::core::topology::TopologyContainer * topology = d_topology_container.get();
    if (!topology)
        return;

    sofa::helper::ReadAccessor<Data<sofa::helper::OptionsGroup>> options = d_element_type;
    auto id = options->getSelectedId();
    if (id == 0 || (id < 0 || id >= options->size())) {
        m_element_type = ElementType::Invalid;

        if (id < 0 || id >= options->size()) {
            std::ostringstream oss;
            oss << "Wrong element type specified. The choices are: [";
            for (size_t i = 0; i < options->size(); ++i) {
                oss << d_element_type.virtualGetValue()[i] << " (" << std::to_string(i) << ")";
                if (i+1 < options->size())
                    oss << ", ";
            }
            oss << "]";
            msg_error() << oss.str();
        }

        // Automatically deduce the element type from the topology container type
        if (dynamic_cast<HexahedronSetTopologyContainer*>(topology))
            m_element_type = ElementType::Hexahedron;
        else if (dynamic_cast<TetrahedronSetTopologyContainer*>(topology))
            m_element_type = ElementType::Tetrahedron;
        else if (dynamic_cast<TriangleSetTopologyContainer*>(topology))
            m_element_type = ElementType::Triangle;
        else if (dynamic_cast<QuadSetTopologyContainer*>(topology))
            m_element_type = ElementType::Quad;
        else if (dynamic_cast<EdgeSetTopologyContainer*>(topology))
            m_element_type = ElementType::Edge;
        else {
            msg_warning() << "Failed to deduced the element type from the topology container type ("
            << topology->getClassName() << "). Trying to deduce from the elements contained inside the container...";
            if (topology->getNbHexahedra() > 0)
                m_element_type = ElementType::Hexahedron;
            else if (topology->getNbTetrahedra() > 0)
                m_element_type = ElementType::Tetrahedron;
            else if (topology->getNbTriangles() > 0)
                m_element_type = ElementType::Triangle;
            else if (topology->getNbQuads() > 0)
                m_element_type = ElementType::Quad;
            else if (topology->getNbEdges() > 0)
                m_element_type = ElementType::Edge;
            else
                msg_error() << "Failed to deduced the element type from the elements contained inside the container.";
        }

        if (m_element_type != ElementType::Invalid)
            msg_info() << "Element type automatically deduced to '" << element_type_string() << "'.";
    } else {
        m_element_type = element_type(id);
    }
}

template<class DataTypes>
void ElasticForce<DataTypes>::addForce(
        const MechanicalParams * /*mparams*/,
        Data<VecDeriv> & d_f,
        const Data<VecCoord> & d_x,
        const Data<VecDeriv> & /*d_v*/)
{
    if (m_element_type == ElementType::Invalid)
        return;

    if (m_element_type == ElementType::Hexahedron) {

    }

}

template<class DataTypes>
void ElasticForce<DataTypes>::computeBBox(const sofa::core::ExecParams* params, bool onlyVisible)
{
    if( !onlyVisible ) return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = this->mstate->read(sofa::core::VecCoordId::position());

    static const Real max_real = std::numeric_limits<Real>::max();
    static const Real min_real = std::numeric_limits<Real>::lowest();
    Real maxBBox[3] = {min_real,min_real,min_real};
    Real minBBox[3] = {max_real,max_real,max_real};
    for (size_t i=0; i<x.size(); i++)
    {
        for (int c=0; c<3; c++)
        {
            if (x[i][c] > maxBBox[c]) maxBBox[c] = (Real)x[i][c];
            else if (x[i][c] < minBBox[c]) minBBox[c] = (Real)x[i][c];
        }
    }

    this->f_bbox.setValue(params,sofa::defaulttype::TBoundingBox<Real>(minBBox,maxBBox));
}

template<class DataTypes>
void ElasticForce<DataTypes>::draw(const sofa::core::visual::VisualParams* vparams)
{
    sofa::core::topology::TopologyContainer * topology = d_topology_container.get();
    if (!topology)
        return;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);

    vparams->drawTool()->disableLighting();

    const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    if (m_element_type == ElementType::Hexahedron) {
        for (auto node_indices : topology->getHexahedra()) {
            std::vector< sofa::defaulttype::Vector3 > points[6];

            auto a = node_indices[0];
            auto b = node_indices[1];
            auto d = node_indices[3];
            auto c = node_indices[2];
            auto e = node_indices[4];
            auto f = node_indices[5];
            auto h = node_indices[7];
            auto g = node_indices[6];


            Coord center = (x[a]+x[b]+x[c]+x[d]+x[e]+x[g]+x[f]+x[h])*0.125;
            Real percentage = 0.15;
            Coord pa = x[a]-(x[a]-center)*percentage;
            Coord pb = x[b]-(x[b]-center)*percentage;
            Coord pc = x[c]-(x[c]-center)*percentage;
            Coord pd = x[d]-(x[d]-center)*percentage;
            Coord pe = x[e]-(x[e]-center)*percentage;
            Coord pf = x[f]-(x[f]-center)*percentage;
            Coord pg = x[g]-(x[g]-center)*percentage;
            Coord ph = x[h]-(x[h]-center)*percentage;



            points[0].push_back(pa);
            points[0].push_back(pb);
            points[0].push_back(pc);
            points[0].push_back(pa);
            points[0].push_back(pc);
            points[0].push_back(pd);

            points[1].push_back(pe);
            points[1].push_back(pf);
            points[1].push_back(pg);
            points[1].push_back(pe);
            points[1].push_back(pg);
            points[1].push_back(ph);

            points[2].push_back(pc);
            points[2].push_back(pd);
            points[2].push_back(ph);
            points[2].push_back(pc);
            points[2].push_back(ph);
            points[2].push_back(pg);

            points[3].push_back(pa);
            points[3].push_back(pb);
            points[3].push_back(pf);
            points[3].push_back(pa);
            points[3].push_back(pf);
            points[3].push_back(pe);

            points[4].push_back(pa);
            points[4].push_back(pd);
            points[4].push_back(ph);
            points[4].push_back(pa);
            points[4].push_back(ph);
            points[4].push_back(pe);

            points[5].push_back(pb);
            points[5].push_back(pc);
            points[5].push_back(pg);
            points[5].push_back(pb);
            points[5].push_back(pg);
            points[5].push_back(pf);


            vparams->drawTool()->setLightingEnabled(false);
            vparams->drawTool()->drawTriangles(points[0], sofa::defaulttype::Vec<4,float>(0.7f,0.7f,0.1f,1.0f));
            vparams->drawTool()->drawTriangles(points[1], sofa::defaulttype::Vec<4,float>(0.7f,0.0f,0.0f,1.0f));
            vparams->drawTool()->drawTriangles(points[2], sofa::defaulttype::Vec<4,float>(0.0f,0.7f,0.0f,1.0f));
            vparams->drawTool()->drawTriangles(points[3], sofa::defaulttype::Vec<4,float>(0.0f,0.0f,0.7f,1.0f));
            vparams->drawTool()->drawTriangles(points[4], sofa::defaulttype::Vec<4,float>(0.1f,0.7f,0.7f,1.0f));
            vparams->drawTool()->drawTriangles(points[5], sofa::defaulttype::Vec<4,float>(0.7f,0.1f,0.7f,1.0f));

        }
    }
}

} // namespace SofaCaribou::GraphComponents::forcefield

SOFA_DECL_CLASS(FEMForcefield)
static int FEMForcefieldClass = sofa::core::RegisterObject("Caribou FEM Forcefield")
                                        .add< SofaCaribou::GraphComponents::forcefield::ElasticForce<sofa::defaulttype::Vec3dTypes> >(true)
;