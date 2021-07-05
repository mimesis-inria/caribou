#include <SofaCaribou/config.h>
#include <SofaCaribou/Forcefield/TetrahedronElasticForce.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/AdvancedTimer.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 200699)
namespace sofa { using Size = unsigned int; }
#endif

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210699)
namespace sofa::type {
using Vector3 = ::sofa::defaulttype::Vector3;
using Mat3x3 = ::sofa::defaulttype::Mat3x3;
template <sofa::Size N, sofa::Size M, typename Real>
using Mat = ::sofa::defaulttype::Mat<N, M, Real>;
template <sofa::Size N, typename Real>
using Vec = ::sofa::defaulttype::Vec<N, Real>;
using RGBAColor = ::sofa::helper::types::RGBAColor;
template <typename Real>
using TBoundingBox = ::sofa::defaulttype::TBoundingBox<Real> ;
}
#endif

#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

namespace SofaCaribou::forcefield {
using namespace caribou::mechanics;

TetrahedronElasticForce::TetrahedronElasticForce()
: d_youngModulus(initData(&d_youngModulus,
    Real(1000), "youngModulus",
    "Young's modulus of the material",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_poissonRatio(initData(&d_poissonRatio,
    Real(0.3),  "poissonRatio",
    "Poisson's ratio of the material",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_corotated(initData(&d_corotated,
    bool(true), "corotated",
    "Whether or not to use corotated elements for the strain computation.",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_topology_container(initLink(
    "topology_container", "Topology that contains the elements on which this force will be computed."))
{
}

void TetrahedronElasticForce::init()
{
    Inherit::init();
    if (not d_topology_container.get()) {
        auto containers = this->getContext()->template getObjects<BaseMeshTopology>(BaseContext::Local);
        auto node = dynamic_cast<const sofa::simulation::Node *> (this->getContext());
        if (containers.empty()) {
            msg_error() << "No topology were found in the context node '" << node->getPathName() << "'.";
        } else if (containers.size() > 1) {
            msg_error() <<
                        "Multiple topology were found in the node '" << node->getPathName() << "'." <<
                        " Please specify which one contains the elements on which this force field will be computed " <<
                        "by explicitly setting the container's path in the 'topology_container' parameter.";
        } else {
            d_topology_container.set(containers[0]);
            msg_info() << "Automatically found the topology '" << d_topology_container.get()->getPathName() << "'.";
        }
    }

    if (!this->mstate.get())
        msg_error() << "No mechanical state set. Please add a mechanical object in the current node or specify a path "
                       "to one in the 'mstate' parameter.";


    reinit();
}

void TetrahedronElasticForce::reinit()
{
    sofa::core::topology::BaseMeshTopology * topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (topology->getNbTetrahedra() == 0) {
        msg_warning() << "The topology container ('" << topology->getPathName() << "') does not contain any tetrahedron.";
        return;
    }

    const sofa::helper::ReadAccessor<Data<VecCoord>> X = state->readRestPositions();

    p_stiffness_matrices.resize(0);
    p_quadrature_nodes.resize(0);

    // Make sure every node of the tetrahedron have its coordinates inside the mechanical state vector
    for (sofa::Index tetrahedron_id = 0; tetrahedron_id < topology->getNbTetrahedra(); ++tetrahedron_id) {
        const auto & node_indices = topology->getTetrahedron(tetrahedron_id);
        for (sofa::Index j = 0; j < NumberOfNodes; ++j) {
            const auto & node_id = node_indices[j];
            if (node_id > X.size()-1) {
                msg_error() << "Some tetrahedrons have node indices outside of the state's position vector. Make sure "
                               "that the mechanical object '" << state->getPathName() << "' contains the position of "
                                                                                         "every tetrahedron nodes.";
                return;
            }
        }
    }

    p_initial_rotation.resize(topology->getNbTetrahedra(), Rotation::Identity());
    p_current_rotation.resize(topology->getNbTetrahedra(), Rotation::Identity());

    // Initialize the initial frame of each tetrahedron
    if (d_corotated.getValue()) {
        for (std::size_t tetrahedron_id = 0; tetrahedron_id < topology->getNbTetrahedra(); ++tetrahedron_id) {
            Tetrahedron tetra = tetrahedron(tetrahedron_id, X);
            p_initial_rotation[tetrahedron_id] = tetra.frame();
        }
    }

    // Gather the integration points for each tetrahedron
    p_quadrature_nodes.resize(topology->getNbTetrahedra());
    for (std::size_t tetrahedron_id = 0; tetrahedron_id < topology->getNbTetrahedra(); ++tetrahedron_id) {
        auto   tetra = tetrahedron(tetrahedron_id, X);
        auto & gauss_node = tetra.gauss_node(0);
        const auto J = tetra.jacobian(gauss_node.position);
        const Mat33 Jinv = J.inverse();
        const auto detJ = J.determinant();
        const  Matrix<NumberOfNodes, 3> dN_dx = (Jinv.transpose() * tetra.dL(gauss_node.position).transpose()).transpose();

        p_quadrature_nodes[tetrahedron_id] = {
            gauss_node.weight,
            detJ,
            dN_dx
        };
    }

    // Initialize the stiffness matrix of every tetrahedrons
    p_stiffness_matrices.resize(topology->getNbTetrahedra());

    // Compute the initial tangent stiffness matrix
    compute_K();
}

void TetrahedronElasticForce::addForce(
        const MechanicalParams* mparams,
        Data<VecDeriv>& d_f,
        const Data<VecCoord>& d_x,
        const Data<VecDeriv>& d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    auto topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbTetrahedra())
        return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> x0 =  state->readRestPositions();
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    const std::vector<Rotation> & initial_rotation = p_initial_rotation;
    std::vector<Rotation> & current_rotation = p_current_rotation;

    bool corotated = d_corotated.getValue();

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addForce");
    const auto number_of_elements = topology->getNbTetrahedra();
    for (std::size_t element_id = 0; element_id < number_of_elements; ++element_id) {
        Tetrahedron tetra = tetrahedron(element_id, x);

        const Rotation &R0 = initial_rotation[element_id];
        const Rotation R0t = R0.transpose();

        Rotation &R = current_rotation[element_id];


        // Extract the tetrahedron's frame
        if (corotated)
            R = tetra.frame();

        const Rotation Rt = R.transpose();

        // Gather the displacement vector
        Vector<12> U;
        size_t i = 0;
        for (const auto &node_id : topology->getTetrahedron(static_cast<sofa::Index>(element_id))) {
            const Vec3 r0{x0[node_id][0], x0[node_id][1], x0[node_id][2]};
            const Vec3 r  {x[node_id][0],  x[node_id][1],  x[node_id][2]};

            const Vec3 u = Rt * r - R0t * r0;

            U[i++] = u[0];
            U[i++] = u[1];
            U[i++] = u[2];
        }

        // Compute the force vector
        const auto &K = p_stiffness_matrices[element_id];
        Vector<NumberOfNodes * 3> F = K.template selfadjointView<Eigen::Upper>() * U;

        // Write the forces into the output vector
        i = 0;
        for (const auto &node_id : topology->getTetrahedron(static_cast<sofa::Index>(element_id))) {
            Vec3 force{F[i * 3 + 0], F[i * 3 + 1], F[i * 3 + 2]};
            force = (R * force).eval();

            f[node_id][0] -= force[0];
            f[node_id][1] -= force[1];
            f[node_id][2] -= force[2];
            ++i;
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addForce");
}

void TetrahedronElasticForce::addDForce(
        const MechanicalParams* mparams,
        Data<VecDeriv>& d_df,
        const Data<VecDeriv>& d_dx)
{
    auto * topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbTetrahedra())
        return;

    auto kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addDForce");
    const auto number_of_elements = topology->getNbTetrahedra();
    for (std::size_t element_id = 0; element_id < number_of_elements; ++element_id) {

        const Rotation & R  = current_rotation[element_id];
        const Rotation   Rt = R.transpose();

        // Gather the displacement vector
        Vector<NumberOfNodes*3> U;
        size_t i = 0;
        for (const auto & node_id : topology->getTetrahedron(static_cast<sofa::Index>(element_id))) {
            const Vec3 v = {dx[node_id][0], dx[node_id][1], dx[node_id][2]};
            const Vec3 u = Rt*v;

            U[i++] = u[0];
            U[i++] = u[1];
            U[i++] = u[2];
        }

        // Compute the force vector
        const auto & K = p_stiffness_matrices[element_id];
        Vector<NumberOfNodes*3> F = K.template selfadjointView<Eigen::Upper>()*U*kFactor;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : topology->getTetrahedron(static_cast<sofa::Index>(element_id))) {
            Vec3 force {F[i*3+0], F[i*3+1], F[i*3+2]};
            force = R*force;

            df[node_id][0] -= force[0];
            df[node_id][1] -= force[1];
            df[node_id][2] -= force[2];

            ++i;
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addDForce");
}

void TetrahedronElasticForce::addKToMatrix(
        sofa::defaulttype::BaseMatrix * matrix,
        SReal kFact,
        unsigned int & offset)
{
    auto * topology = d_topology_container.get();

    if (!topology)
        return;

    std::vector<Rotation> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addKToMatrix");

    const auto number_of_elements = topology->getNbTetrahedra();
    for (std::size_t element_id = 0; element_id < number_of_elements; ++element_id) {
        const auto & node_indices = topology->getTetrahedron(static_cast<sofa::Index>(element_id));
        sofa::type::Mat3x3 R;
        for (unsigned int m = 0; m < 3; ++m) {
            for (unsigned int n = 0; n < 3; ++n) {
                R(m, n) = current_rotation[element_id](m, n);
            }
        }
        const auto   Rt = R.transposed();

        // Since the matrix K is block symmetric, we only kept the DxD blocks on the upper-triangle the matrix.
        // Here we need to accumulate the full matrix into Sofa's BaseMatrix.
        const auto & K = p_stiffness_matrices[element_id];

        // Blocks on the diagonal
        for (sofa::Index i = 0; i < NumberOfNodes; ++i) {
            const auto x = static_cast<Eigen::Index>(i*Dimension);
            sofa::type::Mat<Dimension, Dimension, Real> k;
            for (Eigen::Index m = 0; m < Dimension; ++m) {
                for (Eigen::Index n = 0; n < Dimension; ++n) {
                    k(static_cast<unsigned int>(m), static_cast<unsigned int>(n)) = K(x+m, x+n);
                }
            }

            k = -1. * R*k*Rt *kFact;

            matrix->add(offset+node_indices[i]*Dimension, offset+node_indices[i]*Dimension, k);
        }

        // Blocks on the upper triangle
        for (sofa::Index i = 0; i < NumberOfNodes; ++i) {
            for (sofa::Index j = i+1; j < NumberOfNodes; ++j) {
                const auto x = static_cast<Eigen::Index>(i*Dimension);
                const auto y = static_cast<Eigen::Index>(j*Dimension);

                sofa::type::Mat<Dimension, Dimension, Real> k;
                for (Eigen::Index m = 0; m < Dimension; ++m) {
                    for (Eigen::Index n = 0; n < Dimension; ++n) {
                        k(static_cast<unsigned int>(m), static_cast<unsigned int>(n)) = K(x+m, y+n);
                    }
                }

                k = -1. * R*k*Rt *kFact;

                matrix->add(offset+node_indices[i]*Dimension, offset+node_indices[j]*Dimension, k);
                matrix->add(offset+node_indices[j]*Dimension, offset+node_indices[i]*Dimension, k.transposed());
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addKToMatrix");
}

void TetrahedronElasticForce::computeBBox(const sofa::core::ExecParams*, bool onlyVisible)
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

    this->f_bbox.setValue(sofa::type::TBoundingBox<Real>(minBBox,maxBBox));
}

void TetrahedronElasticForce::draw(const sofa::core::visual::VisualParams* vparams)
{
    auto *topology = d_topology_container.get();
    if (!topology)
        return;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    vparams->drawTool()->saveLastState();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0, true);

    vparams->drawTool()->disableLighting();

    const VecCoord &x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    std::vector< sofa::type::Vec<3, double> > points[4];
    const auto number_of_elements = topology->getNbTetrahedra();
    for(Topology::TetrahedronID i = 0 ; i<number_of_elements;++i)
    {
        const auto t=topology->getTetra(i);

        const auto & a = t[0];
        const auto & b = t[1];
        const auto & c = t[2];
        const auto & d = t[3];
        Coord center = (x[a]+x[b]+x[c]+x[d])*0.125;
        Coord pa = (x[a]+center)*(Real)0.666667;
        Coord pb = (x[b]+center)*(Real)0.666667;
        Coord pc = (x[c]+center)*(Real)0.666667;
        Coord pd = (x[d]+center)*(Real)0.666667;

        points[0].push_back(pa);
        points[0].push_back(pb);
        points[0].push_back(pc);

        points[1].push_back(pb);
        points[1].push_back(pc);
        points[1].push_back(pd);

        points[2].push_back(pc);
        points[2].push_back(pd);
        points[2].push_back(pa);

        points[3].push_back(pd);
        points[3].push_back(pa);
        points[3].push_back(pb);
    }

    sofa::type::RGBAColor face_colors[4] = {
            {1.0, 0.0, 0.0, 1.0},
            {1.0, 0.0, 0.5, 1.0},
            {1.0, 1.0, 0.0, 1.0},
            {1.0, 0.5, 1.0, 1.0}
    };

    vparams->drawTool()->drawTriangles(points[0], face_colors[0]);
    vparams->drawTool()->drawTriangles(points[1], face_colors[1]);
    vparams->drawTool()->drawTriangles(points[2], face_colors[2]);
    vparams->drawTool()->drawTriangles(points[3], face_colors[3]);

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

void TetrahedronElasticForce::compute_K()
{
    auto topology = d_topology_container.get();

    if (!topology)
        return;

    if (p_stiffness_matrices.size() != topology->getNbTetrahedra())
        return;

    const Real youngModulus = d_youngModulus.getValue();
    const Real poissonRatio = d_poissonRatio.getValue();

    const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
    const Real mu = youngModulus / (2 * (1 + poissonRatio));

    Eigen::Matrix<Real, 6, 6> C;
    C <<
      l + 2*mu,    l,          l,       0,  0,  0,
        l,       l + 2*mu,     l,       0,  0,  0,
        l,         l,        l + 2*mu,  0,  0,  0,
        0,         0,          0,      mu,  0,  0,
        0,         0,          0,       0, mu,  0,
        0,         0,          0,       0,  0, mu;

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::compute_k");

    const auto number_of_elements = topology->getNbTetrahedra();
    for (std::size_t element_id = 0; element_id < number_of_elements; ++element_id) {
        auto & K = p_stiffness_matrices[element_id];
        K.fill(0.);

        const auto & gauss_node = p_quadrature_nodes[element_id];
        const auto detJ = gauss_node.jacobian_determinant;
        const auto dN_dx = gauss_node.dN_dx;
        const auto w = gauss_node.weight;

        // Computation of the element tangent-stiffness matrix
        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
            const Vec3 dxi = dN_dx.row(i).transpose();

            Matrix<6,3> Bi;
            Bi <<
                dxi[0],    0  ,    0  ,
                   0  , dxi[1],    0  ,
                   0  ,    0  , dxi[2],
                dxi[1], dxi[0],    0  ,
                   0  , dxi[2], dxi[1],
                dxi[2],    0  , dxi[0];
            for (std::size_t j = i; j < NumberOfNodes; ++j) {
                // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                const Vec3 dxj = dN_dx.row(j).transpose();
                Matrix<6,3> Bj;
                Bj <<
                   dxj[0],    0  ,    0  ,
                      0  , dxj[1],    0  ,
                      0  ,    0  , dxj[2],
                   dxj[1], dxj[0],    0  ,
                      0  , dxj[2], dxj[1],
                   dxj[2],    0  , dxj[0];

                K.template block<3, 3>(i*3, j*3).noalias() += (Bi.transpose()*C*Bj) * detJ * w;
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::compute_k");
}

static int TetrahedronElasticForceClass = RegisterObject("Caribou tetrahedron FEM Forcefield")
    .add< TetrahedronElasticForce >(true)
;

} // namespace SofaCaribou::forcefield
