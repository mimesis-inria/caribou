#include "IBMForcefield.h"
#include "../Helper/MirtichIntegration.h"
#include "../Helper/Hexahedron.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace sofa
{

using namespace core::objectmodel;
using namespace core::behavior;
using namespace core::topology;
using namespace component::topology;
using namespace defaulttype;

namespace caribou
{

namespace forcefield
{

template<class DataTypes>
IBMForcefield<DataTypes>::IBMForcefield()
        : d_youngModulus(initData(&d_youngModulus, Real(1000), "youngModulus", "[IN] Young's modulus of the material"))
        , d_poissonRatio(initData(&d_poissonRatio, Real(0.3),  "poissonRatio", "[IN] Poisson's ratio of the material"))
        , d_hexahedrons(initData(&d_hexahedrons, "hexahedrons",
                "[IN] List of hexahedrons [n1H1, n2H1, n3H1, ..., n8Hm] where niHj is the node id of vertex i (0 to 7) of the hexahedron j (0 to m)."))
        , d_hexahedrons_flags(initData(&d_hexahedrons_flags, "hexahedrons_flags", "[IN] Flags (Outside -1, Boundary 0 or Inside 1) of the hexahedrons"))
        , d_triangle_positions(initData(&d_triangle_positions, "triangle_positions", "[IN] List of triangles vertices positions."))
        , d_triangles(initData(&d_triangles, "triangles",
                "[IN] List of triangles per hexahedron (they should represent both the cut surface and the hexa's subsurfaces under the cut).", false))
{
            f_listening.setValue(true);
}

template<class DataTypes>
void IBMForcefield<DataTypes>::init()
{
    Inherit::init();

    // If no link to a vector of initial positions is set and no initial positions were manually added to the input, try
    // to find a mechanical state that contain those positions
    if (!d_initial_positions.getParent() && d_initial_positions.getValue().empty()) {
        auto state = this->getMState();
        if (!state)
            state = this->getContext()->template get<MechanicalState<DataTypes>>();
        if (state) {
            d_initial_positions.setParent(state->findData("rest_position"));
            msg_info() << "Initial positions automatically linked to the mechanical object '"
                       << state->getPathName() << ".rest_position'";
            if (d_initial_positions.getValue().empty()) {
                msg_warning()
                        << "Initial positions were automatically linked to the mechanical object '"
                        << state->getPathName() << "' but it contains an empty set of positions.";
            }
        } else {
            msg_error() << "No initial positions were provided and no mechanical state was found in the current context.";
        }
    }

    // If no link to an hexahedron list is set and no hexahedron were manually added to the input, try to find a
    // list of hexahedron in the current context.
    if (!d_hexahedrons.getParent() && d_hexahedrons.getValue().empty()) {
        auto topology = this->getContext()->template get<BaseMeshTopology>();
        if (topology && topology->findData("hexahedra")) {
            d_hexahedrons.setParent(topology->findData("hexahedra"));
            msg_info() << "Hexahedron list automatically linked to the mesh topology '"
                << topology->getPathName() << "'";
            if (topology->getNbHexahedra() == 0) {
                msg_warning()
                    << "Hexahedron list was automatically linked to the mesh topology '"
                    << topology->getPathName() << "' but it contains an empty set of hexahedrons.";
            }
        } else {
            msg_error() << "No hexahedrons list was provided and no mesh topology was found in the current context.";
        }
    }

    // If no hexahedron flags set, try to find  CutGridEngine in the current context
    if (!d_hexahedrons_flags.getParent() && d_hexahedrons_flags.getValue().empty()) {
        auto engine = this->getContext()->template get<engine::CutGridEngine>();
        if (engine && engine->findData("hexahedrons_flags")) {
            d_hexahedrons_flags.setParent(engine->findData("hexahedrons_flags"));
            msg_info() << "Hexahedron flag list automatically linked to the cut grid engine '"
                << engine->getPathName() << "'";
        }
    }

    // If no triangles set, try to find a CutGridEngine in the current context
    if (!d_triangle_positions.getParent() || !d_triangles.getParent()) {
        auto engine = this->getContext()->template get<engine::CutGridEngine>();
        if (engine && engine->findData("triangle_positions") && engine->findData("triangles")) {
            d_triangle_positions.setParent(engine->findData("triangle_positions"));
            d_triangles.setParent(engine->findData("triangles"));
            msg_info() << "Hexahedron triangles list automatically linked to the cut grid engine '"
                       << engine->getPathName() << "'";
        }
    }

    reinit();
}

template<class DataTypes>
void IBMForcefield<DataTypes>::reinit()
{
    Inherit::reinit();
    update();
}

template<class DataTypes>
void IBMForcefield<DataTypes>::update()
{

    if (   !d_initial_positions.isDirty()
        && !d_hexahedrons_flags.isDirty()
        && !d_hexahedrons.isDirty()
        && !d_triangle_positions.isDirty()
        && !d_triangles.isDirty()
        && !d_poissonRatio.isDirty()
        && !d_youngModulus.isDirty()) {
        cleanDirty();
        return;
    }

    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Hexahedron>>> hexahedrons = d_hexahedrons;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Coord>>> initial_positions = d_initial_positions;

    cleanDirty();

    m_hexahedrons.clear();

    // Lame's coefficients and the material matrix
    const Real & youngModulus = this->d_youngModulus.getValue();
    const Real & poissonRatio = this->d_poissonRatio.getValue();

    const Real l = youngModulus*poissonRatio / ((1 + poissonRatio)*(1 - 2*poissonRatio));
    const Real m = youngModulus / (2 * (1 + poissonRatio));
    const Real a = l + 2*m;
    Mat66 C;
    C[0][0] = a; C[0][1] = l; C[0][2] = l; C[0][3] = 0; C[0][4] = 0; C[0][5] = 0;
    C[1][0] = l; C[1][1] = a; C[1][2] = l; C[1][3] = 0; C[1][4] = 0; C[1][5] = 0;
    C[2][0] = l; C[2][1] = l; C[2][2] = a; C[2][3] = 0; C[2][4] = 0; C[2][5] = 0;
    C[3][0] = 0; C[3][1] = 0; C[3][2] = 0; C[3][3] = m; C[3][4] = 0; C[3][5] = 0;
    C[4][0] = 0; C[4][1] = 0; C[4][2] = 0; C[4][3] = 0; C[4][4] = m; C[4][5] = 0;
    C[5][0] = 0; C[5][1] = 0; C[5][2] = 0; C[5][3] = 0; C[5][4] = 0; C[5][5] = m;

    // Make sure that the hexahedrons are rectangular and that their nodes are contained in the initial positions vector
    for (const Hexahedron & nodes : hexahedrons) {

        // Step 1 (/2) : Check that all nodes are contain in the initial positions vector
        for (const PointID & node : nodes) {
            if (node > initial_positions.size()) {
                msg_error() << "Some of the hexahedron's nodes don't match the size of the initial positions vector.";
                m_hexahedrons.clear();
                return;
            }
        }

        // Step 2 (/2) : Check that all hexahedrons are rectangular
        Coord lx = initial_positions[nodes[1]]-initial_positions[nodes[0]];
        Coord ly = initial_positions[nodes[3]]-initial_positions[nodes[0]];
        Coord lz = initial_positions[nodes[4]]-initial_positions[nodes[0]];

        if ( not (
                (initial_positions[nodes[3]]+lx-initial_positions[nodes[2]]).norm() < lx.norm()*0.001 &&
                (initial_positions[nodes[0]]+lz-initial_positions[nodes[4]]).norm() < lz.norm()*0.001 &&
                (initial_positions[nodes[1]]+lz-initial_positions[nodes[5]]).norm() < lz.norm()*0.001 &&
                (initial_positions[nodes[2]]+lz-initial_positions[nodes[6]]).norm() < lz.norm()*0.001 &&
                (initial_positions[nodes[3]]+lz-initial_positions[nodes[7]]).norm() < lz.norm()*0.001
                ))
        {
            msg_error() << "Some of the hexahedrons are not rectangular.";
            m_hexahedrons.clear();
            return;
        }
    }

    // Compute the element stiffness matrices
    m_hexahedrons.resize(hexahedrons.size());
    for (size_t hexa_id = 0; hexa_id < hexahedrons.size(); ++hexa_id) {
        Hexa & hexa = m_hexahedrons[hexa_id];

        std::array<Coord, 8> nodes;

        for (size_t node_id = 0; node_id < hexahedrons[hexa_id].size(); ++node_id) {
            hexa.nodes[node_id] = hexahedrons[hexa_id][node_id];
            nodes[node_id] = initial_positions[hexa.nodes[node_id]];
        }

        if (hexa_id == 0) {
            const Real hx = (nodes[1]-nodes[0]).norm();
            const Real hy = (nodes[3]-nodes[0]).norm();
            const Real hz = (nodes[4]-nodes[0]).norm();
            std::cout << "hx, hy, hz = " << hx << ", " << hy << ", " << hz <<std::endl;
        }

        // Jacobian matrix computation at local position (xi, eta, zeta) :
        // If x(xi, eta, zeta) = sum_i N_i (xi, eta, zeta) x_i
        // Then
        //      | dx/dxi    dy/dxi   dz/dxi   |
        //  J = | dx/deta   dy/deta  dz/deta  |
        //      | dx/dzeta  dy/dzeta dz/dzeta |
        //
        // Where dx/dXi = sum_i [dNi/dXi]*x_i
        auto Jacobian = [&nodes] (Real xi, Real eta, Real zeta) -> Mat33 {
            // local position aren't used as the derivatives of the shape function are constants
            // Normally, you would loop on each nodes to compute the jacobian as J = sum [dNi/dXi]*x_i
            SOFA_UNUSED(xi);
            SOFA_UNUSED(eta);
            SOFA_UNUSED(zeta);

            const Real hx = (nodes[1]-nodes[0]).norm();
            const Real hy = (nodes[3]-nodes[0]).norm();
            const Real hz = (nodes[4]-nodes[0]).norm();

            return
            {
                {   hx/2.,     0  ,    0    },
                {     0  ,   hy/2.,    0    },
                {     0  ,     0  ,   hz/2. }
            };
        };

        // Integrand used to compute the stiffness
        auto f = [this, &Jacobian, &C] (Real xi, Real eta, Real zeta) -> Mat24_24 {
            Mat6_24 B;
            const Coord x(xi, eta, zeta);
            Mat33 J = Jacobian(xi, eta, zeta);
            Real detJ = defaulttype::determinant(J);
            if (detJ < 1.0e-100) {
                msg_error() << "Some of the hexahedrons are badly shaped and their jacobian cannot be inverted.";
                m_hexahedrons.clear();
                return Mat24_24();
            }

            Mat33 Jinverted = J.inverted();

            // Building the strain-displacement matrix B(xi, eta, zeta)
            for (unsigned int i = 0; i < 8; ++i) {
                const auto & Xi = sofa::caribou::helper::hexahedron::local_coordinate_of_node[i];
                const Deriv dN_dXi = 1/8. * Deriv(
                              Xi[0]     * (1. + Xi[1]*eta) * (1. + Xi[2]*zeta),
                        (1. + Xi[0]*xi) *       Xi[1]      * (1. + Xi[2]*zeta),
                        (1. + Xi[0]*xi) * (1. + Xi[1]*eta) *       Xi[2]
                );
                const Deriv dN_dx = Jinverted * dN_dXi; // Derivative of the shape function at local position (xi, eta, zeta)
                B[0][i*3+0] = dN_dx[0];  B[0][i*3+1] = 0;         B[0][i*3+2] = 0;
                B[1][i*3+0] = 0;         B[1][i*3+1] = dN_dx[1];  B[1][i*3+2] = 0;
                B[2][i*3+0] = 0;         B[2][i*3+1] = 0;         B[2][i*3+2] = dN_dx[2];
                B[3][i*3+0] = dN_dx[1];  B[3][i*3+1] = dN_dx[0];  B[3][i*3+2] = 0;
                B[4][i*3+0] = 0;         B[4][i*3+1] = dN_dx[2];  B[4][i*3+2] = dN_dx[1];
                B[5][i*3+0] = dN_dx[2];  B[5][i*3+1] = 0;         B[5][i*3+2] = dN_dx[0];
            }

            const Mat24_6 Bt = B.transposed();

            return Bt * C * B * detJ;
        };

        // Gauss quadrature
        Mat24_24 & K = hexa.K;
        K.fill(0.);

        const std::array<Coord, 8> integration_points = {{
                {-1./sqrt(3.0), -1./sqrt(3.0), -1./sqrt(3.0)},
                {-1./sqrt(3.0), -1./sqrt(3.0),  1./sqrt(3.0)},
                {-1./sqrt(3.0),  1./sqrt(3.0), -1./sqrt(3.0)},
                {-1./sqrt(3.0),  1./sqrt(3.0),  1./sqrt(3.0)},
                { 1./sqrt(3.0), -1./sqrt(3.0), -1./sqrt(3.0)},
                { 1./sqrt(3.0), -1./sqrt(3.0),  1./sqrt(3.0)},
                { 1./sqrt(3.0),  1./sqrt(3.0), -1./sqrt(3.0)},
                { 1./sqrt(3.0),  1./sqrt(3.0),  1./sqrt(3.0)}
        }};

        for (const auto & Xi : integration_points) {
            const auto & xi   = Xi[0];
            const auto & eta  = Xi[1];
            const auto & zeta = Xi[2];
            K += f(xi, eta, zeta);
        }

        msg_info_when(hexa_id == 0) << K[0][0] << " " << K[0][1] << " " << K[0][2];
    }
}

template<class DataTypes>
void IBMForcefield<DataTypes>::reset()
{
    Inherit::reset();
}

template<class DataTypes>
void IBMForcefield<DataTypes>::addForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);
    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> x0 = d_initial_positions;
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    for (const auto & hexa : m_hexahedrons) {
        // Gather the displacement vector
        Vec24 U;
        size_t i = 0;
        for (const auto & node_id : hexa.nodes) {
            U[i++] = x[node_id][0] - x0[node_id][0];
            U[i++] = x[node_id][1] - x0[node_id][1];
            U[i++] = x[node_id][2] - x0[node_id][2];
        }

        // Compute the force vector
        const auto & K = hexa.K;
        Vec24 F = K*U;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : hexa.nodes) {
            f[node_id][0] -= F[i++];
            f[node_id][1] -= F[i++];
            f[node_id][2] -= F[i++];
        }
    }
}

template<class DataTypes>
void IBMForcefield<DataTypes>::addDForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_df, const Data<VecDeriv>& d_dx)
{
    auto kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;

    for (const auto & hexa : m_hexahedrons) {
        // Gather the displacement vector
        Vec24 U;
        size_t i = 0;
        for (const auto & node_id : hexa.nodes) {
            U[i++] = dx[node_id][0];
            U[i++] = dx[node_id][1];
            U[i++] = dx[node_id][2];
        }

        // Compute the force vector
        const auto & K = hexa.K;
        Vec24 F = K*U*kFactor;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : hexa.nodes) {
            df[node_id][0] -= F[i++];
            df[node_id][1] -= F[i++];
            df[node_id][2] -= F[i++];
        }
    }
}

template<class DataTypes>
void IBMForcefield<DataTypes>::addKToMatrix(BaseMatrix * matrix, SReal kFact, unsigned int &offset)
{
    SOFA_UNUSED(offset);
    SOFA_UNUSED(matrix);
    SOFA_UNUSED(kFact);
}

template<class DataTypes>
void IBMForcefield<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    SOFA_UNUSED(vparams);
}

SOFA_DECL_CLASS(IBMForcefield)
static int IBMForcefieldClass = core::RegisterObject("Caribou IBM Forcefield")
        .add< IBMForcefield<sofa::defaulttype::Vec3dTypes> >(true)
;

} // namespace forcefield

} // namespace caribou

} // namespace sofa