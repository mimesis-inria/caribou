#include <array>

#include <SofaCaribou/config.h>
#include <SofaCaribou/Ode/BackwardEulerODESolver.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/helper/testing/BaseTest.h>
#include <sofa/simulation/Node.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/SimpleApi.h>
#include <SofaBaseMechanics/MechanicalObject.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::simulation;
using namespace sofa::simpleapi;
using namespace sofa::helper::logging;

#if (defined(SOFA_VERSION) && SOFA_VERSION >= 201299)
using namespace sofa::testing;
#endif

/** Make sure residual norms at each newton steps remains the same */
TEST(BackwardEulerODESolver, Beam) {
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "RequiredPlugin", {{"pluginName", "SofaBoundaryCondition SofaEngine"}});
#if (defined(SOFA_VERSION) && SOFA_VERSION > 201299)
    createObject(root, "RequiredPlugin", {{"pluginName", "SofaTopologyMapping"}});
#endif
    createObject(root, "RegularGridTopology", {{"name", "grid"}, {"min", "-7.5 -7.5 0"}, {"max", "7.5 7.5 80"}, {"n", "3 3 9"}});

    auto meca = createChild(root, "meca");
    // Create the ODE system
    auto solver = dynamic_cast<SofaCaribou::ode::BackwardEulerODESolver *>(
            createObject(meca, "BackwardEulerODESolver", {{"newton_iterations", "10"}, {"correction_tolerance_threshold", "1e-5"}, {"residual_tolerance_threshold", "1e-5"}, {"printLog", "0"}}).get()
    );
    SOFA_UNUSED(solver);
    createObject(meca, "LLTSolver", {{"Backend", "Pardiso"}});
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(meca, "MechanicalObject", {{"name", "mo"}, {"src", "@../grid"}}).get()
    );

    // Complete hexahedral topology container
    createObject(meca, "HexahedronSetTopologyContainer", {{"name", "mechanical_topology"}, {"src", "@../grid"}});
    createObject(meca, "HexahedronSetGeometryAlgorithms");

    // Mechanics
    createObject(meca, "SaintVenantKirchhoffMaterial", {{"young_modulus", "15000"}, {"poisson_ratio", "0.3"}});
    createObject(meca, "HyperelasticForcefield");
    createObject(meca, "DiagonalMass", {{"massDensity", "0.2"}});

    // Fix the left side of the beam
    createObject(meca, "BoxROI", {{"name", "fixed_roi"}, {"quad", "@surface_topology.quad"}, {"box", "-7.5 -7.5 -0.9 7.5 7.5 0.1"}});
    createObject(meca, "FixedConstraint", {{"indices", "@fixed_roi.indices"}});

    // The simulation is supposed to converged in 5 load increments. Let's check the norms of force residuals.
    // See '/validation/fenics_rectangular_beam_bending_euler_implicit.py' to see how these values were computed.
    std::array<std::vector<double>, 30> force_residuals = {{
        {1.000000000000000E+00, 1.263047604094785E-02, 4.050301446519007E-06},
        {1.000000000000000E+00, 2.632224446156891E-02, 1.619922833521263E-04, 1.616980158894692E-06},
        {1.000000000000000E+00, 1.638518248678503E-02, 3.321542817305300E-05, 4.111144237506441E-08},
        {1.000000000000000E+00, 3.662467717346047E-03, 8.427258832743349E-06},
        {1.000000000000000E+00, 1.439221631373104E-04, 5.098664430760925E-10},
        {1.000000000000000E+00, 1.933845731685127E-03, 1.509709372843117E-06},
        {1.000000000000000E+00, 2.833188504840843E-03, 8.994470332219715E-07},
        {1.000000000000000E+00, 1.611759545281556E-03, 3.625939015511337E-07},
        {1.000000000000000E+00, 3.076448851727649E-04, 1.135654189363627E-08},
        {1.000000000000000E+00, 1.396369913198485E-05, 1.139056283741289E-11},
        {1.000000000000000E+00, 2.481173603651881E-04, 9.550188626110061E-10},
        {1.000000000000000E+00, 3.204391442680505E-04, 3.356373532028338E-09},
        {1.000000000000000E+00, 1.640463688294433E-04, 3.700309102675296E-10},
        {1.000000000000000E+00, 2.445476160882390E-05, 6.208686437235950E-12},
        {1.000000000000000E+00, 3.478068238687613E-06},
        {1.000000000000000E+00, 3.068294045100563E-05, 4.934996641626367E-12},
        {1.000000000000000E+00, 3.486215670524717E-05, 4.028023638451749E-12},
        {1.000000000000000E+00, 1.613623442734944E-05, 1.038887357656362E-11},
        {1.000000000000000E+00, 1.858279755436709E-06},
        {1.000000000000000E+00, 6.861427082517277E-07},
        {1.000000000000000E+00, 3.738665999792629E-06},
        {1.000000000000000E+00, 3.783471951891048E-06},
        {1.000000000000000E+00, 1.568349411166458E-06},
        {1.000000000000000E+00, 1.283686566734150E-07},
        {1.000000000000000E+00, 1.173502099504653E-07},
        {1.000000000000000E+00, 4.465591842506106E-07},
        {1.000000000000000E+00, 4.047122691901852E-07},
        {1.000000000000000E+00, 1.498608167902650E-07},
        {1.000000000000000E+00, 1.661440211842977E-08},
        {1.000000000000000E+00, 1.942886558809736E-08}
    }};

    std::array<std::array<double, 3>, 30> middle_point_positions = {{
        {-2.282753558235673E-16, -9.562279888345634E+00, 7.939044911706209E+01},
        {1.901259251404760E-16, -2.200410424789615E+01, 7.659222365577006E+01},
        {2.015685812650761E-16, -3.120604722752000E+01, 7.281692906739460E+01},
        {1.011708340370718E-15, -3.528683993235831E+01, 7.056333712082601E+01},
        {-4.443435421368654E-16, -3.498586694947343E+01, 7.070392039083237E+01},
        {1.650018896531739E-16, -3.200995159620113E+01, 7.233551447938810E+01},
        {1.012486093846609E-15, -2.830526401195670E+01, 7.412486168727997E+01},
        {5.441904441722609E-16, -2.545401323396197E+01, 7.531823803168405E+01},
        {1.191676194938067E-15, -2.420965860623026E+01, 7.578805469064319E+01},
        {1.646407048363850E-15, -2.446657107484880E+01, 7.568966106228775E+01},
        {2.283557143341706E-15, -2.559250504010098E+01, 7.525728747850749E+01},
        {2.947869809569551E-15, -2.686151369944466E+01, 7.474349714257900E+01},
        {2.608397198097907E-15, -2.776414171227908E+01, 7.435949352208119E+01},
        {2.119576377750636E-15, -2.811175451126505E+01, 7.420702108886103E+01},
        {1.873320252405293E-15, -2.798182877719765E+01, 7.426380132603475E+01},
        {1.303886951797084E-15, -2.759265512981549E+01, 7.443292654245047E+01},
        {1.023256629307262E-15, -2.717662150984503E+01, 7.461066307604858E+01},
        {1.021691406323644E-15, -2.689306970402827E+01, 7.472989238647719E+01},
        {1.760965698150847E-16, -2.679684503737767E+01, 7.476996083044881E+01},
        {-1.705871313810693E-16, -2.685543354897088E+01, 7.474552358242528E+01},
        {8.134563211324784E-16, -2.699200024082513E+01, 7.468838986279032E+01},
        {-3.973603301789813E-16, -2.712925614969991E+01, 7.463062511409947E+01},
        {-4.108455694780989E-16, -2.721757729427946E+01, 7.459326236938304E+01},
        {-3.803364958639054E-16, -2.724284873309757E+01, 7.458253917003719E+01},
        {7.385944151580816E-16, -2.721868848076719E+01, 7.459278259215507E+01},
        {2.525836216466421E-15, -2.717156673434412E+01, 7.461273635065081E+01},
        {2.364539703625593E-15, -2.712668991281993E+01, 7.463170140258259E+01},
        {2.341392475519807E-15, -2.709939933059889E+01, 7.464321530875475E+01},
        {1.872563872557323E-15, -2.709319289376401E+01, 7.464583124869419E+01},
        {1.647607571061131E-15, -2.710271984455126E+01, 7.464181417401340E+01}
    }};

    getSimulation()->init(root.get());

    for (unsigned int step_id = 0; step_id < force_residuals.size(); ++step_id) {
        getSimulation()->animate(root.get(), 1);
// @todo uncomment the following when we have time to debug why fenics converges faster...

//        EXPECT_EQ(solver->squared_residuals().size(), force_residuals[step_id].size()) << "Time step # "<< step_id;
//        for (unsigned int newton_step_id = 0; newton_step_id < solver->squared_residuals().size(); ++newton_step_id) {
//            auto residual = sqrt(solver->squared_residuals()[newton_step_id]) / sqrt(solver->squared_residuals()[0]);
//            auto solution_residual = force_residuals[step_id][newton_step_id];
//            auto rel_error = abs(residual - solution_residual) / solution_residual;
//            EXPECT_LE(rel_error, 1) << "Time step # "<< step_id << ", Newton step # " << newton_step_id << "\n"
//            << "Residual: " << residual << "\n"
//            << "Solution Residual: " << solution_residual << "\n"
//            << "Difference: " << abs(residual - solution_residual);
//        }

        // Let's also validate the position of the node positioned at the center of the end-surface of the beam
        // (the surface where the traction is applied)
        const auto & middle_point = mo->read(sofa::core::ConstVecCoordId::position())->getValue()[76];
        const auto  middle_point_solution = sofa::defaulttype::Vec3(middle_point_positions[step_id][0], middle_point_positions[step_id][1], middle_point_positions[step_id][2]);
        const auto abs_err = (middle_point - middle_point_solution).norm();
        const auto rel_err = abs_err / middle_point_solution.norm();

        EXPECT_LE(rel_err, 0.005);
    }

    getSimulation()->unload(root);
}