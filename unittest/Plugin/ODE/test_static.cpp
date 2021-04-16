#include <array>

#include <SofaCaribou/config.h>
#include <SofaCaribou/Ode/StaticODESolver.h>

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

/** Initialization without any linear solver (expecting an error) */
TEST(StaticODESolver, InitWithoutSolver) {
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_EMIT(Error) ;

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "StaticODESolver", {{"printLog", "true"}});
    getSimulation()->init(root.get());
    getSimulation()->unload(root);
}

/** Initialization without any compatible linear solver (expecting an error) */
TEST(StaticODESolver, InitSofaSolver) {
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_EMIT(Error) ;

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "StaticODESolver", {{"printLog", "true"}});
    createObject(root, "CGLinearSolver");
    createObject(root, "CGLinearSolver");
    getSimulation()->init(root.get());
    getSimulation()->unload(root);
}

/** Initialization with both a Caribou linear solver and a SOFA solver (should inform the user of the choice made) */
TEST(StaticODESolver, InitCaribouSolver) {
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Warning, Error);
    EXPECT_MSG_EMIT(Info) ;

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "StaticODESolver", {{"printLog", "true"}});
    createObject(root, "CGLinearSolver");
    createObject(root, "LDLTSolver");
    getSimulation()->init(root.get());
    getSimulation()->unload(root);
}

/** Initialization with two Caribou linear solvers (should warn the user of the choice made) */
TEST(StaticODESolver, InitMultipleCaribouSolver) {
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Info, Error);
    EXPECT_MSG_EMIT(Warning);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "StaticODESolver", {{"printLog", "true"}});
    createObject(root, "LDLTSolver", {{"name", "first_solver"}});
    createObject(root, "LDLTSolver", {{"name", "second_solver"}});
    getSimulation()->init(root.get());
    getSimulation()->unload(root);
}

/** Make sure residual norms at each newton steps remains the same */
TEST(StaticODESolver, Beam) {
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
    auto solver = dynamic_cast<SofaCaribou::ode::StaticODESolver *>(
            createObject(meca, "StaticODESolver", {{"newton_iterations", "10"}, {"correction_tolerance_threshold", "1e-5"}, {"residual_tolerance_threshold", "1e-5"}}).get()
    );
    createObject(meca, "LDLTSolver");
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(meca, "MechanicalObject", {{"name", "mo"}, {"src", "@../grid"}}).get()
    );

    // Complete hexahedral topology container
    createObject(meca, "HexahedronSetTopologyContainer", {{"name", "mechanical_topology"}, {"src", "@../grid"}});

    // Mechanics
    createObject(meca, "SaintVenantKirchhoffMaterial", {{"young_modulus", "3000"}, {"poisson_ratio", "0.499"}});
    createObject(meca, "HyperelasticForcefield");

    // Fix the left side of the beam
    createObject(meca, "BoxROI", {{"name", "fixed_roi"}, {"quad", "@surface_topology.quad"}, {"box", "-7.5 -7.5 -0.9 7.5 7.5 0.1"}});
    createObject(meca, "FixedConstraint", {{"indices", "@fixed_roi.indices"}});

    // Apply traction on the right side of the beam
    createObject(meca, "BoxROI", {{"name", "top_roi"}, {"quad", "@surface_topology.quad"}, {"box", "-7.5 -7.5 79.9 7.5 7.5 80.1"}});
    createObject(meca, "QuadSetTopologyContainer", {{"name", "traction_container"}, {"quads", "@top_roi.quadInROI"}});
    createObject(meca, "TractionForce", {{"traction", "0 -30 0"}, {"slope", "0.2"}, {"quads", "@traction_container.quads"}});

    // The simulation is supposed to converged in 5 load increments. Let's check the norms of force residuals.
    // todo(jnbrunet): These values should be validated against another external software
    std::array<std::vector<double>, 5> force_residuals = {{
            {1.000000000000000e+00, 4.802014099846802e-03, 3.457785823647809e-04, 6.417769508095073e-08},  // Step 1
            {1.000000000000000e+00, 4.700052322118300e-03, 3.771851553383748e-04, 1.289757638852382e-05, 2.125041655083045e-08}, // Step 2
            {1.000000000000000e+00, 4.416064536488316e-03, 4.703733156787435e-04, 3.659722934024732e-05, 7.791507698853802e-08}, // Step 3
            {1.000000000000000e+00, 4.003458215116851e-03, 6.263051713537671e-04, 4.996976075176631e-05, 1.422069948624354e-07}, // Step 4
            {1.000000000000000e+00, 3.526942203674829e-03, 8.307813177405512e-04, 4.667215114394798e-05, 1.646730071539153e-07}  // Step 5
    }};

    getSimulation()->init(root.get());

    for (unsigned int step_id = 0; step_id < force_residuals.size(); ++step_id) {
        getSimulation()->animate(root.get(), 1);
        EXPECT_EQ(solver->squared_residuals().size(), force_residuals[step_id].size());
        for (unsigned int newton_step_id = 0; newton_step_id < solver->squared_residuals().size(); ++newton_step_id) {
            double residual = solver->squared_residuals()[newton_step_id] / solver->squared_residuals()[0];
            EXPECT_NEAR(sqrt(residual), force_residuals[step_id][newton_step_id], 1e-10);
        }
    }

    // Let's also validate the position of the node positioned at the center of the end-surface of the beam
    // (the surface where the traction is applied)
    const auto & middle_point = mo->read(sofa::core::ConstVecCoordId::position())->getValue()[76];
    EXPECT_NEAR(middle_point[0],   0.000, 1e-3); // x
    EXPECT_NEAR(middle_point[1], -21.016, 1e-3); // y
    EXPECT_NEAR(middle_point[2],  76.190, 1e-3); // z

    getSimulation()->unload(root);
}