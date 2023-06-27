#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/testing/BaseTest.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/behavior/DefaultMultiMatrixAccessor.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/simulation/graph/SimpleApi.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/component/mass/MeshMatrixMass.inl>
DISABLE_ALL_WARNINGS_END

#include <SofaCaribou/Mass/CaribouMass[Tetrahedron].h>
#include <SofaCaribou/Mass/CaribouMass[Hexahedron].h>
#include <SofaCaribou/Algebra/EigenMatrix.h>

using sofa::helper::system::PluginManager ;
using namespace sofa::simulation;
using namespace sofa::simpleapi;
using namespace sofa::helper::logging;
using namespace SofaCaribou::mass;
using namespace caribou::geometry;
using namespace caribou;
using namespace sofa::testing;

TEST(CaribouMass, LinearTetrahedron) {
    using namespace sofa::core::objectmodel;
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "RequiredPlugin", {{"pluginName", "Sofa.Component.Engine.Select"}});
    createObject(root, "RequiredPlugin", {{"pluginName", "Sofa.Component.Topology.Mapping"}});

    // Some component to avoid warnings
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");

    createObject(root, "RegularGridTopology", {{"name", "grid"}, {"min", "-7.5 -7.5 0"}, {"max", "7.5 7.5 80"}, {"n", "3 3 9"}});

    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
        createObject(root, "MechanicalObject", {{"name", "mo"}, {"src", "@grid"}}).get()
    );
    createObject(root, "TetrahedronSetTopologyContainer", {{"name", "topology"}});
    createObject(root, "TetrahedronSetTopologyModifier");
    createObject(root, "TetrahedronSetGeometryAlgorithms");
    createObject(root, "Hexa2TetraTopologicalMapping", {{"input", "@grid"}, {"output", "@topology"}});
    auto caribou_mass = dynamic_cast<CaribouMass<Tetrahedron> *> (
        createObject(root, "CaribouMass", {{"name", "caribou_mass"}, {"topology", "@topology"}, {"density", "2"}}).get()
    );
    auto sofa_mass = dynamic_cast<sofa::component::mass::MeshMatrixMass<sofa::defaulttype::Vec3Types> *> (
            createObject(root, "MeshMatrixMass", {{"name", "sofa_mass"}, {"topology", "@topology"}, {"massDensity", "2"}, {"lumping", "false"}}).get()
    );

    auto sofa_mass_diagonal = dynamic_cast<sofa::component::mass::MeshMatrixMass<sofa::defaulttype::Vec3Types> *> (
            createObject(root, "MeshMatrixMass", {{"name", "sofa_mass_diagonal"}, {"topology", "@topology"}, {"massDensity", "2"}, {"lumping", "true"}}).get()
    );
    getSimulation()->init(root.get());

    // Get M from caribou
    const Eigen::SparseMatrix<double> M = caribou_mass->M();
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & M_diag = caribou_mass->M_diag();

    // Get M from caribou using SOFA API
    SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double>> M2;
    sofa::core::behavior::DefaultMultiMatrixAccessor accessor;
    M2.resize((signed) mo->getSize()*3, (signed) mo->getSize()*3);
    accessor.setGlobalMatrix(&M2);
    accessor.addMechanicalState(mo);
    accessor.setupMatrices();
    sofa::core::MechanicalParams mechanical_parameters;
    mechanical_parameters.setMFactor(1);
    caribou_mass->addMToMatrix(&mechanical_parameters, &accessor);
    M2.compress();

    EXPECT_DOUBLE_EQ(M.sum(), M2.matrix().sum());

    // Get lumped M from caribou using SOFA API
    caribou_mass->findData("lumped")->read("true");
    SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double>> M2_diag;
    M2_diag.resize((signed) mo->getSize()*3, (signed) mo->getSize()*3);
    accessor.setGlobalMatrix(&M2_diag);
    caribou_mass->addMToMatrix(&mechanical_parameters, &accessor);
    M2_diag.compress();

    EXPECT_DOUBLE_EQ(M.sum(), M2_diag.matrix().sum());

    // Get M from SOFA
    SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double>> SofaM;
    SofaM.resize((signed) mo->getSize()*3, (signed) mo->getSize()*3);
    accessor.setGlobalMatrix(&SofaM);
    sofa_mass->addMToMatrix(&mechanical_parameters, &accessor);
    SofaM.compress();

    EXPECT_DOUBLE_EQ(M.sum(), SofaM.matrix().sum());

    // Get M diagonal from SOFA
    SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double>> SofaM_diagonal;
    SofaM_diagonal.resize((signed) mo->getSize()*3, (signed) mo->getSize()*3);
    accessor.setGlobalMatrix(&SofaM_diagonal);
    sofa_mass_diagonal->addMToMatrix(&mechanical_parameters, &accessor);
    SofaM_diagonal.compress();

    EXPECT_DOUBLE_EQ(M_diag.diagonal().sum(), SofaM_diagonal.matrix().sum());

    // AddForce
    using Real = CaribouMass<Tetrahedron>::Real;
    using VecDeriv = CaribouMass<Tetrahedron>::VecDeriv;
    using DataVecDeriv = CaribouMass<Tetrahedron>::DataVecDeriv;
    DataVecDeriv d_f_caribou (VecDeriv(mo->getSize(), {0, 0, 0}));
    DataVecDeriv d_f_caribou_dia (VecDeriv(mo->getSize(), {0, 0, 0}));
    DataVecDeriv d_f_sofa (VecDeriv(mo->getSize(), {0, 0, 0}));
    DataVecDeriv d_f_sofa_dia (VecDeriv(mo->getSize(), {0, 0, 0}));

    caribou_mass->findData("lumped")->read("false");
    caribou_mass->addForce(&mechanical_parameters, d_f_caribou, mo->x, mo->v);

    caribou_mass->findData("lumped")->read("true");
    caribou_mass->addForce(&mechanical_parameters, d_f_caribou_dia, mo->x, mo->v);

    sofa_mass->addForce(&mechanical_parameters, d_f_sofa, mo->x, mo->v);
    sofa_mass_diagonal->addForce(&mechanical_parameters, d_f_sofa_dia, mo->x, mo->v);

    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor>> f_caribou ((d_f_caribou.getValue().data()->data()),  mo->getSize(), 3);
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor>> f_caribou_dia ((d_f_caribou_dia.getValue().data()->data()),  mo->getSize(), 3);
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor>> f_sofa ((d_f_sofa.getValue().data()->data()),  mo->getSize(), 3);
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor>> f_sofa_dia ((d_f_sofa_dia.getValue().data()->data()),  mo->getSize(), 3);


    EXPECT_DOUBLE_EQ(f_caribou.norm(), f_caribou_dia.norm());
    EXPECT_DOUBLE_EQ(f_caribou_dia.norm(), f_sofa_dia.norm());
}

TEST(CaribouMass, LinearHexahedron) {
    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "RequiredPlugin", {{"pluginName", "Sofa.Component.Engine.Select"}});
    createObject(root, "RequiredPlugin", {{"pluginName", "Sofa.Component.Topology.Mapping"}});

    // Some component to avoid warnings
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");

    createObject(root, "RegularGridTopology", {{"name", "grid"}, {"min", "-7.5 -7.5 0"}, {"max", "7.5 7.5 80"}, {"n", "3 3 9"}});

    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"src", "@grid"}}).get()
    );
    createObject(root, "HexahedronSetTopologyContainer", {{"name", "topology"}, {"src", "@grid"}});
    createObject(root, "HexahedronSetGeometryAlgorithms");
    auto caribou_mass = dynamic_cast<const CaribouMass<Hexahedron> *> (
            createObject(root, "CaribouMass", {{"name", "caribou_mass"}, {"topology", "@topology"}, {"density", "2"}}).get()
    );

    auto sofa_mass = dynamic_cast<sofa::component::mass::MeshMatrixMass<sofa::defaulttype::Vec3Types> *> (
            createObject(root, "MeshMatrixMass", {{"name", "sofa_mass"}, {"topology", "@topology"}, {"massDensity", "2"}}).get()
    );

    getSimulation()->init(root.get());

    // Get M from caribou
    const Eigen::SparseMatrix<double> M = caribou_mass->M();

    // Get M from SOFA
    SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double>> SofaM;
    sofa::core::behavior::DefaultMultiMatrixAccessor accessor;
    SofaM.resize((signed) mo->getSize()*3, (signed) mo->getSize()*3);
    accessor.setGlobalMatrix(&SofaM);
    accessor.addMechanicalState(mo);
    accessor.setupMatrices();
    sofa::core::MechanicalParams mechanical_parameters;
    mechanical_parameters.setMFactor(1);
    sofa_mass->addMToMatrix(&mechanical_parameters, &accessor);
    SofaM.compress();

    EXPECT_DOUBLE_EQ(M.sum(), SofaM.matrix().sum());
}
