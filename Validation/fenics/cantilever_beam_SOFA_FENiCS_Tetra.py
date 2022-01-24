# Required import for python
import Sofa
import SofaCaribou
import numpy as np


# Choose in your script to activate or not the GUI
USE_GUI = True


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    root = Sofa.Core.Node("root")
    createScene(root)
    Sofa.Simulation.init(root)

    if not USE_GUI:
        for iteration in range(10):
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


def createScene(root):
    root.addObject('DefaultVisualManagerLoop')
    root.addObject('DefaultAnimationLoop')
    root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    root.addObject('RequiredPlugin', pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")


    import meshio
    beam_p1 = meshio.read("../meshes/beam_p1.vtu")

    sofa_node = root.addChild("tetra_node_SOFA")
    sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
    sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sofa_node.addObject('MechanicalObject', name="mo", position=beam_p1.points.tolist())
    sofa_node.addObject('CaribouTopology', name='topology', template='Tetrahedron', indices=beam_p1.cells_dict['tetra'].tolist())
    sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    sofa_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
    sofa_node.addObject('ConstantForceField', force="0 -1000 0", indices="@top_roi.indices")
    sofa_node.addObject('SaintVenantKirchhoffMaterial', young_modulus="3000", poisson_ratio="0.3")
    sofa_node.addObject('HyperelasticForcefield', printLog=True)


    fenics_node = root.addChild("tetra_node_FEniCS")
    fenics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
    fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    fenics_node.addObject('MechanicalObject', name="mo", position=beam_p1.points.tolist())
    fenics_node.addObject('CaribouTopology', name='topology', template='Tetrahedron', indices=beam_p1.cells_dict['tetra'].tolist())
    fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    fenics_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
    fenics_node.addObject('ConstantForceField', force="0 -1000 0", indices="@top_roi.indices")
    fenics_node.addObject('SaintVenantKirchhoffMaterial_FEniCS', young_modulus="3000", poisson_ratio="0.3")
    fenics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)

    # import meshio
    # beam_p1 = meshio.read("./meshes/beam_p2.vtu")

    # sofa_node = root.addChild("tetra_node_SOFA")
    # sofa_node.addObject('StaticSolver', newton_iterations="2", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
    # sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    # sofa_node.addObject('MechanicalObject', name="mo", position=beam_p1.points.tolist())
    # sofa_node.addObject('CaribouTopology', name='topology', template='Tetrahedron10', indices=beam_p1.cells_dict['tetra10'].tolist())
    # sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    # sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    # sofa_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
    # sofa_node.addObject('ConstantForceField', force="0 -10 0", indices="@top_roi.indices")
    # sofa_node.addObject('SaintVenantKirchhoffMaterial', young_modulus="3000", poisson_ratio="0.3")
    # sofa_node.addObject('HyperelasticForcefield', printLog=True)


    # fenics_node = root.addChild("tetra_node_FEniCS")
    # fenics_node.addObject('StaticSolver', newton_iterations="2", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
    # fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    # fenics_node.addObject('MechanicalObject', name="mo", position=beam_p1.points.tolist())
    # fenics_node.addObject('CaribouTopology', name='topology', template='Tetrahedron10', indices=beam_p1.cells_dict['tetra10'].tolist())
    # fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    # fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    # fenics_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
    # fenics_node.addObject('ConstantForceField', force="0 -10 0", indices="@top_roi.indices")
    # fenics_node.addObject('SaintVenantKirchhoffMaterial_FEniCS', young_modulus="3000", poisson_ratio="0.3")
    # fenics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)



    return root


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()