import numpy as np
import meshio

import Sofa
import SofaCaribou


class ControlFrame(Sofa.Core.Controller):


    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):
        
        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin', pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")

        beam_q2 = meshio.read("/home/sidaty/Desktop/external_plugins/caribou/Validation/meshes/beam_q2.vtu")

        # TODO improve the manual permutation for matching the redefinition of the hexahedron
        # TODO redefine the visualization of the hexaedron
        indices = np.empty(beam_q2.cells_dict['hexahedron20'].shape)
        indices = beam_q2.cells_dict['hexahedron20'][:, [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]


        sofa_node = root.addChild("tetra_node_SOFA")
        sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
        sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.sofa_mo = sofa_node.addObject('MechanicalObject', name="mo", position=beam_q2.points.tolist())
        sofa_node.addObject('CaribouTopology', name='topology', template='Hexahedron20', indices=indices.tolist())
        sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
        sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
        sofa_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
        sofa_node.addObject('ConstantForceField', force="0 -100 0", indices="@top_roi.indices")
        sofa_node.addObject('SaintVenantKirchhoffMaterial', young_modulus="3000", poisson_ratio="0.3")
        sofa_node.addObject('HyperelasticForcefield', printLog=True)


        fenics_node = root.addChild("tetra_node_FEniCS")
        fenics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
        fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.fenics_mo = fenics_node.addObject('MechanicalObject', name="mo", position=beam_q2.points.tolist())
        fenics_node.addObject('CaribouTopology', name='topology', template='Hexahedron20', indices=indices.tolist())
        fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
        fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
        fenics_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
        fenics_node.addObject('ConstantForceField', force="0 -100 0", indices="@top_roi.indices")
        fenics_node.addObject('SaintVenantKirchhoffMaterial_FEniCS', young_modulus="3000", poisson_ratio="0.3")
        fenics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        return root

    def onSimulationInitDoneEvent(self, event):
        self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())
        self.fenics_rest_position = np.array(self.fenics_mo.position.value.copy().tolist())
        

    def onAnimateBeginEvent(self, event):

        pass

    def onAnimateEndEvent(self, event): 
        sofa_current_positions = np.array(self.sofa_mo.position.value.copy().tolist())
        fenics_current_positions = np.array(self.fenics_mo.position.value.copy().tolist())

        errors = []
        for sofa_current_point, fenics_current_point, sofa_initial_point in zip(sofa_current_positions, fenics_current_positions, self.sofa_rest_position): 
            if np.linalg.norm(sofa_current_point - sofa_initial_point) != 0:
                errors.append(np.linalg.norm(sofa_current_point - fenics_current_point)/np.linalg.norm(sofa_current_point - sofa_initial_point))

        mean_error = np.mean(np.array(errors))
        
        print(f"Relative Mean Error: {100*mean_error} %")

    
    

def createScene(node):
    node.addObject(ControlFrame(node))




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




# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()