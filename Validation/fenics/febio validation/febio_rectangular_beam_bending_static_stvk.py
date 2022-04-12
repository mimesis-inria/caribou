#!/usr/bin/python3.7

"""
Bending rectangular beam simulated with FEBio.

Dimensions: 15x15x80
Material: St-Venant-Kirchhoff (young modulus 3000, poisson ratio 0.499)
ODE solver: Static Newton-Raphson
Left face: Clamped
Right face: Traction of [0, -30, 0]
Number of load increments: 5

To run with podman:
podman run --rm --userns keep-id \
    -v $PWD:/opt/shared:z \
    -w /opt/shared \
    febio_validation python febio_rectangular_beam_bending_static_stvk.py
"""
import Sofa
import SofaCaribou
import numpy as np
import os, meshio, tempfile, subprocess
from febio_scene import *
from parameters import parameters, meshio_to_febio_element_types




# SOniCS scene and displacement computing
class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node, febio_current_points):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)
        self.febio_current_points = febio_current_points

    def CreateGraph(self, root):

        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin',
                       pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")

        fenics_node = root.addChild("fenics_node")
        fenics_node.addObject('StaticSolver', newton_iterations=parameters['nnewtonsteps'], relative_correction_tolerance_threshold=parameters['displacement_tolerance'],
                              relative_residual_tolerance_threshold=parameters['residual_tolerance'], printLog="1")
        fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.fenics_mo = fenics_node.addObject('MechanicalObject', name="mo", position=parameters['mesh'].points.tolist())
        fenics_node.addObject('CaribouTopology', name='topology', template=parameters['element'],
                              indices=parameters['indices'].tolist())
        fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
        fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
        fenics_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
        fenics_node.addObject('ConstantForceField', force=parameters['traction'].tolist(), indices="@top_roi.indices")
        fenics_node.addObject('FEniCS_Material', template=parameters['element'], young_modulus="3000",
                              poisson_ratio="0.3", C01=parameters['c01'], C10=parameters['c10'], k=parameters['k'], material_name=parameters['sonics_material'], path="/home/..")

        fenics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)


        return root

    def onSimulationInitDoneEvent(self, event):
        self.sonics_initial_points = np.array(self.fenics_mo.position.value.copy().tolist())

    def onAnimateBeginEvent(self, event):

        pass

    def onAnimateEndEvent(self, event):
        sonics_current_positions = np.array(self.fenics_mo.position.value.copy().tolist())
        sonics_displacement = sonics_current_positions - self.sonics_initial_points
        
        febio_displacement = self.febio_current_points - self.sonics_initial_points
        
        errors = []
        
        for febio_current_point, sonics_current_point, sonics_initial_point in zip(self.febio_current_points,
                                                                                sonics_current_positions,
                                                                                self.sonics_initial_points):
            if np.linalg.norm(sonics_current_point - sonics_initial_point) != 0:
                errors.append(np.linalg.norm(febio_current_point - sonics_current_point) / np.linalg.norm(
                    sonics_current_point - sonics_initial_point))
                
                
        mean_error = np.mean(np.array(errors))

        print(f"Relative Mean Error: {100 * mean_error} %")



def createScene(node, febio_current_points):
    node.addObject(ControlFrame(node, febio_current_points))


# Choose in your script to activate or not the GUI
USE_GUI = True

def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    root = Sofa.Core.Node("root")


    # FEBIO displacements computing 
    elements = sort_elements(meshio_to_febio_element_types, parameters)
    febio_current_points = febio_scene(parameters, elements)
    #print(displacements_febio)

    createScene(root, febio_current_points)
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