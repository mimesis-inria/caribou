import Sofa
import SofaCaribou
import numpy as np
import meshio

mesh = meshio.read("../meshes/liver.obj")
indices = mesh.cells_dict["triangle"]


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node, refine=10):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node, refine)

    def CreateGraph(self, root, refine=10):
        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin',
                       pluginName="SofaCaribou CImgPlugin SofaMiscCollision SofaGeneralSimpleFem SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine SofaGeneralDeformable SofaMiscFem SofaMeshCollision SofaGeneralLoader")

        root.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-30",
                       absolute_correction_tolerance_threshold="1e-30",
                       absolute_residual_tolerance_threshold="1e-5",
                       relative_residual_tolerance_threshold="1e-10",
                       printLog="1")
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        root.gravity = [0, 0, 0]
        self.refine = refine
        root.addObject('RegularGridTopology', name='grid', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80],
                       n=[12, 12, self.refine])

        meca = root.addChild("mechanical")
        # - Mechanics
        self.dofs = meca.addObject('MechanicalObject', name='mo', src='@../grid')

        # Complete hexa container (only needed for the BC BoxROI later on, a bug from SOFA's side)
        # meca.addObject('HexahedronSetTopologyContainer', src='@../grid', name='mechanical_topology')
        # meca.addObject('HexahedronSetTopologyContainer', hexahedra='@grid.hexahedra',
        #                name='material_1_hexahedrons')
        meca.addObject('TetrahedronSetTopologyContainer', name='topo')
        meca.addObject('TetrahedronSetTopologyModifier')
        meca.addObject('Hexa2TetraTopologicalMapping', input='@../grid', output='@topo', swapping=True)
        # meca.addObject('NeoHookeanMaterial', young_modulus=3000, poisson_ratio=0.3, name='material_1')
        # meca.addObject('HyperelasticForcefield', topology="@topo", material='@material_1')

        # meca.addObject('FEniCS_Material', template="Tetrahedron", young_modulus=3000,
        #                poisson_ratio=0.3, material_name="NeoHookean")
        meca.addObject('FEniCS_Material', template="Tetrahedron", bulk_modulus=1000000, a=1180, b=8, a_f=18.5 * 10e4,
                       b_f=16, a_s=2.5 * 10e4, b_s=11.1, a_fs=2160, b_fs=11.4, material_name="Ogden")

        meca.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        # Fix the left side of the beam
        meca.addObject('BoxROI', name='fixed_roi', quad='@../grid.quad', box=[-7.5, -7.5, -0.9, 7.5, 7.5, 0.1])
        meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

        # Apply traction on the right side of the beam
        meca.addObject('BoxROI', name='top_roi', box=[-7.5, -7.5, 79.9, 7.5, 7.5, 80.1])
        # meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
        # meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1 / 5, topology='@quad_container')
        # meca.addObject("ConstantForceField", indices="@top_roi.indices", totalForce=[0, -1000, 0])
        meca.addObject("ConstantForceField", indices="@top_roi.indices", totalForce=[0, 0, 100000])

        return root

    def onSimulationInitDoneEvent(self, event):
        pass

    def onAnimateBeginEvent(self, event):
        pass

    def onAnimateEndEvent(self, event):
        current_positions = np.array(self.dofs.position.value.copy().tolist()[-1])
        rest_positions = np.array(self.dofs.rest_position.value.copy().tolist()[-1])
        displacement = []
        for current_point, initial_point in zip(current_positions, rest_positions):
            if np.linalg.norm(current_point - initial_point) != 0:
                displacement.append(np.linalg.norm(current_point - initial_point))
        with open('./convergence_study/displacement_' + str(len(self.dofs.position)) + ".txt", 'w') as f:
            f.write(str(np.mean(displacement)))
        print(np.mean(displacement))


def createScene(node, refine):
    node.addObject(ControlFrame(node, refine))


# Choose in your script to activate or not the GUI
USE_GUI = False


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    if not USE_GUI:
        refinement = np.arange(200, 250, 10)
        for refine in refinement:
            root = Sofa.Core.Node("root")
            createScene(root, refine)
            Sofa.Simulation.init(root)
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        root = Sofa.Core.Node("root")
        createScene(root)
        Sofa.Simulation.init(root)
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()
