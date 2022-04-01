import Sofa
import SofaCaribou
import meshio
import numpy as np
from manufactured_solution import assemble, integrate, compute_solution, ConstantForceField

ELEMENT_TYPE = "Hexahedron"
ELEMENT_APPROXIMATION_DEGREE = 2
MATERIAL_MODEL = "SaintVenantKirchhoff"
mu = 1.0
l  = 1.25
rad = 1
length = 3
poisson_ratio = 0.3
young_modulus = 3000
P, f, u_s = compute_solution(mu, l, rad, length)
# TODO improve the manual permutation for matching the redefinition of the hexahedron
# TODO redefine the visualization of the hexaedron

""" if ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 1:
    element = "Tetrahedron"
    mesh = meshio.read("./meshes/beam_p1.vtu")
    indices = np.empty(mesh.cells_dict['tetra'].shape)
    indices = mesh.cells_dict['tetra']
elif ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
    element = "Tetrahedron10"
    mesh = meshio.read("./meshes/beam_p2.vtu")
    indices = np.empty(mesh.cells_dict['tetra10'].shape)
    indices = mesh.cells_dict['tetra10'][:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
elif ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 1:
    element = "Hexahedron"
    mesh = meshio.read("./meshes/beam_q1.vtu")
    indices = np.empty(mesh.cells_dict['hexahedron'].shape)
    indices = mesh.cells_dict['hexahedron'][:, [4, 5, 0, 1, 7, 6, 3, 2]]
    #indices = mesh.cells_dict['hexahedron'][:, [0, 1, 2, 3, 4, 5, 6, 7]]
elif ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
    element = "Hexahedron20"
    mesh = meshio.read("./meshes/beam_q2.vtu")
    indices = np.empty(mesh.cells_dict['hexahedron20'].shape)
    indices = mesh.cells_dict['hexahedron20'][:,
              [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]

else:
    raise ValueError('The element or the approximation order is not implemented yet.')
 """

if MATERIAL_MODEL == "SaintVenantKirchhoff" or MATERIAL_MODEL == "NeoHookean":
    material = MATERIAL_MODEL + "Material"
else:
    raise ValueError('The material model is not implemented yet.')

element = "Tetrahedron"
mesh = meshio.read("./meshes/cylinder_p1.vtu")
indices = np.empty(mesh.cells_dict['tetra'].shape)
indices = mesh.cells_dict['tetra']

class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):

        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin',
                       pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")

        sofa_node = root.addChild("sofa_node")
        sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                            relative_residual_tolerance_threshold="1e-10", printLog="1")
        sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        sofa_node.addObject('MechanicalObject', name="mo", position=mesh.points.tolist())
        sofa_node.addObject('CaribouTopology', name='topology', template=element, indices=mesh.cells[1].data.tolist())
        sofa_node.addObject('CaribouTopology', name='dirichlet_boundary', template='Triangle', indices=mesh.cells[0].data[np.ma.masked_equal(mesh.cell_data['gmsh:physical'][0], 1).mask].tolist())
        sofa_node.addObject('CaribouTopology', name='neumann_boundary',   template='Triangle', indices=mesh.cells[0].data[np.ma.masked_inside(mesh.cell_data['gmsh:physical'][0], 2, 3).mask].tolist())
   
        sofa_node.addObject(material, young_modulus=young_modulus, poisson_ratio=poisson_ratio)
        sofa_node.addObject('HyperelasticForcefield', printLog=True)

        sofa_node.addObject('FixedConstraint', indices=np.unique(root.sofa_node.dirichlet_boundary.indices.array()).tolist())
        """ sofa_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
        sofa_node.addObject('ConstantForceField', name='external_forces', force="0 -100 0", indices="@top_roi.indices")
         """
        sofa_node.addObject(ConstantForceField(name='external_forces'))
        """ fenics_node = root.addChild("fenics_node")
        fenics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                              relative_residual_tolerance_threshold="1e-10", printLog="1")
        fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.fenics_mo = fenics_node.addObject('MechanicalObject', name="mo", position=mesh.points.tolist())
        fenics_node.addObject('CaribouTopology', name='topology', template=element,
                              indices=indices.tolist())
        fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
        fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
        fenics_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
        fenics_node.addObject('ConstantForceField', force="0 -100 0", indices="@top_roi.indices")
        fenics_node.addObject(material + '_FEniCS', template=element, young_modulus="3000",
                              poisson_ratio="0.3")
        fenics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)
         """
        return root

    def onSimulationInitDoneEvent(self, event):
        #self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())
        #self.fenics_rest_position = np.array(self.fenics_mo.position.value.copy().tolist())
        pass
    def onAnimateBeginEvent(self, event):

        pass

    def onAnimateEndEvent(self, event):
        #sofa_current_positions = np.array(self.sofa_mo.position.value.copy().tolist())
        #fenics_current_positions = np.array(self.fenics_mo.position.value.copy().tolist())
        pass
        """ errors = []
        for sofa_current_point, fenics_current_point, sofa_initial_point in zip(sofa_current_positions,
                                                                                fenics_current_positions,
                                                                                self.sofa_rest_position):
            if np.linalg.norm(sofa_current_point - sofa_initial_point) != 0:
                errors.append(np.linalg.norm(sofa_current_point - fenics_current_point) / np.linalg.norm(
                    sofa_current_point - sofa_initial_point))
                
                
        mean_error = np.mean(np.array(errors))

        print(f"Relative Mean Error: {100 * mean_error} %") """


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
    
    return root


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    root = main()
    print('Assembling the external force vector...', end='', flush=True)
    external_forces = \
        assemble(root.sofa_node.topology.domain(), lambda x, y, z, _: f(x, y, z)) + \
        assemble(root.sofa_node.neumann_boundary.domain(), lambda x, y, z, t: np.dot(P(x, y, z), t.normal()))
    print(' Done.', flush=True)

    exact_error = np.sqrt(integrate(root.sofa_node.topology.domain(), lambda x, y, z, _: np.dot(u_s(x, y, z), u_s(x, y, z))))
    print(f"Exact error is {exact_error}")

    for load in [1e-3, 1e-2, 1e-1, 0.15, 0.5, 1.0]:
        root.sofa_node.external_forces.forces = (external_forces*load)
        Sofa.Simulation.animate(root, 1)
        u_h = (root.sofa_node.mo.position.array() - root.sofa_node.mo.rest_position.array())
        error_L2 = np.sqrt(integrate(
            root.sofa_node.topology.domain(),
            lambda x, y, z, u_g, _: np.dot((u_s(x, y, z) - u_g), (u_s(x, y, z) - u_g)),
            u_h
        ))
        print(f"relative L2 error at {int(load*100)}% load: {error_L2/exact_error}")
