import sys
import unittest
from parameterized import parameterized
from pathlib import Path

import Sofa
import SofaCaribou
import numpy as np
from scipy.sparse import csr_matrix
import meshio

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')
YOUNG_MODULUS = 3000
POISSON_RATIO = 0.3

TEST_CASES = [
    ("Kirchhoff_tetra", "SaintVenantKirchhoff", "Tetrahedron"),
    ("Kirchhoff_tetra10", "SaintVenantKirchhoff", "Tetrahedron10"),
    ("Kirchhoff_hexa", "SaintVenantKirchhoff", "Hexahedron"),
    ("Kirchhoff_hexa20", "SaintVenantKirchhoff", "Hexahedron20"),
    ("NeoHooke_tetra", "NeoHookean", "Tetrahedron"),
    ("NeoHooke_tetra10", "NeoHookean", "Tetrahedron10"),
    ("NeoHooke_hexa", "NeoHookean", "Hexahedron"),
    ("NeoHooke_hexa20", "NeoHookean", "Hexahedron20"),
]


def generate_geometry(element):
    element_sofa, element_fenics, mesh, indices_sofa, indices_fenics, fixed_indices = None, None, None, None, None, None
    if element == "Tetrahedron":
        element_sofa = "Tetrahedron"
        element_fenics = element_sofa
        mesh = meshio.read("../../meshes/fenics/cube_tetra_1.msh")
        indices_sofa = mesh.cells_dict['tetra']
        indices_fenics = indices_sofa
        fixed_indices = get_fixed_indices(mesh, "triangle")

    elif element == "Tetrahedron10":
        element_sofa = "Tetrahedron10"
        element_fenics = element_sofa
        mesh = meshio.read("../../meshes/fenics/cube_tetra_2.msh")
        indices_sofa = mesh.cells_dict['tetra10']
        indices_fenics = indices_sofa[:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
        fixed_indices = get_fixed_indices(mesh, "triangle6")

    elif element == "Hexahedron":
        element_sofa = "Hexahedron"
        element_fenics = element_sofa + "_FEniCS"
        mesh = meshio.read("../../meshes/fenics/cube_hexa_1.msh")
        indices_sofa = mesh.cells_dict['hexahedron']
        indices_fenics = indices_sofa[:, [0, 1, 3, 2, 4, 5, 7, 6]]
        fixed_indices = get_fixed_indices(mesh, "quad")

    elif element == "Hexahedron20":
        element_sofa = "Hexahedron20"
        element_fenics = "Hexahedron_FEniCS20"
        mesh = meshio.read("../../meshes/fenics/cube_hexa_2.msh")
        indices_sofa = mesh.cells_dict['hexahedron20']
        indices_fenics = indices_sofa[:, [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]
        fixed_indices = get_fixed_indices(mesh, "quad8")
    return element_sofa, element_fenics, mesh.points, indices_sofa, indices_fenics, fixed_indices


def get_fixed_indices(mesh, mesh_type):
    a = [i.size for i in mesh.cell_data["gmsh:physical"]]
    indices = np.cumsum(a)
    return mesh.cells_dict[mesh_type][0:indices[0]]


def createScene(node, element_type, material_model):
    node.addObject('DefaultVisualManagerLoop')
    node.addObject('DefaultAnimationLoop')
    node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    node.addObject('RequiredPlugin',
                   pluginName="SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaBoundaryCondition")
    element_sofa, element_fenics, position, indices_sofa, indices_fenics, fixed_indices = generate_geometry(
        element_type)
    node.gravity = [0, -0.1, 0]
    sofa_node = node.addChild("sofa_node")
    sofa_node.addObject('StaticSolver', newton_iterations="15", relative_correction_tolerance_threshold="1e-15",
                        relative_residual_tolerance_threshold="1e-10", printLog="0")
    sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sofa_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    sofa_node.addObject('CaribouTopology', name='topology', template=element_sofa, indices=indices_sofa.tolist())
    sofa_node.addObject('FixedConstraint', indices=fixed_indices.tolist())
    sofa_node.addObject("CaribouMass", density=3)
    sofa_node.addObject(material_model + "Material", young_modulus=YOUNG_MODULUS, poisson_ratio=POISSON_RATIO)
    sofa_node.addObject('HyperelasticForcefield', name="ff", printLog=False)

    fenics_node = node.addChild("fenics_node")
    fenics_node.addObject('StaticSolver', newton_iterations="15", relative_correction_tolerance_threshold="1e-15",
                          relative_residual_tolerance_threshold="1e-10", printLog="0")
    fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    fenics_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    fenics_node.addObject('CaribouTopology', name='topology', template=element_fenics,
                          indices=indices_fenics.tolist())
    fenics_node.addObject('FixedConstraint', indices=fixed_indices.tolist())
    fenics_node.addObject("UniformMass", totalMass=3, separateGravity=True)
    fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    fenics_node.addObject('FEniCS_Material', template=element_fenics, young_modulus=YOUNG_MODULUS,
                          poisson_ratio=POISSON_RATIO, material_name=material_model)
    fenics_node.addObject('HyperelasticForcefield_FEniCS', name="ff", printLog=False)


class TestHyperelasticForcefield(unittest.TestCase):

    def test_stiffness_assembly(self):
        for name, material, element in TEST_CASES:
            with self.subTest(msg=name):
                root = Sofa.Core.Node()
                createScene(root, element, material)
                Sofa.Simulation.init(root)
                Sofa.Simulation.animate(root, root.dt.value)
                position_sofa = root.sofa_node.mo.position.value
                position_fenics = root.fenics_node.mo.position.value

                if element == "Hexahedron20":
                    position_sofa = position_sofa[np.unique(root.sofa_node.topology.indices.value[:, :8])]
                    position_fenics = position_fenics[np.unique(root.fenics_node.topology.indices.value[:, :8])]

                self.assertMatrixQuasiEqual(position_sofa, position_fenics)

    def assertMatrixQuasiEqual(self, A, B):
        """ absolute(a - b) <= (atol + rtol * absolute(b)) """
        return self.assertTrue(np.allclose(A, B),
                               f"Matrices are not equal, with \nA={A} \nand\nB={B}")


if __name__ == '__main__':
    unittest.main()
