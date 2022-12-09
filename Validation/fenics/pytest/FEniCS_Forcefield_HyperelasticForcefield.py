import sys
import unittest
from parameterized import parameterized
from pathlib import Path

import Sofa
import SofaCaribou
import numpy as np
from scipy.sparse import csr_matrix

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
    # ("Kirchhoff_hexa20", "SaintVenantKirchhoff", "Hexahedron20"),
    ("NeoHooke_tetra", "NeoHookean", "Tetrahedron"),
    ("NeoHooke_tetra10", "NeoHookean", "Tetrahedron10"),
    ("NeoHooke_hexa", "NeoHookean", "Hexahedron"),
    # ("NeoHooke_hexa20", "NeoHookean", "Hexahedron20"),
]


def generate_geometry(element):
    hexa_p_0 = [0, 0, 0]
    hexa_p_1 = [1, 0, 0]
    hexa_p_2 = [1, 1, 0]
    hexa_p_3 = [0, 1, 0]
    hexa_p_4 = [0, 0, 1]
    hexa_p_5 = [1, 0, 1]
    hexa_p_6 = [1, 1, 1]
    hexa_p_7 = [0, 1, 1]
    hexa_p_8 = [0.5, 0, 0]
    hexa_p_9 = [1, 0.5, 0]
    hexa_p_10 = [0.5, 1, 0]
    hexa_p_11 = [0, 0.5, 0]
    hexa_p_12 = [0.5, 0, 1]
    hexa_p_13 = [1, 0.5, 1]
    hexa_p_14 = [0.5, 1, 1]
    hexa_p_15 = [0, 0.5, 1]
    hexa_p_16 = [0, 0, 0.5]
    hexa_p_17 = [1, 0, 0.5]
    hexa_p_18 = [1, 1, 0.5]
    hexa_p_19 = [0, 1, 0.5]

    tetra_p_0 = [0, 0, 0]
    tetra_p_1 = [1, 0, 0]
    tetra_p_2 = [0, 1, 0]
    tetra_p_3 = [0, 0, 1]
    tetra_p_4 = [0.5, 0, 0]
    tetra_p_5 = [0.5, 0.5, 0]
    tetra_p_6 = [0, 0.5, 0]
    tetra_p_7 = [0, 0, 0.5]
    tetra_p_8 = [0.5, 0, 0.5]
    tetra_p_9 = [0, 0.5, 0.5]

    if element == "Tetrahedron":
        return element, element, np.array([tetra_p_0, tetra_p_1, tetra_p_2, tetra_p_3]), np.arange(4), np.arange(4)
    elif element == "Tetrahedron10":
        return element, element, np.array(
            [tetra_p_0, tetra_p_1, tetra_p_2, tetra_p_3, tetra_p_4, tetra_p_5, tetra_p_6, tetra_p_7, tetra_p_8,
             tetra_p_9]), np.arange(10), np.array([0, 1, 2, 3, 9, 8, 5, 7, 6, 4])
    elif element == "Hexahedron":
        return element, element + "_FEniCS", np.array(
            [hexa_p_0, hexa_p_1, hexa_p_2, hexa_p_3, hexa_p_4, hexa_p_5, hexa_p_6, hexa_p_7]), np.arange(8), np.array(
            [0, 1, 3, 2, 4, 5, 7, 6])
    elif element == "Hexahedron20":
        return element, "Hexahedron_FEniCS20", np.array(
            [hexa_p_0, hexa_p_1, hexa_p_2, hexa_p_3, hexa_p_4, hexa_p_5, hexa_p_6, hexa_p_7, hexa_p_8, hexa_p_9,
             hexa_p_10, hexa_p_11, hexa_p_12, hexa_p_13, hexa_p_14, hexa_p_15, hexa_p_16, hexa_p_17, hexa_p_18,
             hexa_p_19]), np.arange(20), np.array(
            [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]
        )


def createScene(node, element_type, material_model):
    node.addObject('DefaultVisualManagerLoop')
    node.addObject('DefaultAnimationLoop')
    node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    node.addObject('RequiredPlugin',
                   pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")
    element_sofa, element_fenics, position, indices_sofa, indices_fenics = generate_geometry(element_type)

    sofa_node = node.addChild("sofa_node")
    sofa_node.addObject('StaticSolver', newton_iterations="1", relative_correction_tolerance_threshold="1e-15",
                        relative_residual_tolerance_threshold="1e-10", printLog="1")
    sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sofa_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    sofa_node.addObject('CaribouTopology', name='topology', template=element_sofa, indices=indices_sofa.tolist())
    sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    sofa_node.addObject(material_model + "Material", young_modulus=YOUNG_MODULUS, poisson_ratio=POISSON_RATIO)
    sofa_node.addObject('HyperelasticForcefield', name="ff", printLog=False)

    fenics_node = node.addChild("fenics_node")
    fenics_node.addObject('StaticSolver', newton_iterations="1", relative_correction_tolerance_threshold="1e-15",
                          relative_residual_tolerance_threshold="1e-10", printLog="1")
    fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    fenics_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    fenics_node.addObject('CaribouTopology', name='topology', template=element_fenics,
                          indices=indices_fenics.tolist())
    fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
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
                K_fenics = csr_matrix(root.fenics_node.ff.K(), copy=True)
                K_sofa = csr_matrix(root.sofa_node.ff.K(), copy=True)
                self.assertMatrixQuasiEqual(K_fenics, K_sofa)

    def assertMatrixQuasiEqual(self, A, B):
        """ absolute(a - b) <= (atol + rtol * absolute(b)) """
        return self.assertTrue(np.allclose(A.todense(), B.todense()),
                               f"Matrices are not equal, with \nA={A} \nand\nB={B}")


if __name__ == '__main__':
    unittest.main()
