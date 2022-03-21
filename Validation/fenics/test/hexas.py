import sys
import unittest
from pathlib import Path

import Sofa
import SofaCaribou
import numpy as np
from scipy.sparse import csr_matrix, linalg

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')

number_of_steps = 10
number_of_newton_iterations = 10
threshold = 1e-15
radius = 5
length = 60
cell_size = 1.5
nx = int(2 * radius / cell_size) + 1
nz = int(length / cell_size) + 1
eps = cell_size / 10


# TODO: first draft of the tests to compare the 2 stiffness matrices.
#  Need to be done for each element, interpolation, material law.
#  We need to find a threshold for the matrices differences.


def createScene(node):
    node.addObject('DefaultVisualManagerLoop')
    node.addObject('DefaultAnimationLoop')
    node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    node.addObject('RequiredPlugin',
                   pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")

    p_0 = [0, 0, 0]
    p_1 = [1, 0, 0]
    p_2 = [1, 1, 0]
    p_3 = [0, 1, 0]
    p_4 = [0, 0, 1]
    p_5 = [1, 0, 1]
    p_6 = [1, 1, 1]
    p_7 = [0, 1, 1]

    position = np.array(
        [p_0, p_1, p_2, p_3, p_4, p_5, p_6, p_7])
    indices = np.array([4, 5, 0, 1, 7, 6, 3, 2])

    sofa_node = node.addChild("sofa_node")
    sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                        relative_residual_tolerance_threshold="1e-10", printLog="1")
    sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sofa_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    sofa_node.addObject('CaribouTopology', name='topology', template='Hexahedron', indices=indices.tolist())
    sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    sofa_node.addObject('NeoHookeanMaterial', young_modulus="3000", poisson_ratio="0.3")
    sofa_node.addObject('HyperelasticForcefield', name="ff", printLog=True)

    fenics_node = node.addChild("fenics_node")
    fenics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                          relative_residual_tolerance_threshold="1e-10", printLog="1")
    fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    fenics_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    fenics_node.addObject('CaribouTopology', name='topology', template='Hexahedron', indices=indices.tolist())
    fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    fenics_node.addObject('NeoHookeanMaterial_FEniCS', template="Hexahedron", young_modulus="3000",
                          poisson_ratio="0.3")
    fenics_node.addObject('HyperelasticForcefield_FEniCS', name="ff", printLog=True)


class TestHyperelasticForcefield(unittest.TestCase):
    def assertMatrixEqual(self, A, B):
        return self.assertTrue(np.array_equal(A.todense(), B.todense()),
                               f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixNotEqual(self, A, B):
        return self.assertFalse(np.array_equal(A.todense(), B.todense()), f"Matrix are equal, with \nA{A}\nand\nB={B}")

    def test_stiffness_assembly(self):
        root = Sofa.Core.Node()
        createScene(root)
        Sofa.Simulation.init(root)
        # x_fenics = np.array(root.fenics_node.mo.position.array(), dtype=np.float64, order='C', copy=False)
        # x_sofa = np.array(root.sofa_node.mo.position.array(), dtype=np.float64, order='C', copy=False)
        K_fenics = csr_matrix(root.fenics_node.ff.K(), copy=True)
        K_sofa = csr_matrix(root.sofa_node.ff.K(), copy=True)
        print(linalg.norm(K_fenics - K_sofa))


if __name__ == '__main__':
    unittest.main()
