import sys
import unittest
from pathlib import Path

import Sofa
import SofaCaribou
import meshio
import numpy as np
from scipy.sparse import csr_matrix

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


def createScene(node):
    node.addObject('DefaultVisualManagerLoop')
    node.addObject('DefaultAnimationLoop')
    node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    node.addObject('RequiredPlugin',
                   pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")

    beam_q1 = meshio.read("../meshes/beam_q1.vtu")

    # TODO improve the manual permutation for matching the redefinition of the hexahedron
    # TODO redefine the visualization of the hexaedron
    indices = np.empty(beam_q1.cells_dict['hexahedron'].shape)
    indices = beam_q1.cells_dict['hexahedron'][:, [4, 5, 0, 1, 7, 6, 3, 2]]

    sofa_node = node.addChild("sofa_node")
    sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                        relative_residual_tolerance_threshold="1e-10", printLog="1")
    sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sofa_node.addObject('MechanicalObject', name="mo", position=beam_q1.points.tolist())
    sofa_node.addObject('CaribouTopology', name='topology', template='Hexahedron', indices=indices.tolist())
    sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    sofa_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
    sofa_node.addObject('ConstantForceField', force="0 -1000 0", indices="@top_roi.indices")
    sofa_node.addObject('SaintVenantKirchhoffMaterial', young_modulus="3000", poisson_ratio="0.3")
    sofa_node.addObject('HyperelasticForcefield', name='ff', printLog=True)


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
        x = np.array(root.sofa_node.mo.position.array(), dtype=np.float64, order='C', copy=False)
        K1 = csr_matrix(root.sofa_node.ff.K(), copy=True)

        # Manually trigger the matrix assembly using a different position vector
        x2 = x * 10
        root.sofa_node.ff.assemble_stiffness(x2)
        K2 = csr_matrix(root.sofa_node.ff.K(), copy=False)
        self.assertMatrixNotEqual(K1, K2)

        # Manually trigger the matrix assembly using the same position vector
        root.sofa_node.ff.assemble_stiffness(x)
        K3 = csr_matrix(root.sofa_node.ff.K(), copy=False)
        self.assertMatrixEqual(K1, K3)


if __name__ == '__main__':
    unittest.main()
