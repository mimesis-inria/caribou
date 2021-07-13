import Sofa, SofaRuntime
import sys
import unittest
import numpy as np
from scipy.sparse import csr_matrix
from pathlib import Path

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')

import SofaCaribou

number_of_steps = 10
number_of_newton_iterations = 10
threshold = 1e-15
radius = 5
length = 60
cell_size=1.5
nx = int(2 * radius / cell_size) + 1
nz = int(length / cell_size) + 1
eps = cell_size / 10


def createScene(node):
    node.addObject('RequiredPlugin', pluginName=[
        'SofaBaseMechanics', 'SofaEngine', 'SofaGraphComponent', 'SofaOpenglVisual',
        'SofaPreconditioner', 'SofaBoundaryCondition', 'SofaSparseSolver', 'SofaTopologyMapping'])
    node.addObject('RegularGridTopology', name='grid', min=[-radius, -radius, -length / 2], max=[radius, radius, length / 2], n=[nx, nx, nz])
    node.addObject('StaticODESolver', newton_iterations=number_of_newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=False)
    node.addObject('LLTSolver', backend='Pardiso')
    node.addObject('MechanicalObject', name='mo', position='@grid.position')
    node.addObject('HexahedronSetTopologyContainer', name='mechanical_topology', src='@grid')
    node.addObject('SaintVenantKirchhoffMaterial', young_modulus=3000, poisson_ratio=0)
    node.addObject('HyperelasticForcefield', name='ff', topology='@mechanical_topology')
    node.addObject('BoxROI', name='base_roi', box=[-radius - eps, -radius - eps, -length / 2 - eps, radius + eps, radius + eps, -length / 2 + eps])
    node.addObject('BoxROI', name='top_roi', box=[-radius - eps, -radius - eps, +length / 2 - eps, radius + eps, radius + eps, +length / 2 + eps], quad='@mechanical_topology.quads')
    node.addObject('FixedConstraint', indices='@base_roi.indices')
    node.addObject('QuadSetTopologyContainer', name='neumann_topology', quads='@top_roi.quadInROI')
    node.addObject('TractionForcefield', topology='@neumann_topology', traction=[0, -30, 0], slope=1 / 5)


class TestHyperelasticForcefield(unittest.TestCase):
    def assertMatrixEqual(self, A, B):
        return self.assertTrue(np.array_equal(A.todense(), B.todense()), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixNotEqual(self, A, B):
        return self.assertFalse(np.array_equal(A.todense(), B.todense()), f"Matrix are equal, with \nA{A}\nand\nB={B}")

    def test_stiffness_assembly(self):
        root = Sofa.Core.Node()
        createScene(root)
        Sofa.Simulation.init(root)
        x = np.array(root.mo.position.array(), dtype=np.float64, order='C', copy=False)
        K1 = csr_matrix(root.ff.K(), copy=True)

        # Manually trigger the matrix assembly using a different position vector
        x2 = x*10
        root.ff.assemble_stiffness(x2)
        K2 = csr_matrix(root.ff.K(), copy=False)
        self.assertMatrixNotEqual(K1, K2)

        # Manually trigger the matrix assembly using the same position vector
        root.ff.assemble_stiffness(x)
        K3 = csr_matrix(root.ff.K(), copy=False)
        self.assertMatrixEqual(K1, K3)


if __name__ == '__main__':
    unittest.main()
