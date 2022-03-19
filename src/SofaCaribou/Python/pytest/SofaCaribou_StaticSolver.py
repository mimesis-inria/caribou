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
from SofaCaribou import NewtonRaphsonSolver

E = NewtonRaphsonSolver.Event

number_of_steps = 10
number_of_newton_iterations = 10
threshold = 1e-15
radius = 7.5
length = 80
nx = 3
nz = 9
eps = 1e-5


def createScene(node):
    node.addObject('RequiredPlugin', pluginName=[
        'SofaBaseMechanics', 'SofaEngine', 'SofaGraphComponent', 'SofaOpenglVisual',
        'SofaPreconditioner', 'SofaBoundaryCondition', 'SofaTopologyMapping'])
    node.addObject('RegularGridTopology', name='grid', min=[-radius, -radius, 0], max=[radius, radius, length], n=[nx, nx, nz])
    node.addObject('StaticODESolver', name='solver', newton_iterations=number_of_newton_iterations, correction_tolerance_threshold=1e-5, residual_tolerance_threshold=1e-5, printLog=False)
    node.addObject('LDLTSolver', backend='Pardiso')
    node.addObject('MechanicalObject', name='mo', position='@grid.position')
    node.addObject('HexahedronSetTopologyContainer', name='mechanical_topology', src='@grid')
    node.addObject('SaintVenantKirchhoffMaterial', young_modulus=3000, poisson_ratio=0.499)
    node.addObject('HyperelasticForcefield', name='ff', topology='@mechanical_topology')
    node.addObject('BoxROI', name='base_roi', box=[-radius - eps, -radius - eps,    0   - eps, radius + eps, radius + eps,    0   + eps])
    node.addObject('BoxROI', name='top_roi', box=[ -radius - eps, -radius - eps, length - eps, radius + eps, radius + eps, length + eps], quad='@mechanical_topology.quads')
    node.addObject('FixedConstraint', indices='@base_roi.indices')
    node.addObject('QuadSetTopologyContainer', name='neumann_topology', quads='@top_roi.quadInROI')
    node.addObject('TractionForcefield', topology='@neumann_topology', traction=[0, -30, 0], slope=0.2)


class TestStaticSolver(unittest.TestCase):
    def assertMatrixEqual(self, A, B):
        return self.assertTrue(np.allclose(A.todense(), B.todense()), f"Matrices are not equal, with \nA={A} \nand\nB={B}")
    def assertVectorAlmostEqual(self, A, B):
        return self.assertTrue(np.allclose(A, B), f"Vectors are not equal, with \nA={A} \nand\nB={B}")

    def test_residuals(self):
        root = Sofa.Core.Node()
        createScene(root)
        Sofa.Simulation.init(root)

        # The simulation is supposed to converged in 5 load increments. Let's check the norms of force residuals.
        # todo(jnbrunet): These values should be validated against another external software
        force_residuals = [
            [1.000000000000000e+00, 4.802014099846802e-03, 3.457785823647809e-04, 6.417769508095073e-08],  # Step 1
            [1.000000000000000e+00, 4.700052322118300e-03, 3.771851553383748e-04, 1.289757638852382e-05, 2.125041655083045e-08],  # Step 2
            [1.000000000000000e+00, 4.416064536488316e-03, 4.703733156787435e-04, 3.659722934024732e-05, 7.791507698853802e-08],  # Step 3
            [1.000000000000000e+00, 4.003458215116851e-03, 6.263051713537671e-04, 4.996976075176631e-05, 1.422069948624354e-07],  # Step 4
            [1.000000000000000e+00, 3.526942203674829e-03, 8.307813177405512e-04, 4.667215114394798e-05, 1.646730071539153e-07],  # Step 5
        ]

        def increment(n): n +=1

        step_number = 0
        count = [0, 0, 0, 0, 0, 0, 0, 0]
        root.solver.register_callback(E.ITERATION_BEGIN,   lambda s: increment(count[E.ITERATION_BEGIN]))
        root.solver.register_callback(E.MATRIX_ASSEMBLED,  lambda s: increment(count[E.MATRIX_ASSEMBLED]))
        root.solver.register_callback(E.MATRIX_ANALYZED,   lambda s: increment(count[E.MATRIX_ANALYZED]))
        root.solver.register_callback(E.MATRIX_FACTORIZED, lambda s: increment(count[E.MATRIX_FACTORIZED]))
        root.solver.register_callback(E.INCREMENT_SOLVED,  lambda s: increment(count[E.INCREMENT_SOLVED]))
        root.solver.register_callback(E.INCREMENT_PROPAGATED, lambda s: increment(count[E.INCREMENT_PROPAGATED]))
        root.solver.register_callback(E.RESIDUAL_UPDATED,  lambda s: increment(count[E.RESIDUAL_UPDATED]))
        root.solver.register_callback(E.ITERATION_END,     lambda s: increment(count[E.ITERATION_END]))

        for step_number in range(5):
            Sofa.Simulation.animate(root, 1)
            self.assertEqual(root.solver.current_iteration, len(force_residuals[step_number]), f"for step {step_number}")
            r0 = root.solver.squared_residuals[0]
            residuals = [np.sqrt(r/r0) for r in root.solver.squared_residuals]
            self.assertVectorAlmostEqual(residuals, force_residuals[step_number])

            self.assertEqual(count[E.ITERATION_BEGIN], root.solver.current_iteration)
            count = [0, 0, 0, 0, 0, 0, 0, 0]


if __name__ == '__main__':
    unittest.main()
