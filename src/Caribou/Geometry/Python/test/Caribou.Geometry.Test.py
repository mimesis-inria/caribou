#! python3

import sys
import unittest
import numpy as np

sys.path.insert(0, "@CARIBOU_PYTHON_LIB_PATH@")

import Caribou
from Caribou.Geometry import Segment


def p1(p):
    if len(p) == 1:
        return 5 + 2*p[0]
    elif len(p) == 2:
        return 5 + 2*p[0] + 3*p[1]
    else:
        5 + 2*p[0] + 3*p[1] + 4*p[2]


class TestSegment(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def test_constructor_linear(self):
        s = Segment()
        self.assertTrue((s.nodes() == np.array([-1., 1.])).all())

        # 1D
        s = Segment(-5, 5)
        self.assertTrue((s.nodes() == np.array([-5., 5.])).all())

        s = Segment([-5, 5])
        self.assertTrue((s.nodes() == np.array([-5., 5.])).all())

        s = Segment([[-5], [5]])
        self.assertTrue((s.nodes() == np.array([-5., 5.])).all())

        # 2D
        s = Segment([-5, 4], [5, -4])
        self.assertTrue((s.nodes() == np.array([[-5, 4], [5, -4]])).all())

        s = Segment([[-5, 4], [5, -4]])
        self.assertTrue((s.nodes() == np.array([[-5, 4], [5, -4]])).all())

        # 3D
        s = Segment([-5, 4, 1], [5, -4, 1])
        self.assertTrue((s.nodes() == np.array([[-5, 4, 1], [5, -4, 1]])).all())

        s = Segment([[-5, 4, 1], [5, -4, 1]])
        self.assertTrue((s.nodes() == np.array([[-5, 4, 1], [5, -4, 1]])).all())

    def test_constructor_quadratic(self):
        s = Segment(Caribou.Quadratic)
        self.assertTrue((s.nodes() == np.array([-1., 1., 0.])).all())

        # 1D
        s = Segment(-5, 5, 0)
        self.assertTrue((s.nodes() == np.array([-5., 5., 0])).all())

        s = Segment([-5, 5, 0])
        self.assertTrue((s.nodes() == np.array([-5., 5., 0])).all())

        s = Segment([[-5], [5], [0]])
        self.assertTrue((s.nodes() == np.array([-5., 5., 0])).all())

        # 2D
        s = Segment([-5, 5], [5, -5], [0, 0])
        self.assertTrue((s.nodes() == np.array([[-5, 5], [5, -5], [0, 0]])).all())

        s = Segment([[-5, 5], [5, -5], [0, 0]])
        self.assertTrue((s.nodes() == np.array([[-5, 5], [5, -5], [0, 0]])).all())

        # 3D
        s = Segment([-5, 4, 1], [5, -4, 1], [0, 0, 1])
        self.assertTrue((s.nodes() == np.array([[-5, 4, 1], [5, -4, 1], [0, 0, 1]])).all())

        s = Segment([[-5, 4, 1], [5, -4, 1], [0, 0, 1]])
        self.assertTrue((s.nodes() == np.array([[-5, 4, 1], [5, -4, 1], [0, 0, 1]])).all())

    def test_linear_1D(self):
        # Shape functions
        s = Segment()
        self.assertTrue((s.nodes() == np.array([-1., 1.])).all())
        self.assertTrue((s.L([-1]) == [1, 0]).all())
        self.assertTrue((s.L([+1]) == [0, 1]).all())

        s = Segment(-5.5, 1.1)

        # Center
        self.assertEqual(s.center()[0], -2.2)

        # Inverse transformation
        for gauss_node in s.gauss_nodes():
            self.assertMatrixEqual(gauss_node.position, s.local_coordinates(s.world_coordinates(gauss_node.position)))

        # Contains point
        s.contains_local(s.local_coordinates(s.center()))

        # Integration
        numerical_solution = 0.
        for gauss_node in s.gauss_nodes():
            x = gauss_node.position
            w = gauss_node.weight
            detJ = s.jacobian(x)[0]
            numerical_solution += p1(s.world_coordinates(x)) * w * detJ
        analytic_solution = (5*1.1 + (1.1*1.1)) - (-5.5*5 + (-5.5)*(-5.5))
        self.assertAlmostEqual(numerical_solution, analytic_solution, delta=1e-10)


if __name__ == '__main__':
    unittest.main()
