#! python3

import sys
import unittest
import numpy as np
from pathlib import Path

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')
import Caribou
from Caribou.Geometry import Segment, Segment3
from Caribou.Geometry import Quad, Quad8
from Caribou.Geometry import Triangle, Triangle6
from Caribou.Geometry import Tetrahedron, Tetrahedron10
from Caribou.Geometry import Hexahedron, Hexahedron20


def p1(p):
    if len(p) == 1:
        return 5 + 2*p[0]
    elif len(p) == 2:
        return 5 + 2*p[0] + 3*p[1]
    else:
        5 + 2*p[0] + 3*p[1] + 4*p[2]


def p2(p):
    if len(p) == 1:
        return 5 + 2*p[0]*p[0]
    elif len(p) == 2:
        return 5 + 2*p[0]*p[1] + 3*p[1]*p[1]
    else:
        return 5 + 2*p[0]*p[1] + 3*p[1]*p[2]+ 4*p[2]*p[2]


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
        s = Segment3()
        self.assertTrue((s.nodes() == np.array([-1., 1., 0.])).all())

        # 1D
        s = Segment3(-5, 5, 0)
        self.assertTrue((s.nodes() == np.array([-5., 5., 0])).all())

        s = Segment3([-5, 5, 0])
        self.assertTrue((s.nodes() == np.array([-5., 5., 0])).all())

        s = Segment3([[-5], [5], [0]])
        self.assertTrue((s.nodes() == np.array([-5., 5., 0])).all())

        # 2D
        s = Segment3([-5, 5], [5, -5], [0, 0])
        self.assertTrue((s.nodes() == np.array([[-5, 5], [5, -5], [0, 0]])).all())

        s = Segment3([[-5, 5], [5, -5], [0, 0]])
        self.assertTrue((s.nodes() == np.array([[-5, 5], [5, -5], [0, 0]])).all())

        # 3D
        s = Segment3([-5, 4, 1], [5, -4, 1], [0, 0, 1])
        self.assertTrue((s.nodes() == np.array([[-5, 4, 1], [5, -4, 1], [0, 0, 1]])).all())

        s = Segment3([[-5, 4, 1], [5, -4, 1], [0, 0, 1]])
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


class TestQuad(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(A, B, rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_constructor_linear(self):
        s = Quad()
        self.assertTrue((s.nodes() == np.array([[-1., -1.], [+1, -1], [+1, +1], [-1, +1]])).all())

        # 2D
        s = Quad([-5, -4], [5, -4], [5, 4], [-5, 4])
        self.assertTrue((s.nodes() == np.array([[-5, -4], [5, -4], [5, 4], [-5, 4]])).all())

        s = Quad([[-5, -4], [5, -4], [5, 4], [-5, 4]])
        self.assertTrue((s.nodes() == np.array([[-5, -4], [5, -4], [5, 4], [-5, 4]])).all())

        # 3D
        s = Quad([-5, -4, 1], [5, -4, 1], [5, 4, 1], [-5, 4, 1])
        self.assertTrue((s.nodes() == np.array([[-5, -4, 1], [5, -4, 1], [5, 4, 1], [-5, 4, 1]])).all())

        s = Quad([[-5, -4, 1], [5, -4, 1], [5, 4, 1], [-5, 4, 1]])
        self.assertTrue((s.nodes() == np.array([[-5, -4, 1], [5, -4, 1], [5, 4, 1], [-5, 4, 1]])).all())

    def test_quadratic_3D(self):
        q = Quad8(
            [-5, -53./15, -53./15], [+5, -53./15, -53./15],
            [+10,+53./15, +53./15], [-10,+53./15, +53./15]
        )
        self.assertEqual(q.number_of_boundary_elements(), 4)
        self.assertEqual(q.number_of_nodes(), 8)

        # Center
        self.assertMatrixEqual(q.center(), [0, 0, 0])

        # Inverse transformation
        for gauss_node in q.gauss_nodes():
            self.assertMatrixAlmostEqual(gauss_node.position, q.local_coordinates(q.world_coordinates(gauss_node.position)))

        for node in q.nodes():
            self.assertMatrixAlmostEqual(node, q.world_coordinates(q.local_coordinates(node)))

        # Contains point
        q.contains_local(q.local_coordinates(q.center()))

        # Integration
        numerical_solution = 0.
        for gauss_node in q.gauss_nodes():
            x = gauss_node.position
            w = gauss_node.weight
            J = q.jacobian(x)
            detJ = np.sqrt(np.linalg.det(J.T.dot(J)))
            numerical_solution += p2(q.world_coordinates(x)) * w * detJ
        analytic_solution = 5116.3690626590305
        self.assertAlmostEqual(numerical_solution, analytic_solution, delta=1e-10)


class TestTriangle(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(A, B, rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_constructor_linear(self):
        t = Triangle()
        self.assertTrue((t.nodes() == np.array([[0, 0], [1, 0], [0, 1]])).all())

        # 2D
        t = Triangle([50, 50], [60, 50], [55, 55])
        self.assertTrue((t.nodes() == np.array([[50, 50], [60, 50], [55, 55]])).all())

        t = Triangle([[50, 50], [60, 50], [55, 55]])
        self.assertTrue((t.nodes() == np.array([[50, 50], [60, 50], [55, 55]])).all())

        # 3D
        t = Triangle([50, 50, 5], [60, 50, 5], [55, 55, 5])
        self.assertTrue((t.nodes() == np.array([[50, 50, 5], [60, 50, 5], [55, 55, 5]])).all())

        t = Triangle([[50, 50, 5], [60, 50, 5], [55, 55, 5]])
        self.assertTrue((t.nodes() == np.array([[50, 50, 5], [60, 50, 5], [55, 55, 5]])).all())

    def test_quadratic_3D(self):
        t = Triangle6(
            [50, 50, 5], [60, 50, 5], [55, 55, 5]
        )
        self.assertEqual(t.number_of_boundary_elements(), 3)
        self.assertEqual(t.number_of_nodes(), 6)

        # Center
        self.assertMatrixAlmostEqual(t.center(), (t.node(0) + t.node(1) + t.node(2)) / 3.)

        # Inverse transformation
        for gauss_node in t.gauss_nodes():
            self.assertMatrixAlmostEqual(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)))

        for node in t.nodes():
            self.assertMatrixAlmostEqual(node, t.world_coordinates(t.local_coordinates(node)))

        # Contains point
        self.assertTrue(t.contains_local(t.local_coordinates(t.center())))

        # Integration
        numerical_solution = 0.
        for gauss_node in t.gauss_nodes():
            x = gauss_node.position
            w = gauss_node.weight
            J = t.jacobian(x)
            detJ = np.sqrt(np.linalg.det(J.T.dot(J)))
            numerical_solution += p2(t.world_coordinates(x)) * w * detJ
        analytic_solution = 164083.3333333333
        self.assertAlmostEqual(numerical_solution, analytic_solution, delta=1e-10)


class TestTetrahedron(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(A, B, rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_constructor_linear(self):
        t = Tetrahedron()
        self.assertTrue((t.nodes() == np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])).all())

        # 3D
        t = Tetrahedron([50, 50, 5], [60, 50, 5], [55, 55, 5], [0, 1, 2])
        self.assertTrue((t.nodes() == np.array([[50, 50, 5], [60, 50, 5], [55, 55, 5], [0, 1, 2]])).all())

        t = Tetrahedron([[50, 50, 5], [60, 50, 5], [55, 55, 5], [0, 1, 2]])
        self.assertTrue((t.nodes() == np.array([[50, 50, 5], [60, 50, 5], [55, 55, 5], [0, 1, 2]])).all())

    def test_quadratic_3D(self):
        t = Tetrahedron10(
            [50, 50, 0], [60, 50, 0], [55, 55, 0], [55, 52.5, -5]
        )
        self.assertEqual(t.number_of_boundary_elements(), 4)
        self.assertEqual(t.number_of_nodes(), 10)

        # Center
        self.assertMatrixAlmostEqual(t.center(), (t.node(0) + t.node(1) + t.node(2) + t.node(3)) / 4.)

        # Inverse transformation
        for gauss_node in t.gauss_nodes():
            self.assertMatrixAlmostEqual(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)))

        for node in t.nodes():
            self.assertMatrixAlmostEqual(node, t.world_coordinates(t.local_coordinates(node)))

        # Contains point
        self.assertTrue(t.contains_local(t.local_coordinates(t.center())))

        # Integration
        numerical_solution = 0.
        for gauss_node in t.gauss_nodes():
            x = gauss_node.position
            w = gauss_node.weight
            J = t.jacobian(x)
            detJ = np.linalg.det(J)
            numerical_solution += p2(t.world_coordinates(x)) * w * detJ
        analytic_solution = -230260.4166666666
        self.assertAlmostEqual(numerical_solution, analytic_solution, delta=1e-10)


class TestHexahedron(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(A, B, rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_constructor_linear(self):
        t = Hexahedron()
        self.assertTrue((t.nodes() == np.array([[-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1], [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1]])).all())

        # 3D
        t = Hexahedron([-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1], [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1])
        self.assertTrue((t.nodes() == np.array([[-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1], [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1]])).all())

        t = Hexahedron([[-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1], [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1]])
        self.assertTrue((t.nodes() == np.array([[-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1], [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1]])).all())

    def test_quadratic_3D(self):
        t = Hexahedron20()
        self.assertEqual(t.number_of_boundary_elements(), 6)
        self.assertEqual(t.number_of_nodes(), 20)

        # Center
        self.assertMatrixAlmostEqual(t.center(), [0, 0, 0])

        # Inverse transformation
        for gauss_node in t.gauss_nodes():
            self.assertMatrixAlmostEqual(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)))

        for node in t.nodes():
            self.assertMatrixAlmostEqual(node, t.world_coordinates(t.local_coordinates(node)))

        # Contains point
        self.assertTrue(t.contains_local(t.local_coordinates(t.center())))

        # Integration
        numerical_solution = 0.
        for gauss_node in t.gauss_nodes():
            x = gauss_node.position
            w = gauss_node.weight
            J = t.jacobian(x)
            detJ = np.linalg.det(J)
            numerical_solution += p2(t.world_coordinates(x)) * w * detJ
        analytic_solution = 50.666666666666664
        self.assertAlmostEqual(numerical_solution, analytic_solution, delta=1e-10)


if __name__ == '__main__':
    unittest.main()
