#! python3

import sys, os
import unittest
import numpy as np
import meshio

sys.path.insert(0, "@CARIBOU_PYTHON_LIB_PATH@")

import Caribou
from Caribou.Topology import Mesh
from Caribou.Geometry import Segment
from Caribou.Geometry import Quad
from Caribou.Geometry import Triangle
from Caribou.Geometry import Tetrahedron
from Caribou.Geometry import Hexahedron


class TestMesh(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def test_constructor_1d(self):
        nodes = np.array([1., 2., 3., 4., 5.])
        mesh = Mesh(nodes)
        self.assertMatrixEqual([1., 2., 3., 4., 5.], mesh.positions([0, 1, 2, 3, 4]))

        mesh = Mesh([1., 2., 3., 4., 5.])
        self.assertMatrixEqual([1., 2., 3., 4., 5.], mesh.positions([0, 1, 2, 3, 4]))

    def test_contructor_3d(self):
        mesh = Mesh([[1,2,3], [4,5,6]])
        self.assertMatrixEqual([[1,2,3], [4,5,6]], mesh.positions([0, 1]))


class TestDomain(unittest.TestCase):
    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B):
        return np.allclose(A, B)

    def test_segment_domain(self):
        for order, element_type, meshfile in [(Caribou.Linear, 'line', '1D_linear.vtk'), (Caribou.Quadratic, 'line3', '1D_quadratic.vtk')]:
            m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', meshfile))
            mesh = Mesh(m.points[:, 0])
            domain = mesh.add_domain("segments", Segment(Caribou._1D, order), m.cells_dict[element_type])
            self.assertEqual(domain.number_of_elements(), 10)
            self.assertMatrixEqual(domain.element_indices(0), m.cells_dict[element_type][0])
            self.assertMatrixEqual(domain.element_indices(9), m.cells_dict[element_type][9])

    def test_triangle_domain(self):
        for dimension, order, element_type, meshfile in [
            (Caribou._2D, Caribou.Linear, 'triangle', '2D_triangle_linear.vtk'),
            (Caribou._3D, Caribou.Linear, 'triangle', '3D_triangle_linear.vtk'),
            (Caribou._2D, Caribou.Quadratic, 'triangle6', '2D_triangle_quadratic.vtk'),
            (Caribou._3D, Caribou.Quadratic, 'triangle6', '3D_triangle_quadratic.vtk'),
        ]:
            m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', meshfile))
            if dimension == Caribou._2D:
                mesh = Mesh(m.points[:, 0:2])
            else:
                mesh = Mesh(m.points)
            domain = mesh.add_domain("triangles", Triangle(dimension, order), m.cells_dict[element_type])
            self.assertEqual(domain.number_of_elements(), 32)
            self.assertMatrixEqual(domain.element_indices(0), m.cells_dict[element_type][0])
            self.assertMatrixEqual(domain.element_indices(31), m.cells_dict[element_type][-1])

        m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', '2D_triangle_quadratic.vtk'))
        mesh = Mesh(m.points[:, 0:2])
        self.assertMatrixAlmostEqual(mesh.position(0), m.points[0][0:2])
        self.assertMatrixAlmostEqual(mesh.position(len(m.points)-1), m.points[-1][0:2])
        mesh.add_domain("triangles", Triangle(Caribou._2D, Caribou.Quadratic), m.cells_dict[element_type])
        domain = mesh.domain("triangles")
        numerical_solution = 0.
        for element_id in range(domain.number_of_elements()):
            element = domain.element(element_id)
            for gauss_node in element.gauss_nodes():
                x = gauss_node.position
                w = gauss_node.weight
                J = element.jacobian(x)
                detJ = np.abs(np.linalg.det(J))
                numerical_solution += w * detJ

        analytic_solution = 100
        self.assertAlmostEqual(numerical_solution, analytic_solution, delta=1e-10)


if __name__ == '__main__':
    unittest.main()
