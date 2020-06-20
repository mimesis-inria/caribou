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

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(A, B, rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_segment_domain(self):
        for order, element_type, meshfile in [(Caribou.Linear, 'line', '1D_linear.vtk'), (Caribou.Quadratic, 'line3', '1D_quadratic.vtk')]:
            m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', meshfile))
            mesh = Mesh(m.points[:, 0])

            if type(m.cells) == dict:
                cells = m.cells[element_type]
            else:
                cells = m.cells_dict[element_type]

            domain = mesh.add_domain("segments", Segment(Caribou._1D, order), cells)
            self.assertEqual(domain.number_of_elements(), 10)
            self.assertMatrixEqual(domain.element_indices(0), cells[0])
            self.assertMatrixEqual(domain.element_indices(9), cells[9])

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

            if type(m.cells) == dict:
                cells = m.cells[element_type]
            else:
                cells = m.cells_dict[element_type]

            domain = mesh.add_domain("triangles", Triangle(dimension, order), cells)
            self.assertEqual(domain.number_of_elements(), 32)
            self.assertMatrixEqual(domain.element_indices(0), cells[0])
            self.assertMatrixEqual(domain.element_indices(31), cells[-1])

        m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', '2D_triangle_quadratic.vtk'))
        mesh = Mesh(m.points[:, 0:2])

        if type(m.cells) == dict:
            cells = m.cells[element_type]
        else:
            cells = m.cells_dict[element_type]

        self.assertMatrixAlmostEqual(mesh.position(0), m.points[0][0:2])
        self.assertMatrixAlmostEqual(mesh.position(len(m.points)-1), m.points[-1][0:2])
        mesh.add_domain("triangles", Triangle(Caribou._2D, Caribou.Quadratic), cells)
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

    def test_barycentric_mapping(self):
        m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', '3D_tetrahedron_quadratic.vtk'))
        mesh = Mesh(m.points)

        if type(m.cells) == dict:
            cells = m.cells['tetra10']
        else:
            cells = m.cells_dict['tetra10']

        domain = mesh.add_domain("triangles", Tetrahedron(Caribou.Quadratic), cells)

        gauss_points_global_coordinates = []
        for element_id in range(domain.number_of_elements()):
            element = domain.element(element_id)
            for gauss_node in element.gauss_nodes():
                gauss_points_global_coordinates.append(element.world_coordinates(gauss_node.position))

        interpolated_positions = domain.embed(gauss_points_global_coordinates).interpolate(m.points)
        self.assertMatrixAlmostEqual(gauss_points_global_coordinates, interpolated_positions, rtol=0, atol=1e-3)

    def test_deformed_liver_tetra(self):
        m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', 'deformed_liver_volume_tetrahedrons.vtu'))
        mesh = Mesh(m.points)

        if type(m.cells) == dict:
            cells = m.cells['tetra']
        else:
            cells = m.cells_dict['tetra']

        domain = mesh.add_domain("tetra", Tetrahedron(Caribou.Linear), cells)

        undeformed_markers = np.array([[265.83,  -155.11, -1539.04],
                                       [260.11,  -126.47, -1540.98],
                                       [274.25,  -176.07, -1505.73],
                                       [285.61,  -211.77, -1503.08],
                                       [260.49,  -124.21, -1514.15],
                                       [264.56,  -143.91, -1478.4 ],
                                       [231.42,  -167.12, -1478.51],
                                       [204.36,  -193.61, -1501.44],
                                       [210.5 ,  -149.3 , -1519.34],
                                       [232.38,  -248.64, -1507.46],
                                       [271.51,  -192.76, -1462.05],
                                       [294.9 ,  -162.29, -1478.03],
                                       [312.83,  -191.94, -1482.31],
                                       [302.98,  -269.28, -1493.  ],
                                       [315.78,  -226.09, -1504.54]])

        deformed_markers = np.array([[270.19066814,  -161.84690798, -1552.47975919],
                                     [263.68932938,  -129.28575865, -1551.84152932],
                                     [281.32091292,  -178.37405176, -1513.74672457],
                                     [294.15755424,  -213.88345885, -1521.76090069],
                                     [271.89878464,  -125.95830477, -1525.02450386],
                                     [279.99528928,  -146.00187886, -1492.49149538],
                                     [258.70099261,  -171.47231679, -1504.13209461],
                                     [242.71558686,  -196.94053417, -1534.58006884],
                                     [243.3667501 ,  -155.00591598, -1550.16871752],
                                     [251.82411482,  -253.24382597, -1533.04677523],
                                     [288.71467753,  -196.48153973, -1479.45513226],
                                     [305.6466301 ,  -166.99511463, -1488.22911258],
                                     [320.86269498,  -198.04169886, -1497.0985972 ],
                                     [304.89113271,  -273.96234025, -1512.40119427],
                                     [320.51991513,  -232.28761329, -1523.22695582]])

        interpolated_displacements = domain.embed(undeformed_markers).interpolate(m.point_data['u'])
        self.assertMatrixAlmostEqual(interpolated_displacements, (deformed_markers - undeformed_markers))

    def test_deformed_liver_hexa(self):
        m = meshio.read(os.path.join(os.path.dirname(__file__), 'meshes', 'deformed_liver_volume_hexahedrons.vtu'))
        mesh = Mesh(m.points)

        if type(m.cells) == dict:
            cells = m.cells['hexahedron']
        else:
            cells = m.cells_dict['hexahedron']

        domain = mesh.add_domain("hexa", Hexahedron(Caribou.Linear), cells)

        undeformed_markers = np.array([[265.83,  -155.11, -1539.04],
                                       [260.11,  -126.47, -1540.98],
                                       [274.25,  -176.07, -1505.73],
                                       [285.61,  -211.77, -1503.08],
                                       [260.49,  -124.21, -1514.15],
                                       [264.56,  -143.91, -1478.4 ],
                                       [231.42,  -167.12, -1478.51],
                                       [204.36,  -193.61, -1501.44],
                                       [210.5 ,  -149.3 , -1519.34],
                                       [232.38,  -248.64, -1507.46],
                                       [271.51,  -192.76, -1462.05],
                                       [294.9 ,  -162.29, -1478.03],
                                       [312.83,  -191.94, -1482.31],
                                       [302.98,  -269.28, -1493.  ],
                                       [315.78,  -226.09, -1504.54]])

        deformed_markers = np.array([[270.9119593 , -157.71215294, -1555.40946454],
                                     [262.14820601, -128.06090875, -1552.47007965],
                                     [277.19378591, -177.43088076, -1516.55486109],
                                     [293.97624055, -213.41924186, -1522.43362963],
                                     [265.88433766, -126.62070365, -1524.0396069 ],
                                     [266.88009204, -147.60863833, -1490.24440459],
                                     [262.09361709, -173.36263398, -1499.35976378],
                                     [243.77807901, -199.20167085, -1529.01202647],
                                     [246.71909105, -156.88296646, -1546.59187345],
                                     [252.34989369, -255.47322977, -1531.41912502],
                                     [275.70481717, -198.33630816, -1478.29025749],
                                     [296.66730156, -168.05975732, -1487.89185992],
                                     [316.39825186, -196.90427925, -1493.97027093],
                                     [302.8963988 , -274.06481598, -1510.41305086],
                                     [319.35443171, -231.06013115, -1521.27203137]])

        interpolated_displacements = domain.embed(undeformed_markers).interpolate(m.point_data['u'])
        self.assertMatrixAlmostEqual(interpolated_displacements, (deformed_markers - undeformed_markers))


if __name__ == '__main__':
    unittest.main()
