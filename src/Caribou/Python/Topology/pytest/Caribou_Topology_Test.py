#! python3

import sys, os
import unittest
import numpy as np
import meshio
from pathlib import Path

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')
import Caribou
from Caribou.Topology import Mesh
from Caribou.Topology import Grid3D
from Caribou.Geometry import Segment, Segment3
from Caribou.Geometry import Triangle, Triangle6
from Caribou.Geometry import Tetrahedron, Tetrahedron10
from Caribou.Geometry import Hexahedron, Hexahedron20

meshfolder = os.path.join(os.path.dirname(__file__), '..', 'meshes')

class TestMesh(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def test_constructor_1d(self):
        nodes = np.array([1., 2., 3., 4., 5.])
        mesh = Mesh(nodes)
        self.assertMatrixEqual([1., 2., 3., 4., 5.], mesh.positions([0, 1, 2, 3, 4]))
        self.assertEqual(mesh.dimension(), 1)

        mesh = Mesh([1., 2., 3., 4., 5.])
        self.assertMatrixEqual([1., 2., 3., 4., 5.], mesh.positions([0, 1, 2, 3, 4]))
        self.assertEqual(mesh.dimension(), 1)

    def test_contructor_3d(self):
        mesh = Mesh([[1,2,3], [4,5,6]])
        self.assertMatrixEqual([[1,2,3], [4,5,6]], mesh.positions([0, 1]))
        self.assertEqual(mesh.dimension(), 3)


class TestDomain(unittest.TestCase):
    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(A, B, rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_segment_domain(self):
        for element, element_type, meshfile in [(Segment, 'line', '1D_linear.vtk'), (Segment3, 'line3', '1D_quadratic.vtk')]:
            m = meshio.read(os.path.join(meshfolder, meshfile))
            mesh = Mesh(m.points[:, 0])

            if type(m.cells) == dict:
                cells = m.cells[element_type]
            else:
                cells = m.cells_dict[element_type]

            # Creation with name
            domain = mesh.add_domain("segments", element(Caribou._1D), cells)
            self.assertEqual(domain.number_of_elements(), 10)
            self.assertMatrixEqual(domain.element_indices(0), cells[0])
            self.assertMatrixEqual(domain.element_indices(9), cells[9])

            # Creation without name
            domain = mesh.add_domain(element(Caribou._1D), cells)
            self.assertEqual(domain.number_of_elements(), 10)
            self.assertMatrixEqual(domain.element_indices(0), cells[0])
            self.assertMatrixEqual(domain.element_indices(9), cells[9])

    def test_triangle_domain(self):
        for dimension, element, element_type, meshfile in [
            (Caribou._2D, Triangle,  'triangle',  '2D_triangle_linear.vtk'),
            (Caribou._3D, Triangle,  'triangle',  '3D_triangle_linear.vtk'),
            (Caribou._2D, Triangle6, 'triangle6', '2D_triangle_quadratic.vtk'),
            (Caribou._3D, Triangle6, 'triangle6', '3D_triangle_quadratic.vtk'),
        ]:
            m = meshio.read(os.path.join(meshfolder, meshfile))
            if dimension == Caribou._2D:
                mesh = Mesh(m.points[:, 0:2])
            else:
                mesh = Mesh(m.points)

            if type(m.cells) == dict:
                cells = m.cells[element_type]
            else:
                cells = m.cells_dict[element_type]

            domain = mesh.add_domain("triangles", element(dimension), cells)
            self.assertEqual(domain.number_of_elements(), 32)
            self.assertMatrixEqual(domain.element_indices(0), cells[0])
            self.assertMatrixEqual(domain.element_indices(31), cells[-1])

        m = meshio.read(os.path.join(meshfolder, '2D_triangle_quadratic.vtk'))
        mesh = Mesh(m.points[:, 0:2])

        if type(m.cells) == dict:
            cells = m.cells[element_type]
        else:
            cells = m.cells_dict[element_type]

        self.assertMatrixAlmostEqual(mesh.position(0), m.points[0][0:2])
        self.assertMatrixAlmostEqual(mesh.position(len(m.points)-1), m.points[-1][0:2])
        mesh.add_domain("triangles", Triangle6(Caribou._2D), cells)
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
        m = meshio.read(os.path.join(meshfolder, '3D_tetrahedron_quadratic.vtk'))
        mesh = Mesh(m.points)

        if type(m.cells) == dict:
            cells = m.cells['tetra10']
        else:
            cells = m.cells_dict['tetra10']

        domain = mesh.add_domain("triangles", Tetrahedron10(), cells)

        gauss_points_global_coordinates = []
        for element_id in range(domain.number_of_elements()):
            element = domain.element(element_id)
            for gauss_node in element.gauss_nodes():
                gauss_points_global_coordinates.append(element.world_coordinates(gauss_node.position))

        interpolated_positions = domain.embed(gauss_points_global_coordinates).interpolate(m.points)
        self.assertMatrixAlmostEqual(gauss_points_global_coordinates, interpolated_positions, rtol=0, atol=1e-3)

    def test_deformed_liver_tetra(self):
        m = meshio.read(os.path.join(meshfolder, 'deformed_liver_volume_tetrahedrons.vtu'))
        mesh = Mesh(m.points)

        if type(m.cells) == dict:
            cells = m.cells['tetra']
        else:
            cells = m.cells_dict['tetra']

        domain = mesh.add_domain("tetra", Tetrahedron(), cells)

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
        m = meshio.read(os.path.join(meshfolder, 'deformed_liver_volume_hexahedrons.vtu'))
        mesh = Mesh(m.points)

        if type(m.cells) == dict:
            cells = m.cells['hexahedron']
        else:
            cells = m.cells_dict['hexahedron']

        domain = mesh.add_domain("hexa", Hexahedron(), cells)

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


class TestGrid(unittest.TestCase):
    def assertMatrixEqual(self, A, B):
        return self.assertTrue((np.asarray(A) == np.asarray(B)).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def assertMatrixAlmostEqual(self, A, B, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        return self.assertTrue(np.allclose(np.asarray(A), np.asarray(B), rtol, atol, equal_nan), f"Matrices are not almost equal, with \nA={A} \nand\nB={B}")

    def test_grid(self):
        # Grid creation
        g = Grid3D(
            [0.25, 0.5, 0.75], # Anchor point (first node at the corner of the grid)
            [2, 2, 2], # Subdivisions (number of cells in the x, y and z directions respectively)
            [100, 100, 100]) # Dimensions (size of the grid in the x, y and z directions respectively)

        # General properties
        self.assertEqual(g.number_of_nodes(), 27)

        # Cell numbering test
        self.assertEqual(g.cell_index_at([1, 0, 1]), 5)
        self.assertMatrixEqual(g.cell_coordinates_at(5), [1, 0, 1])

        # Cell positioning
        self.assertFalse(g.contains([ 0.24,    0.50,   0.75]))
        self.assertFalse(g.contains([ 0.25,    0.49,   0.75]))
        self.assertFalse(g.contains([ 0.25,    0.50,   0.74]))
        self.assertTrue(g.contains([ 0.25,    0.50,   0.75]))
        self.assertTrue(g.contains([ 50.00,  50.00,  50.00]))
        self.assertTrue(g.contains([100.25, 100.50, 100.75]))
        self.assertFalse(g.contains([100.26, 100.50, 100.75]))
        self.assertFalse(g.contains([100.25, 100.51, 100.75]))

        # Cells around nodes
        self.assertTrue(0 == len(g.cells_around([-50, -50, -50])))
        self.assertMatrixEqual(g.cells_around([0.25, 0.50, 0.75]), [0])
        self.assertMatrixEqual(g.cells_around([50.25, 0.50, 0.75]), [0,1])
        self.assertMatrixEqual(g.cells_around([100.25, 0.50, 0.75]), [1])
        self.assertMatrixEqual(g.cells_around([50.25, 50.50, 50.75]), [0, 4, 2, 6, 1, 5, 3, 7])

        # Cells around faces
        self.assertMatrixEqual(g.cells_around([0.25, 25.50, 25.75]), [0])
        self.assertMatrixEqual(g.cells_around([25.25, 25.50, 0.75]), [0])
        self.assertMatrixEqual(g.cells_around([25.25, 0.50, 25.75]), [0])
        self.assertMatrixEqual(g.cells_around([50.25, 25.50, 25.75]), [0, 1])
        self.assertMatrixEqual(g.cells_around([25.25, 50.50, 25.75]), [0, 2])
        self.assertMatrixEqual(g.cells_around([25.25, 25.50, 50.75]), [0, 4])
        self.assertMatrixEqual(g.cells_around([75.25, 100.50, 75.75]), [7])
        self.assertMatrixEqual(g.cells_around([75.25, 75.50, 100.75]), [7])
        self.assertMatrixEqual(g.cells_around([100.25, 75.50, 75.75]), [7])

        # Node position queries
        self.assertMatrixAlmostEqual(g.node(0), [  0.25,   0.5, 0.75])
        self.assertMatrixAlmostEqual(g.node(2), [100.25,   0.5, 0.75])
        self.assertMatrixAlmostEqual(g.node(6), [  0.25, 100.5, 0.75])
        self.assertMatrixAlmostEqual(g.node(8), [100.25, 100.5, 0.75])

        self.assertMatrixAlmostEqual(g.node(18), [  0.25,   0.5, 100.75])
        self.assertMatrixAlmostEqual(g.node(20), [100.25,   0.5, 100.75])
        self.assertMatrixAlmostEqual(g.node(24), [  0.25, 100.5, 100.75])
        self.assertMatrixAlmostEqual(g.node(26), [100.25, 100.5, 100.75])

        # Node indexing
        for i in range(g.number_of_nodes()):
            self.assertEqual(i, g.node_index_at(g.node_coordinates_at(i)))

        # Edge queries
        # First slice (2D grid)
        self.assertEqual(g.number_of_edges(), 3*12 + 2*9)
        self.assertMatrixEqual(g.edge(0), [0, 1])
        self.assertMatrixEqual(g.edge(0), [0, 1])
        self.assertMatrixEqual(g.edge(1), [1, 2])
        self.assertMatrixEqual(g.edge(2), [0, 3])
        self.assertMatrixEqual(g.edge(3), [1, 4])
        self.assertMatrixEqual(g.edge(4), [2, 5])
        self.assertMatrixEqual(g.edge(5), [3, 4])
        self.assertMatrixEqual(g.edge(6), [4, 5])
        self.assertMatrixEqual(g.edge(7), [3, 6])
        self.assertMatrixEqual(g.edge(8), [4, 7])
        self.assertMatrixEqual(g.edge(9), [5, 8])
        self.assertMatrixEqual(g.edge(10), [6, 7])
        self.assertMatrixEqual(g.edge(11), [7, 8])

        # Between the slice 1 and slice 2
        self.assertMatrixEqual(g.edge(12), [0, 9])
        self.assertMatrixEqual(g.edge(13), [1, 10])
        self.assertMatrixEqual(g.edge(14), [2, 11])
        self.assertMatrixEqual(g.edge(15), [3, 12])
        self.assertMatrixEqual(g.edge(16), [4, 13])
        self.assertMatrixEqual(g.edge(17), [5, 14])
        self.assertMatrixEqual(g.edge(18), [6, 15])
        self.assertMatrixEqual(g.edge(19), [7, 16])
        self.assertMatrixEqual(g.edge(20), [8, 17])

        # Face queries
        self.assertEqual(g.number_of_faces(), 36)
        # First slice (2D grid)
        self.assertMatrixEqual(g.face(0), [0, 1, 4, 3])
        self.assertMatrixEqual(g.face(1), [1, 2, 5, 4])
        self.assertMatrixEqual(g.face(2), [3, 4, 7, 6])
        self.assertMatrixEqual(g.face(3), [4, 5, 8, 7])
        # Between the slice 1 and slice 2
        # - xz axis
        self.assertMatrixEqual(g.face(4), [0, 1, 10,  9])
        self.assertMatrixEqual(g.face(5), [1, 2, 11, 10])
        self.assertMatrixEqual(g.face(6), [3, 4, 13, 12])
        self.assertMatrixEqual(g.face(7), [4, 5, 14, 13])
        self.assertMatrixEqual(g.face(8), [6, 7, 16, 15])
        self.assertMatrixEqual(g.face(9), [7, 8, 17, 16])
        # - yz axis
        self.assertMatrixEqual(g.face(10), [0, 3, 12,  9])
        self.assertMatrixEqual(g.face(11), [1, 4, 13, 10])
        self.assertMatrixEqual(g.face(12), [2, 5, 14, 11])
        self.assertMatrixEqual(g.face(13), [3, 6, 15, 12])
        self.assertMatrixEqual(g.face(14), [4, 7, 16, 13])
        self.assertMatrixEqual(g.face(15), [5, 8, 17, 14])

        # Cell Node indices
        for i in range(g.number_of_cells()):
            node_indices = g.node_indices_of(i)
            nodes = np.array([
                g.node(node_indices[0]), g.node(node_indices[1]), g.node(node_indices[2]), g.node(node_indices[3]),
                g.node(node_indices[4]), g.node(node_indices[5]), g.node(node_indices[6]), g.node(node_indices[7])
            ])
            self.assertMatrixAlmostEqual(nodes, g.cell_at(i).nodes())

        # Cell queries by position
        for i in range(g.number_of_cells()):
            cell = g.cell_at(i)
            for gauss_node in cell.gauss_nodes():
                p = cell.world_coordinates(gauss_node.position)
                self.assertEqual(g.cell_index_containing(p, False), i)


if __name__ == '__main__':
    unittest.main()
