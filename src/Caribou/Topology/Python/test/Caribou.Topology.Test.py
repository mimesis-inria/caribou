#! python3

import sys
import unittest
import numpy as np

sys.path.insert(0, "@CARIBOU_PYTHON_LIB_PATH@")

import Caribou
from Caribou.Topology import Mesh


class TestMesh(unittest.TestCase):

    def assertMatrixEqual(self, A, B):
        return self.assertTrue((A == B).all(), f"Matrices are not equal, with \nA={A} \nand\nB={B}")

    def test_constructor_1d(self):
        mesh = Mesh(Caribou._1D)
        mesh.add_node(1)
        mesh.add_nodes([2, 3, 4, 5])
        self.assertMatrixEqual([1., 2., 3., 4., 5.], mesh.positions([0, 1, 2, 3, 4]))

    def test_contructor_3d(self):
        mesh = Mesh(Caribou._3D)


if __name__ == '__main__':
    unittest.main()
