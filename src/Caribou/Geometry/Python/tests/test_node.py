import pytest

import sys
sys.path.append("@PACKAGE_INSTALL_PREFIX@")
print sys.path

from Caribou.Geometry import Point1D, Point2D, Point3D

def test_1D():
    p1 = Point1D()
    p1.x = 1
    assert p1.x == 1

    p1 = Point1D(1)
    assert p1.x == 1
    p1.x = 0

    p2 = p1
    assert p2.x == 0

    p2.x = 1
    assert p1.x == 1
    assert p2.x == 1

def test_2D():
    p1 = Point2D()
    p1.x = 1
    p1.y = 2
    assert p1.x == 1
    assert p1.y == 2


def test_3D():
    p1 = Point3D()
    p1.x = 1
    p1.y = 2
    p1.z = 3
    assert p1.x == 1
    assert p1.y == 2
    assert p1.z == 3