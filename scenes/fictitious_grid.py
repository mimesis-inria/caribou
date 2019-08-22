#!/usr/bin/python3

import Sofa
import SofaCaribou
import numpy as np
from numpy import pi
from math import sqrt
import Tools

radius = 5
length = 60

use_implicit = False

m = Tools.cylinder2(radius=radius, length=length, size=0.5, min_size=None, max_size=None, dimension=2, order=1, traction_box=None, verbose=False)

def is_inside(p):
    x, y, z = p
    if z < -length/2.:
        return 1
    elif z > length/2.:
        return 1
    else:
        return x*x + y*y - radius*radius


def createScene(root):
    n = root.addChild('meca')
    n.addObject('Mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist())
    grid = n.addObject('FictitiousGrid', n=[5, 5, 20], min=[-radius, -radius, -length/2], max=[radius, radius, length/2], use_implicit_surface=use_implicit, surface_positions=m.points.tolist(), )
    grid.set_implicit_test_function(is_inside)
