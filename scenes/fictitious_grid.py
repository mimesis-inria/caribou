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

m = Tools.cylinder2(radius=radius, length=length, size=2, min_size=None, max_size=None, dimension=2, order=1, traction_box=None, verbose=False)


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
    n.addObject('Mesh', name='mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist()[50:])
    grid = n.addObject('FictitiousGrid', n=[11, 11, 41], min=[-radius, -radius, -length/2], max=[radius, radius, length/2], use_implicit_surface=use_implicit, printLog=True)
    # grid.set_implicit_test_function(is_inside)

    v = n.addChild('visual')
    v.addObject('OglModel', src='@../mesh')
