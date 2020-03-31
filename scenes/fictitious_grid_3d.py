#!/usr/bin/python3

import Sofa
import SofaCaribou
import numpy as np
from numpy import pi
from math import sqrt
import meshio
import os

radius = 5
length = 60
n = [10, 10, 40]
subdivisions = 2

mx = (radius / ((n[0])*pow(2, subdivisions)))/2
my = (radius / ((n[1])*pow(2, subdivisions)))/2
mz = (length / ((n[2])*pow(2, subdivisions)))/2

use_implicit = False

m = meshio.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cylinder_surface.vtu'))

def is_inside(p):
    x, y, z = p
    if z < -length/2.:
        return 1
    elif z > length/2.:
        return 1
    else:
        return x*x + y*y - radius*radius

def createScene(root):
    root.addObject('APIVersion', level='17.06')

    if isinstance(m.cells, dict):
        root.addObject('Mesh', name='surface_mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist())
    else:
        for cells in m.cells:
            if cells.type == 'triangle':
                root.addObject('Mesh', name='surface_mesh', position=m.points.tolist(), triangles=cells.data.tolist())
                break

    grid = root.addObject('FictitiousGrid',
                          template='Vec3',
                          name='integration_grid',
                          n=n,
                          min=[-radius-mx, -radius-my, -length/2-mz],
                          max=[+radius+mx, +radius+my, +length/2+mz],
                          use_implicit_surface=use_implicit,
                          maximum_number_of_subdivision_levels=subdivisions,
                          printLog=True,
                          draw_boundary_cells=True,
                          draw_outside_cells=True,
                          draw_inside_cells=True,)
    if use_implicit:
        grid.set_implicit_test_function(is_inside)
