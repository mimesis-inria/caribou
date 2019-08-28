#!/usr/bin/python3

import Sofa
import SofaCaribou
import numpy as np
from numpy import pi
from math import sqrt
import Tools

radius = 5
length = 60
n = [10, 10, 40]
subdivisions = 1

mx = (radius / ((n[0])*pow(2, subdivisions)))/2
my = (radius / ((n[1])*pow(2, subdivisions)))/2
mz = (length / ((n[2])*pow(2, subdivisions)))/2

use_implicit = False

m = Tools.cylinder(center1=[0,0,-length/2.], center2=[0,0,length/2.], radius=radius, number_of_cuts=25, dimension=2)


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
    root.addObject('Mesh', name='mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist())
    grid = root.addObject('FictitiousGrid',
                          name='grid',
                          n=n,
                          min=[-radius-mx, -radius-my, -length/2-mz],
                          max=[+radius+mx, +radius+my, +length/2+mz],
                          use_implicit_surface=use_implicit,
                          maximum_number_of_subdivision_levels=subdivisions,
                          printLog=True)
    if use_implicit:
        grid.set_implicit_test_function(is_inside)

    node = root.addChild('meca')
    node.addObject('StaticODESolver',
                   newton_iterations=20,
                   correction_tolerance_threshold=1e-5,
                   residual_tolerance_threshold=1e-5,
                   printLog=True)
    node.addObject('SparseLUSolver')

    node.addObject('MechanicalObject', position='@../grid.positions')
    node.addObject('HexahedronSetTopologyContainer', hexahedra='@../grid.hexahedrons')
    node.addObject('HexahedronElasticForce',
                   youngModulus=5000,
                   poissonRatio=0.49,
                   corotated=False,
                   linearStrain=False,
                   printLog=True)

    node.addObject('BoxROI', name='left_roi',
                   box=[-radius - 1.5*mx, -radius - 1.5*my, -length / 2.0 - 1.5*mz,
                        +radius + 1.5*mx, +radius + 1.5*mx, -length / 2.0 + 1.5*mz],
                   drawBoxes=True)

    node.addObject('FixedConstraint', indices='@left_roi.indices')


    v = node.addChild('visual')
    v.addObject('OglModel', src='@../../mesh', color='red')
    v.addObject('MechanicalObject', src='@../../mesh')
    v.addObject('TriangleSetTopologyContainer', src='@../../mesh')
    v.addObject('BoxROI', box=[4, 0, 0, 5, 1, 1], drawBoxes=True)

