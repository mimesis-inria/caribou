#!/usr/bin/python3

import Sofa
import SofaCaribou
import numpy as np
from numpy import pi
from math import sqrt
import Tools

radius = 5
length = 60
n = [11, 11, 41]
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
    node = root.addChild('meca')
    node.addObject('Mesh', name='mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist())
    grid = node.addObject('FictitiousGrid',
                          n=n,
                          min=[-radius-mx, -radius-my, -length/2-mz],
                          max=[+radius+mx, +radius+my, +length/2+mz],
                          # n=[1, 2, 2],
                          # min=[+radius+mx - dx, -dx/2, -dz/2],
                          # max=[+radius+mx, -dx/2 + 2*dx, -dz/2 + 2*dz],
                          use_implicit_surface=use_implicit,
                          maximum_number_of_subdivision_levels=subdivisions,
                          printLog=True)
    if use_implicit:
        grid.set_implicit_test_function(is_inside)

    v = node.addChild('visual')
    v.addObject('OglModel', src='@../mesh', color='red')
    v.addObject('MechanicalObject', src='@../mesh')
    v.addObject('TriangleSetTopologyContainer', src='@../mesh')
    v.addObject('BoxROI', box=[4, 0, 0, 5, 1, 1], drawBoxes=True)

