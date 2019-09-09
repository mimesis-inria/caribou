#!/usr/bin/python3

import Sofa
import SofaCaribou
import numpy as np
from numpy import pi
from math import sqrt


radius = 5
n = [5, 5]
subdivisions = 10
mx = (radius / ((n[0])*pow(2, subdivisions)))/2
my = (radius / ((n[1])*pow(2, subdivisions)))/2


def is_inside(p):
    x, y = p
    return x*x + y*y - radius*radius


def createScene(root):
    root.addObject('APIVersion', level='17.06')

    grid = root.addObject('FictitiousGrid',
                          template='Vec2d',
                          name='integration_grid',
                          n=n,
                          min=[-radius-mx, -radius-my],
                          max=[+radius+mx, +radius+my],
                          use_implicit_surface=True,
                          maximum_number_of_subdivision_levels=subdivisions,
                          printLog=True,
                          draw_boundary_cells=True,
                          draw_outside_cells=True,
                          draw_inside_cells=True,)
    grid.set_implicit_test_function(is_inside)
