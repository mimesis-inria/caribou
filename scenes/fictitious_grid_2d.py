#!/usr/bin/python3

import Sofa
import SofaCaribou
from math import sqrt, pi, cos

radius = 5
n = [4, 4]
subdivisions = 10

cell_size = 2*radius/(n[0]-1)
subcell_size = cell_size / pow(2, subdivisions)
eps = subcell_size/2

s = 0
circle_corner = cos(pi/4)*radius
mx = ((circle_corner + 2*circle_corner) - radius)*s + eps


def createScene(root):
    root.bbox = [-radius, -radius, -1, radius, radius, 1]
    root.addObject('APIVersion', level='20.06')
    root.addObject('InteractiveCamera', position=[0, 0, 1], lookAt=[0, 0, 0], projectionType=1, printLog=True)

    root.addObject('CircleIsoSurface', radius=radius, center=[0, 0])

    root.addObject('FictitiousGrid',
                   template='Vec2d',
                   name='integration_grid',
                   n=n,
                   min=[-radius - mx, -radius - mx],
                   max=[+radius + mx, +radius + mx],
                   maximum_number_of_subdivision_levels=subdivisions,
                   printLog=True,
                   draw_boundary_cells=True,
                   draw_outside_cells=True,
                   draw_inside_cells=True, )
