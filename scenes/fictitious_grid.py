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
subdivisions = 2

mx = (radius / ((n[0])*pow(2, subdivisions)))/2
my = (radius / ((n[1])*pow(2, subdivisions)))/2
mz = (length / ((n[2])*pow(2, subdivisions)))/2

use_implicit = False

m = Tools.cylinder(center1=[0,0,-length/2.], center2=[0,0,length/2.], radius=radius, number_of_cuts=50, dimension=2)


def is_inside(p):
    x, y, z = p
    if z < -length/2.:
        return 1
    elif z > length/2.:
        return 1
    else:
        return x*x + y*y - radius*radius


def create_mechanical(node, use_fictitious, color):
    node.addObject('Mesh', name='mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist())
    if use_fictitious:
        grid = node.addObject('FictitiousGrid',
                              name='integration_grid',
                              n=n,
                              min=[-radius-mx, -radius-my, -length/2-mz],
                              max=[+radius+mx, +radius+my, +length/2+mz],
                              use_implicit_surface=use_implicit,
                              maximum_number_of_subdivision_levels=subdivisions,
                              printLog=True)
        if use_implicit:
            grid.set_implicit_test_function(is_inside)
    else:
        node.addObject('SparseGridTopology',
                              name='integration_grid',
                              n=[(n[0]*pow(2, subdivisions))+1, (n[1]*pow(2, subdivisions))+1, (n[2]*pow(2, subdivisions))+1],
                              min=[-radius-mx, -radius-my, -length/2-mz],
                              max=[+radius+mx, +radius+my, +length/2+mz],
                              src='@mesh',
                              printLog=True)

    meca = node.addChild('meca')

    meca.addObject('StaticODESolver',
                   newton_iterations=20,
                   correction_tolerance_threshold=1e-5,
                   residual_tolerance_threshold=1e-5,
                   printLog=True)
    meca.addObject('CGLinearSolver', iterations=1000)

    if use_fictitious:
        meca.addObject('MechanicalObject', position='@../integration_grid.positions')
        meca.addObject('HexahedronSetTopologyContainer', name='hexahedrons_contaier', hexahedra='@../integration_grid.hexahedrons')
        meca.addObject('FictitiousGridElasticForce',
                       youngModulus=5000,
                       poissonRatio=0.49,
                       corotated=False,
                       linearStrain=False,
                       fictitious_grid='@../integration_grid',
                       integration_method='SubdividedGauss',
                       printLog=True)
    else:
        meca.addObject('SparseGridTopology',
                       name='grid',
                       n=[n[0]+1, n[1]+1, n[2]+1],
                       min=[-radius-mx, -radius-my, -length/2-mz],
                       max=[+radius+mx, +radius+my, +length/2+mz],
                       src='@../mesh',
                       printLog=True)

        meca.addObject('MechanicalObject', position='@grid.position')
        meca.addObject('HexahedronSetTopologyContainer', name='hexahedrons_contaier', hexahedra='@grid.hexahedra')

        meca.addObject('HexahedronElasticForce',
                       youngModulus=5000,
                       poissonRatio=0.49,
                       corotated=False,
                       linearStrain=False,
                       topology_container='@hexahedrons_contaier',
                       integration_method='SubdividedGauss',
                       integration_grid='@../integration_grid',
                       number_of_subdivisions=subdivisions,
                       printLog=True)

    meca.addObject('BoxROI', name='left_roi',
                   box=[-radius - 1.5*mx, -radius - 1.5*my, -length / 2.0 - mz - mz/2,
                        +radius + 1.5*mx, +radius + 1.5*mx, -length / 2.0 - mz + mz/2],
                   drawBoxes=True)

    meca.addObject('FixedConstraint', indices='@left_roi.indices')

    meca.addObject('BoxROI', name='right_roi',
                   strict=True,
                   box=[-radius - 1.5*mx, -radius - 1.5*my, +length / 2.0 + mz - mz/2,
                        +radius + 1.5*mx, +radius + 1.5*mx, +length / 2.0 + mz + mz/2],
                   drawBoxes=True)
    meca.addObject('QuadSetTopologyContainer', name='quads_container', quads='@right_roi.quadInROI')
    meca.addObject('TriangleSetTopologyContainer', name='triangles_container')
    meca.addObject('TriangleSetTopologyModifier')
    meca.addObject('Quad2TriangleTopologicalMapping', input='@quads_container', output='@triangles_container')
    meca.addObject('TractionForce', traction=[0, -30, 0], slope=1/5., triangles='@triangles_container.triangles')


    v = meca.addChild('visual')
    v.addObject('OglModel', src='@../../mesh', color=color)
    v.addObject('BarycentricMapping')


def createScene(root):
    root.addObject('APIVersion', level='17.06')

    create_mechanical(root.addChild('sparse'), False, 'red')
    create_mechanical(root.addChild('fictitious'), True, 'blue')


