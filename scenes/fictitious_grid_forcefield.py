#!/usr/bin/python3

import Sofa
import SofaCaribou

import meshio
import os

radius = 5
length = 60
n = [10, 10, 40]
subdivisions = 2

mx = (radius / ((n[0])*pow(2, subdivisions)))/2
my = (radius / ((n[1])*pow(2, subdivisions)))/2
mz = (length / ((n[2])*pow(2, subdivisions)))/2

m = meshio.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cylinder_surface.vtu'))


def create_mechanical(node, use_fictitious, color):
    if isinstance(m.cells, dict):
        node.addObject('Mesh', name='surface_mesh', position=m.points.tolist(), triangles=m.cells['triangle'].tolist())
    else:
        for cells in m.cells:
            if cells.type == 'triangle':
                node.addObject('Mesh', name='surface_mesh', position=m.points.tolist(), triangles=cells.data.tolist())
                break

    if use_fictitious:
        grid = node.addObject('FictitiousGrid',
                              name='integration_grid',
                              n=n,
                              min=[-radius-mx, -radius-my, -length/2-mz],
                              max=[+radius+mx, +radius+my, +length/2+mz],
                              maximum_number_of_subdivision_levels=subdivisions,
                              printLog=True)
    else:
        node.addObject('SparseGridTopology',
                              name='integration_grid',
                              n=[(n[0]*pow(2, subdivisions))+1, (n[1]*pow(2, subdivisions))+1, (n[2]*pow(2, subdivisions))+1],
                              min=[-radius-mx, -radius-my, -length/2-mz],
                              max=[+radius+mx, +radius+my, +length/2+mz],
                              src='@surface_mesh',
                              printLog=True)

    meca = node.addChild('meca')

    meca.addObject('StaticODESolver',
                   newton_iterations=20,
                   correction_tolerance_threshold=1e-5,
                   residual_tolerance_threshold=1e-5,
                   printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=1000, residual_tolerance_threshold=1e-5, preconditioning_method='Diagonal')

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
                       src='@../surface_mesh',
                       printLog=True)

        meca.addObject('MechanicalObject', position='@grid.position')
        meca.addObject('HexahedronSetTopologyContainer', name='hexahedrons_container', hexahedra='@grid.hexahedra')

        meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=5000, poisson_ratio=0.49)
        meca.addObject('HyperelasticForcefield', topology='@hexahedrons_container')

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
    meca.addObject('TractionForce', traction=[0, -30, 0], slope=1/5., quads='@quads_container.quads')


    v = meca.addChild('visual')
    v.addObject('OglModel', src='@../../surface_mesh', color=color)
    v.addObject('BarycentricMapping')


def createScene(root):
    root.addObject('APIVersion', level='17.06')

    # create_mechanical(root.addChild('sparse'), False, 'red')
    create_mechanical(root.addChild('fictitious'), True, 'blue')


