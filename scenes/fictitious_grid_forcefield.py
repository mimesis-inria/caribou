#!/usr/bin/python3

import Sofa
import SofaCaribou

import meshio
import os

radius = 5
length = 60
n = [10, 10, 40]
subdivisions = 2

mx = radius / (n[0]-1) / pow(2, subdivisions) / 100
my = radius / (n[1]-1) / pow(2, subdivisions) / 100
mz = length / (n[2]-1) / pow(2, subdivisions) / 100

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
                              name='grid',
                              n=n,
                              min=[-radius-mx, -radius-my, -length/2-mz],
                              max=[+radius+mx, +radius+my, +length/2+mz],
                              maximum_number_of_subdivision_levels=subdivisions,
                              printLog=True)
    else:
        node.addObject('SparseGridTopology',
                              name='grid',
                              n=n,
                              min=[-radius-mx, -radius-my, -length/2-mz],
                              max=[+radius+mx, +radius+my, +length/2+mz],
                              src='@surface_mesh',
                              printLog=True)

    meca = node.addChild('meca')

    meca.addObject('StaticODESolver',
                   newton_iterations=20,
                   correction_tolerance_threshold=1e-3,
                   residual_tolerance_threshold=1e-3,
                   printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=100, residual_tolerance_threshold=1e-3, preconditioning_method='IncompleteCholesky', printLog=False)

    meca.addObject('MechanicalObject', position='@../grid.position')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=5000, poisson_ratio=0.49)
    meca.addObject('HexahedronSetTopologyContainer', name='hexahedrons_container', hexahedra='@../grid.hexahedra')

    if use_fictitious:
        meca.addObject('FictitiousGridHyperelasticForcefield',
                       fictitious_grid='@../grid',
                       integration_method='SubdividedVolume',
                       printLog=True)
    else:
        meca.addObject('HyperelasticForcefield', topology='@hexahedrons_container')

    meca.addObject('BoxROI', name='left_roi',
                   box=[-radius - 1.5*mx, -radius - 1.5*my, -length / 2.0 - mz - mz/2,
                        +radius + 1.5*mx, +radius + 1.5*mx, -length / 2.0 - mz + mz/2],
                   drawBoxes=True)

    meca.addObject('FixedConstraint', indices='@left_roi.indices')

    t = meca.addChild('traction')
    t.addObject('MechanicalObject', position='@../../surface_mesh.position')
    t.addObject('TriangleSetTopologyContainer', triangles='@../../surface_mesh.triangles')
    t.addObject('BoxROI', name='right_roi', strict=True,
                box=[-radius - mx, -radius - my, +length / 2.0 - mz,
                     +radius + mx, +radius + mx, +length / 2.0 + mz],
                drawBoxes=True)
    t.addObject('TractionForce', traction=[0, -30, 0], slope=1/10, triangles='@right_roi.trianglesInROI', printLog=True)
    t.addObject('BarycentricMapping')

    v = meca.addChild('visual')
    v.addObject('OglModel', src='@../../surface_mesh', color=color)
    v.addObject('BarycentricMapping')


def createScene(root):
    root.addObject('APIVersion', level='17.06')

    create_mechanical(root.addChild('sparse'), False, 'red')
    create_mechanical(root.addChild('fictitious'), True, 'blue')


