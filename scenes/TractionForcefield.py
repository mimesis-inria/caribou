#!/usr/bin/python3

import Sofa
import SofaCaribou
import numpy as np
import meshio
import os


m = meshio.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rectangular_beam.vtu'))

base_id = 1
top_id = 2

base_mask = np.array(np.ma.masked_not_equal(m.cell_data['quad']['gmsh:physical'], base_id).mask)
base_faces = m.cells['quad'][~base_mask]
base_nodes = np.unique(base_faces).astype(int)

top_mask = np.array(np.ma.masked_not_equal(m.cell_data['quad']['gmsh:physical'], top_id).mask)
top_faces = m.cells['quad'][~top_mask]


def createScene(root):
    root.addObject('APIVersion', level='17.06')

    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    meca = root.addChild("meca")
    meca.addObject('StaticODESolver', newton_iterations=2, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=100, residual_tolerance_threshold=1e-5, preconditioning_method="IncompleteCholesky", printLog=True)
    meca.addObject('MechanicalObject', position=m.points.tolist())
    meca.addObject('HexahedronSetTopologyContainer', name='topo', hexahedra=m.cells['hexahedron'].tolist())
    meca.addObject('HexahedronElasticForce', topology_container='@topo', youngModulus=3000, poissonRatio=0, corotated=False, linearStrain=False)

    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads=top_faces.tolist())
    meca.addObject('TractionForce', traction=[0, -30, 0], slope=1/5, quads='@quad_container.quads')

    meca.addObject('FixedConstraint', indices=base_nodes.tolist())


