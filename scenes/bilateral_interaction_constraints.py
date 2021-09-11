#!/usr/bin/python3

import Sofa


def createScene(root):
    root.dt = 1

    root.addObject('RequiredPlugin', pluginName=[
        'SofaSparseSolver', 'SofaTopologyMapping', 'SofaBoundaryCondition', 'SofaEngine',
        'SofaImplicitOdeSolver', 'SofaConstraint', 'SofaSimpleFem', 'SofaGeneralLinearSolver'
    ])
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showForceFields')

    root.addObject('FreeMotionAnimationLoop')
    root.addObject('GenericConstraintSolver', tolerance=1e-3, maxIterations=1000)

    root.addObject('RegularGridTopology', name='grid1', min=[0, -7.5, -7.5],  max=[40, 7.5, 7.5], n=[9, 3, 3])
    root.addObject('RegularGridTopology', name='grid2', min=[40, -7.5, -7.5], max=[80, 7.5, 7.5], n=[9, 3, 3])

    ##########
    # BEAM 1 #
    ##########

    beam_1 = root.addChild("beam_1")
    beam_1.addObject('EulerImplicitSolver')
    beam_1.addObject('LDLTSolver', backend="Eigen")
    # beam_1.addObject('BackwardEulerODESolver', newton_iterations=10, rayleigh_stiffness=0, rayleigh_mass=0, residual_tolerance_threshold=1e-5, pattern_analysis_strategy="ALWAYS", printLog=True)
    # beam_1.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    beam_1.addObject('MechanicalObject', name='mo', src='@../grid1')

    # Complete hexa container
    beam_1.addObject('HexahedronSetTopologyContainer', src='@../grid1', name='mechanical_topology')

    # - Mass
    beam_1.addObject('HexahedronSetGeometryAlgorithms')
    beam_1.addObject('DiagonalMass', massDensity=0.5)

    # - Mechanics
    beam_1.addObject('SaintVenantKirchhoffMaterial', young_modulus=7500, poisson_ratio=0.3)
    beam_1.addObject('HyperelasticForcefield')

    # Fix the left side of the beam
    beam_1.addObject('BoxROI', name='fixed_roi', quad='@surface_topology.quad', box=[-0.9, -7.5, -7.5, 0.1, 7.5, 7.5])
    beam_1.addObject('FixedConstraint', indices='@fixed_roi.indices')

    # Get the indices of the right side of the beam
    beam_1.addObject('BoxROI', name='bilateral_roi', box=[40-0.9, -7.5, -7.5, 40+0.1, 7.5, 7.5])

    # Constraints
    beam_1.addObject('LinearSolverConstraintCorrection')

    ##########
    # BEAM 2 #
    ##########

    beam_2 = root.addChild("beam_2")
    beam_2.addObject('EulerImplicitSolver')
    beam_2.addObject('LDLTSolver', backend="Eigen")
    # beam_2.addObject('BackwardEulerODESolver', newton_iterations=10, rayleigh_stiffness=0, rayleigh_mass=0, residual_tolerance_threshold=1e-5, pattern_analysis_strategy="ALWAYS", printLog=True)
    # beam_2.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    beam_2.addObject('MechanicalObject', name='mo', src='@../grid2')

    # Complete hexa container
    beam_2.addObject('HexahedronSetTopologyContainer', src='@../grid2', name='mechanical_topology')

    # - Mass
    beam_2.addObject('HexahedronSetGeometryAlgorithms')
    beam_2.addObject('DiagonalMass', massDensity=0.5)

    # - Mechanics
    beam_2.addObject('SaintVenantKirchhoffMaterial', young_modulus=7500, poisson_ratio=0.3)
    beam_2.addObject('HyperelasticForcefield')

    # Fix the right side of the beam
    beam_2.addObject('BoxROI', name='fixed_roi', quad='@surface_topology.quad', box=[80-0.9, -7.5, -7.5, 80+0.1, 7.5, 7.5])
    beam_2.addObject('FixedConstraint', indices='@fixed_roi.indices')

    # Get the indices of the left side of the beam
    beam_2.addObject('BoxROI', name='bilateral_roi', box=[40-0.9, -7.5, -7.5, 40+0.1, 7.5, 7.5])

    # Constraints
    beam_2.addObject('LinearSolverConstraintCorrection')

    # Add constraints at the interface between the two beams
    root.addObject('BilateralInteractionConstraint', object1='@beam_1/mo', object2='@beam_2/mo', first_point='@beam_1/bilateral_roi.indices', second_point='@beam_2/bilateral_roi.indices')
