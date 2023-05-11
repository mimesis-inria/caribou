import Sofa

number_of_steps = 10
number_of_newton_iterations = 10
threshold = 1e-15
radius = 5
length = 60
cell_size = 3
nx = int(2 * radius / cell_size) + 1
nz = int(length / cell_size) + 1
eps = cell_size / 10


def createScene(root):
    root.dt = 0.001
    root.addObject('RequiredPlugin', pluginName=[
        'Sofa.Component.SceneUtility', # APIVersion
        'Sofa.Component.Constraint.Projective', # FixedConstraint
        'Sofa.Component.Engine.Select', # BoxROI
        'Sofa.Component.ODESolver.Forward', # CentralDifferenceSolver
        'Sofa.Component.StateContainer', # MechanicalObject
        'Sofa.Component.Topology.Container.Grid', # RegularGridTopology
        'Sofa.Component.Visual' # VisualStyle
    ])
    root.addObject('APIVersion', level='23.06.99')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showForceFields')
    root.addObject('RegularGridTopology', name='grid', min=[-length/2, -radius, -radius], max=[length/2, radius, radius], n=[nz, nx, nx])
    root.addObject('DefaultAnimationLoop')

    # Caribou lumped mass
    root.addChild('caribou_lumped')
    root.caribou_lumped.addObject('CentralDifferenceSolver')
    root.caribou_lumped.addObject('MechanicalObject', name='mo', src='@../grid')
    root.caribou_lumped.addObject('CaribouTopology', name='topology', indices='@../grid.hexahedra', template='Hexahedron')
    root.caribou_lumped.addObject('CaribouMass', name='mass', topology='@topology', density=2, lumped=True)
    root.caribou_lumped.addObject('SaintVenantKirchhoffMaterial', young_modulus=5000, poisson_ratio=0)
    root.caribou_lumped.addObject('HyperelasticForcefield', name='ff', topology='@topology')
    root.caribou_lumped.addObject('BoxROI', name='base_roi', box=[-length/2 - eps, -radius - eps, -radius - eps, -length/2+eps, radius + eps, radius + eps])
    root.caribou_lumped.addObject('FixedConstraint', indices='@base_roi.indices')
