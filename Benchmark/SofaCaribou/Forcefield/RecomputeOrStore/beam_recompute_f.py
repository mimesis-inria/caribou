#!/usr/bin/python3

import re
import SofaRuntime
import Sofa
from SofaRuntime import Timer
import SofaCaribou
import numpy as np

number_of_steps = 10
number_of_newton_iterations = 10
threshold = 1e-15
radius = 5
length = 60
cell_size = 2
use_tetrahedron_mesh = False


class Controller(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.use_sofa_profiler_timer = False
        self.steps = []

    def onAnimateBeginEvent(self, e):
        if len(Timer.getRecords('Animate')):
            self.use_sofa_profiler_timer = True
        else:
            Timer.setEnabled("cg_timer", True)
            Timer.begin("cg_timer")

    def onAnimateEndEvent(self, e):
        if self.use_sofa_profiler_timer:
            records = Timer.getRecords("Animate")
        else:
            records = Timer.getRecords("cg_timer")

        if "AnimateVisitor" in records:
            methods = []
            for k, v in zip(records['AnimateVisitor'].keys(), records['AnimateVisitor'].values()):
                match = re.search('Mechanical \((.*)\)', k)
                if match is not None:
                    methods.append({
                        'name': match.group(1),
                        'newton_steps': extract_newton_steps(v)
                    })
                elif k == 'Mechanical':
                    if isinstance(v, list):
                        for i, v in zip(range(len(v)), v):
                            methods.append({
                                'name': 'method_' + str(i),
                                'newton_steps': extract_newton_steps(v)
                            })
                    else:
                        methods.append({
                            'name': 'method',
                            'newton_steps': extract_newton_steps(v)
                        })
            self.steps.append(methods)
        if not self.use_sofa_profiler_timer:
            Timer.end("cg_timer")


def add_test_case(node, forcefield='HyperelasticForcefieldRecomputeF', tetrahedron=False, cell_size=1.5):
    nx = int(2 * radius / cell_size) + 1
    nz = int(length / cell_size) + 1
    eps = cell_size / 10

    node.addObject('RegularGridTopology', name='grid', min=[-radius, -radius, -length / 2], max=[radius, radius, length / 2], n=[nx, nx, nz])

    node.addObject('StaticODESolver', newton_iterations=number_of_newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=False)
    node.addObject('LLTSolver', backend='Pardiso')
    node.addObject('MechanicalObject', name='mo', position='@grid.position')
    if tetrahedron:
        node.addObject('TetrahedronSetGeometryAlgorithms')
        node.addObject('TetrahedronSetTopologyModifier')
        node.addObject('TetrahedronSetTopologyContainer', name='mechanical_topology')
        node.addObject('Hexa2TetraTopologicalMapping', input='@grid', output='@mechanical_topology')
    else:
        node.addObject('HexahedronSetTopologyContainer', name='mechanical_topology', src='@grid')
    node.addObject('SaintVenantKirchhoffMaterial', young_modulus=3000, poisson_ratio=0)
    node.addObject(forcefield, topology='@mechanical_topology')

    node.addObject('BoxROI', name='base_roi', box=[-radius - eps, -radius - eps, -length / 2 - eps, radius + eps, radius + eps, -length / 2 + eps])
    node.addObject('BoxROI', name='top_roi', box=[-radius - eps, -radius - eps, +length / 2 - eps, radius + eps, radius + eps, +length / 2 + eps], quad='@mechanical_topology.quads')
    node.addObject('FixedConstraint', indices='@base_roi.indices')
    if use_tetrahedron_mesh:
        node.addObject('TriangleSetTopologyContainer', name='neumann_topology', triangles='@top_roi.trianglesInROI')
    else:
        node.addObject('QuadSetTopologyContainer', name='neumann_topology', quads='@top_roi.quadInROI')
    node.addObject('TractionForcefield', topology='@neumann_topology', traction=[0, -30, 0], slope=1 / 5)


def createScene(root):

    root.addObject(Controller(name='timing_controller'))

    root.addObject('RequiredPlugin', pluginName=['SofaCaribou', 'SofaCaribou.Benchmark'])
    root.addObject('RequiredPlugin', pluginName=[
        'SofaBaseMechanics', 'SofaEngine', 'SofaGraphComponent', 'SofaOpenglVisual',
        'SofaPreconditioner', 'SofaBoundaryCondition', 'SofaSparseSolver', 'SofaTopologyMapping'])
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    root.addObject('APIVersion', level='17.06')

    add_test_case(root.addChild('RecomputeF'), 'HyperelasticForcefieldRecomputeF', tetrahedron=use_tetrahedron_mesh, cell_size=cell_size)
    add_test_case(root.addChild('StoreF'), 'HyperelasticForcefieldStoreF', tetrahedron=use_tetrahedron_mesh, cell_size=cell_size)
    add_test_case(root.addChild('StoreF&S'), 'HyperelasticForcefieldStoreFAndS', tetrahedron=use_tetrahedron_mesh, cell_size=cell_size)


def main() :
    import Sofa.Simulation
    import Sofa.Core
    import SofaRuntime

    root = Sofa.Core.Node()
    createScene(root)
    Sofa.Simulation.init(root)
    number_of_nodes = len(next(root.children).mo.position)
    if use_tetrahedron_mesh:
        number_of_elements = len(next(root.children).mechanical_topology.tetrahedra.array())
    else:
        number_of_elements = len(next(root.children).mechanical_topology.hexahedra.array())

    print(f"Number of nodes = {number_of_nodes}")
    print(f"Number of elements = {number_of_elements}")

    print("Computing... (this may take a while)")
    for _ in range(number_of_steps):
        Sofa.Simulation.animate(root, 1)

    # import Sofa.Gui
    # Sofa.Gui.GUIManager.Init("simple_beam", "qglviewer")
    # Sofa.Gui.GUIManager.createGUI(root, __file__)
    # Sofa.Gui.GUIManager.SetDimension(1080, 1080)
    # Sofa.Gui.GUIManager.MainLoop(root)
    # Sofa.Gui.GUIManager.closeGUI()

    pretty_print_methods(root.timing_controller.steps, print_only_averages=True)


def pretty_print_methods(steps, number_format='{:.3f}', print_only_averages=False):
    if len(steps) == 0:
        return

    maximum_number_of_newton_steps = 0
    fields = {}
    for step_id, step in enumerate(steps):
        for method in step:
            method_name = method['name']
            maximum_number_of_newton_steps = max(maximum_number_of_newton_steps, len(method['newton_steps']))
            for newton_step_id, newton_step in enumerate(method['newton_steps']):
                for k, v in zip(newton_step.keys(), newton_step.values()):
                    if k not in fields:
                        fields[k] = {'width': len(k), 'methods':{}}
                    if method_name not in fields[k]['methods']:
                        fields[k]['methods'][method_name] = {'width': len(method_name), 'steps':[]}
                    if step_id == len(fields[k]['methods'][method_name]['steps']):
                        fields[k]['methods'][method_name]['steps'].append([])
                    fields[k]['methods'][method_name]['steps'][step_id].append(v)
                    if isinstance(v, str):
                        fields[k]['methods'][method_name]['width'] = max(fields[k]['methods'][method_name]['width'], len(v))
                    else:
                        fields[k]['methods'][method_name]['width'] = max(fields[k]['methods'][method_name]['width'], len(number_format.format(v)))
    for field in fields.values():
        field['width'] = max(field['width'], len(' '.join(['{{: ^{}}}'.format(m['width']).format('') for m in field['methods'].values()]))+2)

    # Print Header
    ni_col_width = len("Newton it. # ") + len(str(maximum_number_of_newton_steps)) + 1
    cols = ["{{: <{}}}".format(ni_col_width).format('')] + [" {{: ^{}}} ".format(v['width']).format(k) for k,v in zip(fields.keys(), fields.values())]
    fields_header = "|".join(cols)
    cols = ["{{: <{}}}".format(ni_col_width).format('')] + [" {{: ^{}}} ".format(v['width']).format(' '.join([
        "{{: ^{}}}".format(m['width']).format(name) for name, m in zip(v['methods'].keys(), v['methods'].values())
    ])) for v in fields.values()]
    methods_header = "|".join(cols)
    line = f'{{:_^{len(methods_header)}}}'.format('_')
    print(line)
    print(fields_header)
    print(methods_header)
    print(line)

    # Print steps
    if not print_only_averages:
        for step_id in range(len(steps)):
            for newton_step_id in range(maximum_number_of_newton_steps):
                cols = ["Newton it. # " + f"{{: >{len(str(maximum_number_of_newton_steps))}}} ".format(newton_step_id)]
                for f in fields.values():
                    method_values = []
                    for m in f['methods'].values():
                        if step_id >= len(m['steps']) or newton_step_id >= len(m['steps'][step_id]):
                            v = '-'
                        else:
                            v = m['steps'][step_id][newton_step_id]
                            v = '-' if v < 0 else number_format.format(v)
                        method_values.append(f'{{: ^{m["width"]}}}'.format(v))
                    cols += [" {{: ^{}}} ".format(f['width']).format(' '.join(method_values))]
                print("|".join(cols))
            print(line)

    # print averages
    cols = ["{{: <{}}}".format(ni_col_width).format('Mean')]
    for field_name, f in fields.items():
        methods_means = []
        for m in f['methods'].values():
            values = np.asarray([ns for s in m['steps'] for ns in s])
            mean = values[values >= 0].mean()
            methods_means.append(f'{{: ^{m["width"]}}}'.format(number_format.format(mean)))
        cols += [" {{: ^{}}} ".format(f['width']).format(' '.join(methods_means))]
    print("|".join(cols))


def extract_newton_steps(record):
    if 'StaticODESolver::Solve' not in record:
        return []
    newton_steps_records = record['StaticODESolver::Solve']['NewtonStep']
    newton_steps = []
    for newton_record in newton_steps_records:
        if not 'MBKBuild' in newton_record or not 'MBKSolve' in newton_record:
            continue
        MBKBuild = newton_record['MBKBuild']
        MBKSolve = newton_record['MBKSolve']
        MBKFactorize = newton_record['MBKFactorize']
        UpdateForce = newton_record['UpdateForce']

        data = {}
        data['Total time'] = newton_record['total_time']
        data['LHS'] = -1
        data['RHS'] = UpdateForce['total_time']
        data['SOL'] = MBKFactorize['total_time']
        data['ANA'] = -1

        if 'MBKAnalyze' in newton_record:
            data['ANA'] = newton_record['MBKAnalyze']['total_time']

        if 'ConjugateGradient::ComputeGlobalMatrix' in MBKBuild:
            if 'BuildMatrix' in MBKBuild['ConjugateGradient::ComputeGlobalMatrix']:
                data['LHS'] = MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['BuildMatrix']['total_time']
            else:
                data['LHS'] = MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['total_time']
        else:
            data['LHS'] = MBKBuild['total_time']

        if 'PCGLinearSolver::solve' in MBKSolve:
            MBKSolve = MBKSolve['PCGLinearSolver::solve']
        elif 'CG-Solve' in MBKSolve:
            MBKSolve = MBKSolve['CG-Solve']

        if 'ConjugateGradient::solve' in MBKSolve:
            if 'HyperelasticForcefield::update_stiffness' in MBKSolve['ConjugateGradient::solve']:
                update_matrix_time = MBKSolve['ConjugateGradient::solve']['HyperelasticForcefield::update_stiffness']['total_time']
                data['LHS'] += update_matrix_time
        elif 'HyperelasticForcefield::addDForce' in MBKSolve:
            if 'HyperelasticForcefield::update_stiffness' in MBKSolve:
                update_matrix_time = MBKSolve['HyperelasticForcefield::update_stiffness']['total_time']
                data['LHS'] += update_matrix_time

        for k, v in zip(data.keys(), data.values()):
            if isinstance(v, str):
                continue
            if -1e-4 < v < 1e-4:
                data[k] = 0

        newton_steps.append(data)
    return newton_steps


if __name__ == "__main__":
    main()
