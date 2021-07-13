#!/usr/bin/python3

import re
import SofaRuntime
import Sofa
from SofaRuntime import Timer
import SofaCaribou

number_of_newton_iterations = 10
threshold = 1e-15
cell_size = 1.35
radius = 5
length = 60

nx = int(2*radius / cell_size)+1
nz = int(length / cell_size) + 1
eps = cell_size/10

direct_solvers = [
    # LLT
    {'name': 'LLTEigen', 'solver': 'LLTSolver', 'arguments': {'backend':'Eigen'}},
    {'name': 'LLTPardi', 'solver': 'LLTSolver', 'arguments': {'backend':'Pardiso'}},
    {'name': 'LLTSofa',  'solver': 'SparseCholeskySolver', 'arguments': {}},

    # # LDLT
    # {'name': 'LDLTEigen', 'solver': 'LDLTSolver', 'arguments': {'backend':'Eigen'}},
    # {'name': 'LDLTPardi', 'solver': 'LDLTSolver', 'arguments': {'backend':'Pardiso'}},
    # {'name': 'LDLTSofa',  'solver': 'SparseLDLSolver', 'arguments': {}},

    # # LU
    # {'name': 'LUEigen', 'solver': 'LUSolver', 'arguments': {'backend':'Eigen'}},
    # {'name': 'LUPardi', 'solver': 'LUSolver', 'arguments': {'backend':'Pardiso'}},
    # {'name': 'LUSofa',  'solver': 'SparseLUSolver', 'arguments': {}},
]

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

        data = {}
        data['Total'] = newton_record['total_time']
        data['Build'] = MBKBuild['total_time']
        data['Solve'] = MBKSolve['total_time']

        for k, v in zip(data.keys(), data.values()):
            if isinstance(v, str):
                continue
            if v < 1e-4:
                data[k] = '0'

        newton_steps.append(data)
    return newton_steps

def pretty_print_methods(methods, number_format='{:.3f}'):
    if len(methods) == 0:
        return
    maximum_number_of_newton_steps = 0
    fields = {}
    for method in methods:
        method_name = method['name']
        maximum_number_of_newton_steps = max(maximum_number_of_newton_steps, len(method['newton_steps']))
        for newton_step in method['newton_steps']:
            for k, v in zip(newton_step.keys(), newton_step.values()):
                if k not in fields:
                    fields[k] = {'width': len(k), 'methods':{}}
                if method_name not in fields[k]['methods']:
                    fields[k]['methods'][method_name] = {'width': len(method_name), 'values':[]}
                fields[k]['methods'][method_name]['values'].append(v)
                if isinstance(v, str):
                    fields[k]['methods'][method_name]['width'] = max(fields[k]['methods'][method_name]['width'], len(v))
                else:
                    fields[k]['methods'][method_name]['width'] = max(fields[k]['methods'][method_name]['width'], len(number_format.format(v)))
    for field in fields.values():
        field['width'] = max(field['width'], len(' '.join(['{{: ^{}}}'.format(m['width']).format('') for m in field['methods'].values()]))+2)

    # Print Header
    ni_col_width = len("Newton it. # ") + len(str(maximum_number_of_newton_steps)) + 1
    cols = ["{{: <{}}}".format(ni_col_width).format('')] + [" {{: ^{}}} ".format(v['width']).format(k) for k,v in zip(fields.keys(), fields.values())]
    print("|".join(cols))
    cols = ["{{: <{}}}".format(ni_col_width).format('')] + [" {{: ^{}}} ".format(v['width']).format(' '.join([
        "{{: ^{}}}".format(m['width']).format(name) for name, m in zip(v['methods'].keys(), v['methods'].values())
    ])) for v in fields.values()]
    print("|".join(cols))

    # Print newton iterations
    for it in range(maximum_number_of_newton_steps):
        cols = ["Newton it. # " + "{{: >{}}} ".format(len(str(maximum_number_of_newton_steps))).format(it)] + \
               [" {{: ^{}}} ".format(f['width']).format(' '.join(["{{: ^{}}}".format(m['width']).format('-' if it>=len(m['values']) else (m['values'][it] if isinstance(m['values'][it], str) else number_format.format(m['values'][it]))) for m in f['methods'].values()])) for f in fields.values()]
        print("|".join(cols))


class Controller(Sofa.Core.Controller):
    def __init__(self):
        super().__init__(self)
        self.use_sofa_profiler_timer = False

    def onAnimateBeginEvent(self, e):
        if len(Timer.getRecords('Animate')):
            self.use_sofa_profiler_timer = True
        else:
            Timer.setEnabled("timer", True)
            Timer.begin("timer")

    def onAnimateEndEvent(self, e):
        print("Done")
        if self.use_sofa_profiler_timer:
            records = Timer.getRecords("Animate")
        else:
            records = Timer.getRecords("timer")

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
                                'name': 'method_'+str(i),
                                'newton_steps': extract_newton_steps(v)
                            })
                    else:
                        methods.append({
                            'name': 'method',
                            'newton_steps': extract_newton_steps(v)
                        })
            print("Here are the results. Copy and paste them in a text editor without word wrap to visualize them.")
            pretty_print_methods(methods)
        if not self.use_sofa_profiler_timer:
            Timer.end("timer")


def createScene(root):
    root.addObject(Controller())

    root.addObject('RequiredPlugin', name='SofaSparseSolver')
    root.addObject('RequiredPlugin', name='SofaGraphComponent')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaPreconditioner')
    root.addObject('RequiredPlugin', name='SofaBoundaryCondition')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    root.addObject('APIVersion', level='17.06')

    root.addObject('RegularGridTopology', name='grid', min=[-radius, -radius, -length/2], max=[radius, radius, length/2], n=[nx, nx, nz])

    for s in direct_solvers:
        name = s['name']
        cg_solver = s['solver']
        arguments = s['arguments']

        meca = root.addChild(name)
        meca.addObject('StaticODESolver', newton_iterations=number_of_newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=False)
        meca.addObject(cg_solver, **arguments)

        meca.addObject('MechanicalObject', name='mo', position='@../grid.position')
        meca.addObject('HexahedronSetTopologyContainer', name='mechanical_topology', src='@../grid')
        meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=3000, poisson_ratio=0)
        meca.addObject('HyperelasticForcefield')

        meca.addObject('BoxROI', name='base_roi', box=[-radius-eps, -radius-eps, -length/2-eps, radius+eps, radius+eps, -length/2+eps])
        meca.addObject('BoxROI', name='top_roi',  box=[-radius-eps, -radius-eps, +length/2-eps, radius+eps, radius+eps, +length/2+eps], quad='@mechanical_topology.quads')

        meca.addObject('FixedConstraint', indices='@base_roi.indices')
        meca.addObject('TractionForce', traction=[0, -30, 0], slope=1/5, quads='@top_roi.quadInROI')


if __name__ == "__main__":
    import Sofa.Simulation
    import Sofa.Core
    import SofaRuntime


    root = Sofa.Core.Node()
    createScene(root)
    Sofa.Simulation.init(root)
    print("Computing... (this may take a while)")

    # Sofa.Simulation.animate(root, 1)

    import Sofa.Gui
    Sofa.Gui.GUIManager.Init("simple_beam", "qglviewer")
    Sofa.Gui.GUIManager.createGUI(root, __file__)
    Sofa.Gui.GUIManager.SetDimension(1080, 1080)
    Sofa.Gui.GUIManager.MainLoop(root)
    Sofa.Gui.GUIManager.closeGUI()
