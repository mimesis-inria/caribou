#!/usr/bin/python3

import numpy as np
import re
import Sofa
from SofaRuntime import Timer
import SofaCaribou

number_of_newton_iterations = 3
number_of_cg_iterations = 1000
cell_size = 1.5
radius = 5
length = 60

nx = int(2*radius / cell_size)+1
nz = int(length / cell_size) + 1
eps = cell_size/10

cg_solvers = [
    {'name':'None', 'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'None',     'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':1e-5}},
    {'name':'Id',   'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'Identity', 'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':1e-5}},
    {'name':'Dia',  'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'Diagonal', 'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':1e-5}},
#    {'name':'LSDia', 'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'LeastSquareDiagonal', 'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':1e-5}},
    {'name':'iChol',  'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'IncompleteCholesky',  'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':1e-5}},
    {'name':'iLU',  'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'IncompleteLU',  'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':1e-5}},
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
        data['Total time'] = newton_record['total_time']
        data['Update global matrix'] = 0
        data['Precond Analysis'] = 0
        data['Precond Factorize'] = 0
        data['Nb of CG iterations'] = 0
        data['Mean CG iteration time'] = 0
        data['Total CG time'] = 0
        if 'ConjugateGradient::ComputeGlobalMatrix' in MBKBuild:
            if 'BuildMatrix' in MBKBuild['ConjugateGradient::ComputeGlobalMatrix']:
                data['Update global matrix'] += MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['BuildMatrix']['total_time']
                data['Precond Factorize'] += MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['PreconditionerFactorization']['total_time']
                if 'PreconditionerAnalysis' in MBKBuild['ConjugateGradient::ComputeGlobalMatrix']:
                    data['Precond Analysis'] += MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['PreconditionerAnalysis']['total_time']
            else:
                data['Update global matrix'] += MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['total_time']

        if 'ConjugateGradient::solve' in MBKSolve:
            update_matrix_time = 0.
            if 'HexahedronElasticForce::compute_k' in MBKSolve['ConjugateGradient::solve']:
                update_matrix_time = MBKSolve['ConjugateGradient::solve']['HexahedronElasticForce::compute_k']['total_time']
                data['Update global matrix'] += update_matrix_time
            data['Nb of CG iterations'] = MBKSolve['ConjugateGradient::solve']['nb_iterations']
            mean_time = 0.
            if 'cg_iteration' in MBKSolve['ConjugateGradient::solve']['cg_iteration']:
                for cg_iteration in MBKSolve['ConjugateGradient::solve']['cg_iteration']:
                    mean_time += cg_iteration['total_time']
                if len(MBKSolve['ConjugateGradient::solve']['cg_iteration']):
                    data['Mean CG iteration time'] = mean_time / len(MBKSolve['ConjugateGradient::solve']['cg_iteration'])
            data['Total CG time'] = MBKSolve['ConjugateGradient::solve']['total_time'] - update_matrix_time
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
                fields[k]['methods'][method_name]['width'] = max(fields[k]['methods'][method_name]['width'], len(number_format.format(v)))
    for field in fields.values():
        field['width'] = max(field['width'], len(' '.join(['{{: ^{}}}'.format(m['width']).format('') for m in field['methods'].values()]))+2)

    # Print Header
    ni_col_width = len("Newton iteration # ") + len(str(maximum_number_of_newton_steps)) + 1
    cols = ["{{: <{}}}".format(ni_col_width).format('')] + [" {{: ^{}}} ".format(v['width']).format(k) for k,v in zip(fields.keys(), fields.values())]
    print("|".join(cols))
    cols = ["{{: <{}}}".format(ni_col_width).format('')] + [" {{: ^{}}} ".format(v['width']).format(' '.join([
        "{{: ^{}}}".format(m['width']).format(name) for name, m in zip(v['methods'].keys(), v['methods'].values())
    ])) for v in fields.values()]
    print("|".join(cols))

    # Print newton iterations
    for it in range(maximum_number_of_newton_steps):
        cols = ["Newton iteration # " + "{{: >{}}} ".format(len(str(maximum_number_of_newton_steps))).format(it)] + \
               [" {{: ^{}}} ".format(f['width']).format(' '.join(["{{: ^{}}}".format(m['width']).format('-' if it>=len(m['values']) else number_format.format(m['values'][it])) for m in f['methods'].values()])) for f in fields.values()]
        print("|".join(cols))


class Controller(Sofa.Core.Controller):
    def __init__(self):
        super().__init__(self)
        self.use_sofa_profiler_timer = False

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
                                'name': 'method_'+str(i),
                                'newton_steps': extract_newton_steps(v)
                            })
                    else:
                        methods.append({
                            'name': 'method',
                            'newton_steps': extract_newton_steps(v)
                        })
            pretty_print_methods(methods)
        if not self.use_sofa_profiler_timer:
            Timer.end("cg_timer")



def createScene(root):
    root.addObject(Controller())
    root.addObject('APIVersion', level='17.06')

    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaPreconditioner')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    root.addObject('RegularGridTopology', name='grid', min=[-radius, -radius, -length/2], max=[radius, radius, length/2], n=[nx, nx, nz])

    for s in cg_solvers:
        name = s['name']
        cg_solver = s['solver']
        arguments = s['arguments']

        meca = root.addChild(name)
        meca.addObject('StaticODESolver', newton_iterations=number_of_newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=False)
        meca.addObject(cg_solver, **arguments)

        meca.addObject('MechanicalObject', name='mo', position='@../grid.position')
        meca.addObject('HexahedronSetTopologyContainer', name='mechanical_topology', src='@../grid')
        meca.addObject('HexahedronElasticForce', topology_container='@mechanical_topology', youngModulus=3000, poissonRatio=0, corotated=False, linearStrain=False)

        meca.addObject('BoxROI', name='base_roi', box=[-radius-eps, -radius-eps, -length/2-eps, radius+eps, radius+eps, -length/2+eps])
        meca.addObject('BoxROI', name='top_roi',  box=[-radius-eps, -radius-eps, +length/2-eps, radius+eps, radius+eps, +length/2+eps], quad='@mechanical_topology.quads')

        meca.addObject('FixedConstraint', indices='@base_roi.indices')
        meca.addObject('TractionForce', traction=[0, -30, 0], slope=1/5, quads='@top_roi.quadInROI')

if __name__ == "__main__":
    import Sofa.Simulation
    import Sofa.Core
    import SofaRuntime
    SofaRuntime.PluginRepository.addFirstPath('/Users/jnbrunet/Sources/sofa/build/lib')
    root = Sofa.Core.Node()
    createScene(root)
    Sofa.Simulation.init(root)
    Sofa.Simulation.animate(root, 1)
