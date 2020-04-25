#!/usr/bin/python3

import numpy as np
import re
import SofaRuntime
import Sofa
from SofaRuntime import Timer
import SofaCaribou

number_of_newton_iterations = 3
number_of_cg_iterations = 1000
threshold = 1e-15
cell_size = 1.5
radius = 5
length = 60

nx = int(2*radius / cell_size)+1
nz = int(length / cell_size) + 1
eps = cell_size/10

cg_solvers = [
    {'name':'None', 'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'None',     'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':threshold}},
    {'name':'Id',   'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'Identity', 'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':threshold}},
    {'name':'Dia',  'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'Diagonal', 'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':threshold}},
#    {'name':'LSDia', 'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'LeastSquareDiagonal', 'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':threshold}},
    {'name':'iChol',  'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'IncompleteCholesky',  'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':threshold}},
    # {'name':'iLU',  'solver':'ConjugateGradientSolver', 'arguments' : {'preconditioning_method':'IncompleteLU',  'maximum_number_of_iterations':number_of_cg_iterations, 'residual_tolerance_threshold':threshold}},

# Sofa solvers
    {'name':'sNone', 'solver':'CGLinearSolver', 'arguments':  {'tolerance':threshold, 'threshold':1e-25, 'iterations':number_of_cg_iterations}},
    {'name':'bJac',  'solver':'PCGLinearSolver', 'arguments': {'tolerance':threshold*threshold, 'iterations':number_of_cg_iterations}, 'precond':'BlockJacobiPreconditioner'},
    {'name':'SSOR',  'solver':'PCGLinearSolver', 'arguments': {'tolerance':threshold*threshold, 'iterations':number_of_cg_iterations}, 'precond':'SSORPreconditioner'},
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
        data['Update global matrix'] = '-'
        data['Precond Analysis'] = '-'
        data['Precond Factorize'] = '-'
        data['Nb of CG iterations'] = '-'
        data['Mean CG iteration time'] = '-'
        data['Total CG time'] = '-'
        if 'ConjugateGradient::ComputeGlobalMatrix' in MBKBuild:
            if 'BuildMatrix' in MBKBuild['ConjugateGradient::ComputeGlobalMatrix']:
                data['Update global matrix'] = MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['BuildMatrix']['total_time']
                data['Precond Factorize'] = MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['PreconditionerFactorization']['total_time']
                if 'PreconditionerAnalysis' in MBKBuild['ConjugateGradient::ComputeGlobalMatrix']:
                    data['Precond Analysis'] = MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['PreconditionerAnalysis']['total_time']
            else:
                data['Update global matrix'] = MBKBuild['ConjugateGradient::ComputeGlobalMatrix']['total_time']
        else:
            data['Update global matrix'] = MBKBuild['total_time']

        if 'PCGLinearSolver::solve' in MBKSolve:
            MBKSolve = MBKSolve['PCGLinearSolver::solve']
            data['Nb of CG iterations'] = str(int(MBKSolve['PCG iterations'] - 1))
        elif 'CG-Solve' in MBKSolve:
            data['Nb of CG iterations'] = str(int(MBKSolve['CG iterations'] - 1))
            MBKSolve = MBKSolve['CG-Solve']

        if 'ConjugateGradient::solve' in MBKSolve:
            update_matrix_time = 0.
            if 'HyperelasticForcefield::update_stiffness' in MBKSolve['ConjugateGradient::solve']:
                update_matrix_time = MBKSolve['ConjugateGradient::solve']['HyperelasticForcefield::update_stiffness']['total_time']
                data['Update global matrix'] += update_matrix_time
            data['Nb of CG iterations'] = str(int(MBKSolve['ConjugateGradient::solve']['nb_iterations']))
            mean_time = 0.
            if 'cg_iteration' in MBKSolve['ConjugateGradient::solve']:
                for cg_iteration in MBKSolve['ConjugateGradient::solve']['cg_iteration']:
                    mean_time += cg_iteration['total_time']
                if len(MBKSolve['ConjugateGradient::solve']['cg_iteration']):
                    data['Mean CG iteration time'] = mean_time / len(MBKSolve['ConjugateGradient::solve']['cg_iteration'])
            data['Total CG time'] = MBKSolve['ConjugateGradient::solve']['total_time'] - update_matrix_time
        elif 'HyperelasticForcefield::addDForce' in MBKSolve:
            update_matrix_time = 0
            if 'HyperelasticForcefield::update_stiffness' in MBKSolve:
                update_matrix_time = MBKSolve['HyperelasticForcefield::update_stiffness']['total_time']
                data['Update global matrix'] += update_matrix_time
            data['Total CG time'] = MBKSolve['total_time'] - update_matrix_time
            mean_time = 0.
            for cg_iteration in MBKSolve['HyperelasticForcefield::addDForce']:
                mean_time += cg_iteration['total_time']
            data['Mean CG iteration time'] = mean_time / len(MBKSolve['HyperelasticForcefield::addDForce'])

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
            Timer.setEnabled("cg_timer", True)
            Timer.begin("cg_timer")

    def onAnimateEndEvent(self, e):
        print("Done")
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
            print("Here are the results. Copy and paste them in a text editor without word wrap to visualize them.")
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

        if 'precond' in s:
            meca.addObject(cg_solver, preconditioners='precond', **arguments)
            meca.addObject(s['precond'], name='precond')
        else:
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

    SofaRuntime.importPlugin('SofaComponentAll')

    root = Sofa.Core.Node()
    createScene(root)
    Sofa.Simulation.init(root)
    print("Computing... (this may take a while)")
    Sofa.Simulation.animate(root, 1)
