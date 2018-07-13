from ..Report import HtmlReport
from ..PDE import NewtonRaphsonSolver, StaticSolver, ImplicitEuler
from ..View import ParaView
from ..Utils import generate_n_colors
from Cylinder import CylinderExperiment
from Experiment import Case

from tempfile import NamedTemporaryFile
from math import pi as PI
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class CylinderExperimentReport(HtmlReport):
    def __init__(self, **kwargs):
        HtmlReport.__init__(self, **kwargs)

        # Parameters
        self.experiment = kwargs.get('experiment')
        self.solution = kwargs.get('solution')

        assert isinstance(self.experiment, CylinderExperiment)

        if not self.name:
            self.name = self.experiment.name

    def add_all_cases(self):
        for c in sorted(self.experiment.cases, key=lambda case: case.id):
            self.add_case(c)

        return self

    def add_case(self, case):
        assert isinstance(case, Case)

        nnodes = case.initial_behavior_mesh.vertices.shape[0]
        ntetra = case.initial_behavior_mesh.volume.tetrahedrons.shape[0]
        nhexas = case.initial_behavior_mesh.volume.hexahedrons.shape[0]

        self.open_section('Case #{} : {}'.format(case.id, case.name))

        ttime = 0
        for step in case.steps:
            ttime = ttime + step.duration

        pressure = np.linalg.norm(self.experiment.pressure)

        ntrian = self.experiment.surface_mesh.surface.triangles.shape[0]
        nquads = self.experiment.surface_mesh.surface.quads.shape[0]

        self.add_list(name="Informations", attributes=[
            ('Time to convergence', ttime),
            ('Run date', case.run_date),
            ('Memory available before execution', '{} MB'.format(case.run_memory)),
            ('Radius', '{} {}'.format(self.experiment.radius, case.unit['length'])),
            ('Length', '{} {}'.format(self.experiment.length, case.unit['length'])),
            ('Volume',
             '{} {}<sup>3</sup>'.format(self.experiment.radius * self.experiment.radius * PI * self.experiment.length,
                                        case.unit['length'])),
            ('Density', '{} {}/{}<sup>3</sup>'.format(case.density, case.unit['mass'],
                                                      case.unit['length'])),
            ('Mass', '{} {}'.format(
                self.experiment.radius * self.experiment.radius * PI * self.experiment.length * case.density,
                case.unit['mass'])),
            ('Pressure', '{} {}'.format(pressure, case.unit['pressure'])),
            ('Load surface', '{} {}<sup>2</sup>'.format(self.experiment.radius * self.experiment.radius * PI,
                                                        case.unit['length'])),
            ('Load force', '{} {}'.format(self.experiment.radius * self.experiment.radius * PI * pressure / 1e6,
                                          case.unit['load'])),
            ('Young modulus', '{} {}'.format(case.young_modulus, case.unit['pressure'])),
            ('Poisson ratio', case.poisson_ratio),
            ('Number of steps', case.number_of_steps),
            ('Number of surface elements', ntrian + nquads),
        ])

        self.add_list(name='Mesh', attributes=[
            ('Number of nodes', nnodes),
            ('Number of tetrahedrons', ntetra),
            ('Number of hexahedrons', nhexas),
        ])

        pde_solver = case.solver
        self.add_list('PDE Solver', [('Type', pde_solver.fullname())] + pde_solver.printable_attributes())

        if isinstance(pde_solver, (StaticSolver, ImplicitEuler)):
            linear_solver = pde_solver.solver
            self.add_list(
                'Linear solver',
                [('Type', linear_solver.fullname())] + linear_solver.printable_attributes()
            )

        self.add_list(
            'Material',
            [('Type', case.material.fullname())] + case.material.printable_attributes()
        )

        self.add_list(
            'Behavior',
            [('Type', case.behavior.fullname())] + case.behavior.printable_attributes()
        )

        meshes = [self.experiment.surface_mesh, case.solution_surface_mesh]
        view_attributes = [
                {'line_width': 0.01, 'color': [0, 0, 0], 'opacity': 0.1, 'representation': ParaView.Representation.Wireframe},
                {'representation': ParaView.Representation.Wireframe}
            ]

        if self.solution is not None and isinstance(self.solution, Case) and not case == self.solution:
            meshes.append(self.solution.solution_surface_mesh)
            view_attributes.append({
                'line_width': 0.01, 'color': [1, 0, 0], 'opacity': 0.4, 'representation': ParaView.Representation.Wireframe
            })

        self.add_image_from_meshes(
            name='Solution',
            meshes=meshes,
            view_attributes=view_attributes
        )

        self.add_convergence_comparison(name="Convergence", cases=[case])

        self.close_section()

        return self

    def add_convergence_comparison(self, name, cases=[], compare_with_actual_load=True, filename=None):
        if not cases:
            return

        if not isinstance(cases, list):
            cases = [cases]

        colors = generate_n_colors(len(cases))
        # pressure = np.linalg.norm(self.experiment.pressure) * self.experiment.radius * self.experiment.radius * PI

        plt.figure(figsize=(20, 10), dpi=300)

        i = 0
        total_nb_of_steps = 0 # only used if the number of cases is 1
        max_load = 0#

        for case in cases:
            if not case.steps:
                continue
            color = colors[i]
            internal_forces = []

            if isinstance(case.solver, NewtonRaphsonSolver):
                nb_newtonsteps = case.solver.maxIt
            elif isinstance(case.solver, StaticSolver):
                nb_newtonsteps = case.solver.newton_iterations
            else:
                raise RuntimeError("Solver type not supported")

            total_nb_of_steps = len(case.steps) * nb_newtonsteps
            completions = [float(p)/total_nb_of_steps for p in range(total_nb_of_steps)]

            j = 0
            pressure = case.steps[len(case.steps)-1].load
            max_load = max(max_load, pressure)

            for step in case.steps:
                lastinternalforce = 0
                for k in range(nb_newtonsteps):
                    if k >= len(step.newtonsteps):
                        internal_force = lastinternalforce
                    else:
                        newtonstep = step.newtonsteps[k]
                        if compare_with_actual_load:
                            internal_force = newtonstep.residual + newtonstep.load
                        else:
                            internal_force = newtonstep.residual
                        lastinternalforce = internal_force

                    internal_forces.append(internal_force)
                j = j+1

            if not internal_forces:
                continue

            # states.insert(0, 0)
            # internal_forces.insert(0, internal_forces[0])
            states = [s*100 for s in completions]

            df = pd.DataFrame({'load (%)': states, case.name: internal_forces})
            plt.semilogy('load (%)', case.name, data=df, color=color, linewidth=1)

            i = i+1

        if compare_with_actual_load:
            if len(cases) == 1:
                completions = [float(p) / total_nb_of_steps * 100 for p in range(total_nb_of_steps)]
                loads = [float(c)*max_load/100 for c in completions]
                dp = pd.DataFrame({'load (%)': completions, 'Applied pressure': loads})
                plt.step('load (%)', 'Applied pressure', data=dp, marker='', color='red', linewidth=1, linestyle='dashed', where='post')
            else:
                completions = range(100)
                loads = [float(c)*max_load/100 for c in completions]
                dp = pd.DataFrame({'load (%)': completions, 'Applied pressure': loads})
                plt.semilogy('load (%)', 'Applied pressure', data=dp, marker='', color='red', linewidth=1, linestyle='dashed')

        plt.xlabel('Load (%)')
        if compare_with_actual_load:
            plt.ylabel('Force')
        else:
            plt.ylabel('Residual')

        plt.legend()
        img = NamedTemporaryFile(suffix='.png')
        plt.savefig(img.name, bbox_inches='tight')
        self.add_image(name=name, path=img.name, binary=True)

        if filename:
            plt.savefig(filename, bbox_inches='tight')
