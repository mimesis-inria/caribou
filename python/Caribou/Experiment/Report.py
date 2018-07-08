from ..Report import HtmlReport
from ..Optimization import NonLinearSolver
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

        self.add_list('PDE Solver', [('Type', case.solver.fullname())] + case.solver.printable_attributes())

        system_solver = case.solver.solver
        if isinstance(system_solver, NonLinearSolver):
            self.add_list(
                'Nonlinear solver',
                [('Type', system_solver.fullname())] + system_solver.printable_attributes()
            )
            linear_solver = system_solver.linearSolver
        else:
            linear_solver = system_solver

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

    def add_convergence_comparison(self, name, cases=[]):
        if not cases:
            return

        if not isinstance(cases, list):
            cases = [cases]

        colors = generate_n_colors(len(cases))
        pressure = np.linalg.norm(self.experiment.pressure) * self.experiment.radius * self.experiment.radius * PI

        plt.figure(figsize=(20, 10), dpi=300)
        maximum_number_of_increment = 0

        i = 0
        for case in cases:
            if not case.steps:
                continue
            color = colors[i]
            pressures = []
            internal_forces = []
            slope = 1. / case.number_of_steps
            maximum_number_of_increment = max(maximum_number_of_increment, len(case.steps))

            j = 1
            for step in case.steps:
                nb_newtonsteps = case.solver.solver.maxIt
                newtonslope = 1. / nb_newtonsteps
                load_completion_at_start = slope*j      # (%)
                load_completion_at_end = slope * (j+1)  # (%)
                for k in range(nb_newtonsteps):
                    newton_completion = newtonslope * k  # Newton completion (%)
                    pressures.append(load_completion_at_start + newton_completion * (
                                load_completion_at_end - load_completion_at_start))
                    if k >= len(step.newtonsteps):
                        internal_force = internal_forces[len(internal_forces)-1]
                    else:
                        newtonstep = step.newtonsteps[k]
                        internal_force = newtonstep.residual + load_completion_at_start*pressure

                    internal_forces.append(internal_force)
                j = j+1

            if not internal_forces:
                continue

            pressures.insert(0, 0)
            internal_forces.insert(0, internal_forces[0])

            df = pd.DataFrame({'load (%)': pressures, case.name: internal_forces})
            plt.semilogy('load (%)', case.name, data=df, color=color, linewidth=1)

            i = i+1

        if maximum_number_of_increment > 0:
            slope = 1. / float(maximum_number_of_increment)
            pressure_states = [0] + [slope * float(i+2) for i in range(maximum_number_of_increment)]
            pressures = [pressure*slope] + [pressure*slope * float(i+1) for i in range(maximum_number_of_increment)]
            dp = pd.DataFrame({'load (%)': pressure_states, 'external force': pressures})
            plt.step('load (%)', 'external force', data=dp, marker='', color='red', linewidth=1, linestyle='dashed')

        plt.legend()
        img = NamedTemporaryFile(suffix='.png')
        plt.savefig(img.name, bbox_inches='tight')
        self.add_image(name=name, path=img.name, binary=True)



