from ..Report import HtmlReport
from ..Optimization import NonLinearSolver
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

        assert isinstance(self.experiment, CylinderExperiment)

        pressure = np.linalg.norm(self.experiment.pressure)

        ntrian = self.experiment.surface_mesh.surface.triangles.shape[0]
        nquads = self.experiment.surface_mesh.surface.quads.shape[0]

        self.add_section(name="Simulation")
        self.add_list(name="Parameters", attributes=[
            ('Radius', '{} {}'.format(self.experiment.radius, self.experiment.unit['length'])),
            ('Length', '{} {}'.format(self.experiment.length, self.experiment.unit['length'])),
            ('Volume', '{} {}<sup>3</sup>'.format(self.experiment.radius * self.experiment.radius * PI * self.experiment.length, self.experiment.unit['length'])),
            ('Density', '{} {}/{}<sup>3</sup>'.format(self.experiment.density, self.experiment.unit['mass'], self.experiment.unit['length'])),
            ('Mass', '{} {}'.format(self.experiment.radius * self.experiment.radius * PI * self.experiment.length * self.experiment.density, self.experiment.unit['mass'])),
            ('Pressure', '{} {}'.format(pressure, self.experiment.unit['pressure'])),
            ('Load surface', '{} {}<sup>2</sup>'.format(self.experiment.radius * self.experiment.radius * PI, self.experiment.unit['length'])),
            ('Load force', '{} {}'.format(self.experiment.radius * self.experiment.radius * PI * pressure / 1e6, self.experiment.unit['load'])),
            ('Young modulus', '{} {}'.format(self.experiment.young_modulus, self.experiment.unit['pressure'])),
            ('Poisson ratio', self.experiment.poisson_ratio),
            ('Number of steps', self.experiment.number_of_steps),
            ('Number of surface elements',
             ntrian + nquads),
        ])

    def add_cases(self, cases=[]):
        if len(cases) == 0:
            for c in self.experiment.cases:
                self.add_case(c)

        else:
            for c in cases:
                self.add_case(c)

        return self

    def add_case(self, case):
        assert isinstance(case, Case)

        pressure = np.linalg.norm(self.experiment.pressure)
        nnodes = case.behavior_mesh.vertices.shape[0]
        ntetra = case.behavior_mesh.volume.tetrahedrons.shape[0]
        nhexas = case.behavior_mesh.volume.hexahedrons.shape[0]

        self.add_section('Case #{} : {}'.format(case.id, case.name))
        self.add_image_from_meshes(
            meshes=[self.experiment.surface_mesh, case.solution_mesh],
            view_attributes=[
                {'line_width': 0.01, 'color': [1, 0, 0], 'opacity': 0.1},
                {}
            ]
        )

        ttime = 0
        for step in case.steps:
            ttime = ttime + step.duration

        self.add_list(name="Informations", attributes=[
            ('Time to convergence', ttime),
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

        p = pressure * 1 / self.experiment.number_of_steps * self.experiment.radius * self.experiment.radius * PI
        steptimes = [0]
        newtonsteptimes = [0]
        pressures = [p]
        if len(case.steps) and len(case.steps[0].newtonsteps):
            forces = [p + case.steps[0].newtonsteps[0].residual]
        else:
            forces = [0]

        lasttime = 0

        for i in range(len(case.steps)):
            step = case.steps[i]
            if not len(step.newtonsteps):
                break
            p = pressure * (i + 1) / self.experiment.number_of_steps * self.experiment.radius * self.experiment.radius * PI
            for newtonstep in step.newtonsteps:
                lasttime = lasttime + newtonstep.duration
                newtonsteptimes.append(lasttime)
                forces.append(p + newtonstep.residual)
            steptimes.append(lasttime)
            pressures.append(p)

        if len(newtonsteptimes) < 2:
            self.add_paragraph(name='Convergence', text='Simulation has diverged')
        else:
            # Convergence graph
            plt.figure(figsize=(20, 10), dpi=300)

            df = pd.DataFrame({'time': newtonsteptimes, 'internal force': forces})
            plt.semilogy('time', 'internal force', data=df, color='skyblue', linewidth=1)

            dp = pd.DataFrame({'time': steptimes, 'external force': pressures})
            plt.step('time', 'external force', data=dp, marker='', color='red', linewidth=1, linestyle='dashed')

            for vline in steptimes:
                plt.axvline(x=vline, color='k', linestyle='--')

            plt.legend()
            img = NamedTemporaryFile(suffix='.png')
            plt.savefig(img.name, bbox_inches='tight')
            self.add_image(name="Convergence", path=img.name, binary=True)

        return self

    def add_convergence_comparison(self, cases=[]):
        pass

