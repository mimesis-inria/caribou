from .Experiment import *
from .. import Mesh
from ..Simulation import Simulation
from ..Boundary import *
from ..Material import Material
from ..Behavior import Behavior, MeshlessGalerkin
from ..Simulation import SofaSceneBuilder
from ..Utils import escape, memory_usage, bbox
from ..Report import HtmlReport
from ..PDE import *
from ..Optimization import *
from ..View import ParaView

import json
import os, sys
from math import pi as PI
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class CylinderExperiment(Experiment):
    def __init__(self, **kwargs):
        Experiment.__init__(self, **kwargs)

        # Parameters
        self.radius = kwargs.get('radius', 7.5)
        self.length = kwargs.get('length', 80)
        self.pressure = kwargs.get('pressure', [0, -1e4, 0])
        self.surface_size = kwargs.get('surface_size')

        # Members
        self.surface_mesh = Mesh.cylinder(center1=[0, 0, 0], center2=[0, 0, self.length], radius=self.radius,
                                          size=self.surface_size, dimension=2, quads=False)

    def add(self, case):
        if not isinstance(case.behavior_mesh, Mesh.Mesh):
            case.behavior_mesh = Mesh.cylinder(
                center1=[0, 0, 0], center2=[0, 0, self.length], radius=self.radius, size=case.element_size,
                dimension=3, quads=False)

        # Material setup
        mat_options = {
            'part': case.behavior_mesh.volume,
            'young_modulus': self.young_modulus,
            'poisson_ratio': self.poisson_ratio,
            'density': self.density,
        }

        if isinstance(case.material, tuple):
            m, options = case.material
            if isinstance(options, dict):
                mat_options.update(options)
            case.material = m(**mat_options)
        elif type(case.material) == type and issubclass(case.material, Material):
            case.material = case.material(**mat_options)

        # Behavior setup
        beh_options = {
            'part': case.behavior_mesh.volume,
        }

        if isinstance(case.behavior, tuple):
            b, options = case.behavior
            if isinstance(options, dict):
                beh_options.update(options)
            case.behavior = b(**beh_options)
        elif type(case.behavior) == type and issubclass(case.behavior, Behavior):
            case.behavior = case.behavior(**beh_options)

        return Experiment.add(self, case)

    def serialize(self):
        return dict(Experiment.serialize(self), **{
            'radius': self.radius,
            'length': self.length,
            'pressure': self.pressure,
            'surface_size': self.surface_size,
            'surface_mesh': self.surface_mesh,
        })

    def save(self, filepath=None):

        # todo(jnbrunet2000@gmail.com): Exporting as vtk file will failed when further import (the field_data will be lost)
        # # INITIAL SURFACE VTK EXPORT
        # export_directory = os.getcwd()
        # export_filename = "initial_surface_{}.vtk".format(escape(self.name))
        #
        # if filepath is not None:
        #     if os.path.isfile(filepath):
        #         export_directory = os.path.dirname(filepath)
        #         export_filename = os.path.basename(filepath)
        #     elif os.path.isdir(filepath):
        #         export_directory = filepath
        #
        # vtkfilepath = os.path.join(export_directory, export_filename)
        #
        # if self.surface_mesh is not None and self.surface_mesh.vertices.size > 0:
        #     self.surface_mesh.save(vtkfilepath)
        #
        # # Behavior cases meshes
        # for case in self.cases:
        #     if case.behavior_mesh is not None:
        #         export_filename = "initial_volume_{}.vtk".format(escape(case.name))
        #         vtkfilepath = os.path.join(export_directory, export_filename)
        #         case.behavior_mesh.save(vtkfilepath)

        Experiment.save(self, filepath=filepath)

    def create_report(self, filepath):
        export_directory = os.getcwd()
        export_filename = "{}.html".format(escape(self.name))

        if filepath is not None:
            if os.path.isfile(filepath):
                export_directory = os.path.dirname(filepath)
                export_filename = os.path.basename(filepath)
            elif os.path.isdir(filepath):
                export_directory = filepath

        filepath = os.path.join(export_directory, export_filename)

        pressure = np.linalg.norm(self.pressure)

        ntrian = self.surface_mesh.surface.triangles.shape[0]
        nquads = self.surface_mesh.surface.quads.shape[0]

        report = HtmlReport(name=self.name)
        report.add_section(name="Simulation")
        report.add_list(name="Parameters", attributes=[
            ('Radius', '{} {}'.format(self.radius, self.unit['length'])),
            ('Length', '{} {}'.format(self.length, self.unit['length'])),
            ('Volume', '{} {}<sup>3</sup>'.format(self.radius * self.radius * PI * self.length, self.unit['length'])),
            ('Density', '{} {}/{}<sup>3</sup>'.format(self.density, self.unit['mass'], self.unit['length'])),
            ('Mass', '{} {}'.format(self.radius * self.radius * PI * self.length * self.density, self.unit['mass'])),
            ('Pressure', '{} {}'.format(pressure, self.unit['pressure'])),
            ('Load surface', '{} {}<sup>2</sup>'.format(self.radius * self.radius * PI, self.unit['length'])),
            ('Load force', '{} {}'.format(self.radius * self.radius * PI * pressure / 1e6, self.unit['load'])),
            ('Young modulus', '{} {}'.format(self.young_modulus, self.unit['pressure'])),
            ('Poisson ratio', self.poisson_ratio),
            ('Number of steps', self.number_of_steps),
            ('Number of surface elements',
             ntrian + nquads),
        ])

        count = 0
        for case in self.cases:
            count = count + 1

            nnodes = case.behavior_mesh.vertices.shape[0]
            ntetra = case.behavior_mesh.volume.tetrahedrons.shape[0]
            nhexas = case.behavior_mesh.volume.hexahedrons.shape[0]

            report.add_section('Experiment {} : {}'.format(count, case.name))
            report.add_image_from_meshes(
                meshes=[self.surface_mesh, case.solution_mesh],
                view_attributes=[
                    {'line_width': 0.01,'color': [1, 0, 0],'opacity': 0.1},
                    {}
                ]
            )
            report.add_list(name='Mesh', attributes=[
                ('Number of nodes', nnodes),
                ('Number of tetrahedrons', ntetra),
                ('Number of hexahedrons', nhexas),
            ])

            report.add_list('PDE Solver', [('Type', case.solver.fullname())] + case.solver.printable_attributes())

            system_solver = case.solver.solver
            if isinstance(system_solver, NonLinearSolver):
                report.add_list(
                    'Nonlinear solver',
                    [('Type', system_solver.fullname())] + system_solver.printable_attributes()
                )
                linear_solver = system_solver.linearSolver
            else:
                linear_solver = system_solver

            report.add_list(
                'Linear solver',
                [('Type', linear_solver.fullname())] + linear_solver.printable_attributes()
            )

            report.add_list(
                'Material',
                [('Type', case.material.fullname())] + case.material.printable_attributes()
            )

            report.add_list(
                'Behavior',
                [('Type', case.behavior.fullname())] + case.behavior.printable_attributes()
            )

            p = pressure * 1 / self.number_of_steps * self.radius * self.radius * PI
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
                p = pressure * (i + 1) / self.number_of_steps * self.radius * self.radius * PI
                for newtonstep in step.newtonsteps:
                    lasttime = lasttime + newtonstep.duration
                    newtonsteptimes.append(lasttime)
                    forces.append(p + newtonstep.residual)
                steptimes.append(lasttime)
                pressures.append(p)

            if len(newtonsteptimes) < 2:
                report.add_paragraph(name='Convergence', text='Simulation has diverged')
                continue

            # Convergence graph
            plt.figure(figsize=(20, 10), dpi=300)

            df = pd.DataFrame({'time': newtonsteptimes, 'internal force': forces})
            plt.semilogy('time', 'internal force', data=df, color='skyblue', linewidth=1)

            dp = pd.DataFrame({'time': steptimes, 'external force': pressures})
            plt.step('time', 'external force', data=dp, marker='', color='red', linewidth=1, linestyle='dashed')

            for vline in steptimes:
                plt.axvline(x=vline, color='k', linestyle='--')

            plt.legend()
            img = os.path.realpath(
                os.path.join(export_directory, 'convergence_graph_{}.png'.format(escape(case.name))))
            plt.savefig(img, bbox_inches='tight')
            report.add_image(name="Convergence", path=img)
            print "Convergence exported at {}".format(img)

        report_path = os.path.join(export_directory, 'report_{}.html'.format(escape(self.name)))
        report.write(filepath)
        print "Report exported at '{}'".format(report_path)

    def run(self):
        sofa_simulation = self.sofa.createSimulation("DAG", self.name)
        self.sofa.setSimulation(sofa_simulation)
        self.sofa.timerSetEnabled(self.name, True)
        self.sofa.timerSetInterval(self.name, 1)
        self.sofa.timerSetOutputType(self.name, 'json')

        count = 0
        for case in self.cases:
            print "======= RUNNING CASE {} =======".format(case.name)
            count = count + 1

            simulation = Simulation()
            simulation.add_meshes([
                self.surface_mesh,
                case.behavior_mesh
            ])
            simulation.set_PDE_solver(case.solver)

            boundaries = [
                FixedBoundary(part=case.behavior_mesh.base, linked_to=case.behavior_mesh.volume),
                PressureBoundary(part=case.behavior_mesh.top, pressure=self.pressure, slope=1. / self.number_of_steps,
                                 linked_to=case.behavior_mesh.volume),
            ]
            watcher = WatcherBoundary(part=self.surface_mesh.surface, linked_to=case.behavior_mesh.volume, link_type=case.link_type)
            simulation.add_boundaries(boundaries + [watcher])

            simulation.add_materials([
                case.material
            ])

            simulation.add_behaviors([
                case.behavior
            ])

            # Launch the sofa's simulation
            print "Memory usage before scene creation : {} MB".format(memory_usage())
            root = self.sofa.createNode("root")
            SofaSceneBuilder(simulation=simulation, node=root)
            sofa_simulation.init(root)
            print "Memory usage after scene creation : {} MB".format(memory_usage())

            if isinstance(case.behavior, MeshlessGalerkin):
                if case.behavior.object is not None:
                    case.behavior.stats_number_of_integration_points = case.behavior.object.number_of_integration_points
                    case.behavior.stats_integration_points_per_particle = case.behavior.object.stats_integration_points_per_particle[0]
                    case.behavior.stats_particles_per_integration_point = case.behavior.object.stats_particles_per_integration_point[0]

            for i in range(self.number_of_steps):
                self.sofa.timerBegin(self.name)
                root.simulationStep(1)
                timer_output = '{' + str(self.sofa.timerEnd(self.name, root)) + '}'

                if timer_output not in ['{None}', '{}']:
                    j = json.loads(timer_output)
                    try:
                        step_timer_output = j[j.keys()[0]]['records'][self.name]
                        mechanical = step_timer_output['Simulation::animate']['AnimateVisitor']['Mechanical']
                    except KeyError as err:
                        print "[ERROR] No timing records in the simulation's output. (missing key '{}')"\
                            .format(err.message)
                        break
                    try:
                        newtonraphson = mechanical['NewtonRaphsonSolver::Solve']
                        try:
                            newton_nb_iterations = newtonraphson['nb_iterations']
                        except KeyError:
                            print "[ERROR] No newton iterations at step {}".format(i)
                            break

                        newton_steps = []
                        for ii in range(int(newton_nb_iterations)):
                            newton_step = mechanical['NewtonRaphsonSolver::Solve']['step_{}'.format(ii)]
                            newton_steps.append(NewtonStep(
                                duration=newton_step['end_time'] - newton_step['start_time'],
                                residual=newton_step['residual'],
                                correction=newton_step['correction']
                            ))
                        case.add_step(Step(
                            duration=mechanical['end_time'] - mechanical['start_time'],
                            newtonsteps=newton_steps
                        ))
                    except KeyError:
                        case.add_step(Step(
                            duration=mechanical['end_time'] - mechanical['start_time'],
                        ))

            case.solution_mesh = Mesh.Mesh(
                vertices=np.array(watcher.state.position),
                parts=self.surface_mesh.parts,
                gmsh=self.surface_mesh.gmsh
            )

            # End the sofa's simulation
            sofa_simulation.unload(root)
