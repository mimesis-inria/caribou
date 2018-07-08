from .Experiment import *
from .. import Mesh
from ..Simulation import Simulation
from ..Boundary import *
from ..Material import Material
from ..Behavior import Behavior, MeshlessGalerkin
from ..Simulation import SofaSceneBuilder
from ..Utils import escape, memory_usage, bbox, memory_available
from ..Report import HtmlReport
from ..PDE import *
from ..Optimization import *
from ..View import ParaView
from .Case import Case, NewtonStep, Step

import json
import os, sys
from math import pi as PI
import math
import numpy as np
from datetime import datetime


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

    def __eq__(self, other):
        """
        Compare the CylinderExperiment with another one
        :param other: The other CylinderExperiment
        :type other: CylinderExperiment
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        return self.radius == other.radius and \
               self.length == other.length and \
               self.pressure == other.pressure and \
               Experiment.__eq__(self, other)

    def add(self, case):
        assert isinstance(case, Case)

        if not isinstance(case.initial_behavior_mesh, Mesh.Mesh):
            case.initial_behavior_mesh = Mesh.cylinder(
                center1=[0, 0, 0], center2=[0, 0, self.length], radius=self.radius, size=case.element_size,
                dimension=3, quads=False)

        # Material setup
        mat_options = {
            'part': case.initial_behavior_mesh.volume,
            'young_modulus': case.young_modulus,
            'poisson_ratio': case.poisson_ratio,
            'density': case.density,
        }

        if isinstance(case.material, tuple):
            m, options = case.material
            if isinstance(options, dict):
                mat_options.update(options)
            case.material = m(**mat_options)
        elif type(case.material) == type and issubclass(case.material, Material):
            case.material = case.material(**mat_options)
        elif isinstance(case.material, Material):
            if not case.material.part == case.initial_behavior_mesh.volume:
                case.material.part = case.initial_behavior_mesh.volume
        else:
            raise AttributeError("The case must have a material class that derives from Material")

        # Behavior setup
        beh_options = {
            'part': case.initial_behavior_mesh.volume,
        }

        if isinstance(case.behavior, tuple):
            b, options = case.behavior
            if isinstance(options, dict):
                beh_options.update(options)
            case.behavior = b(**beh_options)
        elif isinstance(case.behavior, type) and issubclass(case.behavior, Behavior):
            case.behavior = case.behavior(**beh_options)
        elif isinstance(case.behavior, Behavior):
            if not case.behavior.part == case.initial_behavior_mesh.volume:
                case.behavior.part = case.initial_behavior_mesh.volume
        else:
            raise AttributeError("The case must have a behavior class that derives from Behavior")

        return Experiment.add(self, case)

    def serialize(self):
        return dict(Experiment.serialize(self), **{
            'radius': self.radius,
            'length': self.length,
            'pressure': self.pressure,
            'surface_size': self.surface_size,
            'surface_mesh': self.surface_mesh,
        })

    def save(self):
        Experiment.save(self)

    def run(self):
        self.save()
        sofa_simulation = self.sofa.createSimulation("DAG", self.name)
        self.sofa.setSimulation(sofa_simulation)
        self.sofa.timerSetEnabled(self.name, True)
        self.sofa.timerSetInterval(self.name, 1)
        self.sofa.timerSetOutputType(self.name, 'json')

        count = 0
        for case in self.cases:
            count = count + 1

            print "======= RUNNING CASE {} =======".format(case.name)
            if case.solution_behavior_mesh is not None:
                print "-> Solution of case #{} ({}) previously computed.".format(count, case.name)
                continue

            case.run_date = datetime.now()
            case.run_memory = memory_available()

            simulation = Simulation()
            simulation.add_meshes([
                self.surface_mesh,
                case.initial_behavior_mesh
            ])
            simulation.set_PDE_solver(case.solver)

            boundaries = [
                FixedBoundary(part=case.initial_behavior_mesh.base, linked_to=case.initial_behavior_mesh.volume),
                PressureBoundary(part=case.initial_behavior_mesh.top, pressure=self.pressure, slope=1. / case.number_of_steps,
                                 linked_to=case.initial_behavior_mesh.volume),
            ]
            watcher = WatcherBoundary(part=self.surface_mesh.surface, linked_to=case.initial_behavior_mesh.volume, link_type=case.link_type)
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

            for i in range(case.number_of_steps):
                self.sofa.timerBegin(self.name)
                root.simulationStep(1)
                timer_output = '{' + str(self.sofa.timerEnd(self.name, root)) + '}'

                if timer_output not in ['{None}', '{}']:
                    j = json.loads(timer_output)
                    try:
                        step_timer_output = j[j.keys()[0]]['records'][self.name]
                        mechanical = step_timer_output['Simulation::animate']['AnimateVisitor']['Mechanical']
                    except KeyError as err:
                        print "[ERROR] No timing records in the simulation's output. (missing key '{}')" \
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

            case.solution_surface_mesh = Mesh.Mesh(
                vertices=np.array(watcher.state.position),
                parts=self.surface_mesh.parts,
                gmsh=self.surface_mesh.gmsh
            )

            case.solution_behavior_mesh = Mesh.Mesh(
                vertices=np.array(case.behavior.state.position),
                parts=case.initial_behavior_mesh.parts,
                gmsh=case.initial_behavior_mesh.gmsh
            )

            case.save()

            # End the sofa's simulation
            sofa_simulation.unload(root)
