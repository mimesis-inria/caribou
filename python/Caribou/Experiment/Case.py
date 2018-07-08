from ..Base import BaseObject, serialize
from ..Behavior import Behavior
from ..Material import Material
from ..Utils import escape

import os
import json


class Case(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.name = kwargs.get('name')
        self.element_size = kwargs.get('element_size')
        self.material = kwargs.get('material')  # material OR material_type OR (material_type, options)
        self.behavior = kwargs.get('behavior')  # behavior OR behavior_type OR (behavior_type, options)
        self.solver = kwargs.get('solver')
        self.number_of_steps = kwargs.get('number_of_steps', 10)
        self.young_modulus = kwargs.get('young_modulus', 5e7)
        self.poisson_ratio = kwargs.get('poisson_ratio', 0.45)
        self.density = kwargs.get('density', 2.3)
        self.unit = kwargs.get('unit', {'length': 'mm', 'mass': 'mg', 'pressure': 'Pa', 'load': 'N'})
        self.link_type = kwargs.get('link_type', None)
        self.solution_surface_mesh = kwargs.get('solution_surface_mesh', None)
        self.solution_behavior_mesh = kwargs.get('solution_behavior_mesh', None)
        self.initial_behavior_mesh = kwargs.get('initial_behavior_mesh', None)
        self.id = kwargs.get('id', None)
        self.experiment = kwargs.get('experiment')
        self.run_date = kwargs.get('run_date')
        self.run_memory = kwargs.get('run_memory')
        self.filepath = kwargs.get('filepath')

        # Members
        self.steps = []

        assert isinstance(self.material, Material) or isinstance(self.material, tuple) or issubclass(self.material, Material)
        assert isinstance(self.behavior, Behavior) or isinstance(self.behavior, tuple) or issubclass(self.behavior, Behavior)

    def __eq__(self, other):
        """
        Compare the case with another one
        :param other: The other case
        :type other: Case
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = self.element_size == other.element_size
        res = res and (self.material == other.material)
        res = res and (self.behavior == other.behavior)
        res = res and (self.solver == other.solver)
        res = res and (self.link_type == other.link_type)
        res = res and (self.initial_behavior_mesh == other.initial_behavior_mesh)
        res = res and (self.number_of_steps == other.number_of_steps)
        res = res and (self.young_modulus == other.young_modulus)
        res = res and (self.poisson_ratio == other.poisson_ratio)
        res = res and (self.density == other.density)
        return res

    def add_step(self, step):
        assert isinstance(step, Step)
        self.steps.append(step)

    def save(self):
        if self.experiment is None:
            raise RuntimeError('Trying to save a case without an experiment set.')

        if self.solution_surface_mesh is not None:
            # todo(jnbrunet2000@gmail.com): Exporting as vtk file will failed when further import (the field_data will be lost)
            filename = 'solution_surface_{}_{}.vtk'.format(escape(self.experiment.name), escape(self.name))
            solution_mesh_filepath = os.path.join(self.experiment.directory, filename)
            self.solution_surface_mesh.save(solution_mesh_filepath)

        if self.solution_behavior_mesh is not None:
            # todo(jnbrunet2000@gmail.com): Exporting as vtk file will failed when further import (the field_data will be lost)
            filename = 'solution_behavior_{}_{}.vtk'.format(escape(self.experiment.name), escape(self.name))
            solution_mesh_filepath = os.path.join(self.experiment.directory, filename)
            self.solution_behavior_mesh.save(solution_mesh_filepath)

        filepath = os.path.join(self.experiment.directory, 'case_{}_{}.json'.format(escape(self.experiment.name), escape(self.name)))
        with open(filepath, 'w') as f:
            json.dump(serialize(self), f)
            print "Case exported at '{}'".format(filepath)

        self.filepath = filepath

    def serialize(self):
        return dict(BaseObject.serialize(self), **{
            'name': self.name,
            'steps': self.steps,
            'material': self.material,
            'behavior': self.behavior,
            'number_of_steps': self.number_of_steps,
            'young_modulus': self.young_modulus,
            'poisson_ratio': self.poisson_ratio,
            'density': self.density,
            'unit': self.unit,
            'solution_surface_mesh': self.solution_surface_mesh,
            'solution_behavior_mesh': self.solution_behavior_mesh,
            'initial_behavior_mesh': self.initial_behavior_mesh,
            'experiment': self.experiment,
            'run_date': self.run_date,
            'run_memory': self.run_memory,
        })

    @classmethod
    def deserialize(cls, **kwargs):
        case = cls(**kwargs)
        steps = kwargs.get('steps', [])
        for step in steps:
            case.add_step(step)
        return case


class Step(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.duration = kwargs.get('duration', 0)
        self.newtonsteps = kwargs.get('newtonsteps', [])

    def serialize(self):
        return dict(BaseObject.serialize(self), **{
            'newtonsteps': self.newtonsteps
        })


class NewtonStep(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.duration = kwargs.get('duration', 0)
        self.residual = kwargs.get('residual', 0)
        self.correction = kwargs.get('correction', 0)