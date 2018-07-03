from ..Base import BaseObject, serialize, deserialize
from ..Behavior import Behavior
from ..Material import Material
from ..Utils import escape

import json, os


class Experiment(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.name = kwargs.get('name')
        self.number_of_steps = kwargs.get('number_of_steps', 10)
        self.young_modulus = kwargs.get('young_modulus', 5e7)
        self.poisson_ratio = kwargs.get('poisson_ratio', 0.45)
        self.density = kwargs.get('density', 2.3)
        self.unit = kwargs.get('unit', {'length': 'mm', 'mass': 'mg', 'pressure': 'Pa', 'load': 'N'})
        self.sofa = kwargs.get('sofa')

        # Members
        self.cases = []

    @classmethod
    def deserialize(cls, **kwargs):
        e = cls(**kwargs)
        cases = kwargs.get('cases', [])
        for c in cases:
            e.add(c)

        return e

    def serialize(self):
        att = BaseObject.serialize(self)
        if 'sofa' in att:
            att.pop('sofa')

        return dict(att, **{
            'name': self.name,
            'number_of_steps': self.number_of_steps,
            'young_modulus': self.young_modulus,
            'poisson_ratio': self.poisson_ratio,
            'density': self.density,
            'unit': self.unit,
            'cases': self.cases
        })

    def save(self, filepath=None):
        export_directory = os.getcwd()
        export_filename = "{}.json".format(escape(self.name))

        if filepath is not None:
            if os.path.isfile(filepath):
                export_directory = os.path.dirname(filepath)
                export_filename = os.path.basename(filepath)
            elif os.path.isdir(filepath):
                export_directory = filepath

        for case in self.cases:
            case.save(export_directory)

        filepath = os.path.join(export_directory, export_filename)
        with open(filepath, 'w') as f:
            json.dump(serialize(self), f)
            print "Experiment exported at '{}'".format(filepath)

    @classmethod
    def load(cls, filepath):
        """

        :param filepath: File path of the json exported experiment
        :return: The experiment object fully loaded
        :rtype: Experiment
        """
        with open(filepath, 'r') as f:
            return deserialize(json.load(f))

    def add(self, case):
        assert isinstance(case, Case)
        self.cases.append(case)
        case.setExperiment(self)

        return case

    def create_report(self, filepath):
        raise NotImplementedError(
            "The tojson functionality isn't implemented for experience of type {}".format(self.classname))

    def run(self):
        raise NotImplementedError(
            "The run functionality isn't implemented for experience of type {}".format(self.classname))


class Case(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.name = kwargs.get('name')
        self.element_size = kwargs.get('element_size')
        self.material = kwargs.get('material')  # material OR material_type OR (material_type, options)
        self.behavior = kwargs.get('behavior')  # behavior OR behavior_type OR (behavior_type, options)
        self.solver = kwargs.get('solver')
        self.link_type = kwargs.get('link_type', None)
        self.solution_mesh = kwargs.get('solution_mesh', None)
        self.behavior_mesh = kwargs.get('behavior_mesh', None)

        # Members
        self.steps = []

        # Private members
        self.__experiment = None

        assert isinstance(self.material, Material) or isinstance(self.material, tuple) or issubclass(self.material, Material)
        assert isinstance(self.behavior, Behavior) or isinstance(self.behavior, tuple) or issubclass(self.behavior, Behavior)

    def setExperiment(self, e):
        self.__experiment = e

    def add_step(self, step):
        assert isinstance(step, Step)
        self.steps.append(step)

    def save(self, directory):
        if self.solution_mesh is None:
            return

        # todo(jnbrunet2000@gmail.com): Exporting as vtk file will failed when further import (the field_data will be lost)
        # filename = 'solution_surface_{}.vtk'.format(escape(self.name))
        # solution_mesh_filepath = os.path.join(directory, filename)
        # self.solution_mesh.save(solution_mesh_filepath)

    def serialize(self):
        return dict(BaseObject.serialize(self), **{
            'steps': self.steps,
            'material': self.material,
            'behavior': self.behavior,
            'solution_mesh': self.solution_mesh,
            'behavior_mesh': self.behavior_mesh,
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