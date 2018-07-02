from ..Base import BaseObject
from ..Behavior import Behavior
from ..Material import Material

import json


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
        self.export_directory = kwargs.get('export_directory')

        # Members
        self.cases = []

    def serialize(self, keys=[], recursive=True):
        keys = keys + ['name', 'number_of_steps', 'young_modulus', 'poisson_ratio', 'density', 'unit']
        return BaseObject.serialize(self, keys)

    def save(self, filepath=None):
        if filepath is None:
            filename = self.export_directory

    def add(self, case):
        assert isinstance(case, Case)
        self.cases.append(case)
        return case

    def create_report(self):
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

        # Members
        self.solution = None
        self.solution_mesh_filename = ""
        self.solution_image_filename = ""

        assert isinstance(self.material, Material) or isinstance(self.material, tuple) or issubclass(self.material, Material)
        assert isinstance(self.behavior, Behavior) or isinstance(self.behavior, tuple) or issubclass(self.behavior, Behavior)

        # Members
        self.steps = []
        self.behavior_mesh = None

    def add_step(self, step):
        assert isinstance(step, Step)
        self.steps.append(step)


class Step(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.duration = kwargs.get('duration', 0)
        self.newtonsteps = kwargs.get('newtonsteps', [])


class NewtonStep(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.duration = kwargs.get('duration', 0)
        self.residual = kwargs.get('residual', 0)
        self.correction = kwargs.get('correction', 0)