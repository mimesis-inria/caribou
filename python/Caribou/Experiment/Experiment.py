from ..Base import BaseObject, serialize, deserialize
from ..Behavior import Behavior
from ..Material import Material
from ..Utils import escape, mkdir

import json, os, glob


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
        self.directory = kwargs.get('directory')

        # Members
        self.cases = []
        self.__curid = 0

        # Initialization
        if self.directory is None:
            self.directory = os.path.abspath(os.getcwd())
        elif not os.path.isdir(self.directory):
            mkdir(self.directory)

    def __eq__(self, other):
        """
        Compare the experiment with another one
        :param other: The other experiment
        :type other: Experiment
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        return self.number_of_steps == other.number_of_steps and \
               self.young_modulus == other.young_modulus and \
               self.poisson_ratio == other.poisson_ratio and \
               self.density == other.density

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
            'directory': self.directory
        })

    def case(self, id):
        if isinstance(id, (str, unicode)):
            for c in self.cases:
                if c.name == id:
                    return c

        if isinstance(id, str) and str.isdigit(id) and int(id)-1 < len(self.cases):
            return self.cases[int(id)-1]
        elif isinstance(id, unicode) and unicode.isdigit(id) and int(id)-1 < len(self.cases):
            return self.cases[int(id)-1]
        elif isinstance(id, int) and id-1 < len(self.cases):
            return self.cases[id-1]

        raise LookupError('Failed to find a case with index "{}"'.format(id))

    def save(self):
        export_filename = "{}.json".format(escape(self.name))

        filepath = os.path.join(self.directory, export_filename)
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

    def load_cases_from_directory(self, directory=None):
        if directory is None or not os.path.isdir(directory):
            directory = self.directory

        assert os.path.isdir(directory)

        for filename in glob.iglob(os.path.join(directory, 'case_*.json')):
            try:
                with open(filename, 'r') as f:
                    case = deserialize(json.load(f))
            except Exception as e:
                print "File '{}' is found but can't be loaded and will be ignored.\n\tError: {}".format(filename, e)
                continue

            if not isinstance(case, Case):
                print "File '{}' found but does not contain a case object".format(filename)
                continue

            if case.experiment is None or not case.experiment == self:
                print "Case '{}' was found in file '{}' but was run with a different experiment.".format(case.name, filename)
                continue

            case.filepath = filename
            self.add(case)
            print "Case '{}' automatically loaded from file '{}'".format(case.name, filename)

        self.__curid = 0 # Reset id counter

    def add(self, case):
        assert isinstance(case, Case)
        self.__curid = self.__curid + 1
        case.id = self.__curid

        for c in self.cases:
            if c == case:
                if not (c.name == case.name):
                    print "Renaming case '{}' for '{}'".format(c.name, case.name)
                    c.name = case.name
                c.id = case.id
                return c

        self.cases.append(case)
        case.experiment = self

        return case

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
        res = res and (self.behavior_mesh == other.behavior_mesh)
        return res

    def add_step(self, step):
        assert isinstance(step, Step)
        self.steps.append(step)

    def save(self):
        if self.experiment is None:
            raise RuntimeError('Trying to save a case without an experiment set.')

        if self.solution_mesh is not None:
            # todo(jnbrunet2000@gmail.com): Exporting as vtk file will failed when further import (the field_data will be lost)
            filename = 'solution_surface_{}_{}.vtk'.format(escape(self.experiment.name), escape(self.name))
            solution_mesh_filepath = os.path.join(self.experiment.directory, filename)
            self.solution_mesh.save(solution_mesh_filepath)

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
            'solution_mesh': self.solution_mesh,
            'behavior_mesh': self.behavior_mesh,
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