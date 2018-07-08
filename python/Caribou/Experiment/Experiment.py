from ..Base import BaseObject, serialize, deserialize
from ..Behavior import Behavior
from ..Material import Material
from ..Utils import escape, mkdir
from .Case import Case

import json, os, glob


class Experiment(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.name = kwargs.get('name')

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

    def serialize(self):
        att = BaseObject.serialize(self)
        if 'sofa' in att:
            att.pop('sofa')

        return dict(att, **{
            'name': self.name,
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
