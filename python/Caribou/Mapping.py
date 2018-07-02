from .Mesh import Part
from .Base import BaseObject


class Mapping(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.input = kwargs.get('input', None)
        self.output = kwargs.get('output', None)

        assert isinstance(self.input, Part)
        assert isinstance(self.output, Part)


class IdentityMapping(Mapping):
    def __init__(self, **kwargs):
        Mapping.__init__(self, **kwargs)


class LinearMapping(Mapping):
    def __init__(self, **kwargs):
        Mapping.__init__(self, **kwargs)


class BarycentricMapping(Mapping):
    def __init__(self, **kwargs):
        Mapping.__init__(self, **kwargs)


class SubsetMapping(Mapping):
    def __init__(self, **kwargs):
        Mapping.__init__(self, **kwargs)


class ParticleMapping(Mapping):
    def __init__(self, **kwargs):
        Mapping.__init__(self, **kwargs)

        # Parameters
        self.numberOfNeighbors = kwargs.get('number_of_neighbors', 12)
        self.dilatation = kwargs.get('dilatation', 3)
