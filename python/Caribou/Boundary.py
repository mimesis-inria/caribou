from .Mapping import *
from .Mesh import SurfacePart
from .Base import BaseObject


class Boundary(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.part = kwargs.get('part', None)
        self.printLog = kwargs.get('print_log', False)

        # Implicit mapping (automatically created)
        self.linked_to = kwargs.get('linked_to', None)
        self.link_type = kwargs.get('link_type', None)
        # Or Explicit mapping
        self.mapping = kwargs.get('mapping', None)

        assert isinstance(self.part, SurfacePart)

        if self.mapping:
            assert isinstance(self.mapping, Mapping)
        elif self.linked_to:
            assert isinstance(self.linked_to, Part)

            if not self.link_type:
                if self.part.mesh == self.linked_to.mesh:
                    self.link_type = IdentityMapping
                else:
                    self.link_type = BarycentricMapping
            else:
                assert issubclass(self.link_type, Mapping)

            self.mapping = self.link_type(input=self.linked_to, output=self.part)


class FixedBoundary(Boundary):
    def __init__(self, **kwargs):
        Boundary.__init__(self, **kwargs)

        # Parameters
        self.velocities_fixed = kwargs.get('velocities_fixed', True)


class PressureBoundary(Boundary):
    def __init__(self, **kwargs):
        Boundary.__init__(self, **kwargs)

        # Parameters
        self.pressure = kwargs.get('pressure', [0, -1, 0])
        self.slope = kwargs.get('slope', 0)
        self.number_of_steps_before_increment = kwargs.get('number_of_steps_before_increment', 1)

        assert len(self.pressure) == 3


class ConstantForceBoundary(Boundary):
    def __init__(self, **kwargs):
        Boundary.__init__(self, **kwargs)

        # Parameters
        self.force = kwargs.get('force', [0, -1, 0])

        assert len(self.force) == 3


class VisualBoundary(Boundary):
    def __init__(self, **kwargs):
        Boundary.__init__(self, **kwargs)

        # Parameters
        self.pressure = kwargs.get('color', [1, 0, 0])


class WatcherBoundary(Boundary):
    def __init__(self, **kwargs):
        Boundary.__init__(self, **kwargs)

        # Private members
        self.__mo = None

    def setState(self, state):
        self.__mo = state

    @property
    def state(self):
        return self.__mo
