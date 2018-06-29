from .Mapping import Mapping, BarycentricMapping, IdentityMapping
from .Mesh import Part, SurfacePart


class Boundary(object):
    def __init__(self, **kwargs):
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
                assert isinstance(self.link_type, Mapping)

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
