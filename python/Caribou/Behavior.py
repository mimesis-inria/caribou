from .Mesh import Part


class Behavior(object):
    def __init__(self, **kwargs):
        # Parameters
        self.part = kwargs.get('part', None)

        assert isinstance(self.part, Part)


class FEMForceField(Behavior):
    def __init__(self, **kwargs):
        Behavior.__init__(self, **kwargs)


class GravityForceField(Behavior):
    def __init__(self, **kwargs):
        Behavior.__init__(self, **kwargs)
