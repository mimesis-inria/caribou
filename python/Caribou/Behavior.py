from .Mesh import Part


class Behavior(object):
    def __init__(self, **kwargs):
        # Parameters
        self.part = kwargs.get('part', None)
        self.printLog = kwargs.get('print_log', False)

        assert isinstance(self.part, Part)


class FEMForceField(Behavior):
    def __init__(self, **kwargs):
        Behavior.__init__(self, **kwargs)


class GravityForceField(Behavior):
    def __init__(self, **kwargs):
        Behavior.__init__(self, **kwargs)


class MeshlessGalerkin(Behavior):
    def __init__(self, **kwargs):
        Behavior.__init__(self, **kwargs)

        # Parameters
        self.verbose = kwargs.get('verbose', False)
        self.grid = kwargs.get('grid', ((0,0,0), (10,10,10), (9, 9, 9)))
        self.surface = kwargs.get('surface', None)
        self.number_of_neighbors = kwargs.get('number_of_neighbors', 8)
        self.dilatation = kwargs.get('dilatation', 1)
