from .Mesh import Part


class Material(object):
    def __init__(self, **kwargs):
        # Parameters
        self.part = kwargs.get('part', None)
        self.young_modulus = kwargs.get('young_modulus', 5000)
        self.poisson_ratio = kwargs.get('poisson_ratio', 0.49)
        self.density = kwargs.get('density', 1)

        assert isinstance(self.part, Part)
        assert self.poisson_ratio < 0.5


class StVenantKirchhoff(Material):
    def __init__(self, **kwargs):
        Material.__init__(self, **kwargs)


class CorotatedStVenantKirchhoff(StVenantKirchhoff):
    def __init__(self, **kwargs):
        StVenantKirchhoff.__init__(self, **kwargs)