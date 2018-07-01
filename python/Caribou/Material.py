from .Mesh import Part
from .Base import BaseObject


class Material(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.part = kwargs.get('part', None)
        self.young_modulus = kwargs.get('young_modulus', 5000)
        self.poisson_ratio = kwargs.get('poisson_ratio', 0.49)
        self.density = kwargs.get('density', 1)

        assert isinstance(self.part, Part)
        assert self.poisson_ratio < 0.5

    def printable_attributes(self):
        return [
                   ('Young modulus', self.young_modulus),
                   ('Poisson ratio', self.poisson_ratio),
                   ('Density', self.density),
        ] + BaseObject.printable_attributes(self)


class LinearElastic(Material):
    def __init__(self, **kwargs):
        Material.__init__(self, **kwargs)

        # Parameters
        self.corotated = kwargs.get('corotated', False)

    def printable_attributes(self):
        return [
                   ('Corotated', self.corotated),
        ] + Material.printable_attributes(self)


class StVenantKirchhoff(Material):
    def __init__(self, **kwargs):
        Material.__init__(self, **kwargs)