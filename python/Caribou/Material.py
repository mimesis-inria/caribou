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

    def __eq__(self, other):
        """
        Compare the material with another one
        :param other: The other material
        :type other: Material
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        return BaseObject.__eq__(self, other) and \
               self.young_modulus == other.young_modulus and \
               self.poisson_ratio == other.poisson_ratio and \
               self.density == other.density

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

    def __eq__(self, other):
        """
        Compare the LinearElastic material with another one
        :param other: The other material
        :type other: LinearElastic
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        return Material.__eq__(self, other) and \
               self.corotated == other.corotated and \
               Material.__eq__(self, other)

    def printable_attributes(self):
        return [
                   ('Corotated', self.corotated),
        ] + Material.printable_attributes(self)


class StVenantKirchhoff(Material):
    def __init__(self, **kwargs):
        Material.__init__(self, **kwargs)