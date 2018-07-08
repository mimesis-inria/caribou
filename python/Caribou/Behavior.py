from .Mesh import Part
from .Base import BaseObject


class Behavior(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.part = kwargs.get('part', None)
        self.printLog = kwargs.get('print_log', False)

        assert isinstance(self.part, Part)

    def __eq__(self, other):
        """
        Compare the Behavior with another one
        :param other: The other Behavior
        :type other: Behavior
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = BaseObject.__eq__(self, other)
        res = res and self.part == other.part
        return res


class FEMForceField(Behavior):
    class Corotational:
        LARGE = "LARGE"
        POLAR = "POLAR"
        SVD = "SVD"

    def __init__(self, **kwargs):
        Behavior.__init__(self, **kwargs)

        # Parameters
        self.corotational_method = kwargs.get('corotational_method', FEMForceField.Corotational.LARGE)

    def __eq__(self, other):
        return Behavior.__eq__(self, other)


class CaribouFEMForceField(FEMForceField):
    class Corotational:
        LARGE = 1
        POLAR = 2
        SVD = 3

    def __init__(self, **kwargs):
        FEMForceField.__init__(self, **kwargs)

    def __eq__(self, other):
        return FEMForceField.__eq__(self, other)


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
        self.incremental_rotation = kwargs.get('incremental_rotation', True)

        self.stats_number_of_integration_points = kwargs.get('stats_number_of_integration_points')
        self.stats_integration_points_per_particle = kwargs.get('stats_integration_points_per_particle')
        self.stats_particles_per_integration_point = kwargs.get('stats_particles_per_integration_point')

        # Private members
        self.__object = None

    def __eq__(self, other):
        """
        Compare the Behavior with another one
        :param other: The other Behavior
        :type other: MeshlessGalerkin
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = Behavior.__eq__(self, other)
        res = res and self.grid == other.grid
        res = res and (self.number_of_neighbors == other.number_of_neighbors)
        res = res and (self.dilatation == other.dilatation)
        res = res and (self.surface == other.surface)
        res = res and (self.incremental_rotation == other.incremental_rotation)
        return res

    def setObject(self, object):
        """
        Set the relative sofa object once it is created by the builder
        :param object: The sofa's object of the meshless forcefield
        """
        self.__object = object

    def serialize(self):
        return dict(Behavior.serialize(self), **{
            'stats_number_of_integration_points': self.stats_number_of_integration_points,
            'stats_integration_points_per_particle': self.stats_integration_points_per_particle,
            'stats_particles_per_integration_point': self.stats_particles_per_integration_point,
        })

    @property
    def object(self):
        """
        Get the relative sofa object
        :return: The relative sofa object
        """
        return self.__object

    def printable_attributes(self):
        atts = [
            ('Grid', self.grid),
            ('Dilatation', self.dilatation),
            ('Number of neighbors', self.number_of_neighbors),
            ('Incremental rotation', self.incremental_rotation),
        ]

        atts = atts + [
            ('Number of integration points', self.stats_number_of_integration_points),
            ('Number of integration points per particle', 'min: {}, max: {}, avg: {}'.format(
                self.stats_integration_points_per_particle[0],
                self.stats_integration_points_per_particle[1],
                self.stats_integration_points_per_particle[2]
            )),
            ('Number of particles per integration point', 'min: {}, max: {}, avg: {}'.format(
                self.stats_particles_per_integration_point[0],
                self.stats_particles_per_integration_point[1],
                self.stats_particles_per_integration_point[2]
            )),
        ]

        return atts + BaseObject.printable_attributes(self)