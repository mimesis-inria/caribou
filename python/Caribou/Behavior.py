from .Mesh import Part
from .Base import BaseObject


class Behavior(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

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

        self.stats_number_of_integration_points = kwargs.get('stats_number_of_integration_points')
        self.stats_integration_points_per_particle = kwargs.get('stats_integration_points_per_particle')
        self.stats_particles_per_integration_point = kwargs.get('stats_particles_per_integration_point')

        # Private members
        self.__object = None

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