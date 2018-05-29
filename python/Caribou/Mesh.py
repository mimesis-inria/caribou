from .Utils import Struct as S
from .Utils import bbox, translate, rotate

import math
import pygmsh
import meshio
import numpy as np
from numpy import linalg as LA
from numpy import pi as PI


class Part(object):
    __count = 0

    def __init__(self, **kwargs):
        Part.__count += 1

        # Parameters
        self.mesh = kwargs.get('mesh', None)
        self.name = kwargs.get('name', 'unnamed_part_{}'.format(Part.__count))
        assert isinstance(self.mesh, Mesh)
        self.points = kwargs.get('points', np.array([], dtype=int))
        self.edges = kwargs.get('edges', np.array([]))
        self.triangles = kwargs.get('triangles', np.array([]))
        self.quads = kwargs.get('quads', np.array([]))
        self.tetrahedrons = kwargs.get('tetrahedrons', np.array([]))
        self.hexahedrons = kwargs.get('hexahedrons', np.array([]))


class SurfacePart(Part):
    def __init__(self, **kwargs):
        Part.__init__(self, **kwargs)

        # Parameters
        self.points = kwargs.get('points', np.array([]))
        self.edges = kwargs.get('edges', np.array([]))
        self.triangles = kwargs.get('triangles', np.array([]))
        self.quads = kwargs.get('quads', np.array([]))


class VolumePart(SurfacePart):
    def __init__(self, **kwargs):
        SurfacePart.__init__(self, **kwargs)

        # Parameters
        self.tetrahedrons = kwargs.get('tetrahedrons', np.array([]))
        self.hexahedrons = kwargs.get('hexahedrons', np.array([]))


class Mesh(object):
    __count = 0

    def __init__(self, **kwargs):
        Mesh.__count += 1
        # Parameters
        self._params = kwargs.copy()
        self.name = kwargs.get('name', 'unnamed_mesh_{}'.format(Mesh.__count))
        self.vertices = kwargs.get('vertices', np.array([]))
        self.bbox = S({
            'min': [0, 0, 0],
            'max': [0, 0, 0]
        })

        self.__parts = []
        parts = kwargs.get('parts', [])
        for p in parts:
            self.add_part(p)

        if len(self.vertices):
            xmin, xmax, ymin, ymax, zmin, zmax = bbox(self.vertices)
            self.bbox = S({
                'min': [xmin, ymin, zmin],
                'max': [xmax, ymax, zmax]
            })

    @property
    def parts(self):
        return self.__parts

    def add_part(self, part):
        """
        Add a part to the mesh. The part will be accessible with a dynamic attribute `mesh.part_name`.
        :param Part part: The part to add
        """
        assert isinstance(part, Part)

        names = [p.name for p in self.__parts]
        assert (part.name not in names)

        try:
            getattr(self, part.name)
        except AttributeError:
            setattr(self, part.name, part)

        self.__parts.append(part)

    def get_part(self, name):
        """
        Get the mesh part from the name
        :param str name: The name of the part
        :return: The part, or None if it doesn't exists
        :rtype: Part
        """
        for p in self.__parts:
            if p.name == name:
                return p
        return None


def fromGmsh(points, cells, point_data, cell_data, field_data, dimension=2):
    parts = []

    mesh = Mesh(vertices=points)

    for field_name in field_data:
        id = field_data[field_name][0]

        if dimension == 2:
            part = SurfacePart(name=field_name, mesh=mesh)
        else:
            part = VolumePart(name=field_name, mesh=mesh)

        if 'vertex' in cell_data:
            mask = np.array(np.ma.masked_not_equal(cell_data['vertex']['gmsh:physical'], id).mask)
            part.points = cells['vertex'][~mask]

        if 'line' in cell_data:
            mask = np.array(np.ma.masked_not_equal(cell_data['line']['gmsh:physical'], id).mask)
            part.edges = cells['line'][~mask]

        if 'triangle' in cell_data:
            mask = np.array(np.ma.masked_not_equal(cell_data['triangle']['gmsh:physical'], id).mask)
            part.triangles = cells['triangle'][~mask]

        if 'quad' in cell_data:
            mask = np.array(np.ma.masked_not_equal(cell_data['quad']['gmsh:physical'], id).mask)
            part.quads = cells['quad'][~mask]

        if 'tetra' in cell_data:
            mask = np.array(np.ma.masked_not_equal(cell_data['tetra']['gmsh:physical'], id).mask)
            part.tetrahedrons = cells['tetra'][~mask]

        if 'hexahedron' in cell_data:
            mask = np.array(np.ma.masked_not_equal(cell_data['hexahedron']['gmsh:physical'], id).mask)
            part.hexahedrons = cells['hexahedron'][~mask]

        if not part.points.size:
            part.points = np.unique(
                np.concatenate((
                    np.unique(part.edges),
                    np.unique(part.triangles),
                    np.unique(part.quads),
                    np.unique(part.tetrahedrons),
                    np.unique(part.hexahedrons),
                ))
            ).astype(int)

        mesh.add_part(part)

    return mesh


def fromGmshFile(filename):
    points, cells, point_data, cell_data, field_data = meshio.read(filename)
    fromGmsh(points, cells, point_data, cell_data, field_data)


def cylinder(center1, center2, radius, size=10, dimension=2, quads=False):
    geo = pygmsh.built_in.Geometry()
    c1 = np.asarray(center1)
    c2 = np.asarray(center2)

    r = radius
    axis = (c2 - c1)
    l = LA.norm(axis)

    # SURFACE
    n = size
    lcar = 2 * PI * r / n
    points = [[r * math.cos(t), r * math.sin(t), 0] for t in [2 * np.pi / n * a for a in range(n)]]
    polygon = geo.add_polygon(
        points,
        lcar=lcar
    )
    top, extruded, lat = geo.extrude(
        polygon.surface,
        translation_axis=[0, 0, l],
        point_on_axis=[0, 0, 0],
        recombine=quads)

    geo._GMSH_CODE.append("Reverse Surface {" + polygon.surface.id + "};")

    geo.add_physical_surface(polygon.surface, 'base')
    geo.add_physical_surface(top, 'top')
    geo.add_physical_surface(lat, 'side')
    geo.add_physical_surface([polygon.surface, top] + lat, 'surface')

    if dimension == 3:
        geo.add_physical_volume(extruded, 'volume')

    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geo, dim=dimension, verbose=False)

    points = rotate(points, [0, 0, l], axis / l)
    points = translate(points, c1)

    return fromGmsh(points, cells, point_data, cell_data, field_data, dimension)
