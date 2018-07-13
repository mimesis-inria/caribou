from ..Mesh import Mesh


class VolumetricMesh(Mesh):
    def __init__(self, **kwargs):
        Mesh.__init__(self, **kwargs)

        # Members
        self.__o = None  # Related sofa object
        self.__init_callbacks = []  # Callbacks to be called once the sofa's object is created

    def create(self, simulation, node):
        node.createObject(
            'MeshGenerationFromPolyhedron',
            inputPoints=None,  # Rest position coordinates of the degrees of freedom
            inputTriangles=None,  # List of triangles
            inputQuads=None,  # List of quads (if no triangles)
            facetAngle=25.0,  # Lower bound for the angle in degrees of the surface mesh facets
            facetSize=0.15,  # Uniform upper bound for the radius of the surface Delaunay balls
            facetApproximation=0.008,  # facetApproximation", "Upper bound for the center-center distances of the surface mesh facets
            cellRatio=4.0,  # Upper bound for the radius-edge ratio of the tetrahedra
            cellSize=0.2,  # Uniform upper bound for the circumradii of the tetrahedra in the mesh
            sharpEdgeAngle=120,  # Threshold angle to detect sharp edges in input surface (activated if sharpEdgeSize > 0)
            sharpEdgeSize=0.0,  # Meshing size for sharp feature edges
        )


def from_surface_mesh(surface_mesh):
    assert isinstance(surface_mesh, Mesh)

    mesh = VolumetricMesh()

