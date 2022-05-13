from Caribou import _1D, _2D, _3D, Linear, Quadratic
from Caribou.Topology import Mesh
from Caribou.Geometry import Segment, Triangle, Tetrahedron, Hexahedron
import SofaCaribou
import numpy as np


def assemble(domain, f):
    nnodes = domain.mesh().number_of_nodes()
    dimension = domain.mesh().dimension()
    canonical_dimension = domain.canonical_dimension()
    forces = np.zeros((nnodes, dimension))

    for element_id in range(domain.number_of_elements()):
        element = domain.element(element_id)
        node_indices = domain.element_indices(element_id)
        for gauss_node in element.gauss_nodes():
            xi = gauss_node.position
            w = gauss_node.weight
            N = element.L(xi)
            J = element.jacobian(xi)
            X = element.world_coordinates(xi)

            if dimension == canonical_dimension:
                detJ = np.abs(np.linalg.det(J))
            elif canonical_dimension == 1:
                detJ = np.linalg.norm(J)
            else:
                detJ = np.linalg.norm(np.cross(J[:, 0], J[:, 1]))

            for i, node_index in enumerate(node_indices):
                forces[node_index] += np.ravel(f(*X, element)) * N[i] * w * detJ

    return forces


def integrate(domain, f, u_h=None):
    dimension = domain.mesh().dimension()
    canonical_dimension = domain.canonical_dimension()

    v = None

    for element_id in range(domain.number_of_elements()):
        element = domain.element(element_id)
        node_indices = domain.element_indices(element_id)
        if u_h is not None:
            nodes_u = u_h[node_indices]

        for gauss_node in element.gauss_nodes():
            xi = gauss_node.position
            w = gauss_node.weight
            N = element.L(xi)
            J = element.jacobian(xi)
            X = element.world_coordinates(xi)
            if u_h is not None:
                # Interpolate the displacement at the Gauss node
                U = np.sum(np.multiply(nodes_u, N[:, np.newaxis]), axis=0)

            if dimension == canonical_dimension:
                detJ = np.abs(np.linalg.det(J))
            elif canonical_dimension == 1:
                detJ = np.linalg.norm(J)
            else:
                detJ = np.linalg.norm(np.cross(J[:, 0], J[:, 1]))

            if u_h is not None:
                v_ = f(*X, U, element) * w * detJ
            else:
                v_ = f(*X, element) * w * detJ

            if v is None:
                v = v_
            else:
                v += v_

    return v
