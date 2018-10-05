import sys

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.pyplot as plt
from distutils.util import strtobool

from Helpers import Hexahedron, integrate

sx, sy, sz = (10, 10, 10)
nodes = np.array([
    [0, 0, 0],
    [sx, 0, 0],
    [sx, 0, sz],
    [0, 0, sz],
    [0, sy, 0],
    [sx, sy, 0],
    [sx, sy, sz],
    [0, sy, sz]
])

faces = [
    [0, 1, 2, 3],
    [3, 7, 4, 0],
    [7, 6, 5, 4],
    [5, 6, 2, 1],
    [5, 1, 0, 4],
    [2, 6, 7, 3]
]

edges = [
    [0, 1], # 0
    [1, 2], # 1
    [2, 3], # 2
    [3, 0], # 3
    [4, 5], # 4
    [5, 6], # 5
    [6, 7], # 6
    [7, 4], # 7
    [0, 4], # 8
    [1, 5], # 9
    [2, 6], # 10
    [3, 7] # 11
]

# faces_to_integrate = [
#     [nodes[node_id].tolist() for node_id in face] for face in faces
# ]

# faces_to_integrate = [
#     [[0, 0, 0], [5, 0, 0], [5, 0, 5], [0, 0, 5]],
#     [[0, 0, 5], [0, 5, 5], [0, 5, 0], [0, 0, 0]],
#     [[0, 5, 5], [5, 5, 5], [5, 5, 0], [0, 5, 0]],
#     [[5, 5, 0], [5, 5, 5], [5, 0, 5], [5, 0, 0]],
#     [[5, 5, 0], [5, 0, 0], [0, 0, 0], [0, 5, 0]],
#     [[5, 0, 5], [5, 5, 5], [0, 5, 5], [0, 0, 5]]
# ]

# faces_to_integrate = [
#     # Face 0
#     [[0, 0, 0], [5, 0, 0], [5, 0, 5]],
#     [[5, 0, 5], [0, 0, 5], [0, 0, 0]],
#
#     # Face 1
#     [[0, 0, 5], [0, 5, 5], [0, 5, 0]],
#     [[0, 5, 0], [0, 0, 0], [0, 0, 5]],
#
#     # Face 2
#     [[0, 1.25, 5], [5, 1.25, 5], [5, 1.25, 0]],
#     [[5, 1.25, 0], [0, 1.25, 0], [0, 1.25, 5]],
#
#     # Face 3
#     [[5, 5, 0], [5, 5, 5], [5, 0, 5]],
#     [[5, 0, 5], [5, 0, 0], [5, 5, 0]],
#
#     # Face 4
#     [[5, 5, 0], [5, 0, 0], [0, 0, 0]],
#     [[0, 0, 0], [0, 5, 0], [5, 5, 0]],
#
#     # Face 5
#     [[5, 0, 5], [5, 5, 5], [0, 5, 5]],
#     [[0, 5, 5], [0, 0, 5], [5, 0, 5]]
# ]


# print integrate(faces_to_integrate)[0]

def bbox(vertices):
    points_array = np.asarray(vertices)
    m = np.min(points_array, axis=0)
    xmin, ymin, zmin = m[0], m[1], m[2]

    m = np.max(points_array, axis=0)
    xmax, ymax, zmax = m[0], m[1], m[2]

    return xmin, xmax, ymin, ymax, zmin, zmax

def bruteforce_integrate(f, triangles, number_of_subdivisions = 10):
    vertices = []
    triangles = np.array(triangles)
    for t in triangles:
        vertices.append(t[0])
        vertices.append(t[1])
        vertices.append(t[2])
    xmin, xmax, ymin, ymax, zmin, zmax = bbox(vertices)
    cell_size = max((xmax-xmin), (ymax-ymin), (zmax-zmin)) / number_of_subdivisions
    cell_volume = cell_size * cell_size * cell_size
    sum = 0.
    for i in range(number_of_subdivisions-1):
        x = xmin + i*cell_size
        for j in range(number_of_subdivisions-1):
            y = ymin + j*cell_size
            for k in range(number_of_subdivisions-1):
                z = zmin + k*cell_size
                p = np.array([x,y,z])
                isInside = True
                for t in triangles:
                    c = (t[0] + t[1] + t[2]) / 3.
                    n = np.cross(t[1]-t[0],t[2]-t[0])
                    n = n / np.linalg.norm(n)
                    if np.dot((c - p), n) < 0:
                        isInside = False
                        break
                if isInside:
                    sum = sum + (f(x,y,z)*cell_volume)
    return sum


def displayHexa(ax, isInside):
    outsideNodes = []
    insideNodes = []
    for i in range(8):
        if not isInside[i]:
            outsideNodes.append(nodes[i].tolist())
        else:
            insideNodes.append(nodes[i].tolist())


    insideNodes = np.array(insideNodes)
    outsideNodes = np.array(outsideNodes)

    ax.scatter(insideNodes[:, 0], insideNodes[:, 1], insideNodes[:, 2], c='r')
    # ax.scatter(outsideNodes[:, 0], outsideNodes[:, 1], outsideNodes[:, 2], c='r')

    verts = [([nodes[n] for n in f]) for f in faces]

    # plot sides
    pc = Poly3DCollection(verts, facecolor='cyan', linewidths=1, edgecolors='r', alpha=.05)
    ax.add_collection3d(pc)

def displayTriangles(ax, triangles):
    trianglesNodes = [[], [], []]
    trianglesFaces = []
    centers = []
    normals = []
    for t in triangles:
        t = np.array(t)
        for p in t:
            trianglesNodes[0].append(p[0])
            trianglesNodes[1].append(p[1])
            trianglesNodes[2].append(p[2])

        c = (t[0] + t[1] + t[2]) / 3.
        n = np.cross(t[1]-t[0],t[2]-t[0])
        n = n / np.linalg.norm(n)
        centers.append(c)
        normals.append([c, c + n*0.1])

    pc = Poly3DCollection(triangles, facecolor='red', linewidths=1, edgecolors='b')
    ax.add_collection3d(pc)
    pc = Line3DCollection(normals, linewidths=2, colors='g')
    ax.add_collection3d(pc)
    ax.scatter(trianglesNodes[0], trianglesNodes[1], trianglesNodes[2], c='r')
    ax.scatter(np.array(centers)[:, 0], np.array(centers)[:, 1], np.array(centers)[:, 2], c='black')

def polynome(x, y, z):
    return x


def alphatriangle_to_triangle(alphatriangle):
    triangle = []
    for n1, n2, alpha in alphatriangle:
        p = nodes[n1] + (nodes[n2]-nodes[n1])*alpha
        triangle.append(p.tolist())
    return triangle


def main(argv):
    if len(argv) != 8:
        print 'cut_hexahedron.py b b b b b b b b'
        sys.exit(2)

    isInside = [strtobool(i) for i in argv]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    displayHexa(ax, isInside)

    # Populate cut triangles
    edgesIntersection = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    for i in range(12):
        n1 = edges[i][0]
        n2 = edges[i][1]
        if isInside[n1] and not isInside[n2]:
            edgesIntersection[i] = 0.25
        elif not isInside[n1] and isInside[n2]:
            edgesIntersection[i] = 0.75
        else:
            edgesIntersection[i] = -1

    faces_alphatriangles, cut_alphatriangles = Hexahedron.triangulate(isInside, edgesIntersection)

    faces_triangles = []
    for face in faces_alphatriangles:
        for alphatriangle in face:
            faces_triangles.append(alphatriangle_to_triangle(alphatriangle))

    cut_triangles = [
        alphatriangle_to_triangle(alphatriangle) for alphatriangle in cut_alphatriangles
    ]

    triangles = faces_triangles + cut_triangles

    print integrate(triangles)[1]


    # print "Brute force = {}".format(bruteforce_integrate(
    #     f=lambda x, y, z : x,
    #     triangles=triangles,
    #     number_of_subdivisions=100
    # ))
    displayTriangles(ax, triangles)


    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim3d(-.5, sx+.5)
    ax.set_ylim3d(-.5, sy+.5)
    ax.set_zlim3d(-.5, sz+.5)

    # ax.set_aspect(aspect='auto', adjustable='box')

    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])