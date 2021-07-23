import numpy as np
import meshio
import tempfile
import subprocess
import os
import argparse

from meshio.xdmf.common import meshio_to_xdmf_type
meshio_to_xdmf_type['tetra10'] = ['Tet_10']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--type', type=str, choices=('p1', 'p2', 'q1', 'q2'), default='p1', help='Type of elements')
    parser.add_argument('-n', nargs='+', required=True, type=int, help='Number of nodes in the x, y [and z] directions, respectively')
    parser.add_argument('--p0', nargs='+', required=True, type=float, help='First corner of the beam')
    parser.add_argument('--p1', nargs='+', required=True, type=float, help='Second corner of the beam')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--gmsh', type=str, default='/usr/bin/gmsh', help='Path to gmsh executable')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output mesh filename.')
    args = parser.parse_args()

    quads = (args.type in ('q1', 'Q1', 'q2', 'Q2'))
    order = 1 if args.type in ('q1', 'Q1', 'p1', 'P1') else 2
    dimension = len(args.n)

    mesh = rectangle(args.p0, args.p1, args.n, dimension, order, quads, args.verbose)
    options = {}
    meshio.write(args.output, mesh, **options)


def geo_to_msh(geo: str, dimension=2, gmsh='/usr/bin/gmsh', verbose=False, geo_filename=None, msh_filename=None):
    """
    Convert a geo (gmsh) data format to a meshio instance.

    This method will first write down the content of the geo inside the geo_filename parameter, then call gmsh using the gmsh
    parameter, and write down the resulting msh inside the msh_filename parameter. Note that if geo_filename and/or
    msh_filename are not specified, temporary files will be used to generate the mesh.

    :param geo: The input content of the geometry (using gmsh's geo format)
    :param dimension: The dimension of the mesh (default to 2)
    :param gmsh: The path to gmsh executable (default to /usr/bin/gmsh)
    :param verbose: Print gmsh log (default to false)
    :param geo_filename: The path to a file where the content of geo will writen. [Optional]
    :param msh_filename: The path to a file where the content of the msh will be writen. [Optional]
    :return: A meshio Mesh instance containing the newly created mesh
    """
    if geo_filename is None:
        with tempfile.NamedTemporaryFile(suffix=".geo") as f:
            geo_filename = f.name
        del_geo = True
    else:
        del_geo = False

    with open(geo_filename, "w") as f:
        f.write(geo)

    if msh_filename is None:
        with tempfile.NamedTemporaryFile(suffix=".msh") as handle:
            msh_filename = handle.name
        del_msh = True
    else:
        del_msh = False

    args = [
        "-{}".format(dimension),
        geo_filename,
        "-format",
        'msh',
        # "-bin",
        "-o",
        msh_filename,
    ]
    p = subprocess.Popen(
        [gmsh] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    while True:
        line = p.stdout.readline()
        if not line:
            break
        if verbose:
            print(line.decode("utf-8"), end="")

    p.communicate()
    assert p.returncode == 0, "Gmsh exited with error (return code {}).".format(
        p.returncode
    )

    mesh = meshio.read(msh_filename)
    if del_geo:
        os.remove(geo_filename)
    if del_msh:
        os.remove(msh_filename)
    return mesh


def rectangle(p0, p1, n=[2,2], dimension=2, order=1, quads=True, verbose=False):
    if len(p0) == 2:
        p0 = [p0[0], p0[1], 0]
        p1 = [p1[0], p1[1], 0]

    if dimension == 3:
        assert len(n) == 3
        assert n[2] > 1

    assert len(n) in [2,3]
    assert n[0] > 1
    assert n[1] > 1

    width, height, length = np.abs(np.asarray(p1) - np.asarray(p0))
    ex = np.array([1, 0, 0])
    ey = np.array([0, 1, 0])

    # Compute the four corners of the rectangular (or top face in 3D)
    p0 = np.asarray(p0)
    p1 = p0 + ey*height
    p2 = p0 + ex*width + ey*height
    p3 = p0 + ex*width

    p = [[str(a) for a in p0] + ['1.0'], [str(a) for a in p1]+ ['1.0'], [str(a) for a in p2]+ ['1.0'], [str(a) for a in p3]+ ['1.0']]


    s = """
    Point(1) = {{{}}};
    Point(2) = {{{}}};
    Point(3) = {{{}}};
    Point(4) = {{{}}};
    
    Line(1) = {{1, 2}};
    Line(2) = {{2, 3}};
    Line(3) = {{3, 4}};
    Line(4) = {{4, 1}};
    
    Line Loop(5) = {{2, 3, 4, 1}};
    Plane Surface(6) = {{5}};
    
    
    Transfinite Surface {{6}} = {{2, 3, 4, 1}};
    Transfinite Line {{1, 3}} = {} Using Progression 1;
    Transfinite Line {{2, 4}} = {} Using Progression 1;
    Reverse Surface {{6}};
    """.format(",".join(p[0]), ",".join(p[1]), ",".join(p[2]), ",".join(p[3]), n[0], n[1])

    if quads:
        s += "Recombine Surface {6};"

    # 3D CASE:
    if length > 0:
        assert len(n) == 3
        assert n[2] > 1

        s += """
        ex[] = Extrude {{0,0,{}}} {{Surface{{6}}; Layers{{{}}}; {}}};
        """.format(length, n[2]-1, 'Recombine;' if quads else '')

        s += 'Physical Surface("base") = {6};'
        s += 'Physical Surface("top") = {ex[0]};'
        s += 'Physical Surface ("surface") = {6, ex[0], ex[2], ex[3], ex[4], ex[5]};'
        s += 'Physical Volume ("volume") = {ex[1]};'

    s += """
        Mesh.ElementOrder = {};
        Mesh.SecondOrderLinear = 0;
        Mesh.SecondOrderIncomplete = 1;
        Mesh {};
        """.format(order, dimension)

    return geo_to_msh(s, dimension=dimension, verbose=verbose)


if __name__ == "__main__":
    main()
