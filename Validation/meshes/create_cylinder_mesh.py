import numpy as np
import meshio
import tempfile
import subprocess
import os
import argparse

from create_beam_mesh import geo_to_msh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--type', type=str, choices=('p1', 'p2', 'q1', 'q2'), default='p1', help='Type of elements')
    parser.add_argument('-s', '--size', required=True, type=float, help='Size of the mesh elements')
    parser.add_argument('-r', '--radius', required=True, type=float, help='Radius of the cylinder')
    parser.add_argument('-l', '--length', required=True, type=float, help='Length of the cylinder')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--gmsh', type=str, default='/usr/bin/gmsh', help='Path to gmsh executable')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output mesh filename.')
    args = parser.parse_args()

    quads = (args.type in ('q1', 'Q1', 'q2', 'Q2'))
    order = 1 if args.type in ('q1', 'Q1', 'p1', 'P1') else 2

    mesh = cylinder(args.radius, args.length, args.size, 3, order, quads, args.verbose)
    print(mesh)
    options = {}
    meshio.write(args.output, mesh, **options)


def cylinder(radius, length, size=0.5, dimension=2, order=1, verbose=False, geo_filename=None, msh_filename=None):
    s = """
    SetFactory("OpenCASCADE");
    """

    s += f"Mesh.CharacteristicLengthMin = {size};"
    s += f"Mesh.CharacteristicLengthMax = {size};"

    s += f"""
    Circle(1) = {{0, 0, 0, {radius}, 0, 2*Pi}};
    Line Loop(1) = {{1}};
    Surface(1) = {{1}};
    """

    s+= f"""
    ex1[] = Extrude {{0,0,{length}}} {{Surface{{1}};}};
    
    Physical Surface("base") = {{1}};
    Physical Surface("top") = {{ex1[0]}};
    Physical Surface ("middle") = {{ex1[2]}};
    """

    if dimension == 2:
        s += "Physical Surface (\"surface\") = {ex1[0], ex1[1], ex1[2]};"
    else:
        s += "Physical Volume (\"volume\") = {ex1[1]};"

    s += f"""
    Mesh.CharacteristicLengthFromPoints = 1;
    Mesh.ElementOrder = {order};
    Mesh.SecondOrderLinear = 1;
    Mesh {dimension};
    """

    print(s)

    return geo_to_msh(s, dimension=dimension, verbose=verbose)


if __name__ == "__main__":
    main()
