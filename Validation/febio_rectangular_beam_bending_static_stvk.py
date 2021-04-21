#!/usr/bin/python3.7

"""
Bending rectangular beam simulated with FEBio.

Dimensions: 15x15x80
Material: St-Venant-Kirchhoff (young modulus 3000, poisson ratio 0.499)
ODE solver: Static Newton-Raphson
Left face: Clamped
Right face: Traction of [0, -30, 0]
Number of load increments: 5

To run with podman:
podman run --rm --userns keep-id \
    -v $PWD:/opt/shared:z \
    -w /opt/shared \
    febio_validation python febio_rectangular_beam_bending_static_stvk.py
"""

import numpy as np
import os, meshio, tempfile, subprocess

# Parameters
residual_tolerance = "1e-10"
displacement_tolerance = "1e-10"
nsteps = 5
nnewtonsteps = 30
nBFGSsteps = 0
young_modulus = 3000
poisson_ratio = 0.499
traction = np.array([0., -30, 0.])

# Mesh
in_msh_file = "meshes/beam_q2.vtu"
sol_vtk_file = "meshes/beam_q2_solution.vtu"
bas_id = 1  # Id of surface elements forming the base of the beam. Must be defined in the gmsh:physical Cell Data attribute.
top_id = 2  # Id of surface elements forming the top of the beam. Must be defined in the gmsh:physical Cell Data attribute.
vol_id = 4  # Id of volumentric elements forming the top of the beam. Must be defined in the gmsh:physical Cell Data attribute.
mesh = meshio.read(in_msh_file)

meshio_to_febio_element_types = {
    'triangle': 'tri3',
    'triangle6': 'tri6',
    'quad': 'quad4',
    'quad8': 'quad8',
    'tetra': 'tet4',
    'tetra10': 'tet10',
    'hexahedron': 'hex8',
    'hexahedron20': 'hex20'
}

# Sort element indices per tag and type
elements = {}
for block, tags in zip(mesh.cells, mesh.cell_data['gmsh:physical']):
    type, indices = block
    type = meshio_to_febio_element_types[type]
    for t in np.unique(tags):
        if t not in elements:
            elements[t] = {}
        mask = np.array(np.ma.masked_not_equal(tags, t).mask)
        if len(indices[~mask]):
            if type not in elements[t]:
                elements[t][type] = indices[~mask]
            else:
                elements[t][type] = np.concatenate((indices[~mask], elements[t][type]), axis=0)

for tag, types in elements.items():
    if len(types) > 1:
        print(f'Tag "{tag}" has "')
        for type, indices in types.items():
            print(f'\t{len(indices)} elements of type "{type}"')
    else:
        type, indices = list(types.items())[0]
        print(f'Tag "{tag}" has {len(indices)} elements of type "{type}"')


print("Mesh has {} nodes".format(len(mesh.points)))

with tempfile.NamedTemporaryFile(suffix=".txt") as f:
    displacement_filename = f.name

# Create the FEBio spec file
c = []
c.append(f'<?xml version="1.0" encoding="UTF-8"?>')
c.append(f'<febio_spec version="2.0">')
c.append(f'  <Module type="solid"/>')
c.append(f'  <Control>')
c.append(f'    <title>Rectangular beam bending</title>')
c.append(f'    <time_steps>{nsteps}</time_steps>')
c.append(f'    <step_size>1</step_size>')
c.append(f'    <rtol>{residual_tolerance}</rtol>')
c.append(f'    <dtol>{displacement_tolerance}</dtol>')
c.append(f'    <max_ups>{nBFGSsteps}</max_ups>')
c.append(f'    <max_refs>{nnewtonsteps}</max_refs>')
c.append(f'    <lstol>0</lstol>')
c.append(f'    <min_residual>1e-15</min_residual>')
c.append(f'    <output_level>OUTPUT_FINAL</output_level>')
c.append(f'  </Control>')
c.append(f'  <Material>')
c.append(f'    <material id="1" type="isotropic elastic">')
c.append(f'      <E>{young_modulus}</E>')
c.append(f'      <v>{poisson_ratio}</v>')
c.append(f'    </material>')
c.append(f'  </Material>')

c.append(f'  <Geometry>')
c.append(f'    <Nodes name="nodes">')
for id, n in enumerate(mesh.points):
    c.append(f'      <node id="{id+1}">{",".join([str(i) for i in n])}</node>')
c.append(f'    </Nodes>')
c.append(f'    <NodeSet name="base_nodes">')
base_nodes = np.unique(np.concatenate([elements[bas_id][type] for type in elements[bas_id].keys()])).astype(int)
for n in base_nodes:
    c.append(f'      <node id="{n+1}"/>')
c.append(f'    </NodeSet>')
for type, indices in elements[vol_id].items():
    c.append(f'    <Elements type="{type}" mat="1">')
    for id, t in enumerate(indices[0].tolist()):
        c.append(f'      <elem id="{id+1}">{",".join([str(i+1) for i in t])}</elem>')
    c.append(f'    </Elements>')
c.append(f'    <Surface name="base">')
id = 0
for type, indices in elements[bas_id].items():
    for t in indices.tolist():
        c.append(f'      <{type} id="{id+1}">{",".join([str(i+1) for i in t])}</{type}>')
        id += 1
c.append(f'    </Surface>')
c.append(f'    <Surface name="top">')
id = 0
for type, indices in elements[top_id].items():
    for t in indices.tolist():
        c.append(f'      <{type} id="{id+1}">{",".join([str(i+1) for i in t])}</{type}>')
        id += 1
c.append(f'    </Surface>')
c.append(f'  </Geometry>')

c.append(f'  <Boundary>')
c.append(f'    <fix bc="xyz" set="base_nodes"/>')
c.append(f'  </Boundary>')

c.append(f'  <Loads>')
c.append(f'    <surface_load type="traction" extend="constant">')
c.append(f'      <surface set="top"/>')
c.append(f'      <scale lc="1">{np.linalg.norm(traction)}</scale>')
c.append(f'      <traction>{",".join([str(i) for i in traction/np.linalg.norm(traction)])}</traction>')
c.append(f'    </surface_load>')
c.append(f'  </Loads>')

c.append(f'  <LoadData>')
c.append(f'    <loadcurve id="1" type="linear">')
c.append(f'      <loadpoint>0,0</loadpoint>')
c.append(f'      <loadpoint>{nsteps},1</loadpoint>')
c.append(f'    </loadcurve>')
c.append(f'  </LoadData>')

c.append(f'  <Output>')
c.append(f'    <logfile file="{displacement_filename}">')
c.append(f'      <node_data data="ux;uy;uz"></node_data>')
c.append(f'    </logfile>')
c.append(f'  </Output>')

c.append(f'</febio_spec>')

with tempfile.NamedTemporaryFile(suffix=".feb") as f:
    feb_filename = f.name
with open(feb_filename, 'w') as f:
    f.write('\n'.join(c) + '\n')

print('Running febio...', end='')
args = [
    "-i",
    feb_filename
]
p = subprocess.Popen(
    ['febio3'] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
)

while True:
    line = p.stdout.readline()
    if not line:
        break
    # line = line.decode("utf-8")
    # print(line, end = '')

p.communicate()
assert p.returncode == 0, f"febio exited with error (return code {p.returncode})."
print('Done')
with open(displacement_filename, 'r') as f:
    lines = f.readlines()
    last_record_position = 0
    for id, l in zip(range(len(lines)),lines):
        if l.startswith("Data Record"):
            last_record_position = id

    start = last_record_position+5
    end = start + len(mesh.points)
    displacement_lines = lines[start:end]

    displacements = np.zeros((len(displacement_lines), 3))
    for id, l in zip(range(len(displacement_lines)),displacement_lines):
        l = l[:-1].split(' ')[1:]
        displacements[id] = np.array([l[0], l[1], l[2]]).astype(float)

assert len(displacements) == len(mesh.points)
out_mesh = meshio.Mesh(points=mesh.points, cells=mesh.cells, point_data={'u': displacements})
meshio.write(sol_vtk_file, out_mesh)
print(f'Solution mesh has been written to {sol_vtk_file}')

os.remove(feb_filename)
os.remove(displacement_filename)
