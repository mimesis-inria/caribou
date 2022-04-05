import meshio
import numpy as np


output = 0
while output != '1' and output != '2':  
    print('Choose one of the bellow materials: ')
    print('1 - SaintVenantKirchhoff')
    print('2 - Neo-Hookean')
    output = input()

if output == '1': 
    sonics_material = 'SaintVenantKirchhoffMaterial_FEniCS'
    febio_material = 'isotropic elastic'
else:
    sonics_material = 'NeoHookeanMaterial_FEniCS'
    febio_material = 'neo-Hookean'

output = 0
while int(output) < 1 or int(output) > 4:  
    print('Choose one of the bellow elements: ')
    print('1 - Tetrahedron')
    print('2 - Tetrahedron10')
    print('3 - Hexahedron')
    print('4 - Hexahedron20')
    output = input()

if output == '1': 
    element = 'Tetrahedron'
    meshfile = "../../meshes/beam_p1.vtu"
    mesh = meshio.read(meshfile)
    indices = np.empty(mesh.cells_dict['tetra'].shape)
    indices = mesh.cells_dict['tetra']
elif output == '2': 
    element = 'Tetrahedron10'
    meshfile = "../../meshes/beam_p2.vtu"
    mesh = meshio.read(meshfile)
    indices = np.empty(mesh.cells_dict['tetra10'].shape)
    indices = mesh.cells_dict['tetra10'][:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
elif output == '3': 
    element = 'Hexahedron_FEniCS'
    meshfile = "../../meshes/beam_q1.vtu"
    mesh = meshio.read(meshfile)
    indices = np.empty(mesh.cells_dict['hexahedron'].shape)
    indices = mesh.cells_dict['hexahedron'][:, [4, 5, 0, 1, 7, 6, 3, 2]]
else:
    element = 'Hexahedron_FEniCS20'
    meshfile = "../../meshes/beam_q2.vtu"
    mesh = meshio.read(meshfile)
    indices = np.empty(mesh.cells_dict['hexahedron20'].shape)
    indices = mesh.cells_dict['hexahedron20'][:,
              [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]




parameters = {
    'element': element,
    'indices': indices,
    'sonics_material': sonics_material, 
    'febio_material': febio_material,
    'meshfile': meshfile,
    'mesh': mesh,
    'young_modulus': 3000, 
    'poisson_ratio': 0.3, 
    'traction': np.array([0., -100, 0.]), 
    'residual_tolerance':"1e-10",
    'displacement_tolerance': "1e-10",
    'nnewtonsteps': 25,
    'nBFGSsteps': 0,
    'nsteps': 10,
    'bas_id': 1,  # Id of surface elements forming the base of the beam. Must be defined in the gmsh:physical Cell Data attribute.
    'top_id': 2,  # Id of surface elements forming the top of the beam. Must be defined in the gmsh:physical Cell Data attribute.
    'vol_id': 4,  # Id of volumentric elements forming the top of the beam. Must be defined in the gmsh:physical Cell Data attribute.
    'displacement_filename': './febio_displacements.txt',
}

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

