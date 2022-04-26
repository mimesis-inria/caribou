#!/usr/bin/python3

"""
Manufactured solution by explicitly setting the solution displacement functions to:

                | 1e-3 * z *exp(x) |
  u (x, y, z) = | 1e-3 * z *exp(y) |
                | 1e-3 * z *exp(z) |

with (x, y, z) material coordinates (undeformed, at rest). The external forces applied to the domain becomes:

  ext = integrate( Div(P), volume_domain) + integrate( P . N , neummann_domain)

where P is the first Piola-Kirchhoff stress tensor, N is the surface normal at rest. The geometry of the volumetric
domain is a cylinder having a radius of 1 and a length of 3.

To run with docker :
  docker run --rm -v $PWD:/w -w /w jnbrunet/caribou_validation python manufactured_solution.py
"""

import Sofa, SofaRuntime, SofaCaribou
import meshio, numpy as np
from manufactured_solution import assemble, integrate, compute_solution, ConstantForceField

order = "quadratic"
material = "SaintVenantKirchhoff"

if order == "linear":
    mesh = meshio.read('meshes/new_cylinder_p1.vtu')
    surface_template = "Triangle"
    volume_template = "Tetrahedron"
    fenics_volume_indices = mesh.cells[1].data
    length = 80
else:
    mesh = meshio.read('meshes/new_cylinder_p2.vtu')
    surface_template = "Triangle6"
    volume_template = "Tetrahedron10"
    fenics_volume_indices = mesh.cells[1].data[:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
    length = 80

mu = 1.0
l = 1.25
rad = 1
poisson_ratio = 1. / (2 * ((mu / l) + 1))
young_modulus = 2 * mu * (1 + poisson_ratio)
P, f, u_s = compute_solution(mu, l, rad, length)


def create_scene(root):
    root.addObject('RequiredPlugin',
                   pluginName='SofaBaseMechanics SofaBoundaryCondition SofaImplicitOdeSolver SofaSparseSolver')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showBehaviorModels showForceFields')

    sofa_node = root.addChild("sofa_node")
    sofa_node.addObject('StaticODESolver', newton_iterations=10, residual_tolerance_threshold=1e-10, printLog=True)
    sofa_node.addObject('LDLTSolver', backend='Pardiso')
    sofa_node.addObject('MechanicalObject', name='mo', position=mesh.points.tolist())
    sofa_node.addObject('CaribouTopology', name='volume', template=volume_template, indices=mesh.cells[1].data.tolist())
    sofa_node.addObject('CaribouTopology', name='dirichlet_boundary', template=surface_template,
                        indices=mesh.cells[0].data[
                            np.ma.masked_equal(mesh.cell_data['gmsh:physical'][0], 1).mask].tolist())
    sofa_node.addObject('CaribouTopology', name='neumann_boundary', template=surface_template,
                        indices=mesh.cells[0].data[
                            np.ma.masked_inside(mesh.cell_data['gmsh:physical'][0], 2, 3).mask].tolist())

    sofa_node.addObject(material + "Material", young_modulus=young_modulus, poisson_ratio=poisson_ratio)
    sofa_node.addObject('HyperelasticForcefield', topology='@volume', printLog=True)

    sofa_node.addObject('FixedConstraint',
                        indices=np.unique(root.sofa_node.dirichlet_boundary.indices.array()).tolist())
    sofa_node.addObject(ConstantForceField(name='external_forces_sofa'))

    fenics_node = root.addChild("fenics_node")
    fenics_node.addObject('StaticODESolver', newton_iterations=10, residual_tolerance_threshold=1e-10, printLog=True)
    fenics_node.addObject('LDLTSolver', backend='Pardiso')
    fenics_node.addObject('MechanicalObject', name='mo', position=mesh.points.tolist())
    fenics_node.addObject('CaribouTopology', name='volume', template=volume_template,
                          indices=fenics_volume_indices.tolist())
    fenics_node.addObject('CaribouTopology', name='dirichlet_boundary', template=surface_template,
                          indices=mesh.cells[0].data[
                              np.ma.masked_equal(mesh.cell_data['gmsh:physical'][0], 1).mask].tolist())
    fenics_node.addObject('CaribouTopology', name='neumann_boundary', template=surface_template,
                          indices=mesh.cells[0].data[
                              np.ma.masked_inside(mesh.cell_data['gmsh:physical'][0], 2, 3).mask].tolist())
    fenics_node.addObject('FEniCS_Material', template=volume_template, young_modulus=young_modulus,
                          poisson_ratio=poisson_ratio, material_name=material)
    fenics_node.addObject('HyperelasticForcefield_FEniCS', topology='@volume', printLog=True)

    fenics_node.addObject('FixedConstraint',
                          indices=np.unique(root.fenics_node.dirichlet_boundary.indices.array()).tolist())
    fenics_node.addObject(ConstantForceField(name='external_forces_fenics'))


if __name__ == '__main__':
    root = Sofa.Core.Node()
    create_scene(root)
    Sofa.Simulation.init(root)

    print('Assembling the external force vector...', end='', flush=True)
    external_forces_sofa = \
        assemble(root.sofa_node.volume.domain(), lambda x, y, z, _: f(x, y, z)) + \
        assemble(root.sofa_node.neumann_boundary.domain(), lambda x, y, z, t: np.dot(P(x, y, z), t.normal()))
    external_forces_fenics = \
        assemble(root.fenics_node.volume.domain(), lambda x, y, z, _: f(x, y, z)) + \
        assemble(root.fenics_node.neumann_boundary.domain(), lambda x, y, z, t: np.dot(P(x, y, z), t.normal()))
    print(' Done.', flush=True)

    exact_error_sofa = np.sqrt(
        integrate(root.sofa_node.volume.domain(), lambda x, y, z, _: np.dot(u_s(x, y, z), u_s(x, y, z))))
    exact_error_fenics = np.sqrt(
        integrate(root.fenics_node.volume.domain(), lambda x, y, z, _: np.dot(u_s(x, y, z), u_s(x, y, z))))

    print(f"Difference for the exact error is {np.abs(exact_error_sofa - exact_error_fenics)}")

    for load in [1e-3, 1e-2, 1e-1, 0.15, 0.5, 1.0]:
        root.sofa_node.external_forces_sofa.forces = (external_forces_sofa * load)
        root.fenics_node.external_forces_fenics.forces = (external_forces_fenics * load)
        Sofa.Simulation.animate(root, 1)
        u_h_sofa = (root.sofa_node.mo.position.array() - root.sofa_node.mo.rest_position.array())
        error_L2_sofa = np.sqrt(integrate(
            root.sofa_node.volume.domain(),
            lambda x, y, z, u_g, _: np.dot((u_s(x, y, z) - u_g), (u_s(x, y, z) - u_g)),
            u_h_sofa
        ))
        print(f"SOFA relative L2 error at {int(load * 100)}% load: {error_L2_sofa / exact_error_sofa}")

        u_h_fenics = (root.fenics_node.mo.position.array() - root.fenics_node.mo.rest_position.array())
        error_L2_fenics = np.sqrt(integrate(
            root.fenics_node.volume.domain(),
            lambda x, y, z, u_g, _: np.dot((u_s(x, y, z) - u_g), (u_s(x, y, z) - u_g)),
            u_h_fenics
        ))
        print(f"FEniCS relative L2 error at {int(load * 100)}% load: {error_L2_fenics / exact_error_fenics}")
        print(
            f"Difference in relative L2 error at {int(load * 100)}% load: {np.abs((error_L2_fenics / exact_error_fenics) - (error_L2_sofa / exact_error_sofa))}")
