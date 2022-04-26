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

if order == "linear":
    mesh = meshio.read('meshes/new_cylinder_p1.vtu')
    surface_template = "Triangle"
    volume_template = "Tetrahedron"
else:
    mesh = meshio.read('meshes/new_cylinder_p2.vtu')
    surface_template = "Triangle6"
    volume_template = "Tetrahedron10"

mu = 1.0
l = 1.25
length = 80
rad = 7.5
poisson_ratio = 1. / (2 * ((mu / l) + 1))
young_modulus = 2 * mu * (1 + poisson_ratio)
P, f, u_s = compute_solution(mu, l, rad, length)


def create_scene(root):
    root.addObject('RequiredPlugin',
                   pluginName='SofaBaseMechanics SofaBoundaryCondition SofaImplicitOdeSolver SofaSparseSolver')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showBehaviorModels showForceFields')
    root.addObject('StaticODESolver', newton_iterations=10, residual_tolerance_threshold=1e-10, printLog=True)
    root.addObject('LDLTSolver', backend='Pardiso')
    root.addObject('MechanicalObject', name='mo', position=mesh.points.tolist())
    root.addObject('CaribouTopology', name='volume', template=volume_template, indices=mesh.cells[1].data.tolist())
    root.addObject('CaribouTopology', name='dirichlet_boundary', template=surface_template,
                   indices=mesh.cells[0].data[np.ma.masked_equal(mesh.cell_data['gmsh:physical'][0], 1).mask].tolist())
    root.addObject('CaribouTopology', name='neumann_boundary', template=surface_template,
                   indices=mesh.cells[0].data[
                       np.ma.masked_inside(mesh.cell_data['gmsh:physical'][0], 2, 3).mask].tolist())

    root.addObject('SaintVenantKirchhoffMaterial', young_modulus=young_modulus, poisson_ratio=poisson_ratio)
    root.addObject('HyperelasticForcefield', topology='@volume', printLog=True)

    root.addObject('FixedConstraint', indices=np.unique(root.dirichlet_boundary.indices.array()).tolist())
    root.addObject(ConstantForceField(name='external_forces'))


if __name__ == '__main__':
    root = Sofa.Core.Node()
    create_scene(root)
    Sofa.Simulation.init(root)

    print('Assembling the external force vector...', end='', flush=True)
    external_forces = \
        assemble(root.volume.domain(), lambda x, y, z, _: f(x, y, z)) + \
        assemble(root.neumann_boundary.domain(), lambda x, y, z, t: np.dot(P(x, y, z), t.normal()))
    print(' Done.', flush=True)

    exact_error = np.sqrt(integrate(root.volume.domain(), lambda x, y, z, _: np.dot(u_s(x, y, z), u_s(x, y, z))))
    print(f"Exact error is {exact_error}")

    for load in [1e-3, 1e-2, 1e-1, 0.15, 0.5, 1.0]:
        root.external_forces.forces = (external_forces * load)
        Sofa.Simulation.animate(root, 1)
        u_h = (root.mo.position.array() - root.mo.rest_position.array())
        error_L2 = np.sqrt(integrate(
            root.volume.domain(),
            lambda x, y, z, u_g, _: np.dot((u_s(x, y, z) - u_g), (u_s(x, y, z) - u_g)),
            u_h
        ))
        print(f"relative L2 error at {int(load * 100)}% load: {error_L2 / exact_error}")
