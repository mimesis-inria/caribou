import numpy as np
import ufl

from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh, plot

# from gmsh_helpers import read_from_msh
#
# import gmsh
#
# gmsh.initialize()
#
# mesh, cell_tags, facet_tags = read_from_msh("../meshes/p1_from_gmsh.msh", cell_data=True, facet_data=True, gdim=3)

L, W, H = 80.0, 7.5, 7.5
domain = mesh.create_box(MPI.COMM_WORLD, [[-H, -W, 0.0], [H, W, L]], [2, 2, 8], mesh.CellType.tetrahedron)
V = fem.VectorFunctionSpace(domain, ("Lagrange", 1))


def left(x):
    return np.isclose(x[2], 0)


def right(x):
    return np.isclose(x[2], L)


fdim = domain.topology.dim - 1
left_facets = mesh.locate_entities_boundary(domain, fdim, left)
right_facets = mesh.locate_entities_boundary(domain, fdim, right)

# Concatenate and sort the arrays based on facet indices. Left facets marked with 1, right facets with two
marked_facets = np.hstack([left_facets, right_facets])
marked_values = np.hstack([np.full(len(left_facets), 1, dtype=np.int32), np.full(len(right_facets), 2, dtype=np.int32)])
sorted_facets = np.argsort(marked_facets)
facet_tag = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

u_bc = np.array((0,) * domain.geometry.dim, dtype=PETSc.ScalarType)

left_dofs = fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.indices[facet_tag.values == 1])
bcs = [fem.dirichletbc(u_bc, left_dofs, V)]

B = fem.Constant(domain, PETSc.ScalarType((0, 0, 0)))
T = fem.Constant(domain, PETSc.ScalarType((0, 0, 0)))
v = ufl.TestFunction(V)
u = fem.Function(V)

# Spatial dimension
d = len(u)

# Identity tensor
I = ufl.variable(ufl.Identity(d))

# Deformation gradient
F = ufl.variable(I + ufl.grad(u))

# Right Cauchy-Green tensor
C = ufl.variable(F.T * F)

# Invariants of deformation tensors
Ic = ufl.variable(ufl.tr(C))
J = ufl.variable(ufl.det(F))

# Elasticity parameters
E = PETSc.ScalarType(3000)
nu = PETSc.ScalarType(0.3)
mu = fem.Constant(domain, E / (2 * (1 + nu)))
lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
# Stored strain energy density (compressible neo-Hookean model)
psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
# Stress
# Hyper-elasticity
P = ufl.diff(psi, F)

metadata = {"quadrature_degree": 1}
ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
dx = ufl.Measure("dx", domain=domain, metadata=metadata)
F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)

problem = fem.petsc.NonlinearProblem(F, u, bcs)

from dolfinx import nls

solver = nls.petsc.NewtonSolver(domain.comm, problem)

# Set Newton solver options
solver.atol = 1e-05
solver.rtol = 1e-10
solver.convergence_criterion = "incremental"

from dolfinx import log

log.set_log_level(log.LogLevel.INFO)
tval0 = -4000 / (H * W)
# T.value[1] = tval0
for n in range(1, 10):
    T.value[1] = n / 10 * tval0
    num_its, converged = solver.solve(u)
    assert (converged)
    u.x.scatter_forward()

import dolfinx.io

with dolfinx.io.VTKFile(MPI.COMM_WORLD, "output.pvd", "w") as vtk:
    vtk.write_function(u, 0.)

displacement = u.x.array.reshape(int(u.x.array.shape[0] / domain.geometry.dim), domain.geometry.dim)
print(displacement[:, 1].min())
