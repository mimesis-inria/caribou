import numpy as np
import ufl

from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh, plot
import dolfinx.io

from gmsh_helpers import read_from_msh

import gmsh

gmsh.initialize()

L, W, H = 80.0, 7.5, 7.5
element = "hexa"
order = 1
filename = "../meshes/fenics/cube_" + element + "_1"
domain, cell_tags, facet_tag = read_from_msh(filename + ".msh", cell_data=True, facet_data=True, gdim=3)
# domain = mesh.create_box(MPI.COMM_WORLD, [[-H, -W, 0.0], [H, W, L]], [2, 2, 8], mesh.CellType.hexahedron)
# with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "./beam_q1.xdmf", "r") as xdmf:
#     domain = xdmf.read_mesh(name="Grid")
#     ct = xdmf.read_meshtags(mesh, name="Grid")

if element == "hexa":
    V = fem.VectorFunctionSpace(domain, ("S", order))
else:
    V = fem.VectorFunctionSpace(domain, ("Lagrange", order))


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

left_dofs = fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.indices[facet_tag.values == 1])
if element == "hexa" and order == 2:
    bcs = [fem.dirichletbc(fem.Function(V), left_dofs)]
else:
    u_bc = np.array((0,) * domain.geometry.dim, dtype=PETSc.ScalarType)
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
I_C = ufl.variable(ufl.tr(C))
II_C = ufl.variable(0.5*(ufl.tr(C)**2 - ufl.tr(C*C)))

J = ufl.variable(ufl.det(F))

# Elasticity parameters
C_01 = PETSc.ScalarType(200)
C_10 = PETSc.ScalarType(100)
K = PETSc.ScalarType(1000)

# stored strain energy density (nearly incompressible Mooney-Rivlin model)
psi = C_01*(J**(-2.0/3.0)*I_C - 3) + C_10*(J**(-4.0/3.0)*II_C - 3) + 0.5 * K * (ufl.ln(J))**2
P = ufl.diff(psi, F)

metadata = {"quadrature_degree": 2}
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

# single loading
tval0 = -9
T.value[1] = tval0
num_its, converged = solver.solve(u)
assert (converged)
u.x.scatter_forward()

# # incremental loading
# tval0 = -4000 / (H * W)
# for n in range(1, 10):
#     T.value[1] = n / 10 * tval0
#     num_its, converged = solver.solve(u)
#     assert (converged)
#     u.x.scatter_forward()

# with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "./test.xdmf", "w") as xdmf:
#     xdmf.write_mesh(domain)

with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "./cube_" + element + "_" + str(order) + "_solution.xdmf",
                         "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(u, 0)


with dolfinx.io.VTKFile(MPI.COMM_WORLD, "./test.pvd", "w") as vtk:
    vtk.write_function(u, 0.)
