import numpy as np
import ufl

from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh, plot
from gmsh_helpers import read_from_msh

import gmsh

gmsh.initialize()

L, W, H = 80.0, 7.5, 7.5
element = "hexa"
order = 2
# filename = "../meshes/fenics/cube_" + element + "_1"
# domain, cell_tags, facet_tag = read_from_msh(filename + ".msh", cell_data=True, facet_data=True, gdim=3)
# domain = mesh.create_box(MPI.COMM_WORLD, [[0, 0, 0.0], [1, 1, 1]], [2, 2, 2], mesh.CellType.hexahedron)
domain = mesh.create_unit_cube(MPI.COMM_WORLD, 1, 1, 1, mesh.CellType.hexahedron)

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
Ic = ufl.variable(ufl.tr(C))
J = ufl.variable(ufl.det(F))

# Elasticity parameters
E = PETSc.ScalarType(3000)
nu = PETSc.ScalarType(0.3)
mu = fem.Constant(domain, E / (2 * (1 + nu)))
lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
# Stored strain energy density (compressible neo-Hookean model)
psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
P = ufl.diff(psi, F)

metadata = {"quadrature_degree": 2}
ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
dx = ufl.Measure("dx", domain=domain, metadata=metadata)

F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)

problem = fem.petsc.NonlinearProblem(F, u, bcs)

# from dolfinx.jit import ffcx_jit
# import cffi
# import numba
# import numba.core.typing.cffi_utils as cffi_support
# from dolfinx.cpp.fem import Form_complex128, Form_float64
# from dolfinx.fem import (Function, FunctionSpace, IntegralType, dirichletbc,
#                          form, locate_dofs_topological)
#
# # JIT compile individual blocks tabulation kernels
# nptype = "complex128" if np.issubdtype(PETSc.ScalarType, np.complexfloating) else "float64"
# ffcxtype = "double _Complex" if np.issubdtype(PETSc.ScalarType, np.complexfloating) else "double"
# ufcx_form00, _, _ = ffcx_jit(domain.comm, F, form_compiler_params={"scalar_type": ffcxtype})
# kernel00 = getattr(ufcx_form00.integrals(0)[0], f"tabulate_tensor_{nptype}")
# Usize = V.element.space_dimension
#
# ffi = cffi.FFI()
# cffi_support.register_type(ffi.typeof('double _Complex'), numba.types.complex128)
# c_signature = numba.types.void(
#     numba.types.CPointer(numba.typeof(PETSc.ScalarType())),
#     numba.types.CPointer(numba.typeof(PETSc.ScalarType())),
#     numba.types.CPointer(numba.typeof(PETSc.ScalarType())),
#     numba.types.CPointer(numba.types.double),
#     numba.types.CPointer(numba.types.int32),
#     numba.types.CPointer(numba.types.uint8))
#
#
# @numba.jit
# def my_print(a):
#     print(a)
#
# @numba.cfunc(c_signature, nopython=True)
# def tabulate_condensed_tensor_A(A_, w_, c_, coords_, entity_local_index, permutation=ffi.NULL):
#     # Tabulate all sub blocks locally
#     A00 = np.zeros((Usize, Usize), dtype=PETSc.ScalarType)
#     coords = numba.carray(coords_, (Usize,), dtype=PETSc.ScalarType)
#     my_print(coords.shape)
#     # my_print(permutation[0])
#     # print(A00)
#     kernel00(ffi.from_buffer(A00), w_, c_, coords_, entity_local_index, permutation)
#
# # @numba.cfunc(c_signature, nopython=True)
# # def tabulate_condensed_tensor_A(A_, w_, c_, coords_, entity_local_index, permutation=ffi.NULL):
# #     # Tabulate all sub blocks locally
# #     A00 = np.zeros((Usize, Usize), dtype=PETSc.ScalarType)
# #     print(A00)
# #     kernel00(ffi.from_buffer(A00), w_, c_, coords_, entity_local_index, permutation)
#
#
#
# # Prepare a Form with a condensed tabulation kernel
# Form = Form_float64 if PETSc.ScalarType == np.float64 else Form_complex128
#
# integrals = {IntegralType.cell: ([(-1, tabulate_condensed_tensor_A.address)], None)}
# a_cond = Form([V._cpp_object, V._cpp_object], integrals, [], [], False, None)
# fem.petsc.assemble_matrix(a_cond)


from dolfinx import nls

solver = nls.petsc.NewtonSolver(domain.comm, problem)

# Set Newton solver options
solver.atol = 1e-05
solver.rtol = 1e-10
solver.convergence_criterion = "incremental"




#
# from dolfinx import log
#
# log.set_log_level(log.LogLevel.INFO)
#
# # single loading
# tval0 = -10
# T.value[1] = tval0
# num_its, converged = solver.solve(u)
# assert (converged)
# u.x.scatter_forward()
#
# # # incremental loading
# # tval0 = -4000 / (H * W)
# # for n in range(1, 10):
# #     T.value[1] = n / 10 * tval0
# #     num_its, converged = solver.solve(u)
# #     assert (converged)
# #     u.x.scatter_forward()
#
# import dolfinx.io
#
# # with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "initial.xdmf", "w") as xdmf:
# #     xdmf.write_mesh(domain)
#
# with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "../meshes/fenics/cube_" + element + "_" + str(order) + "_solution.xdmf",
#                          "w") as xdmf:
#     xdmf.write_mesh(domain)
#     xdmf.write_function(u, 0)
#
# # with dolfinx.io.VTKFile(MPI.COMM_WORLD, "output.pvd", "w") as vtk:
# #     vtk.write_function(u, 0.)
