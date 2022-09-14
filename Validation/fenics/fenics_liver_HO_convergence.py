import numpy as np
import ufl

from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh, plot
from gmsh_helpers import read_from_msh

import gmsh

gmsh.initialize()
files = np.arange(1, 17)
for file in files:
    element = "tetra"
    order = 1
    filename = "./liver-fenics-"+str(file)
    domain, cell_tags, facet_tag = read_from_msh(filename + ".msh", cell_data=True, facet_data=True, gdim=3)

    if element == "hexa":
        V = fem.VectorFunctionSpace(domain, ("S", order))
    else:
        V = fem.VectorFunctionSpace(domain, ("Lagrange", order))


    def left(x):
        return np.greater(x[0], -1)


    def right(x):
        return np.less(x[0], -3)


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
    #
    # # Identity tensor
    # I = ufl.variable(ufl.Identity(d))
    #
    # # Deformation gradient
    # F = ufl.variable(I + ufl.grad(u))
    #
    # # Right Cauchy-Green tensor
    # C = ufl.variable(F.T * F)
    #
    # # Invariants of deformation tensors
    # Ic = ufl.variable(ufl.tr(C))
    # J = ufl.variable(ufl.det(F))
    #
    # # Elasticity parameters
    # E = PETSc.ScalarType(3000)
    # nu = PETSc.ScalarType(0.3)
    # mu = fem.Constant(domain, E / (2 * (1 + nu)))
    # lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
    # # Stored strain energy density (compressible neo-Hookean model)
    # psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    # P = ufl.diff(psi, F)

    # Kinematics
    I = ufl.Identity(d)  # Identity tensor
    F = ufl.variable(I + ufl.grad(u))  # Deformation gradient
    C = ufl.variable(F.T * F)  # Right Cauchy-Green tensor
    J = ufl.det(F)
    I1 = ufl.tr(C)

    # Elasticity parameters
    bulk_modulus = fem.Constant(domain, PETSc.ScalarType(10e8))
    # a = fem.Constant(domain, PETSc.ScalarType(1180))
    # b = fem.Constant(domain, PETSc.ScalarType(8))
    # a_f = fem.Constant(domain, PETSc.ScalarType(18.5 * 10e4))
    # b_f = fem.Constant(domain, PETSc.ScalarType(16))
    # a_s = fem.Constant(domain, PETSc.ScalarType(2.5 * 10e4))
    # b_s = fem.Constant(domain, PETSc.ScalarType(11.1))
    # a_fs = fem.Constant(domain, PETSc.ScalarType(2160))
    # b_fs = fem.Constant(domain, PETSc.ScalarType(11.4))
    a = fem.Constant(domain, PETSc.ScalarType(1e6))
    b = fem.Constant(domain, PETSc.ScalarType(5))
    a_f = fem.Constant(domain, PETSc.ScalarType(16e4))
    b_f = fem.Constant(domain, PETSc.ScalarType(12.8))
    a_s = fem.Constant(domain, PETSc.ScalarType(18e4))
    b_s = fem.Constant(domain, PETSc.ScalarType(10))
    a_fs = fem.Constant(domain, PETSc.ScalarType(9e3))
    b_fs = fem.Constant(domain, PETSc.ScalarType(12))

    f_0 = ufl.as_vector([0.0, 1.0 / ufl.sqrt(2), 1.0 / ufl.sqrt(2)])
    s_0 = ufl.as_vector([0.0, 1.0 / ufl.sqrt(2), -1.0 / ufl.sqrt(2)])

    I_4_f_0 = ufl.dot(C * f_0, f_0)
    I_4_s_0 = ufl.dot(C * s_0, s_0)
    I_8_f_0_s_0 = ufl.dot(C * s_0, f_0)

    W_vol = bulk_modulus * (J ** 2 - 1.0 - 2.0 * ufl.ln(J)) / 4.0
    W_1 = a * ufl.exp((I1 - 3.0) * b) / (2.0 * b)

    I_4_f_0 = I_4_f_0 - 1.0
    I_4_s_0 = I_4_s_0 - 1.0

    W_4f = a_f * (ufl.exp(((I_4_f_0) ** 2) * b_f) - 1.0) / (2.0 * b_f)
    W_4s = a_s * (ufl.exp(((I_4_s_0) ** 2) * b_s) - 1.0) / (2.0 * b_s)

    W_8_fs = a_fs * (ufl.exp((I_8_f_0_s_0 ** 2) * b_fs) - 1.0) / (2.0 * b_fs)

    # Stored strain energy density (Ogden)
    psi = W_vol + W_1 + W_4f + W_4s + W_8_fs

    P = ufl.diff(psi, F)

    metadata = {"quadrature_degree": 2}
    ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)

    F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)
    #
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
    tval0 = -1000
    T.value[1] = tval0
    num_its, converged = solver.solve(u)
    assert (converged)
    u.x.scatter_forward()
    a = u.x.array.shape[0]
    print(np.min(u.x.array.reshape((int(a/3), 3))[:, 2]))

    with open('./convergence_study/liver/displacement_' + str(u.x.array.shape[0]) + ".txt", 'w') as f:
        f.write(str(np.min(u.x.array.reshape((int(a/3), 3))[:, 2])))

    # print(np.max(u.x.array))

    # with open('./convergence_study/displacement_' + str(u.x.array.shape[0]) + ".txt", 'w') as f:
    #     f.write(str(np.min(u.x.array.reshape(((nx + 1) * (ny + 1) * (nz + 1), 3))[:, 2])))
    # import dolfinx.io
    #
    # #
    # # # with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "initial.xdmf", "w") as xdmf:
    # # #     xdmf.write_mesh(domain)
    # #
    # with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "./liver_solution.xdmf",
    #                          "w") as xdmf:
    #     xdmf.write_mesh(domain)
    #     xdmf.write_function(u, 0)

    # with dolfinx.io.VTKFile(MPI.COMM_WORLD, "output.pvd", "w") as vtk:
    #     vtk.write_function(u, 0.)
