import numpy as np
import ufl

from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh
from dolfinx import nls
from dolfinx.io import gmshio, XDMFFile

L, W, H = 80.0, 7.5, 7.5
elements = ["tetra", "hexa"]
orders = [1, 2]
materials = ["NeoHooke", "SaintVenantKirchhoff"]
YOUNG_MODULUS = 3000
POISSON_RATIO = 0.3
DENSITY = 3
GRAVITY = -0.1

for element in elements:
    for order in orders:
        for material in materials:
            filename = "./meshes/cube_" + element + "_1"
            domain, cell_tags, facet_tag = gmshio.read_from_msh(filename + ".msh", MPI.COMM_WORLD)

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

            # Elasticity parameters
            young = PETSc.ScalarType(YOUNG_MODULUS)
            poisson = PETSc.ScalarType(POISSON_RATIO)
            mu = fem.Constant(domain, young / (2 * (1 + poisson)))
            lmbda = fem.Constant(domain, young * poisson / ((1 + poisson) * (1 - 2 * poisson)))

            # Stored strain energy density (compressible neo-Hookean model)
            if material=="NeoHooke":
                Ic = ufl.variable(ufl.tr(C))
                J = ufl.variable(ufl.det(F))
                psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
            else:
                E = ufl.variable(0.5 * (C - I))
                psi = (lmbda / 2) * ufl.tr(E) ** 2 + mu * ufl.tr(E * E)

            P = ufl.diff(psi, F)

            metadata = {"quadrature_degree": 2}
            ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
            dx = ufl.Measure("dx", domain=domain, metadata=metadata)

            F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)

            problem = fem.petsc.NonlinearProblem(F, u, bcs)

            solver = nls.petsc.NewtonSolver(domain.comm, problem)

            # Set Newton solver options
            solver.atol = 1e-05
            solver.rtol = 1e-10
            solver.convergence_criterion = "incremental"

            # single loading
            # tval0 = -10
            # T.value[1] = tval0
            bval0 = GRAVITY * DENSITY
            B.value[1] = bval0
            num_its, converged = solver.solve(u)
            assert (converged)
            u.x.scatter_forward()


            with XDMFFile(MPI.COMM_WORLD, "./meshes/fenics_"+ material + element + str(order) + "_solution.xdmf",
                                     "w") as xdmf:
                xdmf.write_mesh(domain)
                xdmf.write_function(u, 0)

