#!/usr/bin/python3.7

"""
Bending rectangular beam simulated with fenics.

Dimensions: 15x15x80
Material: St-Venant-Kirchhoff (young modulus 3000, poisson ratio 0.4999
ODE solver: Static Newton-Raphson
"""

from dolfin import *
import os

if os.getenv('FIX_DIJITSO', 0):
    import dijitso
    # Fix bug on podman
    default_params = dijitso.params.default_cache_params()
    default_params['temp_dir_root'] = "/home/fenics/.cache/fenics"
    def _params() :
        return default_params
    dijitso.params.default_cache_params = _params

# Parameters
radius = 7.5
length = 80
nsteps = 5
young_modulus = 3000
poisson_ratio = 0.499
traction = 30.
ncuts = 4
a_tol, r_tol = 1e-10, 1e-10
max_nb_of_newton_iterations = 10

# Load mesh
mesh = Mesh()
with XDMFFile('rectangular_beam_q1.xdmf') as infile:
    infile.read(mesh)

V = VectorFunctionSpace(mesh, 'P', degree=1)


# Sub domain for clamp at left end
class Left(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and near(x[2], 0, tol)


# Sub domain for rotation at right end
class Right(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and near(x[2], length, tol)


# Set up boundary condition at left end
bc = DirichletBC(V,  Constant((0.0, 0.0, 0.0)), Left())

# Create mesh function over the cell facets
boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
Right().mark(boundary_markers, 1)

# Define measure for boundary condition integral
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)


# Constant
mu    = Constant(young_modulus / (2.0*(1.0 + poisson_ratio)))
lmbda = Constant(young_modulus*poisson_ratio / ((1.0 + poisson_ratio)*(1.0 - 2.0*poisson_ratio)))

# Define variational problem
u = Function(V)
v = TrialFunction(V)
w = TestFunction(V)      # Test function
d = u.geometric_dimension()  # space dimension
f = Constant((0, 0, 0))
I = Identity(d)

for i in range(nsteps):
    print(f'Applying a traction of {-traction/nsteps*(i+1)}')
    T = Constant((0, -traction/nsteps*(i+1), 0))
    F = I + nabla_grad(u)
    C = F.T*F
    E = 0.5*(C-I)
    S = lmbda*tr(E)*I + 2*mu*E

    a = inner(S, derivative(E, u, w))*dx
    L = dot(f, w)*dx + dot(T, w)*ds(1)

    # Residual
    R = a - L

    # Jacobian of R
    J = derivative(R, u, v)

    # Newton-Raphson iterations
    du = Function(V)  # displacement increment
    U = Function(V)   # Total displacement since the beginning of the Newton iterations

    for newton_it in range(max_nb_of_newton_iterations):
        A, B = assemble_system(J, -R, [bc])
        solver = LUSolver(A, "mumps")
        solver.parameters["symmetric"] = True

        solver.solve(du.vector(), B)

        # Apply the solution increment
        u.vector()[:] += du.vector()

        # Update the residual
        B = assemble(R)
        bc.apply(B)
        f_norm = B.norm("l2")
        du_norm = du.vector().norm("l2")

        # Update norms
        if newton_it == 0:
            U.vector()[:] = du.vector()
            f0_norm = f_norm
        else:
            U.vector()[:] += du.vector()

        U_norm = U.vector().norm("l2")
        print(f'Newton iteration #{newton_it}: |R|/|R0| = {f_norm/f0_norm:.15E}   |du|/|U| = {du_norm/U_norm:.15E}')

        # Check for convergence
        if f_norm < a_tol:
            print(f'[CONVERGED] |R| = {f_norm:.15E}  < a_tol = {a_tol}')
            break
        if f_norm/f0_norm < r_tol:
            print(f'[CONVERGED] |R|/|R0| = {f_norm/f0_norm:.15E}  < r_tol = {r_tol}')
            break
        if du_norm/U_norm < r_tol:
            print(f'[CONVERGED] |du|/|U| = {du_norm/U_norm:.15E}  < r_tol = {r_tol}')
            break

# Save solution
vtkfile = File('rectangular_beam/solution.pvd')
vtkfile << u
