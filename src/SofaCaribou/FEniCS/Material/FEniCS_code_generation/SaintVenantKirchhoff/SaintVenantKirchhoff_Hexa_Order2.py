from ufl import (Coefficient, Constant, Identity, Mesh, FunctionSpace,
                 TestFunction, TrialFunction, det, inner,
                 VectorElement, derivative, dx, ds, grad,
                 hexahedron, tr, variable)

# Function spaces
cell = hexahedron
d = cell.geometric_dimension()
coordinate_element = VectorElement("Q", cell, 1)
element = VectorElement("S", cell, 2)

mesh = Mesh(coordinate_element)
V = FunctionSpace(mesh, element)

# Trial and test functions
du = TrialFunction(V)  # Incremental displacement
v = TestFunction(V)  # Test function

# Functions
u = Coefficient(V)  # Displacement from previous iteration
B = Coefficient(V)  # Body forces
T = Coefficient(V)  # Traction forces

# Kinematics
I = Identity(d)  # Identity tensor
F = variable(I + grad(u))  # Deformation gradient
C = variable(F.T * F)  # Right Cauchy-Green tensor
E = variable(0.5 * (C - I))

# Elasticity parameters
young = Constant(mesh)
poisson = Constant(mesh)
mu = young / (2 * (1 + poisson))
lmbda = young * poisson / ((1 + poisson) * (1 - 2 * poisson))

# Stored strain energy density (compressible neo-Hookean model)
psi = (lmbda / 2) * tr(E) ** 2 + mu * tr(E * E)

# Total potential energy
Pi = psi * dx(degree=2) - inner(B, u)*dx(degree=2) - inner(T, u)*ds(degree=2)

# First variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Export forms
forms = [F, J, Pi]