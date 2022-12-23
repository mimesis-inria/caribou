from ufl import (Coefficient, Constant, Identity,
                 TestFunction, TrialFunction, det, inner, 
                 VectorElement, derivative, dx, ds, grad,
                 hexahedron, tr, variable, ln)

# Function spaces
cell = hexahedron
d = cell.geometric_dimension()
element = VectorElement("S", cell, 1)

# Trial and test functions
du = TrialFunction(element)  # Incremental displacement
v = TestFunction(element)  # Test function

# Functions
u = Coefficient(element)  # Displacement from previous iteration
B = Coefficient(element)  # Body forces
T = Coefficient(element)  # Traction forces

# Kinematics
I = Identity(d)  # Identity tensor
F = variable(I + grad(u))  # Deformation gradient
C = variable(F.T * F)  # Right Cauchy-Green tensor
E = variable(0.5 * (C - I))
Ic = tr(C)
J = det(F)

# Elasticity parameters
young = Constant(cell)
poisson = Constant(cell)
mu = young / (2 * (1 + poisson))
lmbda = young * poisson / ((1 + poisson) * (1 - 2 * poisson))

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2

# Total potential energy
Pi = psi * dx(degree=2) - inner(B, u)*dx(degree=2) - inner(T, u)*ds(degree=2)

# First variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Export forms
forms = [F, J, Pi]