from ufl import (Coefficient, Constant, Identity,
                 TestFunction, TrialFunction, det, inner,
                 VectorElement, derivative, dx, ds, grad,
                 tetrahedron, tr, variable, ln)

# Function spaces
cell = tetrahedron
d = cell.geometric_dimension()
element = VectorElement("Lagrange", cell, 1)

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
I_C = tr(C)
II_C = (1.0/2.0)*(tr(C)**2 - tr(C*C))
III_C = det(C)
J = III_C**(1.0/2.0)

# Elasticity parameters
C_01 = Constant(cell)
C_10 = Constant(cell)
K = Constant(cell)

# stored strain energy density (nearly incompressible Mooney-Rivlin model)
psi = C_01*(J**(-2.0/3.0)*I_C - 3) + C_10*(J**(-4.0/3.0)*II_C - 3) + 0.5 * K * (ln(J))**2

# Total potential energy
Pi = psi * dx(degree=1) - inner(B, u)*dx(degree=1) - inner(T, u)*ds(degree=1)

# First variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Export forms
forms = [F, J, Pi]
