from ufl import (Coefficient, Constant, Identity, Mesh, FunctionSpace,
                 TestFunction, TrialFunction, det, inner,
                 VectorElement, derivative, dx, ds, grad,
                 hexahedron, tr, variable, ln)

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
I_C = tr(C)
II_C = (1.0/2.0)*(tr(C)**2 - tr(C*C))
III_C = det(C)
J = III_C**(1.0/2.0)

# Elasticity parameters
C_01 = Constant(mesh)
C_10 = Constant(mesh)
K = Constant(mesh)

# stored strain energy density (nearly incompressible Mooney-Rivlin model)
psi = C_01*(J**(-2.0/3.0)*I_C - 3) + C_10*(J**(-4.0/3.0)*II_C - 3) + 0.5 * K * (ln(J))**2

# Total potential energy
Pi = psi * dx(degree=2) - inner(B, u)*dx(degree=2) - inner(T, u)*ds(degree=2)

# First variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Export forms
forms = [F, J, Pi]
