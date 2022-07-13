from ufl import (Coefficient, Constant, Identity, Mesh, FunctionSpace,
                 TestFunction, TrialFunction, inner, ds,
                 VectorElement, derivative, dx, grad,
                 hexahedron, tr, variable, det, as_vector,
                 sqrt, dot, ln, exp)

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
J  = det(F)
I1 = tr(C)

# Elasticity parameters
bulk_modulus = Constant(mesh) 
a =    Constant(mesh)   
b =    Constant(mesh)  
a_f =  Constant(mesh)
b_f =  Constant(mesh)
a_s =  Constant(mesh)
b_s =  Constant(mesh)
a_fs = Constant(mesh)
b_fs = Constant(mesh)

f_0 = as_vector([0.0, 1.0/sqrt(2), 1.0/sqrt(2)])
s_0 = as_vector([0.0, 1.0/sqrt(2), -1.0/sqrt(2)])

I_4_f_0 = dot(C*f_0, f_0)
I_4_s_0 = dot(C*s_0, s_0)
I_8_f_0_s_0 = dot(C*s_0, f_0)

W_vol = bulk_modulus*(J**2-1.0-2.0*ln(J))/4.0
W_1 = a*exp((I1-3.0)*b)/(2.0*b)

I_4_f_0 = I_4_f_0-1.0
I_4_s_0 = I_4_s_0-1.0

W_4f = a_f*(exp(((I_4_f_0)**2)*b_f)-1.0)/(2.0*b_f)
W_4s = a_s*(exp(((I_4_s_0)**2)*b_s)-1.0)/(2.0*b_s)

W_8_fs = a_fs*(exp((I_8_f_0_s_0**2)*b_fs)-1.0)/(2.0*b_fs)

# Stored strain energy density (Ogden)
psi = W_vol + W_1 + W_4f + W_4s + W_8_fs

# Total potential energy
Pi = psi * dx(degree=2) - inner(B, u)*dx(degree=2) - inner(T, u)*ds(degree=2)

# First variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Export forms
forms = [F, J, Pi]