def compute(mu, lmbda, rad, length):
    import sympy as smp
    from sympy import Matrix
    import numpy as np

    Tr  = lambda A: A[0, 0] + A[1, 1] + A[2, 2]
    Div = lambda A, b: Matrix(3, 3, lambda i, j: smp.diff(A[i,j], b[j], 1)) * smp.ones(A.shape[1], 1)

    l, m    = 1.25, 1.
    Id      = smp.eye(3)
    x, y, z = smp.symbols('x y z')
    u, v, w = [1e-3*z*smp.exp(v) for v in [x, y, z]]
    grad_u  = Matrix([u, v, w]).jacobian(Matrix([x, y, z]))

    F       = smp.MatrixSymbol('F', 3, 3)
    C       = F.T*F
    E       = (C-Id)/2
    W       = 0.5*l*Tr(E)**2 + m*Tr(E*E)
    P       = Matrix(3, 3, lambda i, j: smp.diff(W, F[i, j], 1)).subs(F, grad_u + Id)
    f       = - Div(P, [x, y, z])

    u, v, w = smp.lambdify([x, y, z], u), smp.lambdify([x, y, z], v), smp.lambdify([x, y, z], w)

    return smp.lambdify([x, y, z], P), \
           smp.lambdify([x, y, z], f), \
           (lambda x_, y_, z_: np.array([u(x_, y_, z_), v(x_, y_, z_), w(x_, y_, z_)]).reshape(3,))
