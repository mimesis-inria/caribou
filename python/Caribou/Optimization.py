from .Base import BaseObject


class SystemSolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)


class LinearSolver(SystemSolver):
    def __init__(self, **kwargs):
        SystemSolver.__init__(self, **kwargs)

        # Parameters
        self.printLog = kwargs.get('print_log', False)


class CGLinearSolver(LinearSolver):
    def __init__(self, **kwargs) :
        LinearSolver.__init__(self, **kwargs)

        # Parameters
        self.iterations = kwargs.get('maximum_iterations', 2500)
        self.tolerance  = kwargs.get('tolerance', 1e-8)
        self.threshold  = kwargs.get('threshold', 1e-8)


class PardisoSolver(LinearSolver):
    NONSYMETRIC = 0
    SYMETRIC = 1
    SYMETRIC_STRUCTURALLY = -1
    SYMETRIC_POSITIVE_DEFINITE = 2

    def __init__(self, **kwargs) :
        LinearSolver.__init__(self, **kwargs)

        # Parameters
        self.symmetric                = kwargs.get('symmetric', PardisoSolver.SYMETRIC)
        self.iterativeSolverNumbering = kwargs.get('iterativeSolverNumbering', True)
        self.verbose = kwargs.get('verbose', False)


class NonLinearSolver(SystemSolver):
    def __init__(self, **kwargs):
        SystemSolver.__init__(self, **kwargs)

        # Parameters
        self.linearSolver = kwargs.get('linear_solver', None)

        assert isinstance(self.linearSolver, LinearSolver)


class NewtonRaphsonSolver(NonLinearSolver):
    def __init__(self, **kwargs):
        NonLinearSolver.__init__(self, **kwargs)

        # Parameters
        self.maxIt = kwargs.get('maximum_iterations', 20)
        self.correctionTolerance = kwargs.get('correction_tolerance', 1e-6)
        self.convergeOnResidual = kwargs.get('converge_on_residual', False)
        self.residualTolerance = kwargs.get('residual_tolerance', 1e-6)
        self.printLog = kwargs.get('print_log', False)
