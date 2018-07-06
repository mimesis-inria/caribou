from .Base import BaseObject


class SystemSolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.printLog = kwargs.get('print_log', False)


class LinearSolver(SystemSolver):
    def __init__(self, **kwargs):
        SystemSolver.__init__(self, **kwargs)


class CGLinearSolver(LinearSolver):
    def __init__(self, **kwargs) :
        LinearSolver.__init__(self, **kwargs)

        # Parameters
        self.iterations = kwargs.get('maximum_iterations', 2500)
        self.tolerance  = kwargs.get('tolerance', 1e-8)
        self.threshold  = kwargs.get('threshold', 1e-8)

    def __eq__(self, other):
        """
        Compare the CGLinearSolver with another one
        :param other: The other CGLinearSolver
        :type other: CGLinearSolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = LinearSolver.__eq__(self, other)
        res = res and self.iterations == other.iterations
        res = res and self.tolerance == other.tolerance
        res = res and self.threshold == other.threshold
        return res

    def printable_attributes(self):
        return [
                   ('Maximum number of iterations', self.iterations),
                   ('Error tolerance', self.tolerance),
                   ('Denominator threshold', self.threshold)
        ] + LinearSolver.printable_attributes(self)


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

    def __eq__(self, other):
        """
        Compare the PardisoSolver with another one
        :param other: The other PardisoSolver
        :type other: PardisoSolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = LinearSolver.__eq__(self, other)
        res = res and self.symmetric == other.symmetric
        res = res and self.iterativeSolverNumbering == other.iterativeSolverNumbering
        return res

    def printable_attributes(self):
        m_type = "Unknown"
        if self.symmetric == PardisoSolver.NONSYMETRIC:
            m_type = "Non symmetric"
        elif self.symmetric == PardisoSolver.SYMETRIC:
            m_type = 'Symmetric'
        elif self.symmetric == PardisoSolver.SYMETRIC_STRUCTURALLY:
            m_type = 'Structurally symmetric'
        elif self.symmetric == PardisoSolver.SYMETRIC_POSITIVE_DEFINITE:
            m_type = 'Symmetric positive definite'
        return [
                   ('Matrix type', m_type),
                   ('Iterative solver numbering', str(self.iterativeSolverNumbering)),
        ] + LinearSolver.printable_attributes(self)


class NonLinearSolver(SystemSolver):
    def __init__(self, **kwargs):
        SystemSolver.__init__(self, **kwargs)

        # Parameters
        self.linearSolver = kwargs.get('linear_solver', None)

        assert isinstance(self.linearSolver, LinearSolver)

    def __eq__(self, other):
        """
        Compare the NonLinearSolver with another one
        :param other: The other NonLinearSolver
        :type other: NonLinearSolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = SystemSolver.__eq__(self, other)
        res = res and self.linearSolver == other.linearSolver
        return res

    def printable_attributes(self):
        return [
                   ('Linear solver', self.linearSolver.fullname()),
        ] + SystemSolver.printable_attributes(self)


class NewtonRaphsonSolver(NonLinearSolver):
    def __init__(self, **kwargs):
        NonLinearSolver.__init__(self, **kwargs)

        # Parameters
        self.maxIt = kwargs.get('maximum_iterations', 20)
        self.correctionTolerance = kwargs.get('correction_tolerance', 1e-6)
        self.convergeOnResidual = kwargs.get('converge_on_residual', False)
        self.residualTolerance = kwargs.get('residual_tolerance', 1e-6)
        self.cutoff = kwargs.get('cutoff', False)

    def __eq__(self, other):
        """
        Compare the NewtonRaphsonSolver with another one
        :param other: The other NewtonRaphsonSolver
        :type other: NewtonRaphsonSolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = NonLinearSolver.__eq__(self, other)
        res = res and self.maxIt == other.maxIt
        res = res and self.correctionTolerance == other.correctionTolerance
        res = res and self.convergeOnResidual == other.convergeOnResidual
        res = res and self.residualTolerance == other.residualTolerance
        res = res and self.cutoff == other.cutoff
        return res

    def printable_attributes(self):
        return [
                   ('Maximum number of iterations', self.maxIt),
                   ('Correction tolerance (|dx|)', self.correctionTolerance),
                   ('Residual tolerance (|f - K(x0 + dx)|)', self.residualTolerance),
                   ('Converge on residual tolerance', self.convergeOnResidual),
                   ('Stop the iterations when the residual goes up', self.cutoff),
        ] + NonLinearSolver.printable_attributes(self)
