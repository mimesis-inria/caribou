from .Base import BaseObject


class LinearSolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)

        # Parameters
        self.printLog = kwargs.get('print_log', False)


class LUSolver(LinearSolver):
    def __init__(self, **kwargs):
        LinearSolver.__init__(self, **kwargs)

        # Parameters
        self.tolerance = kwargs.get('tolerance', 0.001)
        self.verbose = kwargs.get('verbose', False)

    def __eq__(self, other):
        """
        Compare the CGLinearSolver with another one
        :param other: The other CGLinearSolver
        :type other: CGLinearSolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = LinearSolver.__eq__(self, other)
        res = res and (self.tolerance == other.tolerance)
        return res

    def printable_attributes(self):
        return [
                   ('Tolerance', self.tolerance),
        ] + LinearSolver.printable_attributes(self)


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
