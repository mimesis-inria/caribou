from .Optimization import LinearSolver
from .Base import BaseObject


class PDESolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)
        # Parameters
        self.printLog = kwargs.get('print_log', False)


class StaticSolver(PDESolver):
    def __init__(self, **kwargs):
        PDESolver.__init__(self, **kwargs)

        # Parameters
        self.solver = kwargs.get('solver', None)
        self.newton_iterations = kwargs.get('newton_iterations', 1)
        self.correction_tolerance = kwargs.get('correction_tolerance', 1e-6)
        self.residual_tolerance = kwargs.get('residual_tolerance', 1e-6)

        assert isinstance(self.solver, LinearSolver)

    def __eq__(self, other):
        """
        Compare the PDESolver with another one
        :param other: The other PDESolver
        :type other: StaticSolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = PDESolver.__eq__(self, other)
        res = res and (self.solver == other.solver)
        res = res and (self.newton_iterations == other.newton_iterations)
        res = res and (self.correction_tolerance == other.correction_tolerance)
        res = res and (self.residual_tolerance == other.residual_tolerance)
        return res

    def printable_attributes(self):
        return [
            ('Linear solver', self.solver.fullname()),
            ('Maximum number of iterations', self.newton_iterations),
            ('Correction tolerance (|dx|)', self.correction_tolerance),
            ('Residual tolerance (|f - K(x0 + dx)|)', self.residual_tolerance),
        ] + PDESolver.printable_attributes(self)


class NewtonRaphsonSolver(StaticSolver):
    def __init__(self, **kwargs):
        StaticSolver.__init__(self, **kwargs)

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
        res = StaticSolver.__eq__(self, other)
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
        ] + StaticSolver.printable_attributes(self)


class DynamicSolver(PDESolver):
    def __init__(self, **kwargs):
        PDESolver.__init__(self, **kwargs)
        # Parameters
        self.solver = kwargs.get('solver', None)
        self.mass = kwargs.get('mass', None)


class ImplicitEuler(DynamicSolver):
    def __init__(self, **kwargs):
        DynamicSolver.__init__(self, **kwargs)

        # Parameters
        self.rayleigh_stiffness = kwargs.get('rayleigh_stiffness', 0)
        self.rayleigh_mass = kwargs.get('rayleigh_mass', 0)
        self.velocity_damping = kwargs.get('velocity_damping', 0)
