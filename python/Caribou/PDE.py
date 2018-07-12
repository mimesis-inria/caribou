from .Optimization import SystemSolver
from .Base import BaseObject


class PDESolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)
        # Parameters
        self.solver = kwargs.get('solver', None)
        self.printLog = kwargs.get('print_log', False)

        assert isinstance(self.solver, SystemSolver)

    def __eq__(self, other):
        """
        Compare the PDESolver with another one
        :param other: The other PDESolver
        :type other: PDESolver
        :return: True if their parameters are equal, false otherwise
        ":rtype: bool
        """
        res = BaseObject.__eq__(self, other)
        res = res and self.solver == other.solver
        return res

    def printable_attributes(self):
        return [
            ('System solver', self.solver.fullname())
        ] + BaseObject.printable_attributes(self)


class StaticSolver(PDESolver):
    def __init__(self, **kwargs):
        PDESolver.__init__(self, **kwargs)


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
