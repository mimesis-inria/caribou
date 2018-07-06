from .Optimization import SystemSolver
from .Base import BaseObject


class PDESolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)
        # Parameters
        self.solver = kwargs.get('solver', None)

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
    None