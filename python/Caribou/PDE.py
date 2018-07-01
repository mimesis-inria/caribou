from .Optimization import SystemSolver
from .Base import BaseObject


class PDESolver(BaseObject):
    def __init__(self, **kwargs):
        BaseObject.__init__(self, **kwargs)
        # Parameters
        None


class StaticSolver(PDESolver):
    def __init__(self, **kwargs):
        PDESolver.__init__(self, **kwargs)

        # Parameters
        self.solver = kwargs.get('solver', None)

        assert isinstance(self.solver, SystemSolver)


class DynamicSolver(PDESolver):
    None