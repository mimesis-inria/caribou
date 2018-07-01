
class BaseObject(object):
    def __init__(self, **kwargs):
        self._params = kwargs.copy()

    @property
    def classname(self):
        return type(self).__name__