
class BaseObject(object):
    def __init__(self, **kwargs):
        self._params = kwargs.copy()

    @classmethod
    def get_class_name(cls):
        return cls.__name__

    @property
    def classname(self):
        return type(self).__name__

    def serialize(self, keys=[], recursive=True):
        if not keys or not isinstance(keys, list) or not len(keys):
            keys = self._params.keys()

        s = {
            'Class': self.fullname(),
            'attributes': {}
        }

        for key in keys:
            s['attributes'][key] = _serialize(o=getattr(self, key), recursive=recursive)

        return s

    def fullname(self):
        module = self.__class__.__module__
        if module is None or module == str.__class__.__module__:
            return self.__class__.__name__
        return module + '.' + self.__class__.__name__

    def printable_attributes(self):
        return []


def _serialize(o, recursive=True):
    if isinstance(o, BaseObject):
        if recursive:
            return o.serialize()
        else:
            return o.fullname()
    elif type(o) == type:
        if issubclass(o, BaseObject):
            return o.get_class_name()
        else:
            return "Unknown class"
    elif isinstance(o, tuple):
        if recursive:
            return {
                'Class': 'tuple',
                'attributes': [_serialize(att) for att in o]
            }
        else:
            return 'tuple'
    elif isinstance(o, list):
        if recursive:
            return [_serialize(att) for att in o]
        else:
            return 'list'
    elif isinstance(o, dict):
        if recursive:
            s = {}
            for k in o:
                s[k] = _serialize(o[k])
            return s
        else:
            return 'dict'
    elif type(o) in [str, int, float, bool]:
        return o
    else:
        return "Unknown type to serialize"




