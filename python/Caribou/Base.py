import importlib
import numpy as np

from .Utils import string_to_object


class BaseObject(object):
    def __init__(self, **kwargs):
        self._params = kwargs.copy()

    @classmethod
    def get_class_name(cls):
        return cls.__name__

    @classmethod
    def get_class_fullname(cls):
        module = cls.__module__
        if module is None or module == str.__class__.__module__:
            return cls.__name__
        return module + '.' + cls.__name__

    @classmethod
    def deserialize(cls, **kwargs):
        return cls(**kwargs)

    @property
    def classname(self):
        return type(self).__name__

    def serialize(self):
        return self._params

    def fullname(self):
        module = self.__class__.__module__
        if module is None or module == str.__class__.__module__:
            return self.__class__.__name__
        return module + '.' + self.__class__.__name__

    def printable_attributes(self):
        return []


def deserialize(attributes):
    if isinstance(attributes, dict):
        if 'Class' in attributes:
            classname = attributes['Class']
            if classname == 'type':
                o = string_to_object(attributes['name'])
                return o
            elif classname == 'tuple':
                return tuple([deserialize(o) for o in attributes['attributes']])
            elif classname == 'ndarray':
                return np.array(attributes['value'])
            else:
                o = string_to_object(classname)
                assert issubclass(o, BaseObject)
                if 'attributes' in attributes:
                    att = deserialize(attributes['attributes'])
                    return o.deserialize(**att)
                else:
                    return o.deserialize()
        else:
            ret = {}
            for k in attributes:
                ret[k] = deserialize(attributes[k])
            return ret
    elif isinstance(attributes, list):
        return [deserialize(o) for o in attributes]
    elif isinstance(attributes, str) and attributes == 'None':
        return None
    elif type(attributes) in [str, int, float, bool, unicode]:
        return attributes
    else:
        print 'Unkown attribute type "{}" of {}'.format(type(attributes), attributes)
        return None


def serialize(o):
    if isinstance(o, BaseObject):
        att = o.serialize()
        for k in att:
            try:
                att[k] = serialize(att[k])
            except AttributeError as e:
                raise AttributeError("Unable to serialize the parameter '{}' (type '{}') of the object '{}' :\n\t {}".format(k, type(att[k]), o.fullname(), e))
        s = {
            'Class': o.fullname(),
            'attributes': att
        }
        return s
    elif type(o) == type:
        if issubclass(o, BaseObject):
            return {
                'Class': 'type',
                'name': o.get_class_fullname()
            }
        else:
            raise AttributeError("Unknown class to serialize '{}'".format(type(o)))
    elif isinstance(o, tuple):
        return {
            'Class': 'tuple',
            'attributes': [serialize(att) for att in o]
        }
    elif isinstance(o, list):
        return [serialize(att) for att in o]
    elif isinstance(o, dict):
        s = {}
        for k in o:
            s[k] = serialize(o[k])
        return s
    elif isinstance(o, np.ndarray):
        return {
            'Class': 'ndarray',
            'value': o.tolist()
        }
    elif o is None:
        return 'None'
    elif type(o) in [str, int, float, bool, unicode]:
        return o
    else:
        raise AttributeError("Unknown type to serialize '{}'".format(type(o)))
