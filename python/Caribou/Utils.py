import importlib
import math
import sys
import errno
import os
import seaborn as sns
import numpy as np
from numpy import linalg as LA
from collections import namedtuple


class Struct(object):
    def __new__(cls, data):
        if isinstance(data, dict):
            return namedtuple(
                'Struct', data.iterkeys()
            )(
                *(Struct(val) for val in data.values())
            )
        elif isinstance(data, (tuple, list, set, frozenset)):
            return type(data)(Struct(_) for _ in data)
        else:
            return data


def fullname(o):
    module = o.__class__.__module__
    if module is None or module == str.__class__.__module__:
        return o.__class__.__name__
    return module + '.' + o.__class__.__name__


def string_to_object(dotted_path):
    """
    Import a dotted module path and return the attribute/class designated by the
    last name in the path. Raise ImportError if the import failed.
    """
    try:
        module_path, class_name = dotted_path.rsplit('.', 1)
    except ValueError:
        msg = "%s doesn't look like a module path" % dotted_path
        raise ValueError(msg)

    try:
        module = importlib.import_module(module_path)
    except Exception:
        try:
            module = importlib.import_module("Benchmark." + module_path)
        except Exception:
            msg = 'Module "%s" not found' % (
                module_path)
            raise AttributeError(msg)


    try:
        return getattr(module, class_name)
    except AttributeError:
        msg = 'Module "%s" does not define a "%s" attribute/class' % (
            module_path, class_name)
        raise AttributeError(msg)


def JsonExport(o):
    try:
        return {
            'object' : fullname(o),
            'params' : o.params
        }
    except AttributeError:
        return {
            'object' : fullname(o)
        }


def JsonImport(dict):
    if dict.has_key('object'):
        class_ = string_to_object(dict.get('object'))
        params = dict.get('params', {})

        object = class_(**params)

        return object
    return dict


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def time_string(seconds):
    ms, s = math.modf(seconds)
    m, s = divmod(s, 60)
    h, m = divmod(m, 60)
    ms = ms * 1000
    if h > 0:
        return "%dh %02dm %02ds %03dms" % (h, m, s, ms)
    elif m > 0:
        return "%02dm %02ds %03dms" % (m, s, ms)
    elif s > 0:
        return "%02ds %03dms" % (s, ms)
    else:
        return "%03dms" % ms


def dump(obj):
    for attr in dir(obj):
        print("obj.%s = %s" % (attr, getattr(obj, attr)))


def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def memory_usage():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem


def memory_available():
    # return the memory usage in MB
    import psutil
    mem = psutil.virtual_memory().available / float(2 ** 20)
    return mem


def escape(string):
    if isinstance(string, unicode):
        o = unicode
    else:
        o = str

    s = ""
    for c in string.lower():
        if o.isalpha(c) or o.isdigit(c):
            s = s + c
        elif o.isspace(c):
            s  = s + '_'
    return s


def generate_n_colors(n):
    """
    Generate N different color sets
    (from https://stackoverflow.com/a/17684501)
    """
    if n < 7:
        return sns.color_palette()
    else:
        sns.color_palette("hls", n)


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def lame(young_modulus, poisson_ratio):
    mu = young_modulus / (2. * (1. + poisson_ratio))
    l = young_modulus * poisson_ratio / ((1. + poisson_ratio) * (1. - 2.*poisson_ratio))
    return mu, l


def from_lame(mu, l):
    poisson_ratio = 1. / (2 * ((mu / l) + 1))
    young_modulus = 2*mu*(1 + poisson_ratio)
    return poisson_ratio, young_modulus


def rotate(points, source_direction, target_direction):
    s = np.asarray(source_direction)
    t = np.asarray(target_direction)

    if LA.norm(s) == 0 or LA.norm(t) == 0:
        return points

    s = s / LA.norm(s) # source normalized
    t = t / LA.norm(t) # target normalized

    vs = np.cross(s, t) # axis of rotation * sin
    if LA.norm(vs) == 0:
        return points

    v  = vs / LA.norm(vs) # axis of rotation
    a  =  np.dot(s, t)    # cos angle
    vt = v * (1 - a)

    # rotation matrix
    R = np.asarray([
        [vt[0]*v[0] + a,      vt[0]*v[1] - vs[2],      vt[2]*v[0] + vs[1]],
        [vt[0]*v[1] + vs[2],  vt[1]*v[1] + a,          vt[1]*v[2] - vs[0]],
        [vt[2]*v[0] - vs[1],  vt[1]*v[2] + vs[0],      vt[2]*v[2] + a]
    ])

    result = [R.dot(p) for p in points]
    return np.asarray(result)


def translate(points, translation):
    points_array = np.asarray(points)
    translation_vector = np.asarray(translation)

    result = points_array + translation_vector

    if isinstance(points, list):
        return result.tolist()
    else:
        return result


def bbox(vertices):
    points_array = np.asarray(vertices)
    m = np.min(points_array, axis=0)
    xmin, ymin, zmin = m[0], m[1], m[2]

    m = np.max(points_array, axis=0)
    xmax, ymax, zmax = m[0], m[1], m[2]

    return xmin, xmax, ymin, ymax, zmin, zmax


class Obj(object):
    def __init__(self, filename=None):
        self.filename = filename
        self.vertices = []
        self.faces = []
        self.bbox = Struct({
            'min': [0, 0, 0],
            'max': [0, 0, 0]
        })

        if filename:
            self.load(filename)

    def load(self, filename):
        # https://inareous.github.io/posts/opening-obj-using-py
        try:
            f = open(filename)
            for line in f:
                if line[:2] == "v ":
                    index1 = line.find(" ") + 1
                    index2 = line.find(" ", index1 + 1)
                    index3 = line.find(" ", index2 + 1)

                    vertex = [float(line[index1:index2]), float(line[index2:index3]), float(line[index3:-1])]
                    vertex = [round(vertex[0], 2), round(vertex[1], 2), round(vertex[2], 2)]
                    self.vertices.append(vertex)

                elif line[0] == "f":
                    string = line.replace("//", "/")
                    ##
                    i = string.find(" ") + 1
                    face = []
                    for item in range(string.count(" ")):
                        if string.find(" ", i) == -1:
                            face.append(string[i:-1])
                            break
                        face.append(string[i:string.find(" ", i)])
                        i = string.find(" ", i) + 1
                    ##
                    self.faces.append(face)

            f.close()

            xmin, xmax, ymin, ymax, zmin, zmax = bbox(self.vertices)
            self.bbox = Struct({
                'min': [xmin, ymin, zmin],
                'max': [xmax, ymax, zmax]
            })

        except IOError:
            print("{} file not found.".format(filename))