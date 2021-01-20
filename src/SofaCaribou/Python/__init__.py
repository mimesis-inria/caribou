# We first import Sofa.Core since without this, bindings will be unable
# to have sofa::core::objectmodel::BaseObject as a base class
from Sofa import Core

# Then, we import the actual bindings
from .${MODULENAME} import *