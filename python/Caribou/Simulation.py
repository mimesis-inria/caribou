from .PDE import *
from .Boundary import *
from .Material import *
from .Behavior import *
from .Optimization import *
from .Mapping import *
from .Mesh import *

import Utils

import numpy as np
from numpy import linalg as LA


class Simulation(object):
    def __init__(self, **kwargs):
        # Members
        self.pde_solver = None
        self.behaviors = []
        self.boundaries = []
        self.materials = []
        self.meshes = []
        self.mappings = []

    def set_PDE_solver(self, solver=None):
        assert isinstance(solver, PDESolver)
        self.pde_solver = solver

    def add_meshes(self, meshes):
        for m in meshes:
            assert isinstance(m, Mesh)
            self.meshes.append(m)

    def add_boundaries(self, boundaries):
        for b in boundaries:
            assert isinstance(b, Boundary)
            self.boundaries.append(b)

    def add_materials(self, materials):
        for m in materials:
            assert isinstance(m, Material)
            self.materials.append(m)

    def add_behaviors(self, behaviors):
        for b in behaviors:
            assert isinstance(b, Behavior)
            self.behaviors.append(b)

    def add_mappings(self, mappings):
        for m in mappings:
            assert isinstance(m, Mapping)
            self.mappings.append(m)
