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


class SofaSceneBuilder(object):
    def __init__(self, **kwargs):
        self.simulation = kwargs.get('simulation', None)
        self.node = kwargs.get('node', None)

        assert isinstance(self.simulation, Simulation)
        assert self.node is not None

        # Meshes map set up
        self.meshes_map = {}
        self.meshes_node = self.node.createChild("meshes")

        # Boundary mapping set up
        for b in self.simulation.boundaries:
            if b.mapping is not None:
                self.simulation.add_mappings([
                    b.mapping
                ])

        # Mappings set up
        self.link_map = {}
        for mesh in self.simulation.meshes:
            for part in mesh.parts:
                self.link_map[part] = []

        dependent_parts = []  # List of mesh parts that depend (via a mapping) to another part
        for mapping in self.simulation.mappings:
            self.link_map[mapping.input].append(mapping.output)
            if mapping.output not in dependent_parts:
                dependent_parts.append(mapping.output)

        self.base_part = []  # List of base mesh parts that do not depend (via a mapping) to any other part
        for mesh in self.simulation.meshes:
            for part in mesh.parts:
                if part not in dependent_parts:
                    # only add the part if it is contains a behavior (volume) or boundary (surface)
                    if isinstance(part, VolumePart):
                        if len(self.get_part_behaviors(part)):
                            self.base_part.append(part)
                    else:
                        if len(self.get_part_boundaries(part)):
                            self.base_part.append(part)

        # Set up the base parts
        for part in self.base_part:
            n = self.node.createChild('{}_{}'.format(part.mesh.name, part.name))
            if isinstance(part, VolumePart):
                self.set_volume_part(part, n)
            else:
                self.set_boundary_part(part, n)

    def get_mesh_object(self, mesh):
        try:
            return self.meshes_map[mesh]
        except KeyError:
            if mesh.get_part('volume'):
                self.meshes_map[mesh] = self.meshes_node.createObject(
                    'Mesh',
                    name=mesh.name,
                    position=mesh.vertices.tolist(),
                    edges=mesh.surface.edges.tolist(),
                    triangles=mesh.surface.triangles.tolist(),
                    quads=mesh.surface.quads.tolist(),
                    tetrahedra=mesh.volume.tetrahedrons.tolist(),
                    hexahedra=mesh.volume.hexahedrons.tolist()
                )
            else:
                self.meshes_map[mesh] = self.meshes_node.createObject(
                    'Mesh',
                    name=mesh.name,
                    position=mesh.vertices.tolist(),
                    edges=mesh.surface.edges.tolist(),
                    triangles=mesh.surface.triangles.tolist(),
                    quads=mesh.surface.quads.tolist(),
                )

            return self.meshes_map[mesh]

    def get_part_behaviors(self, part):
        behaviors = []
        for b in self.simulation.behaviors:
            if b.part == part:
                behaviors.append(b)
        return behaviors

    def get_part_materials(self, part):
        materials = []
        for b in self.simulation.materials:
            if b.part == part:
                materials.append(b)
        return materials

    def get_part_boundaries(self, part):
        boundaries = []
        for b in self.simulation.boundaries:
            if b.part == part:
                boundaries.append(b)
        return boundaries

    def get_linked_parts(self, part):
        parts = []
        for m in self.simulation.mappings:
            if m.input == part and m.output not in parts:
                parts.append(m.output)
        return parts

    def get_mappings(self, part):
        mappings = []
        for m in self.simulation.mappings:
            if m.output == part and m not in mappings:
                mappings.append(m)
        return mappings

    def add_topology(self, part, node):

        if part.triangles.size:
            node.createObject(
                'TriangleSetTopologyContainer',
                points=part.points.tolist(),
                edges=part.edges.tolist(),
                triangles=part.triangles.tolist(),
            )
            node.createObject('TriangleSetTopologyModifier')
            node.createObject('TriangleSetTopologyAlgorithms')
            node.createObject('TriangleSetGeometryAlgorithms')

        if part.quads.size:
            node.createObject(
                'QuadSetTopologyContainer',
                points=part.points.tolist(),
                edges=part.edges.tolist(),
                quads=part.quads.tolist(),
            )
            node.createObject('QuadSetTopologyModifier')
            node.createObject('QuadSetTopologyAlgorithms')
            node.createObject('QuadSetGeometryAlgorithms')

        if part.tetrahedrons.size:
            node.createObject(
                'TetrahedronSetTopologyContainer',
                points=part.points.tolist(),
                edges=part.edges.tolist(),
                triangles=part.triangles.tolist(),
                tetrahedra=part.tetrahedrons.tolist(),
            )
            node.createObject('TetrahedronSetTopologyModifier')
            node.createObject('TetrahedronSetTopologyAlgorithms')
            node.createObject('TetrahedronSetGeometryAlgorithms')

        if part.hexahedrons.size:
            node.createObject(
                'HexahedronSetTopologyContainer',
                points=part.points.tolist(),
                edges=part.edges.tolist(),
                quads=part.quads.tolist(),
                hexahedra=part.hexahedrons.tolist(),
            )
            node.createObject('HexahedronSetTopologyModifier')
            node.createObject('HexahedronSetTopologyAlgorithms')
            node.createObject('HexahedronSetGeometryAlgorithms')

    def add_linear_solver(self, node, solver):
        if isinstance(solver, CGLinearSolver):
            node.createObject(
                'CGLinearSolver',
                iterations=solver.iterations,
                tolerance=solver.tolerance,
                threshold=solver.threshold,
                printLog=solver.printLog,
            )
        elif isinstance(solver, PardisoSolver):
            node.createObject(
                'PardisoSolver',
                symmetric=solver.symmetric,
                iterativeSolverNumbering=solver.iterativeSolverNumbering,
                verbose=solver.verbose,
                printLog=solver.printLog,
            )
        else:
            raise NotImplementedError(
                "The linear solver `{}` isn't compatible for the sofa scene builder.".format(
                    solver.__class__.__name__)
            )

    def add_pde_solver(self, node):
        if isinstance(self.simulation.pde_solver, StaticSolver):
            system_solver = self.simulation.pde_solver.solver

            if isinstance(system_solver, NewtonRaphsonSolver):
                node.createObject(
                    'NewtonRaphsonSolver',
                    maxIt=system_solver.maxIt,
                    correctionTolerance=system_solver.correctionTolerance,
                    residualTolerance=system_solver.residualTolerance,
                    convergeOnResidual=system_solver.convergeOnResidual,
                    printLog=system_solver.printLog
                )

                linear_solver = system_solver.linearSolver
                self.add_linear_solver(node, linear_solver)
            elif isinstance(system_solver, LinearSolver):
                node.createObject('StaticSolver')
                self.add_linear_solver(node, system_solver)
            else:
                raise NotImplementedError(
                    "The system solver `{}` isn't compatible for the sofa scene builder.".format(
                        system_solver.__class__.__name__)
                )
        else:
            raise NotImplementedError(
                "The PDE solver `{}` isn't compatible for the sofa scene builder.".format(
                    self.simulation.pde_solver.__class__.__name__)
            )

    def set_volume_part(self, part, node):
        # Check if the part is influenced by at least a behavior in order to add a PDE solver
        dynamic = False
        for b in self.simulation.behaviors:
            if b.part == part:
                dynamic = True
                break

        if dynamic:
            self.add_pde_solver(node)

        node.createObject('MechanicalObject',
                          src=self.get_mesh_object(part.mesh).getLinkPath(),
                          showObject=True,
                          showObjectScale=5,
                          )

        # Topologies creation
        self.add_topology(part, node)

        materials = self.get_part_materials(part)

        # Behaviors creation
        for behavior in self.get_part_behaviors(part):
            if isinstance(behavior, FEMForceField):
                assert part.tetrahedrons.size + part.hexahedrons.size, "FEM Forcefields need either a tetra or hexa topology"
                assert len(materials), "FEM Forcefields need a material"
                assert len(materials) == 1, "FEM Forcefields can only deal with 1 material a the moment"

                material = materials[0]
                if isinstance(material, LinearElastic):
                    method = 'large' if material.corotated else 'small'
                    if part.tetrahedrons.size:
                        node.createObject(
                            'TetrahedronFEMForceField',
                            method=method,
                            youngModulus=material.young_modulus,
                            poissonRatio=material.poisson_ratio
                        )
                    if part.hexahedrons.size:
                        node.createObject(
                            'HexahedronFEMForceField',
                            method=method,
                            youngModulus=material.young_modulus,
                            poissonRatio=material.poisson_ratio
                        )
                elif isinstance(material, StVenantKirchhoff):
                    if part.tetrahedrons.size:
                        mu, l = Utils.lame(young_modulus=material.young_modulus, poisson_ratio=material.poisson_ratio)
                        node.createObject(
                            'TetrahedronHyperelasticityFEMForceField',
                            materialName="StVenantKirchhoff",
                            ParameterSet="{} {}".format(mu, l),
                        )
                    else:
                        raise NotImplementedError("Nonlinear elastic is only supported with tetrahedons.")
                else:
                    raise NotImplementedError(
                        "The material `{}` isn't compatible for the sofa scene builder.".format(
                            material.__class__.__name__)
                    )
            elif isinstance(behavior, MeshlessGalerkin):
                assert len(materials), "Meshless Forcefield need a material"
                assert len(materials) == 1, "Meshless Forcefield can only deal with 1 material a the moment"

                if behavior.surface is not None:
                    assert isinstance(behavior.surface, SurfacePart)
                    node.createObject('DisplacedMeshTopology',
                                      points=part.points.tolist(),
                                      edges=part.edges.tolist(),
                                      triangles=part.triangles.tolist(),
                                      quads=part.quads.tolist(),
                                      )

                material = materials[0]

                (g_min_x, g_min_y, g_min_z), (g_max_x, g_max_y, g_max_z), (nx, ny, nz) = behavior.grid
                grid = node.createObject('MeshlessGridTopology',
                                         n=[nx, ny, nz],
                                         min=[g_min_x, g_min_y, g_min_z],
                                         max=[g_max_x, g_max_y, g_max_z],
                                         grid_type=2,
                                         poissonRatio=material.poisson_ratio,
                                         youngModulus=material.young_modulus,
                                         density=material.density
                                         )

                node.createObject('QuarticSplineKernel', dilatation=behavior.dilatation)

                if isinstance(material, LinearElastic):
                    node.createObject('GalerkinForcefield',
                                      printLog=behavior.printLog,
                                      number_of_neighbors=behavior.number_of_neighbors,
                                      corotated=material.corotated,
                                      gauss_positions=grid.getLinkPath() + '.gauss_positions',
                                      gauss_weights=grid.getLinkPath() + '.gauss_weights',
                                      poissonRatio=material.poisson_ratio,
                                      youngModulus=material.young_modulus,
                                      verbose=behavior.verbose,
                                      )

            elif isinstance(behavior, GravityForceField):
                assert part.tetrahedrons.size + part.hexahedrons.size, "Gravity need either a tetra or hexa topology"
                assert len(materials), "Gravity need a material"
                assert len(materials) == 1, "Gravity can only deal with 1 material a the moment"

                material = materials[0]
                node.createObject(
                    'DiagonalMass',
                    massDensity=material.density,
                )
            else:
                raise NotImplementedError(
                    "The behavior `{}` isn't compatible for the sofa scene builder.".format(
                        behavior.__class__.__name__)
                )

        # Mappings creation
        mappings = self.get_mappings(part)

        # Boundaries creation
        linked_parts = self.get_linked_parts(part)
        for linked_part in linked_parts:
            if isinstance(linked_part, SurfacePart):
                self.set_boundary_part(linked_part, node)
            else:
                self.set_volume_part(linked_part, node)

    def set_boundary_part(self, part, parent_node):
        boundaries = self.get_part_boundaries(part)
        for boundary in boundaries:
            if isinstance(boundary, FixedBoundary):
                # node = parent_node.createChild(
                #     '{}_{}_{}'.format(parent_node.name, part.name, boundary.__class__.__name__))
                # node.createObject('MechanicalObject', src=self.get_mesh_object(part.mesh).getLinkPath())
                parent_node.createObject(
                    'FixedConstraint',
                    indices=part.points.tolist(),
                    activate_projectVelocity=boundary.velocities_fixed
                )
            elif isinstance(boundary, PressureBoundary):
                node = parent_node.createChild(
                    '{}_{}_{}'.format(parent_node.name, part.name, boundary.__class__.__name__))
                node.createObject('MechanicalObject', src=self.get_mesh_object(part.mesh).getLinkPath())

                # Topologies creation
                self.add_topology(part, node)

                node.createObject(
                    'PressureForcefield',
                    pressure=boundary.pressure,
                    slope=boundary.slope,
                    printLog=boundary.printLog,
                    triangles=part.triangles.tolist(),
                )

                if boundary.linked_to is not None:
                    self.set_mapping(parent_node, node, boundary.mapping)

            elif isinstance(boundary, ConstantForceBoundary):
                node = parent_node.createChild(
                    '{}_{}_{}'.format(parent_node.name, part.name, boundary.__class__.__name__))
                node.createObject('MechanicalObject', src=self.get_mesh_object(part.mesh).getLinkPath())

                node.createObject(
                    'ConstantForceField',
                    totalForce=boundary.force,
                    indices=boundary.part.points.tolist(),
                )

                if boundary.linked_to is not None:
                    self.set_mapping(parent_node, node, boundary.mapping)
            elif isinstance(boundary, VisualBoundary):
                node = parent_node.createChild(
                    '{}_{}_{}'.format(parent_node.name, part.name, boundary.__class__.__name__))

                node.createObject('OglModel', src=self.get_mesh_object(part.mesh).getLinkPath())

                if boundary.linked_to is not None:
                    self.set_mapping(parent_node, node, boundary.mapping)
            elif isinstance(boundary, WatcherBoundary):
                if boundary.linked_to is None:
                    raise AttributeError("The Watcher boundary needs a linked part.")

                node = parent_node.createChild(
                    '{}_{}_{}'.format(parent_node.name, part.name, boundary.__class__.__name__))

                mo = node.createObject('MechanicalObject', src=self.get_mesh_object(part.mesh).getLinkPath())
                boundary.setState(mo)

                self.set_mapping(parent_node, node, boundary.mapping)
            else:
                raise NotImplementedError(
                    "The boundary `{}` isn't compatible for the sofa scene builder.".format(
                        boundary.__class__.__name__)
                )

    def set_mapping(self, input_node, output_node, mapping):
        if isinstance(mapping, IdentityMapping):
            output_node.createObject(
                'IdentityMapping',
                input=input_node.getLinkPath(),
                output=output_node.getLinkPath()
            )
        elif isinstance(mapping, LinearMapping):
            output_node.createObject(
                'LinearMapping',
                input=input_node.getLinkPath(),
                output=output_node.getLinkPath()
            )
        elif isinstance(mapping, BarycentricMapping):
            output_node.createObject(
                'BarycentricMapping',
                input=input_node.getLinkPath(),
                output=output_node.getLinkPath()
            )
        elif isinstance(mapping, ParticleMapping):
            output_node.createObject(
                'ParticleMapping',
                input=input_node.getLinkPath(),
                output=output_node.getLinkPath(),
                number_of_neighbors=mapping.numberOfNeighbors,
            )
        else:
            raise NotImplementedError(
                "The mapping `{}` isn't compatible for the sofa scene builder.".format(
                    mapping.__class__.__name__)
            )

