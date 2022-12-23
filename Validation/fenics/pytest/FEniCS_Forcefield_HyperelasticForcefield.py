import sys
import unittest
from parameterized import parameterized
from pathlib import Path

import Sofa
import SofaCaribou
import numpy as np
from scipy.sparse import csr_matrix
import meshio
import docker

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')
YOUNG_MODULUS = 3000
POISSON_RATIO = 0.3

TEST_CASES = [
    ("SaintVenantKirchhofftetra1", "SaintVenantKirchhoff", "Tetrahedron"),
    ("SaintVenantKirchhofftetra2", "SaintVenantKirchhoff", "Tetrahedron10"),
    ("SaintVenantKirchhoffhexa1", "SaintVenantKirchhoff", "Hexahedron"),
    ("SaintVenantKirchhoffhexa2", "SaintVenantKirchhoff", "Hexahedron20"),
    ("NeoHooketetra1", "NeoHookean", "Tetrahedron"),
    ("NeoHooketetra2", "NeoHookean", "Tetrahedron10"),
    ("NeoHookehexa1", "NeoHookean", "Hexahedron"),
    ("NeoHookehexa2", "NeoHookean", "Hexahedron20"),
]


def generate_geometry(element):
    element_sofa, element_sonics, mesh, indices_sofa, indices_sonics, fixed_indices = None, None, None, None, None, None
    if element == "Tetrahedron":
        element_sofa = "Tetrahedron"
        element_sonics = element_sofa
        mesh = meshio.read(current_dir / 'meshes' / 'cube_tetra_1.msh')
        indices_sofa = mesh.cells_dict['tetra']
        indices_sonics = indices_sofa
        fixed_indices = get_fixed_indices(mesh, "triangle")

    elif element == "Tetrahedron10":
        element_sofa = "Tetrahedron10"
        element_sonics = element_sofa
        mesh = meshio.read(current_dir / 'meshes' / 'cube_tetra_2.msh')
        indices_sofa = mesh.cells_dict['tetra10']
        indices_sonics = indices_sofa[:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
        fixed_indices = get_fixed_indices(mesh, "triangle6")

    elif element == "Hexahedron":
        element_sofa = "Hexahedron"
        element_sonics = element_sofa + "_FEniCS"
        mesh = meshio.read(current_dir / 'meshes' / 'cube_hexa_1.msh')
        indices_sofa = mesh.cells_dict['hexahedron']
        indices_sonics = indices_sofa[:, [0, 1, 3, 2, 4, 5, 7, 6]]
        fixed_indices = get_fixed_indices(mesh, "quad")

    elif element == "Hexahedron20":
        element_sofa = "Hexahedron20"
        element_sonics = "Hexahedron_FEniCS20"
        mesh = meshio.read(current_dir / 'meshes' / 'cube_hexa_2.msh')
        indices_sofa = mesh.cells_dict['hexahedron20']
        indices_sonics = indices_sofa[:, [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]
        fixed_indices = get_fixed_indices(mesh, "quad8")
    return element_sofa, element_sonics, mesh.points, indices_sofa, indices_sonics, fixed_indices


def get_fixed_indices(mesh, mesh_type):
    a = [i.size for i in mesh.cell_data["gmsh:physical"]]
    indices = np.cumsum(a)
    return mesh.cells_dict[mesh_type][0:indices[0]]


def createScene(node, element_type, material_model):
    node.addObject('DefaultVisualManagerLoop')
    node.addObject('DefaultAnimationLoop')
    node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    node.addObject('RequiredPlugin',
                   pluginName="SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaBoundaryCondition")
    element_sofa, element_sonics, position, indices_sofa, indices_sonics, fixed_indices = generate_geometry(
        element_type)
    node.gravity = [0, -0.1, 0]
    sofa_node = node.addChild("sofa_node")
    sofa_node.addObject('StaticSolver', newton_iterations="15", relative_correction_tolerance_threshold="1e-15",
                        relative_residual_tolerance_threshold="1e-10", printLog="0")
    sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sofa_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    sofa_node.addObject('CaribouTopology', name='topology', template=element_sofa, indices=indices_sofa.tolist())
    sofa_node.addObject('FixedConstraint', indices=fixed_indices.tolist())
    sofa_node.addObject("CaribouMass", density=3)
    sofa_node.addObject(material_model + "Material", young_modulus=YOUNG_MODULUS, poisson_ratio=POISSON_RATIO)
    sofa_node.addObject('HyperelasticForcefield', name="ff", printLog=False)

    sonics_node = node.addChild("sonics_node")
    sonics_node.addObject('StaticSolver', newton_iterations="15", relative_correction_tolerance_threshold="1e-15",
                          relative_residual_tolerance_threshold="1e-10", printLog="0")
    sonics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
    sonics_node.addObject('MechanicalObject', name="mo", position=position.tolist())
    sonics_node.addObject('CaribouTopology', name='topology', template=element_sonics,
                          indices=indices_sonics.tolist())
    sonics_node.addObject('FixedConstraint', indices=fixed_indices.tolist())
    sonics_node.addObject("UniformMass", totalMass=3, separateGravity=True)
    sonics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    sonics_node.addObject('FEniCS_Material', template=element_sonics, young_modulus=YOUNG_MODULUS,
                          poisson_ratio=POISSON_RATIO, material_name=material_model)
    sonics_node.addObject('HyperelasticForcefield_FEniCS', name="ff", printLog=False)


class TestHyperelasticForcefield(unittest.TestCase):

    def test_stiffness_assembly(self):
        for name, material, element in TEST_CASES:
            with self.subTest(msg=name):
                root = Sofa.Core.Node()
                createScene(root, element, material)
                Sofa.Simulation.init(root)
                K_sonics = csr_matrix(root.sonics_node.ff.K(), copy=True)
                K_sofa = csr_matrix(root.sofa_node.ff.K(), copy=True)
                if element == "Hexahedron20":
                    self.skipTest(
                        "FEniCS employs different " + str(
                            element) + " elements. Hence stiffness matrices cannot be equals")

                self.assertMatrixQuasiEqual(K_sonics.todense(), K_sofa.todense())

    def test_displacement(self):
        for name, material, element in TEST_CASES:
            with self.subTest(msg=name):
                root = Sofa.Core.Node()
                createScene(root, element, material)
                Sofa.Simulation.init(root)
                Sofa.Simulation.animate(root, root.dt.value)
                position_sofa = root.sofa_node.mo.position.value
                position_sonics = root.sonics_node.mo.position.value

                if element == "Hexahedron20":
                    position_sofa = position_sofa[np.unique(root.sofa_node.topology.indices.value[:, :8])]
                    position_sonics = position_sonics[np.unique(root.sonics_node.topology.indices.value[:, :8])]

                self.assertMatrixQuasiEqual(position_sofa, position_sonics)

    def test_fenics(self):
        # TODO is there a better way to test compatibility with fenics?
        if not any(Path(f"{current_dir}/meshes/").glob("*_solution.xdmf")):
            client = docker.from_env()
            image = client.images.pull("dolfinx/dolfinx:v0.5.1")
            container = client.containers.run(image, "python3 fenics_validation.py", volumes=[f'{current_dir}:/root'],
                                              remove=True, detach=True)
            container.stop()

        for name, material, element in TEST_CASES:
            with self.subTest(msg=name):
                root = Sofa.Core.Node()
                createScene(root, element, material)
                Sofa.Simulation.init(root)
                Sofa.Simulation.animate(root, root.dt.value)

                position_sonics = root.sonics_node.mo.position.value
                rest_position_sonics = root.sonics_node.mo.rest_position.value

                if element == "Hexahedron20":
                    position_sonics = position_sonics[np.unique(root.sonics_node.topology.indices.value[:, :8])]
                    rest_position_sonics = rest_position_sonics[
                        np.unique(root.sonics_node.topology.indices.value[:, :8])]

                if element == "Tetrahedron10":
                    position_sonics = position_sonics[np.unique(root.sonics_node.topology.indices.value[:, :4])]
                    rest_position_sonics = rest_position_sonics[
                        np.unique(root.sonics_node.topology.indices.value[:, :4])]

                with meshio.xdmf.TimeSeriesReader(f"{current_dir}/meshes/fenics_" + name + "_solution.xdmf") as reader:
                    points, _ = reader.read_points_cells()
                    _, point_data, _ = reader.read_data(0)

                fenics_rest_position = points
                position_fenics = np.empty(fenics_rest_position.shape)
                # TODO FEniCS and SOFA do not have the same dofs numbering, have to look for the closest point
                from scipy import spatial

                tree = spatial.KDTree(rest_position_sonics)
                for positions, displacement in zip(fenics_rest_position, point_data["f"]):
                    position_fenics[tree.query(positions)[1]] = positions + displacement

                self.assertMatrixQuasiEqual(position_fenics, position_sonics)

    def assertMatrixQuasiEqual(self, A, B):
        """ absolute(a - b) <= (atol + rtol * absolute(b)) """
        return self.assertTrue(np.allclose(A, B),
                               f"Matrices are not equal, with \nA={A} \nand\nB={B}")


if __name__ == '__main__':
    unittest.main()
