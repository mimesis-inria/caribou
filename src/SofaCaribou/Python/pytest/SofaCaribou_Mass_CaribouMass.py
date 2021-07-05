import Sofa, SofaRuntime
import sys
import unittest
import numpy as np
from scipy.sparse import csr_matrix
from pathlib import Path

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')

import SofaCaribou


class TestCaribouMass(unittest.TestCase):
    def test_mass_assembly(self):
        root = Sofa.Core.Node()
        root.addObject('RequiredPlugin', pluginName=[
            'SofaBaseMechanics', 'SofaEngine', 'SofaTopologyMapping'])
        root.addObject('RegularGridTopology', name='grid', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[3, 3, 9])
        root.addObject('MechanicalObject', name='mo', src='@grid')
        root.addObject('HexahedronSetTopologyContainer', name='topology', src='@grid')
        root.addObject('CaribouMass', name='mass', topology='@topology', density=2)
        Sofa.Simulation.init(root)
        M = root.mass.M()
        M_diag = root.mass.M_diag()
        self.assertAlmostEqual(M.sum(), M_diag.sum())
        self.assertAlmostEqual(M.sum(), 108000)
