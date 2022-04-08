try:
    from Sofa.Core import ForceFieldVec3d as ForceField
except:
    from Sofa.Core import ForceField


class ConstantForceField(ForceField):
    def __init__(self, *args, **kwargs):
        ForceField.__init__(self, *args, **kwargs)
        self.forces = None

    def addForce(self, m, out_force, pos, vel):
        with out_force.writeable() as wa:
            wa[:] += self.forces
            
    def addDForce(self, df, dx, params):
        pass

    def addKToMatrix(self, mparams, nNodes, nDofs):
        pass