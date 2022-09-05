import Sofa
import SofaCaribou

CONTACT = False


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):
        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin',
                       pluginName="SofaHaptics Geomagic SofaCaribou CImgPlugin SofaMiscCollision SofaGeneralSimpleFem SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine SofaGeneralDeformable SofaMiscFem SofaMeshCollision SofaGeneralLoader")

        root.gravity = [0, -9.81, 0]
        root.dt = 0.005
        root.addObject("DefaultPipeline", verbose=0, depth=6, draw=0)
        root.addObject("BruteForceBroadPhase")
        root.addObject("BVHNarrowPhase")
        root.addObject("DefaultContactManager", name="Response", response="FrictionContact")
        root.addObject("LocalMinDistance", name="Proximity", alarmDistance=0.15, contactDistance=0.05,
                       angleCone=0.5)
        root.addObject("FreeMotionAnimationLoop")
        root.addObject("LCPConstraintSolver", tolerance=0.001, maxIt=1000)
        root.addObject("GeomagicDriver", name="GeomagicDevice", deviceName="Default Device", scale=1,
                       drawDeviceFrame=True, positionBase=[0, 0, 0], drawDevice=False,
                       orientationBase=[0, 0.707, 0, -0.707])

        omni_node = root.addChild("Omni")
        omni_node.addObject("MechanicalObject", template="Rigid3d", name="DOFs",
                            position="@GeomagicDevice.positionDevice")
        omni_node.addObject("MechanicalStateController", template="Rigid3d", listening=True, mainDirection=[-1, 0, 0],
                            handleEventTriggersUpdate=True)

        omni_visual = omni_node.addChild("VisuAvatar", activated=False)
        omni_visual.addObject("MeshObjLoader", name="meshLoader_0", filename="mesh/sphere.obj", scale=0.1,
                              handleSeams=True)
        omni_visual.addObject("OglModel", name="Visual", src="@meshLoader_0", color="gray")
        omni_visual.addObject("RigidMapping", input="@..", output="@Visual", index=[0])

        omni_ref_model = omni_node.addChild("RefModel")
        omni_ref_model.addObject("MeshObjLoader", name="loader",
                                 filename="Demos/Dentistry/data/mesh/dental_instrument_centerline.obj")
        omni_ref_model.addObject("MeshTopology", src="@loader")
        omni_ref_model.addObject("MechanicalObject", src="@loader", name="instrumentRefState1", ry=-180, rz=-90, dz=3.5,
                                 dx=-0.3)
        omni_ref_model.addObject("RigidMapping")

        omni_ref_model_right = omni_node.addChild("RefModelRight")
        omni_ref_model_right.addObject("MeshObjLoader", name="loader",
                                       filename="Demos/Dentistry/data/mesh/dental_instrument_centerline.obj")
        omni_ref_model_right.addObject("MeshTopology", src="@loader")
        omni_ref_model_right.addObject("MechanicalObject", src="@loader", name="instrumentRefState2", ry=-180, rz=-90,
                                       dz=3.5, dx=-0.3, dy=0.5)
        omni_ref_model_right.addObject("RigidMapping")

        omni_ref_model_left = omni_node.addChild("RefModelLeft")
        omni_ref_model_left.addObject("MeshObjLoader", name="loader",
                                      filename="Demos/Dentistry/data/mesh/dental_instrument_centerline.obj")
        omni_ref_model_left.addObject("MeshTopology", src="@loader")
        omni_ref_model_left.addObject("MechanicalObject", src="@loader", name="instrumentRefState3", ry=-180, rz=-90,
                                      dz=3.5, dx=-0.3, dy=-0.5)
        omni_ref_model_left.addObject("RigidMapping")

        instrument_node = root.addChild("Instrument")
        instrument_node.addObject("EulerImplicitSolver", name="ODE solver", rayleighStiffness=0.05, rayleighMass=1.0)
        instrument_node.addObject("CGLinearSolver", iterations=25, tolerance=1e-10, threshold=10e-10)
        instrument_node.addObject("MechanicalObject", template="Rigid3d", name="instrumentState")
        instrument_node.addObject("UniformMass", name="mass", totalMass=0.5, showAxisSizeFactor=False)
        instrument_node.addObject("LCPForceFeedback", activate=True, forceCoef=0.1)
        instrument_node.addObject("UncoupledConstraintCorrection")

        instrument_visual = instrument_node.addChild("VisuModel")
        instrument_visual.addObject("MeshObjLoader", name="meshLoader_1",
                                    filename="Demos/Dentistry/data/mesh/dental_instrument.obj",
                                    handleSeams=True)
        instrument_visual.addObject("OglModel", name="InstrumentVisualModel", src="@meshLoader_1",
                                    color=[1.0, 0.2, 0.2, 1.0], ry=-180, rz=-90, dz=3.5, dx=-0.3)
        instrument_visual.addObject("RigidMapping", name="MM->VM mapping", input="@instrumentState",
                                    output="@InstrumentVisualModel")

        instrument_collision = instrument_node.addChild("CollisionModel")
        instrument_collision.addObject("MeshObjLoader", name="loader",
                                       filename="Demos/Dentistry/data/mesh/dental_instrument.obj")
        instrument_collision.addObject("MeshTopology", src="@loader")
        instrument_collision.addObject("MechanicalObject", name="instrumentCollisionState1", src="@loader", ry=-180,
                                       rz=-90, dz=3.5, dx=-0.3)
        instrument_collision.addObject("LineCollisionModel")
        instrument_collision.addObject("PointCollisionModel")
        instrument_collision.addObject("RigidMapping", name="MM->CM mapping", input="@instrumentState",
                                       output="@instrumentCollisionState1")

        instrument_ref_model_right = instrument_node.addChild("RefModelRight")
        instrument_ref_model_right.addObject("MeshObjLoader", name="loader",
                                             filename="Demos/Dentistry/data/mesh/dental_instrument_centerline.obj")
        instrument_ref_model_right.addObject("MeshTopology", src="@loader")
        instrument_ref_model_right.addObject("MechanicalObject", src="@loader", name="instrumentCollisionState2",
                                             ry=-180, rz=-90,
                                             dz=3.5, dx=-0.3, dy=0.5)
        instrument_ref_model_right.addObject("RigidMapping", name="MM->CM mapping", input="@instrumentState",
                                             output="@instrumentCollisionState2")

        instrument_ref_model_left = instrument_node.addChild("RefModelLeft")
        instrument_ref_model_left.addObject("MeshObjLoader", name="loader",
                                            filename="Demos/Dentistry/data/mesh/dental_instrument_centerline.obj")
        instrument_ref_model_left.addObject("MeshTopology", src="@loader")
        instrument_ref_model_left.addObject("MechanicalObject", src="@loader", name="instrumentCollisionState3",
                                            ry=-180, rz=-90,
                                            dz=3.5, dx=-0.3, dy=-0.5)
        instrument_ref_model_left.addObject("RigidMapping", name="MM->CM mapping", input="@instrumentState",
                                            output="@instrumentCollisionState3")

        instrument_node.addObject("VectorSpringForceField", template="Vec3d",
                                  object1="@Omni/RefModel/instrumentRefState1",
                                  object2="@Instrument/CollisionModel/instrumentCollisionState1", stiffness=10,
                                  viscosity=0)
        instrument_node.addObject("VectorSpringForceField", template="Vec3d",
                                  object1="@Omni/RefModelRight/instrumentRefState2",
                                  object2="@Instrument/RefModelRight/instrumentCollisionState2", stiffness=10,
                                  viscosity=0)
        instrument_node.addObject("VectorSpringForceField", template="Vec3d",
                                  object1="@Omni/RefModelLeft/instrumentRefState3",
                                  object2="@Instrument/RefModelLeft/instrumentCollisionState3", stiffness=10,
                                  viscosity=0)

        liver = root.addChild('liver')
        liver.addObject("EulerImplicitSolver", name="cg_odesolver")
        liver.addObject("CGLinearSolver", iterations=25, name="linear solver", tolerance=1.0e-9, threshold=1.0e-9)
        liver.addObject("MeshGmshLoader", name="loader", filename="mesh/liver.msh")
        liver.addObject('MechanicalObject', name="mo", src="@loader")
        liver.addObject('CaribouTopology', name='topology', template="Tetrahedron", position="@loader.position",
                        indices="@loader.tetrahedra")
        liver.addObject('CaribouMass', name='mass', topology='@topology', density=1, lumped=True)
        liver.addObject('FEniCS_Material', template="Tetrahedron", young_modulus=8000,
                        poisson_ratio=0.3, material_name="NeoHookean")
        # liver.addObject('FEniCS_Material', template="Tetrahedron", bulk_modulus=10e5, a=1180, b=8, a_f=18.5 * 10e4,
        #                 b_f=16, a_s=2.5 * 10e4, b_s=11.1, a_fs=2160, b_fs=11.4, material_name="Ogden")
        liver.addObject('HyperelasticForcefield_FEniCS', printLog=True)
        liver.addObject("PrecomputedConstraintCorrection", recompute=True)

        liver_collision = liver.addChild('collision')
        liver_collision.addObject("MeshTopology", src="@../loader")
        liver_collision.addObject("TriangleCollisionModel")
        liver_collision.addObject("LineCollisionModel")
        liver_collision.addObject("PointCollisionModel")
        liver.addObject('BoxROI', name="FixedIndices", box="-2 3 0 -1 4 1", drawBoxes=False)
        liver.addObject('FixedConstraint', name="FixedConstraint", indices="@FixedIndices.indices",
                        showObject=False)

        liver_visu = liver.addChild('visu')
        liver_visu.addObject("OglModel", name="liver_ogl", src="@../loader", texturename="textures/liver-texture-square.png")
        liver_visu.addObject("CaribouBarycentricMapping")


        return root


def onSimulationInitDoneEvent(self, event):
    pass


def onAnimateBeginEvent(self, event):
    pass


def onAnimateEndEvent(self, event):
    pass


def createScene(node):
    node.addObject(ControlFrame(node))


# Choose in your script to activate or not the GUI
USE_GUI = True


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    root = Sofa.Core.Node("root")
    createScene(root)
    Sofa.Simulation.init(root)

    if not USE_GUI:
        for iteration in range(10):
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()
