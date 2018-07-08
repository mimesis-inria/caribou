from paraview.simple import *


class ParaView(object):
    class Camera(object):
        def __init__(self, **kwargs):
            # Parameters
            self.position = kwargs.get('position', [-150, 0, 40])
            self.focal_point = kwargs.get('focal_point', [0, 0, 40])
            self.up = kwargs.get('up', [0, 1, 0])
            self.orthogonal = kwargs.get('orthogonal', False)
            self.orthogonal_scale = kwargs.get('orthogonal_scale', 90)
            self.angle = kwargs.get('angle', 30)

    class Representation:
        Glyphs = '3D Glyphs'
        Outline = 'Outline'
        PointGaussian = 'Point Gaussian'
        Points = 'Points'
        Surface = 'Surface'
        SurfaceWithEdges = 'Surface With Edges'
        Volume = 'Volume'
        Wireframe = 'Wireframe'

    class MeshQuality:
        AspectBeta = 'Aspect Beta'
        AspectGamma = 'Aspect Gamma'
        AspectFrobenius = 'Aspect Frobenius'
        AspectRatio = 'Aspect Ratio'
        CollapseRatio = 'Collapse Ratio'
        Condition = 'Condition'
        Distortion = 'Distortion'
        Jacobian = 'Jacobian'
        MinimumDihedralAngle = 'Minimum Dihedral Angle'
        RadiusRatio = 'Radius Ratio'

    class View(object):
        def __init__(self, **kwargs):
            self.vtk_file = kwargs.get('vtk_file', None)

    class ElementsView(View):
        def __init__(self, **kwargs):
            ParaView.View.__init__(self, **kwargs)

            # Parameters
            self.intersect_with = kwargs.get('intersect_with', 'Plane')
            self.intersect_origin = kwargs.get('intersect_origin', [0, 0, 40])
            self.intersect_normal = kwargs.get('intersect_normal', [-1, 0, 0])
            self.display_quality = kwargs.get('display_quality', True)

        def render(self, view):
            reader = OpenDataFile(self.vtk_file)
            SetActiveSource(reader)
            vtkdisplay = Show(reader, view)
            view.ResetCamera()
            view.Update()
            filter = ExtractCellsByRegion(Input=reader)
            filter.IntersectWith = self.intersect_with
            filter.IntersectWith.Origin = self.intersect_origin
            filter.IntersectWith.Normal = self.intersect_normal
            filter.Extractintersected = True
            display = Show(filter, view)
            display.ScaleFactor = 8.0
            display.OSPRayScaleArray = 'gmsh:geometrical'
            Hide(reader, view)
            view.Update()

            if self.display_quality:
                meshQuality1 = MeshQuality(Input=filter)
                qualityLUT = GetColorTransferFunction('Quality')
                qualityPWF = GetOpacityTransferFunction('Quality')
                meshQuality1Display = Show(meshQuality1, view)
                meshQuality1Display.OSPRayScaleArray = 'Quality'
                meshQuality1Display.LookupTable = qualityLUT
                meshQuality1Display.ScalarOpacityFunction = qualityPWF
                meshQuality1.TetQualityMeasure = ParaView.MeshQuality.RadiusRatio
                Hide(filter, view)
                meshQuality1Display.SetScalarBarVisibility(view, False)
                view.Update()
                # qualityLUT.RescaleTransferFunction(1.0, 3.1449686647)
                # qualityPWF.RescaleTransferFunction(1.0, 3.1449686647)

    class SurfaceView(View):
        def __init__(self, **kwargs):
            ParaView.View.__init__(self, **kwargs)

            # Parameters
            self.representation = kwargs.get('representation', ParaView.Representation.Surface)
            self.color = kwargs.get('color', [0.0, 0.0, 0.0])
            self.opacity = kwargs.get('opacity', 1)
            self.line_width = kwargs.get('line_width', 1)

        def render(self, view):
            reader = OpenDataFile(self.vtk_file)
            SetActiveSource(reader)
            surface = Show(reader, view)
            surface.Representation = 'Surface'
            surface.SetRepresentationType(self.representation)
            surface.AmbientColor = self.color
            surface.Opacity = self.opacity
            surface.LineWidth = self.line_width

    def __init__(self, **kwargs):
        # Parameters
        self.size = kwargs.get('size', (1000, 400))
        self.background_color = kwargs.get('background_color', [1, 1, 1])
        self.camera = kwargs.get('camera', ParaView.Camera())
        self.transparent_background = kwargs.get('transparent_background', True)
        self.views = kwargs.get('views', [])

    def addView(self, view):
        assert isinstance(view, ParaView.View)

    def save(self, filename):
        width, height = self.size

        view = CreateRenderView()

        view.ViewSize = [width, height]

        for v in self.views:
            v.render(view)
            ResetCamera(view)

        view.Background = self.background_color
        view.CameraPosition = self.camera.position
        view.CameraFocalPoint = self.camera.focal_point
        view.CameraViewUp = self.camera.up
        view.CameraParallelProjection = self.camera.orthogonal
        view.CameraViewAngle = self.camera.angle
        view.OrientationAxesVisibility = False

        if self.camera.orthogonal:
            view.CameraParallelScale = self.camera.orthogonal_scale

        view.Update()

        # Render(view)
        SaveScreenshot(filename, TransparentBackground=self.transparent_background, view=view)
        Delete(view)
        RemoveViewsAndLayouts()


