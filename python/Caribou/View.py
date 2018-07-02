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

    class View(object):
        def __init__(self, **kwargs):
            # Parameters
            self.vtk_file = kwargs.get('vtk_file', None)
            self.representation = kwargs.get('representation', ParaView.Representation.Surface)
            self.color = kwargs.get('color', [0.0, 0.0, 0.0])
            self.opacity = kwargs.get('opacity', 1)
            self.line_width = kwargs.get('line_width', 1)

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
            reader = OpenDataFile(v.vtk_file)
            SetActiveSource(reader)
            surface = Show(reader, view)
            surface.Representation = 'Surface'
            surface.SetRepresentationType(v.representation)
            surface.AmbientColor = v.color
            surface.Opacity = v.opacity
            surface.LineWidth = v.line_width
            ResetCamera(view)
            # Render(view)

        view.Background = self.background_color
        view.CameraPosition = self.camera.position
        view.CameraFocalPoint = self.camera.focal_point
        view.CameraViewUp = self.camera.up
        view.CameraParallelProjection = self.camera.orthogonal
        view.CameraViewAngle = self.camera.angle
        view.OrientationAxesVisibility = False

        if self.camera.orthogonal:
            view.CameraParallelScale = self.camera.orthogonal_scale

        Render(view)
        SaveScreenshot(filename, TransparentBackground=self.transparent_background, view=view)


