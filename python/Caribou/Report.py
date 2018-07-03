import os
from math import pi as PI
import math
import tempfile
import base64

from .View import ParaView
from .Utils import bbox
from .Mesh import Mesh


class HtmlReport(object):
    def __init__(self, **kwargs):
        # Parameters
        self.name = kwargs.get('name')
        assert not self.name == ""

        # Members
        self.lines = []
        self.lines.append('<html>')
        self.lines.append('<head>')
        self.lines.append('<title>{}</title>'.format(self.name))

        self.lines.append("""
        <style>
            table {
                font-family: arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
            }
            
            td, th {
                border: 1px solid #dddddd;
                text-align: left;
                padding: 8px;
            }
            
            tr:nth-child(even) {
                background-color: #dddddd;
            }
        </style>
        """)

        self.lines.append('</head>')
        self.lines.append('<body>')
        self.lines.append('<h1>{}</h1>'.format(self.name))

    def add_section(self, name):
        self.lines.append('<h2>{}</h2>'.format(name))

    def add_list(self, name=None, attributes=[]):
        assert len(attributes)

        if name:
            self.lines.append('<h3>{}</h3>'.format(name))

        self.lines.append('<table>')
        self.lines.append("""
          <tr>
            <th>Name</th>
            <th>Value</th>
          </tr>
        """)
        for key, value in attributes:
            self.lines.append('<tr>')
            self.lines.append('  <td>{}</td>'.format(key))
            self.lines.append('  <td>{}</td>'.format(value))
            self.lines.append('</tr>')

        self.lines.append('</table>')

    def add_image(self, name=None, path=None, binary=False):
        if name is not None:
            self.lines.append('<h3>{}</h3>'.format(name))

        if not binary:
            self.lines.append('<img src="{}" alt="{}" width="100%"/>'.format(path, name))
        else:
            with open(path, "rb") as image_file:
                encoded_string = base64.b64encode(image_file.read())
                self.lines.append('<img src="data:image/gif;base64,{}" alt="{}" width="100%"/>'.format(encoded_string, name))

    def add_paragraph(self, name=None, text=None):
        if name is not None:
            self.lines.append('<h3>{}</h3>'.format(name))
        if text is not None:
            self.lines.append('<p>{}</p>'.format(text))

    def write(self, filepath):
        lines = list(self.lines)

        lines.append('</body>')
        lines.append('</html>')

        with open(filepath, 'w') as f:
            f.writelines(lines)
            f.flush()
            f.close()

    def add_image_from_meshes(self, name=None, meshes=[], view_attributes=[], image_width=1000):
        assert len(meshes) == len(view_attributes)
        tempdir = tempfile.gettempdir()
        tempfiles = []
        views = []

        xmin, xmax, ymin, ymax, zmin, zmax = (0, 0, 0, 0, 0, 0)
        i = 0
        for mesh in meshes:
            assert isinstance(mesh, Mesh)
            txmin, txmax, tymin, tymax, tzmin, tzmax = bbox(mesh.vertices)
            xmin, xmax, ymin, ymax, zmin, zmax = min(xmin, txmin), min(ymin, tymin), min(zmin, tzmin), max(xmax, txmax), max(ymax, tymax), max(zmax, tzmax)
            if mesh.filepath is None or not os.path.isfile(mesh.filepath):
                temp = os.path.join(tempdir, 'temp_{}.vtk'.format(i))
                mesh.save(temp)
                tempfiles.append(temp)

            views.append(ParaView.View(**dict({
                'vtk_file': mesh.filepath
            }, **view_attributes[i])))
            i = i + 1

        width, height, length = xmax - xmin, ymax - ymin, zmax - zmin
        image_height = int(height / width * image_width)
        camera_angle = 20
        camera_x = xmax + (length / 2. / math.tan(math.radians(camera_angle / 2.))) * 1.15
        camera_y = ymin + (height / 2.)
        camera_z = zmin + (length / 2.)

        temp_image = os.path.join(tempdir, 'temp_image.png')
        tempfiles.append(temp_image)

        ParaView(
            size=(image_width, image_height),
            camera=ParaView.Camera(
                angle=camera_angle,
                position=[-camera_x, camera_y, camera_z],
                focal_point=[xmax, camera_y, camera_z],
            ),
            views=views
        ).save(temp_image)

        self.add_image(name=name, path=temp_image, binary=True)

        for t in tempfiles:
            os.remove(t)
