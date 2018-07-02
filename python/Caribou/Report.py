import os


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

    def add_image(self, name=None, path=None):
        if name is not None:
            self.lines.append('<h3>{}</h3>'.format(name))
        self.lines.append('<img src="{}" alt="{}" width="100%"/>'.format(path, name))

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
