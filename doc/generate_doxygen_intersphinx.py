#! python3

# https://sphobjinv.readthedocs.io/en/latest/syntax.html

import xml.etree.ElementTree as ET
import os
import zlib
import re

xml_output_directory = os.path.join(os.path.dirname(__file__), "Doxygen", "xml")
html_output_directory = os.path.join(os.path.dirname(__file__), "Doxygen", "html")

files = sorted(os.listdir(xml_output_directory))
lines = ""
for file in files:
    try:
        tree = ET.parse(os.path.join(xml_output_directory, file))
    except:
        continue

    root = tree.getroot()

    if not root.tag == "doxygen":
        continue

    for child in root.findall('compounddef'):
        id = child.get('id')
        kind = child.get('kind')
        language = child.get('language')
        dispname = child.find('compoundname').text
        name = dispname.replace(' ', '')
        if dispname == name:
            dispname = u'-'

        if kind in ['class', 'struct']:
            entry = (u'%s %s:%s %s %s %s\n' %
                     (name, language.lower().replace('+', 'p'), 'class', 1, f"{id}.html", dispname))
            lines += entry


def escape(string):
    return re.sub("\\s+", " ", string)


with open(os.path.join(html_output_directory, 'objects.inv'), 'wb') as f:
    # header
    f.write((u'# Sphinx inventory version 2\n'
             u'# Project: %s\n'
             u'# Version: %s\n'
             u'# The remainder of this file is compressed using zlib.\n' %
             (escape('caribou'),
              escape('1.0'))).encode('utf-8'))
    compressor = zlib.compressobj(9)
    # f.write(lines.encode('utf-8'))
    f.write(compressor.compress(lines.encode('utf-8')))
    f.write(compressor.flush())



