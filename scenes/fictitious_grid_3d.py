#!/usr/bin/python3

import Sofa

radius = 5
n = [10, 10, 10]
subdivisions = 4

mx = (radius / ((n[0])*pow(2, subdivisions)))/2
my = (radius / ((n[1])*pow(2, subdivisions)))/2
mz = (radius / ((n[2])*pow(2, subdivisions)))/2


def createScene(root):
    root.addObject('RequiredPlugin', pluginName=[
        'Sofa.Component.SceneUtility', # APIVersion
    ])
    root.addObject('APIVersion', level='23.06.99')
    root.addObject('DefaultAnimationLoop')
    root.addObject('SphereIsoSurface', radius=radius, center=[0, 0, 0])

    root.addObject('FictitiousGrid',
                   template='Vec3',
                   name='integration_grid',
                   n=n,
                   min=[-radius - mx, -radius - my, -radius - mz],
                   max=[+radius + mx, +radius + my, +radius + mz],
                   maximum_number_of_subdivision_levels=subdivisions,
                   printLog=True,
                   draw_boundary_cells=True,
                   draw_outside_cells=False,
                   draw_inside_cells=True
                   )
