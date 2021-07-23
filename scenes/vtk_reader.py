#!/usr/bin/python3

import Sofa


def createScene(root):
    
    root.addObject('APIVersion', level='21.06')
    
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaLoader')

    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    root.addObject('MeshObjLoader', name='LiverSurface', filename='/media/Sidaty/Data/sidaty/external_plugins/caribou/scenes/liver_tom.obj')

    visu = root.addChild('Visu', tags="Visual", gravity="0 -9.81 0")
    visu.addObject('OglModel', name="VisualModel", color="red", src="@../../LiverSurface")
    
