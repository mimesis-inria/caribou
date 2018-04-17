#!/usr/bin/env python2

import sys
from libSofaPython import Sofa as sofa


def on_animate_begin_event():
    print "Beginning the step"


def create_scene(root):
    root.createObject('APIVersion', name=17.12)
    root.createObject('VisualStyle', displayFlags="showBehaviorModels showCollisionModels")
    root.bbox = "0 0 0 10 10 10"


if __name__ == '__main__':
    from Benchmark import Gui
    from Benchmark.Simulation import Simulation

    simulation = Simulation()

    simulation.load()
    Gui.init()

    sofa.loadPlugin("caribou")

    create_scene(simulation.root)
    binder = simulation.root.createObject('PythonEventBinder')
    simulation.init()
    simulation.root.initVisual()

    binder.bind('AnimateBeginEvent', on_animate_begin_event)

    #simulation.step()
    #simulation.step()

    Gui.launch(simulation.root, filename=sys.argv[0])
