import os
import warnings

try:
    from libSofaPython import Sofa
except ImportError as v:
    warnings.warn("The python module Sofa is required.")


def init(prefix=None, config_dir=None, screenshot_dir=None):
    if not prefix:
        prefix = os.environ.get('SOFA_ROOT', os.path.curdir)
    if not config_dir:
        config_dir = os.path.join(prefix, "config")
    if not screenshot_dir:
        screenshot_dir = os.path.join(prefix, "screenshots")

    Sofa.GUIManager.setSofaPrefix(prefix)
    Sofa.GUIManager.setConfigDirectoryPath(config_dir)
    Sofa.GUIManager.setScreenshotDirectoryPath(screenshot_dir)

    Sofa.GUIManager.init("qglviewer")
    Sofa.GUIManager.createGUI()


def launch(root, filename=None, width=800, height=600):
    Sofa.GUIManager.setDimension(width, height)
    Sofa.GUIManager.setScene(root, filename)
    Sofa.GUIManager.MainLoop(root)
    Sofa.GUIManager.closeGUI()
