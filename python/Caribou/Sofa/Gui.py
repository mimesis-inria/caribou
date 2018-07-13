import os
from libSofaPython import Sofa as sofa


def init(prefix=None, config_dir=None, screenshot_dir=None):
    if not prefix:
        prefix = os.environ.get('SOFA_ROOT', os.path.curdir)
    if not config_dir:
        config_dir = os.path.join(prefix, "config")
    if not screenshot_dir:
        screenshot_dir = os.path.join(prefix, "screenshots")

    sofa.GUIManager.setSofaPrefix(prefix)
    sofa.GUIManager.setConfigDirectoryPath(config_dir)
    sofa.GUIManager.setScreenshotDirectoryPath(screenshot_dir)

    sofa.GUIManager.init("qglviewer")
    sofa.GUIManager.createGUI()


def launch(root, filename=None, width=800, height=600):
    sofa.GUIManager.setDimension(width, height)
    sofa.GUIManager.setScene(root, filename)
    sofa.GUIManager.MainLoop(root)
    sofa.GUIManager.closeGUI()
