import warnings

try:
    import SimpleITK as sitk
except ImportError as v:
    warnings.warn("The python module SimpleITK is required.")
    exit(1)

try:
    import PIL.Image
except ImportError as v:
    warnings.warn("The python module Pillow is required.")
    exit(2)

try:
    import numpy as np
except ImportError as v:
    warnings.warn("The python module Numpy is required.")
    exit(3)


class MedicalImage(object):
    def __init__(self, filename):
        self.filename = filename
        self.image = sitk.ReadImage(filename)

    def takeslice(self, axe='x', index=0, x_reversed=False, y_reversed=False):
        s = None
        if axe == 'x':
            s = self.image[index, :, :]
        elif axe == 'y':
            s = self.image[:, index, :]
        else:
            s = self.image[:, :, index]

        if x_reversed:
            s = s[::-1, :]
        if y_reversed:
            s = s[:, ::-1]

        return s

    def points(self, img=None, min_value=None, max_value=None):
        if not img:
            img = self.image

        nda = sitk.GetArrayViewFromImage(img)

        spacing = np.array(img.GetSpacing())[::-1]
        origin = np.array(img.GetOrigin())[::-1]

        if min_value is not None and max_value is not None:
            indices = np.array(np.where((nda >= min_value) & (nda <= max_value)))
        elif min_value is not None:
            indices = np.array(np.where(nda >= min_value))
        elif max_value is not None:
            indices = np.array(np.where(nda <= max_value))
        else:
            indices = nda

        return indices.T * spacing + origin

    def dimensions(self, img=None):
        if not img:
            img = self.image

        nda = sitk.GetArrayViewFromImage(img)
        spacing = img.GetSpacing()
        shape = nda.shape[::-1]
        if nda.ndim is 2:
            n = (nx, ny) = (shape[0], shape[1])
            d = (dx, dy) = (spacing[0], spacing[1])
            size = (int(round(nx * dx)), int(round(ny * dy)))
        else:
            n = (nx, ny, nz) = (shape[0], shape[1], shape[2])
            d = (dx, dy, dz) = (spacing[0], spacing[1], spacing[2])
            size = (int(round(nx * dx)), int(round(ny * dy)), int(round(nz * dz)))

        return n, d, size

    def snapshot(self, img, min_value=-1024, max_value=1024):
        nda = sitk.GetArrayViewFromImage(img)
        if nda.ndim is not 2:
            raise RuntimeError('Can only show 2D images')

        n, d , size = self.dimensions(img)
        image = np.array(nda, copy=True)
        image.clip(min_value, max_value, out=image)
        image -= min_value
        image //= (max_value - min_value) / 255
        image = image.astype(np.uint8)
        im = PIL.Image.fromarray(obj=image, mode='L').resize(size=size)
        return im
