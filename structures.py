import ctypes
from ctypes import POINTER, Structure, c_int, c_double, c_float, c_uint

Pixel = c_double

class Sharpnesses(ctypes.Structure):
    _fields_ = [
        ("N", ctypes.c_int),
        ("sharpness", ctypes.POINTER(Pixel))
    ]

class Blur_Vector(ctypes.Structure):
    _fields_ = [
        ("angle", ctypes.c_int),
        ("magnitude", ctypes.c_float),
    ]

class Blur_Vector_Group(ctypes.Structure):
    _fields_ = [
        ("len_vectors", ctypes.c_int),
        # Assuming blur_vectors_rgb is a pointer to a pointer of Blur_Vector
        ("blur_vectors", ctypes.POINTER(Blur_Vector)),
    ]

class Pixel_HSV(Structure):
    _fields_ = [
        ("parent_id", c_int),
        ("h", c_double),
        ("s", c_double),
        ("v", c_double),
    ]

# Define the C structure in Python
class Image_RGB(ctypes.Structure):
    _fields_ = [
        ("height", ctypes.c_uint),
        ("width", ctypes.c_uint),
        ("r", ctypes.POINTER(ctypes.c_double)),
        ("g", ctypes.POINTER(ctypes.c_double)),
        ("b", ctypes.POINTER(ctypes.c_double)),
    ]


class Color_Palette(Structure):
    _fields_ = [
        ("N", c_int),
        ("averages", POINTER(Pixel_HSV)),
        ("percentages", POINTER(c_double)),
    ]


class RGB_Statistics(ctypes.Structure):
    _fields_ = [
        ("Br", ctypes.c_double),
        ("Bg", ctypes.c_double),
        ("Bb", ctypes.c_double),
        ("Cr", ctypes.c_double),
        ("Cg", ctypes.c_double),
        ("Cb", ctypes.c_double),
    ]

class Crop_Boundaries(ctypes.Structure):
    _fields_ = [
        ("N", ctypes.c_int),
        ("top", ctypes.POINTER(ctypes.c_int)),
        ("bottom", ctypes.POINTER(ctypes.c_int)),
        ("left", ctypes.POINTER(ctypes.c_int)),
        ("right", ctypes.POINTER(ctypes.c_int)),
    ]

class Blur_Profile(ctypes.Structure):
    _fields_ = [
        ("num_angle_bins", ctypes.c_int),
        ("num_radius_bins", ctypes.c_int),
        ("angle_bin_size", ctypes.c_int),
        ("radius_bin_size", ctypes.c_int),
        ("bins", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ]

    def get_bin_values(self):
        angle_bins = []
        for i in range(self.num_angle_bins):
            radius_array = self.bins[i]
            angle_bin = [radius_array[j] for j in range(self.num_radius_bins)]
            angle_bins.append(angle_bin)
        return angle_bins
    
    
class Full_Report_Data(ctypes.Structure):
    _fields_ = [
        ("rgb_stats", ctypes.POINTER(RGB_Statistics)),
        ("color_palette", ctypes.POINTER(Color_Palette)),
        ("blur_profile", ctypes.POINTER(Blur_Profile)),
        ("blur_vectors", ctypes.POINTER(Blur_Vector_Group)),
        ("average_saturation", ctypes.c_double),
        ("sharpness", ctypes.POINTER(Sharpnesses)),
    ]