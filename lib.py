# Import utility functions
import os

import ctypes

from .structures import (
    Sharpnesses,
    Blur_Vector,
    Blur_Vector_Group,
    Pixel_HSV,
    Image_RGB,
    Color_Palette,
    RGB_Statistics,
    Crop_Boundaries,
    Blur_Profile,
    Full_Report_Data
)

directory = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(directory, 'libreport_data.so')
lib = ctypes.CDLL(lib_path)

# Set restype and argtypes for loaded functions
lib.get_full_report_data.restype = ctypes.POINTER(Full_Report_Data)
lib.get_full_report_data.argtypes = [
    ctypes.POINTER(Image_RGB), ctypes.POINTER(Crop_Boundaries),
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_int,
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_float, ctypes.c_float,
    ctypes.c_double, ctypes.c_double, ctypes.c_int
]

lib.get_blur_profile_visual.restype = ctypes.POINTER(Image_RGB)
lib.get_blur_profile_visual.argtypes = [ctypes.POINTER(Blur_Profile), ctypes.c_int, ctypes.c_int]
