import ctypes
import numpy as np
from PIL import Image
from .structures import Image_RGB, Crop_Boundaries, Image_PGM


def hsv_to_rgb(h, s, v):
    h, s, v = h, s, v
    c = v * s
    x = c * (1 - abs((h / 60) % 2 - 1))
    m = v - c

    if h < 60:
        r, g, b = c, x, 0
    elif h < 120:
        r, g, b = x, c, 0
    elif h < 180:
        r, g, b = 0, c, x
    elif h < 240:
        r, g, b = 0, x, c
    elif h < 300:
        r, g, b = x, 0, c
    else:
        r, g, b = c, 0, x

    r, g, b = (r + m) * 255, (g + m) * 255, (b + m) * 255
    return int(r), int(g), int(b)

# Function to convert PIL image to Image_RGB
def pil_image_to_image_rgb(pil_image):
    width, height = pil_image.size
    img_array = np.array(pil_image) / 255.0
    r, g, b = img_array[:,:,0], img_array[:,:,1], img_array[:,:,2]

    r_flat = r.flatten().astype(np.double)
    g_flat = g.flatten().astype(np.double)
    b_flat = b.flatten().astype(np.double)

    r_ctypes = r_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    g_ctypes = g_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    b_ctypes = b_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    pil_image.r_ctypes = r_ctypes
    pil_image.g_ctypes = g_ctypes
    pil_image.b_ctypes = b_ctypes

    return Image_RGB(height=height, width=width, r=r_ctypes, g=g_ctypes, b=b_ctypes)


def image_rgb_to_pillow(image_rgb_ptr, width, height):
    # Assuming image_rgb is an instance of Image_RGB and already filled with data
    # Convert ctypes pointers to numpy arrays
    image_rgb = image_rgb_ptr.contents

    DoubleArray = ctypes.c_double * (width * height)

    r_array = ctypes.cast(image_rgb.r, ctypes.POINTER(DoubleArray)).contents
    g_array = ctypes.cast(image_rgb.g, ctypes.POINTER(DoubleArray)).contents
    b_array = ctypes.cast(image_rgb.b, ctypes.POINTER(DoubleArray)).contents

    r_np = np.ctypeslib.as_array(r_array).reshape(height, width)
    g_np = np.ctypeslib.as_array(g_array).reshape(height, width)
    b_np = np.ctypeslib.as_array(b_array).reshape(height, width)

    img_np = np.stack((r_np, g_np, b_np), axis=-1) * 255
    img_np = np.clip(img_np, 0, 255).astype(np.uint8)

    # Convert the numpy array to a PIL Image and return
    return Image.fromarray(img_np, 'RGB')


def image_pgm_to_pillow(image_pgm_ptr, width, height):
    # Assuming image_pgm is an instance of Image_PGM and already filled with data
    # Convert ctypes pointer to numpy array
    image_pgm = image_pgm_ptr.contents

    DoubleArray = ctypes.c_double * (width * height)
    data_array = ctypes.cast(image_pgm.data, ctypes.POINTER(DoubleArray)).contents

    # Convert to a numpy array and reshape it to match image dimensions
    data_np = np.ctypeslib.as_array(data_array).reshape(height, width)

    # Scale the pixel values to the 0-255 range and convert to uint8 for image representation
    img_np = np.clip(data_np * 255, 0, 255).astype(np.uint8)

    # Convert the numpy array to a PIL Image and return
    return Image.fromarray(img_np, 'L')  # 'L' mode for grayscale image
