import ctypes
from types import SimpleNamespace
import json
import sys
from math import cos, sin, radians
from PIL import Image, ImageDraw, ImageFont
import tkinter as tk
from PIL import ImageTk
import numpy as np
import time
from ctypes import POINTER, Structure, c_int, c_double

# Load the shared library
lib = ctypes.CDLL('./build/libreport_data.so')

Pixel = ctypes.c_double

import ctypes

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


lib.get_full_report_data.restype = ctypes.POINTER(Full_Report_Data)
lib.get_full_report_data.argtypes = [ctypes.POINTER(Image_RGB), ctypes.POINTER(Crop_Boundaries),
                                     ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                     ctypes.c_double, ctypes.c_double,
                                     ctypes.c_double, ctypes.c_int,
                                     ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                     ctypes.c_float, ctypes.c_float,
                                     ctypes.c_double, ctypes.c_double, ctypes.c_int
                                     ]
lib.get_blur_profile_visual.restype = ctypes.POINTER(Image_RGB)
lib.get_blur_profile_visual.argtypes = [ctypes.POINTER(Blur_Profile), ctypes.c_int, ctypes.c_int]

class Report:
    def __init__(self, report_ptr, height, width):
        report_data = report_ptr.contents
        self.data_ptr = report_ptr

        self.rgb_stats = self._convert_rgb_statistics(report_data.rgb_stats)
        self.rgb_stats.height = height
        self.rgb_stats.width = width
        self.color_palette = self._convert_color_palette(report_data.color_palette)
        self.blur_profile = self._convert_blur_profile(report_data.blur_profile)
        self.blur_vectors = self._convert_blur_vectors(report_data.blur_vectors)
        self.average_saturation = report_data.average_saturation
        self.sharpnesses = self._convert_sharpnesses(report_data.sharpness)


    def _convert_sharpnesses(self, sharpnesses_ptr):
        # Ensure the pointer is valid
        if not sharpnesses_ptr:
            return []

        # Dereference the pointer to get the Sharpnesses structure
        c_sharpnesses = sharpnesses_ptr.contents

        # Get the number of sharpness values
        n = c_sharpnesses.N

        # Extract sharpness array into a Python list
        sharpnesses = [c_sharpnesses.sharpness[i] for i in range(n)]

        return sharpnesses


    def _convert_blur_vectors(self, blur_vector_pointer):
        c_blur_vector_group = blur_vector_pointer.contents
        blur_vectors = []

        if c_blur_vector_group.len_vectors > 0:
                c_blur_vectors = c_blur_vector_group.blur_vectors
                for i in range(c_blur_vector_group.len_vectors):
                    c_blur_vector = c_blur_vectors[i]
                    vector = SimpleNamespace()
                    vector.angle = c_blur_vector.angle
                    vector.magnitude = c_blur_vector.magnitude
                    blur_vectors.append(vector)
        return blur_vectors


    def _convert_rgb_statistics(self, rgb_stats_ptr):
        rgb_stats = rgb_stats_ptr.contents
        return rgb_stats
    

    def _convert_color_palette(self, color_palette_ptr):
        color_palette = color_palette_ptr.contents

        # Since averages is a pointer to Pixel_HSV, cast it to the right type
        averages = ctypes.cast(color_palette.averages, POINTER(Pixel_HSV * color_palette.N)).contents
        
        # Convert HSV to RGB and store in a list
        rgb_averages = [hsv_to_rgb(pixel.h, pixel.s, pixel.v) for pixel in averages]
        
        percentages = [color_palette.percentages[i] for i in range(color_palette.N)]

        color_palette.colors = rgb_averages

        color_palette.quantities = percentages

        return color_palette
    

    def _convert_blur_profile(self, blur_profile_ptr):
        c_blur_profile = blur_profile_ptr.contents
        num_angle_bins = c_blur_profile.num_angle_bins
        num_radius_bins = c_blur_profile.num_radius_bins
        blur_profile = SimpleNamespace()

        def _2d_array_from_pointer(ptr, rows, cols):
            array_2d = []
            for i in range(rows):
                row_ptr = ctypes.cast(ptr[i], ctypes.POINTER(ctypes.c_double * cols)).contents
                array_2d.append(list(row_ptr))
            return array_2d

        bins = _2d_array_from_pointer(c_blur_profile.bins, num_angle_bins, num_radius_bins)

        def check_and_correct_nan(bins):
            for angle_index, row in enumerate(bins):
                for radius_index, val in enumerate(row):
                    if np.isnan(val):
                        print(f"NaN found at angle {angle_index}, radius {radius_index}. Correcting to 0.")
                        bins[angle_index][radius_index] = 0.0  # Correct NaN to 0
            return bins

        blur_profile.bins = check_and_correct_nan(bins)

        return blur_profile
    

    def generate_color_palette_image(self):
        num_colors = len(self.color_palette.colors)
        block_size = 50  # Size of each color block

        # Calculate the number of blocks per row to determine the image size
        num_blocks_per_row = int(np.ceil(np.sqrt(num_colors)))
        img_width = num_blocks_per_row * block_size
        img_height = ((num_colors + num_blocks_per_row - 1) // num_blocks_per_row) * block_size

        # Create a new image with a white background
        img = Image.new('RGB', (img_width, img_height), 'black')
        draw = ImageDraw.Draw(img)
        font = ImageFont.load_default()  # Load the default font

        for i, (color, quantity) in enumerate(zip(self.color_palette.colors, self.color_palette.quantities)):
            # Calculate the position of the block
            row = i // num_blocks_per_row
            col = i % num_blocks_per_row
            x1 = col * block_size
            y1 = row * block_size
            x2 = x1 + block_size
            y2 = y1 + block_size

            # Convert color to 0-255 range and draw the block
            display_color = tuple(int(c) for c in color)
            draw.rectangle([x1, y1, x2, y2], fill=display_color)

            # Draw the quantity text in the center of the block
            text = f"{quantity:.1%}"  # Convert to percentage representation
            text_width, text_height = draw.textbbox((0, 0), text, font=font)[2:]
            text_x = x1 + (block_size - text_width) / 2
            text_y = y1 + (block_size - text_height) / 2
            draw.text((text_x, text_y), text, fill='black', font=font)
        self.color_palette_image = img
    

    def display_color_palette_image(self):
        # Display image in window
        window = tk.Tk()
        tk_image = ImageTk.PhotoImage(self.color_palette_image)
        label = tk.Label(window, image=tk_image)
        label.pack()
        window.mainloop()


    def generate_blur_profile_image(self):
        bp_ptr = ctypes.byref(self.blur_profile)
        height = self.rgb_stats.height
        width = self.rgb_stats.width
        image_rgb_ptr = lib.get_blur_profile_visual(bp_ptr, height, width)

        img = image_rgb_to_pillow(image_rgb_ptr, width, height)

        # Crop the image to remove the blank right half
        self.blur_profile_image = img.crop((0, 0, width // 2, height))

    
    def display_blur_profile(self, path):
        # Display the image
        window = tk.Tk()
        screen_width = int(window.winfo_screenwidth() * 0.9)
        screen_height = int(window.winfo_screenheight() * 0.9)

        # Calculate the new size to fit the screen while maintaining aspect ratio
        img_ratio = self.blur_profile_image.width / self.blur_profile_image.height
        screen_ratio = screen_width / screen_height
        if img_ratio > screen_ratio:
            new_width = screen_width
            new_height = int(screen_width / img_ratio)
        else:
            new_height = screen_height
            new_width = int(screen_height * img_ratio)

        # Resize the image
        img_resized = self.blur_profile_image.resize((new_width, new_height), Image.Resampling.LANCZOS)
        self.blur_profile_image.save(path)  # Save the cropped image
        # Display the image
        tk_image = ImageTk.PhotoImage(img_resized)
        label = tk.Label(window, image=tk_image)
        label.pack()

        window.mainloop()


    def display_all(self):
        """
            Creates a tkinter window, displaying all report statistics, the input image, and
                a visual color palette.
            REQUIRED:
                -self.image must be manually set to the input image by the programmer.
                -self.generate_color_palette_image() MUST be called prior to self.display_all().
                -self.bounding_boxes must be manually set IF the programmer wishes to see 
                    the sharpness of bounding boxes in the display window.
        """
        window = tk.Tk()
        window.title("Image Analysis Report")

        # Set up the main image
        max_width = window.winfo_screenwidth() * 0.8
        max_height = window.winfo_screenheight() * 0.9
        scale_width = max_width / self.image.width
        scale_height = max_height / self.image.height
        scale_factor = min(scale_width, scale_height)
        resized_image = self.image.resize((int(self.image.width * scale_factor), int(self.image.height * scale_factor)), Image.Resampling.LANCZOS)
        image_photo = ImageTk.PhotoImage(resized_image)

        image_center_x = image_photo.width() // 2
        image_center_y = image_photo.height() // 2

        # Create a canvas for the main image and draw on it
        canvas = tk.Canvas(window, width=image_photo.width(), height=image_photo.height())
        canvas.pack(side='left', padx=10)
        canvas.create_image(0, 0, anchor='nw', image=image_photo)

        # Draw arrows representing blur vectors
        length_scale_factor = min(image_photo.width()/2, image_photo.height()/2)
        for vector in self.blur_vectors:
            arrow_angle = vector.angle
            arrow_length = (vector.magnitude * length_scale_factor)
            end_x = image_center_x + arrow_length * cos(radians(arrow_angle))
            end_y = image_center_y - arrow_length * sin(radians(arrow_angle))
            canvas.create_line(image_center_x, image_center_y, end_x, end_y, arrow='last', fill='red', width=2)

        # Add bounding box lines and sharpness text
        if "bounding_boxes" in dir(self):
            for i in range(self.bounding_boxes.N):
                x0 = int(self.bounding_boxes.left[i] * scale_factor)
                y0 = int(self.bounding_boxes.top[i] * scale_factor)
                x1 = int(self.bounding_boxes.right[i] * scale_factor)
                y1 = int(self.bounding_boxes.bottom[i] * scale_factor)

                # Draw bounding box lines
                canvas.create_line(x0, y0, x1, y0, fill='red', width=2)  # Top line
                canvas.create_line(x0, y1, x1, y1, fill='red', width=2)  # Bottom line
                canvas.create_line(x0, y0, x0, y1, fill='red', width=2)  # Left line
                canvas.create_line(x1, y0, x1, y1, fill='red', width=2)  # Right line

                # Calculate the position for the text
                text_x = (x0 + x1) / 2
                text_y = y0 - 10  # Adjust this value to position the text slightly above the top line

                # Add text showing the sharpness
                sharpness = self.sharpnesses[i]
                scale = 10**4
                scaled_value = round(sharpness * scale)
                sharpness_value = f"Sharpness: {scaled_value / scale:.4f}"
                canvas.create_text(text_x, text_y, text=sharpness_value, fill='red', anchor="s")

        canvas.image = image_photo  # Keep a reference

        # Display RGB stats, average saturation, and sharpness
        stats_text = f"\nRed Brightness: {self.rgb_stats.Br}\n"
        stats_text += f"Green Brightness: {self.rgb_stats.Bg}\n"
        stats_text += f"Blue Brightness: {self.rgb_stats.Bb}\n"
        stats_text += f"Red Contrast: {self.rgb_stats.Cr}\n"
        stats_text += f"Green Contrast: {self.rgb_stats.Cg}\n"
        stats_text += f"Blue Contrast: {self.rgb_stats.Cb}\n"
        stats_text += f"Saturation: {self.average_saturation}\n"
        stats_label = tk.Label(window, text=stats_text, justify=tk.LEFT)
        stats_label.pack(side='top', padx=10)

        # Display color palette and blur profile images
        color_palette_photo = ImageTk.PhotoImage(self.color_palette_image)
        color_palette_label = tk.Label(window, image=color_palette_photo)
        color_palette_label.pack(side='top', padx=10)
        color_palette_label.image = color_palette_photo  # Keep a reference

        window.mainloop()


    def to_json(self):
        # Define the maximum number of color entries
        max_color_entries = 100
        max_vector_entries = 10
        max_sharpnesses = 10

        # Start with basic statistics and blur profile data
        report_data = {
            'Height': self.rgb_stats.height,
            'Width': self.rgb_stats.width,
            'Average Saturation': self.average_saturation,
            'Red Brightness': self.rgb_stats.Br,
            'Green Brightness': self.rgb_stats.Bg,
            'Blue Brightness': self.rgb_stats.Bb,
            'Red Contrast': self.rgb_stats.Cr,
            'Green Contrast': self.rgb_stats.Cg,
            'Blue Contrast': self.rgb_stats.Cb,
        }
        for i in range(max_vector_entries):
            angle = self.blur_vectors[i].angle
            magnitude = self.blur_vectors[i].magnitude
            report_data[f'Blur Vector {i+1} Angle'] = angle
            report_data[f'Blur Vector {i+1} Magnitude'] = magnitude


        # Add color palette data, padding with zeros if necessary
        for i in range(max_color_entries):
            if i < len(self.color_palette.colors):
                h, s, v = self.color_palette.colors[i]
                percentage = self.color_palette.quantities[i]
            else:
                h, s, v, percentage = 0, 0, 0, 0

            report_data[f'Color {i+1} H'] = h
            report_data[f'Color {i+1} S'] = s
            report_data[f'Color {i+1} V'] = v
            report_data[f'Color {i+1} Percentage'] = percentage

        # Add sharpnesses list, padding with zeros if necessary
        for i in range(max_sharpnesses):
            if i < len(self.sharpnesses):
                sharpness = self.sharpnesses[i]
            else:
                sharpness = 0.0
            report_data[f'Sharpness {i+1}:'] = sharpness

        # Convert the report data to a JSON string
        json_data = json.dumps(report_data, indent=4)
        return json_data
    
    def __del__(self):
        lib.free_full_report(ctypes.byref(self.data_ptr))


def set_bounding_boxes(bounding_boxes):
    """
    bounding_boxes should be a list of tuples or a list of dictionaries,
    where each tuple or dictionary represents one bounding box
    with 'top', 'bottom', 'left', 'right' values.
    """
    n = len(bounding_boxes)
    top_array = (ctypes.c_int * n)()
    bottom_array = (ctypes.c_int * n)()
    left_array = (ctypes.c_int * n)()
    right_array = (ctypes.c_int * n)()

    for i, bbox in enumerate(bounding_boxes):
        top_array[i] = bbox['top']
        bottom_array[i] = bbox['bottom']
        left_array[i] = bbox['left']
        right_array[i] = bbox['right']

    crop_boundaries = Crop_Boundaries(
        N=n,
        top=top_array,
        bottom=bottom_array,
        left=left_array,
        right=right_array
    )

    return crop_boundaries


def get_report(pil_image: Image, salient_characters=ctypes.POINTER(Crop_Boundaries)(),
               h_partitions=18, s_partitions=2, v_partitions=3,
               black_thresh=0.1, gray_thresh=0.1,
               coverage_thresh=0.95, linked_list_size=1000, downsample_rate=1,
               radius_partitions=40, angle_partitions=72,
               quantity_weight=0.1, saturation_value_weight=0.9,
               fft_streak_thresh=1.15, magnitude_thresh=0.3, blur_cutoff_ratio_denom=2):
    # Convert PIL image to Image_RGB
    width = pil_image.width
    height = pil_image.height
    crop_boundaries = salient_characters
    image_rgb = pil_image_to_image_rgb(pil_image)

    image_ctype = ctypes.byref(image_rgb)
    start_time = time.time()  # Start timing

    report_data_ptr = lib.get_full_report_data(image_ctype, crop_boundaries,
                                               h_partitions, s_partitions, v_partitions,
                                               black_thresh, gray_thresh,
                                               coverage_thresh, linked_list_size, downsample_rate,
                                               radius_partitions, angle_partitions,
                                               quantity_weight, saturation_value_weight,
                                               fft_streak_thresh, magnitude_thresh, 
                                               blur_cutoff_ratio_denom
                                               )
    
    end_time = time.time()  # End timing
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    # Check if the report data pointer is not null or an unexpected return type
    if report_data_ptr is None or isinstance(report_data_ptr, int):
        print("Failed to get report data")
        return None

    # Convert the returned Full_Report_Data pointer to a Python object
    report = Report(report_data_ptr, height, width)

    return report


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


def run_demonstration():
    global image_name
    image_path = "images/original/"
    # Check if a custom image path was provided as a command-line argument
    if len(sys.argv) > 1:
        image_name = sys.argv[1]
    else:
        image_name = "10.png"

    # Open image
    image = Image.open(image_path+image_name)

    bounds = [
        {
            "left": 61,
            "right": 383,
            "top": 212,
            "bottom": 897,
        },
        {
            "left": 363,
            "right": 591,
            "top": 130,
            "bottom": 805,
        },
        {
            "left": 467,
            "right": 944,
            "top": 94,
            "bottom": 996,
        },
    ]
    bounding_boxes = set_bounding_boxes(bounds)


    # Set up report variables, and Generate report
    report = get_report(image, salient_characters=bounding_boxes)
    # report = get_report(image)
    report.image = image

    # Display images that represent image: blur profile and color palette
    report.generate_color_palette_image()
    report.bounding_boxes = bounding_boxes
    report.display_all()
    json = report.to_json()
    print(json)


if __name__ == "__main__":
    run_demonstration()
