import csv
import io
import ctypes
import json
import sys
from math import atan2, pi, sqrt, cos, sin, radians
from PIL import Image, ImageDraw, ImageFont
import tkinter as tk
from PIL import ImageTk
from statistics import stdev
import numpy as np
import time
from ctypes import POINTER, Structure, c_int, c_double

# Load the shared library
lib = ctypes.CDLL('./build/libreport_data.so')

Pixel = ctypes.c_double

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

class Blur_Profile_RGB(ctypes.Structure):
    _fields_ = [
        ("num_angle_bins", ctypes.c_int),
        ("num_radius_bins", ctypes.c_int),
        ("angle_bin_size", ctypes.c_int),
        ("radius_bin_size", ctypes.c_int),
        ("r_bins", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("g_bins", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("b_bins", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ]

    def get_bin_values(self):
        angle_bins_r = []
        angle_bins_g = []
        angle_bins_b = []

        for i in range(self.num_angle_bins):
            # Red channel
            radius_array_r = self.r_bins[i]
            angle_bin_r = [radius_array_r[j] for j in range(self.num_radius_bins)]
            angle_bins_r.append(angle_bin_r)

            # Green channel
            radius_array_g = self.g_bins[i]
            angle_bin_g = [radius_array_g[j] for j in range(self.num_radius_bins)]
            angle_bins_g.append(angle_bin_g)

            # Blue channel
            radius_array_b = self.b_bins[i]
            angle_bin_b = [radius_array_b[j] for j in range(self.num_radius_bins)]
            angle_bins_b.append(angle_bin_b)
        
        return {
            "r_bins": angle_bins_r,
            "g_bins": angle_bins_g,
            "b_bins": angle_bins_b
        }
    
class Full_Report_Data(ctypes.Structure):
    _fields_ = [
        ("rgb_stats", ctypes.POINTER(RGB_Statistics)),
        ("color_palette", ctypes.POINTER(Color_Palette)),
        ("blur_profile", ctypes.POINTER(Blur_Profile_RGB)),
        ("average_saturation", ctypes.c_double),
        ("sharpness", ctypes.c_double),
    ]


lib.get_full_report_data.restype = ctypes.POINTER(Full_Report_Data)
lib.get_full_report_data.argtypes = [ctypes.POINTER(Image_RGB), 
                                     ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                     ctypes.c_double, ctypes.c_double,
                                     ctypes.c_double, ctypes.c_int,
                                     ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                     ctypes.c_float, ctypes.c_float]
lib.get_blur_profile_visual.restype = ctypes.POINTER(Image_RGB)
lib.get_blur_profile_visual.argtypes = [ctypes.POINTER(Blur_Profile_RGB), ctypes.c_int, ctypes.c_int]

class Report:
    def __init__(self, report_ptr, height, width):
        # Use contents to get the actual data if report_ptr is a pointer
        # Use the [] operator if report_ptr is an array of structures
        report_data = report_ptr.contents

        self.rgb_stats = self._convert_rgb_statistics(report_data.rgb_stats)
        self.rgb_stats.height = height
        self.rgb_stats.width = width
        self.color_palette = self._convert_color_palette(report_data.color_palette)
        self.blur_profile = self._convert_blur_profile(report_data.blur_profile)
        self.average_saturation = report_data.average_saturation
        self.sharpness = report_data.sharpness


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
        blur_profile = blur_profile_ptr.contents
        num_angle_bins = blur_profile.num_angle_bins
        num_radius_bins = blur_profile.num_radius_bins

        def _2d_array_from_pointer(ptr, rows, cols):
            array_2d = []
            for i in range(rows):
                row_ptr = ctypes.cast(ptr[i], ctypes.POINTER(ctypes.c_double * cols)).contents
                array_2d.append(list(row_ptr))
            return array_2d

        r_bins = _2d_array_from_pointer(blur_profile.r_bins, num_angle_bins, num_radius_bins)
        g_bins = _2d_array_from_pointer(blur_profile.g_bins, num_angle_bins, num_radius_bins)
        b_bins = _2d_array_from_pointer(blur_profile.b_bins, num_angle_bins, num_radius_bins)

        def check_and_correct_nan(bins, color):
            for angle_index, row in enumerate(bins):
                for radius_index, val in enumerate(row):
                    if np.isnan(val):
                        print(f"NaN found in {color} bins at angle {angle_index}, radius {radius_index}. Correcting to 0.")
                        bins[angle_index][radius_index] = 0.0  # Correct NaN to 0
            return bins

        blur_profile.r = check_and_correct_nan(r_bins, 'R')
        blur_profile.g = check_and_correct_nan(g_bins, 'G')
        blur_profile.b = check_and_correct_nan(b_bins, 'B')

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

    
    def display_blur_profile(self):
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
        save_path = "images/output/blur_profile_" + image_name
        self.blur_profile_image.save(save_path)  # Save the cropped image
        # Display the image
        tk_image = ImageTk.PhotoImage(img_resized)
        label = tk.Label(window, image=tk_image)
        label.pack()

        window.mainloop()


    def vectorize_blur_profile(self):
        blur_profile = self.blur_profile
        # blur_profile.r[angle_bin][r_bin]
        bins = blur_profile.r
        angles = []
        avg = 0
        for angle_bin in bins:
            tot = 0
            for r_i in range(len(angle_bin)//2):
                tot += angle_bin[r_i]
            angles.append(tot)
            avg += tot
        avg /= (blur_profile.num_angle_bins)
        # Find maxima of smoothed out image
        maxima, lesser_maxima = find_relative_maxima(angles,avg, smooth_window=5)

        transition_bands = []
        for max_angle in maxima:
            # Adjust angle_index to be within the valid range
            angle_index = (max_angle[0] + blur_profile.num_angle_bins // 2) % blur_profile.num_angle_bins

            # Get the signal for the current angle
            cur_sig = blur_profile.r[angle_index]

            # Calculate the average and standard deviation of the current signal
            average = np.mean(cur_sig)
            standard_deviation = np.std(cur_sig)

            # Define a threshold for identifying the transition band
            # threshold = average + standard_deviation * 0  # Adjust multiplier as needed
            threshold = 0.3

            # Find the index where the signal falls below the threshold
            cur_max_i = len(cur_sig) - 1
            for i, value in enumerate(cur_sig):
                if value < threshold:
                    cur_max_i = i
                    break


            # Store the angle index, radius index, and value at the transition band
            transition_band = [angle_index, cur_max_i, cur_sig[cur_max_i]]
            transition_bands.append(transition_band)
        self.transition_bands = transition_bands
        self.blur_vectors = []
        for band in self.transition_bands:
            angle_index, radius_index, _ = band
            arrow_length = (radius_index / self.blur_profile.num_radius_bins)
            arrow_angle = int(180 * (angle_index / self.blur_profile.num_angle_bins) - 90)
            blur_vector = [arrow_angle, arrow_length]
            self.blur_vectors.append(blur_vector)


    def display_all(self):
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

        # Draw arrows representing transition bands
        for band in self.transition_bands:
            angle_index, radius_index, _ = band
            arrow_length = (radius_index / self.blur_profile.num_radius_bins * image_photo.width())
            arrow_angle = int(180 * (angle_index / self.blur_profile.num_angle_bins) - 90)
            end_x = image_center_x + arrow_length * cos(radians(arrow_angle))
            end_y = image_center_y - arrow_length * sin(radians(arrow_angle))
            canvas.create_line(image_center_x, image_center_y, end_x, end_y, arrow='last', fill='red', width=2)

        canvas.image = image_photo  # Keep a reference

        # Display RGB stats, average saturation, and sharpness
        stats_text = f"\nRed Brightness: {self.rgb_stats.Br}\n"
        stats_text += f"Green Brightness: {self.rgb_stats.Bg}\n"
        stats_text += f"Blue Brightness: {self.rgb_stats.Bb}\n"
        stats_text += f"Red Contrast: {self.rgb_stats.Cr}\n"
        stats_text += f"Green Contrast: {self.rgb_stats.Cg}\n"
        stats_text += f"Blue Contrast: {self.rgb_stats.Cb}\n"
        stats_text += f"Saturation: {self.average_saturation}\n"
        stats_text += f"Sharpness: {self.sharpness}\n"
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

        # Start with basic statistics and blur profile data
        report_data = {
            'Height': self.rgb_stats.height,
            'Width': self.rgb_stats.width,
            'Average Saturation': self.average_saturation,
            'Sharpness': self.sharpness,
            'Red Brightness': self.rgb_stats.Br,
            'Green Brightness': self.rgb_stats.Bg,
            'Blue Brightness': self.rgb_stats.Bb,
            'Red Contrast': self.rgb_stats.Cr,
            'Green Contrast': self.rgb_stats.Cg,
            'Blue Contrast': self.rgb_stats.Cb,
        }
        for i in range(max_vector_entries):
            if i < len(self.blur_vectors):
                angle, magnitude = self.blur_vectors[i]
            else:
                angle, magnitude = 0, 0
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

        # Convert the report data to a JSON string
        json_data = json.dumps(report_data, indent=4)
        return json_data


def find_relative_maxima(data, avg, smooth_window=1):
    # Smooth the data to reduce noise
    smoothed_data = smooth_data(data, window_size=smooth_window)
    ERROR_THRESH = 1.1
    
    # Find maxima in the smoothed data
    maxima = []
    lesser_maxima = []
    for i in range(1, len(smoothed_data) - 1):
        if smoothed_data[i] > smoothed_data[i - 1] and smoothed_data[i] > smoothed_data[i + 1]:
            if smoothed_data[i] > avg * ERROR_THRESH:
                maxima.append((i, data[i]))  # Use original data value for maxima
            else:
                lesser_maxima.append((i, data[i]))  # Use original data value for maxima
    return maxima, lesser_maxima


def smooth_data(data, window_size=3):
    # Simple moving average for smoothing
    return np.convolve(data, np.ones(window_size) / window_size, mode='same')


def get_report(pil_image: Image,
               h_partitions=18, s_partitions=2, v_partitions=3,
               black_thresh=0.1, gray_thresh=0.1,
               coverage_thresh=0.90, linked_list_size=1000, downsample_rate=2,
               radius_partitions=4, angle_partitions=18,
               quantity_weight=0.1, saturation_value_weight=0.9):
    # Convert PIL image to Image_RGB
    width = pil_image.width
    height = pil_image.height
    image_rgb = pil_image_to_image_rgb(pil_image)

    image_ctype = ctypes.byref(image_rgb)
    start_time = time.time()  # Start timing

    report_data_ptr = lib.get_full_report_data(image_ctype, 
                                               h_partitions, s_partitions, v_partitions,
                                               black_thresh, gray_thresh,
                                               coverage_thresh, linked_list_size, downsample_rate,
                                               radius_partitions, angle_partitions,
                                               quantity_weight, saturation_value_weight
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
    report.data_ptr = report_data_ptr

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


image_name = "12.png"
def run_demonstration():
    global image_name
    image_path = "images/original/"
    # Check if a custom image path was provided as a command-line argument
    if len(sys.argv) > 1:
        image_name = sys.argv[1]

    # Open image
    image = Image.open(image_path+image_name)

    # Generate report
    report = get_report(image,
                        coverage_thresh=.96,
                        downsample_rate=3,
                        radius_partitions=40,
                        angle_partitions=72,
                        quantity_weight=0.9,
                        saturation_value_weight=0.1)
    report.image = image

    # Display images that represent image: blur profile and color palette
    report.vectorize_blur_profile()
    report.generate_color_palette_image()
    report.generate_blur_profile_image()
    # Display full user interface
    report.display_all()
    json = report.to_json()
    print(json)

    lib.free_full_report(ctypes.byref(report.data_ptr))


if __name__ == "__main__":
    run_demonstration()
