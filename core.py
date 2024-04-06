import ctypes
import io
from types import SimpleNamespace
from math import cos, sin, radians
from ctypes import POINTER
import json
import time
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import tkinter as tk
import numpy as np
from PIL import ImageTk
from .structures import (
    Image_RGB, Crop_Boundaries, Full_Report_Data, RGB_Statistics, Color_Palette,
    Blur_Profile, Blur_Vector_Group, Sharpnesses, Pixel_HSV
)
from .utils import pil_image_to_image_rgb, image_rgb_to_pillow, hsv_to_rgb, image_pgm_to_pillow

# Assuming 'lib' is your ctypes.CDLL loaded library
from .lib import lib


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
        self.bp_ptr = c_blur_profile
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


    def generate_blur_direction_frequency_response(self):
        blur_vectors = self.blur_vectors
        blur_profile = self.blur_profile
        plt.figure(figsize=(10, 6))

        # Create a normalized radius index scale from 0 to 1
        normalized_radius_indices = np.linspace(0, 1, len(blur_profile.bins[0]))

        # Initialize a list to store all frequency responses for averaging
        all_frequency_responses = []

        for bv in blur_vectors:
            if bv.magnitude == 0.0:
                continue
            # Quantize the angle to match the blur_profile bins
            angle_partitions = len(blur_profile.bins)
            q_ang = int(bv.angle / (361 / angle_partitions) + angle_partitions / 2) % angle_partitions

            # Fetch the frequency response for the quantized angle
            frequency_response = blur_profile.bins[q_ang]
            all_frequency_responses.append(frequency_response)

            # Plot the frequency response for the current angle
            plt.plot(normalized_radius_indices, frequency_response, label=f'Directional Angle: {bv.angle} degrees')

            # Compute and plot the frequency response at 90 degrees from the current angle
            if bv.angle > 0.0:
                perpendicular_angle = (bv.angle - 90)
            else:
                perpendicular_angle = (bv.angle + 90)
            q_perp_ang = int(perpendicular_angle / (361 / angle_partitions) + angle_partitions / 2) % angle_partitions
            perp_frequency_response = blur_profile.bins[q_perp_ang]
            plt.plot(normalized_radius_indices, perp_frequency_response, label=f'Streak at {perpendicular_angle} degrees')

        plt.axhline(y=self.magnitude_threshold, color='r', linestyle='-', label='Blur magnitude threshold')

        # Compute the average for the first half of the radii across all angle bins
        half_radii = len(blur_profile.bins[0]) // self.blur_cutoff_ratio_denom
        fft_thresh_line = np.mean([row[:half_radii] for row in blur_profile.bins]) * self.fft_streak_threshold
        plt.axhline(y=fft_thresh_line, color='b', linestyle='-', label='FFT Streak threshold')

        # Calculate the average frequency response across all angles
        average_frequency_response = np.mean(blur_profile.bins, axis=0)
        plt.plot(normalized_radius_indices, average_frequency_response, label='Average Response', linewidth=2, linestyle='--')

        plt.title('Frequency Response by Angle')
        plt.xlabel('Radius Index')
        plt.ylabel('Magnitude')
        plt.legend()
        plt.grid(True)

        # Save the plot to a BytesIO buffer and then to a PIL.Image
        buffer = io.BytesIO()
        plt.savefig(buffer, format='png')
        buffer.seek(0)
        image = Image.open(buffer)
        self.blur_vector_plot = image
        plt.close()
    

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
        font = ImageFont.truetype("DejaVuSans.ttf", 12)  # You can adjust the font size as needed

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
    

    def generate_blur_profile_image(self):
        bp_ptr = ctypes.byref(self.bp_ptr)
        height = self.rgb_stats.height
        width = self.rgb_stats.width
        image_rgb_ptr = lib.get_blur_profile_visual(bp_ptr, height, width)

        img = image_pgm_to_pillow(image_rgb_ptr, width, height)

        # Crop the image to remove the blank right half
        self.blur_profile_image = img.crop((0, 0, width // 2, height))

    
    def display_color_palette_image(self):
        # Display image in window
        window = tk.Tk()
        tk_image = ImageTk.PhotoImage(self.color_palette_image)
        label = tk.Label(window, image=tk_image)
        label.pack()
        window.mainloop()


    def display_blur_profile(self):
        # Display the image
        window = tk.Tk()

        # Resize the image
        max_width = window.winfo_screenwidth() * 0.8
        max_height = window.winfo_screenheight() * 0.9
        if hasattr(Image, 'Resampling'):
            resampling_filter = Image.Resampling.LANCZOS
        else:
            # Fallback for older Pillow versions
            resampling_filter = Image.LANCZOS
        scale_width = max_width / self.blur_profile_image.width
        scale_height = max_height / self.blur_profile_image.height
        scale_factor = min(scale_width, scale_height)
        resized_image = self.blur_profile_image.resize(
            (int(self.blur_profile_image.width * scale_factor), int(self.blur_profile_image.height * scale_factor)),
            resampling_filter
        )
        # Display the image
        tk_image = ImageTk.PhotoImage(resized_image)
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
        max_width = window.winfo_screenwidth() * 0.6
        max_height = window.winfo_screenheight() * 0.9
        scale_width = max_width / self.image.width
        scale_height = max_height / self.image.height
        scale_factor = min(scale_width, scale_height)

        if hasattr(Image, 'Resampling'):
            resampling_filter = Image.Resampling.LANCZOS
        else:
            # Fallback for older Pillow versions
            resampling_filter = Image.LANCZOS

        resized_image = self.image.resize(
            (int(self.image.width * scale_factor), int(self.image.height * scale_factor)),
            resampling_filter
        )
        image_photo = ImageTk.PhotoImage(resized_image)

        image_center_x = image_photo.width() // 2
        image_center_y = image_photo.height() // 2

        # Create a canvas for the main image and draw on it
        canvas = tk.Canvas(window, width=image_photo.width(), height=image_photo.height())
        canvas.pack(side='left', padx=10)
        canvas.create_image(0, 0, anchor='nw', image=image_photo)
        canvas.image = image_photo

        # Draw arrows representing blur vectors
        length_scale_factor = min(image_photo.width()/2, image_photo.height()/2)
        for vector in self.blur_vectors:
            arrow_angle = vector.angle
            arrow_length = (vector.magnitude * length_scale_factor)
            end_x = image_center_x + arrow_length * cos(radians(arrow_angle))
            end_y = image_center_y - arrow_length * sin(radians(arrow_angle))
            canvas.create_line(image_center_x, image_center_y, end_x, end_y, arrow='last', fill='red', width=2)

        # Add bounding box lines and sharpness text
        if "bounding_boxes" in dir(self) and self.bounding_boxes is not None:
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

        # Setup the right side container
        right_side_frame = tk.Frame(window)
        right_side_frame.pack(side='right', fill='both')

        # Display RGB stats, average saturation, and sharpness
        stats_text = f"\nRed Brightness: {self.rgb_stats.Br}\n" + \
                    f"Green Brightness: {self.rgb_stats.Bg}\n" + \
                    f"Blue Brightness: {self.rgb_stats.Bb}\n" + \
                    f"Red Contrast: {self.rgb_stats.Cr}\n" + \
                    f"Green Contrast: {self.rgb_stats.Cg}\n" + \
                    f"Blue Contrast: {self.rgb_stats.Cb}\n" + \
                    f"Saturation: {self.average_saturation}\n"
        stats_label = tk.Label(right_side_frame, text=stats_text, justify=tk.LEFT)
        stats_label.pack(side='top', padx=10, anchor='n')

        # Define the maximum size for images
        max_image_width = window.winfo_screenwidth() * 0.3  # 30% of screen width
        max_image_height = window.winfo_screenheight() * 0.3  # 30% of screen height

        # Function to resize image within a maximum width and height
        def resize_image(image, max_width, max_height):
            orig_width, orig_height = image.size
            ratio = min(max_width / orig_width, max_height / orig_height)
            new_size = (int(orig_width * ratio), int(orig_height * ratio))
            return image.resize(new_size, resampling_filter)

        # Display color palette image
        color_palette_image_resized = resize_image(self.color_palette_image, max_image_width, max_image_height)
        color_palette_photo = ImageTk.PhotoImage(color_palette_image_resized)
        color_palette_label = tk.Label(right_side_frame, image=color_palette_photo)
        color_palette_label.pack(side='top', anchor='n')
        color_palette_label.image = color_palette_photo  # Keep a reference

        # Try to display the blur_vector_plot if it exists
        if hasattr(self, 'blur_vector_plot') and self.blur_vector_plot is not None:
            blur_vector_image_resized = resize_image(self.blur_vector_plot, max_image_width, max_image_height)
            blur_vector_photo = ImageTk.PhotoImage(blur_vector_image_resized)
            blur_vector_label = tk.Label(right_side_frame, image=blur_vector_photo)
            blur_vector_label.pack(side='top', anchor='n')
            blur_vector_label.image = blur_vector_photo  # Keep a reference

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


def get_report(pil_image: Image, salient_characters=None,
               h_partitions=18, s_partitions=2, v_partitions=3,
               black_thresh=0.1, gray_thresh=0.1,
               coverage_thresh=0.95, linked_list_size=1000, downsample_rate=1,
               radius_partitions=40, angle_partitions=72,
               quantity_weight=0.1, saturation_value_weight=0.9,
               fft_streak_thresh=1.20, magnitude_thresh=0.3, blur_cutoff_ratio_denom=2):
    if salient_characters is None:
        salient_characters = ctypes.POINTER(Crop_Boundaries)()

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
    report.magnitude_threshold = magnitude_thresh
    report.fft_streak_threshold = fft_streak_thresh
    report.blur_cutoff_ratio_denom = blur_cutoff_ratio_denom

    return report


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
