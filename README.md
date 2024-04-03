# PhotoHive_DSP README.md

## Introduction

This guide will help you download and integrate your project with the PhotoHive_DSP library, a digital signal processing tool for pre-Computer Vision feature extraction for images of mountain bike and running races. This tool is for Machine-Learning and Computer Vision engineers with experience in Python and the Pillow library.

## Safety First

(This may seem silly… but code security…) This library is not yet considered stable. Memory issues are possible, and system crashes, null pointer return values, and the like are completely possible. Library failure could cause vulnerabilities and crashes in the code that calls it. Exercise caution when using this library on production environments and make sure to handle errors well when calling the library.

## What does PhotoHive_DSP measure?

This library is built to measure the following variables quickly, precisely, and accurately:
- Average saturation
- Sharpness of salient characters
- Brightness of the red, green, and blue channels
- Color palette
- Contrast of the red, green, and blue channels
- Angle and Magnitude of blurring

## Downloading and Integrating the Project

This section provides instructions for downloading and setting up the library.

### A. Download and build the library

1. Go to the PhotoHive_DSP git repository: https://github.com/Joseph-93/PhotoHive_DSP.
2. Download or clone the full repository within your root folder, enter the repository, and run `cmake` on the library as follows (tutorial uses Linux, but Windows and MacOS versions should work the same):
    ```
    git clone https://github.com/Joseph-93/PhotoHive_DSP.git
    cd PhotoHive_DSP
    rm -rf build && mkdir build
    cmake -S src -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON -DDEBUGGING=ON
    cmake --build build --config Debug
    ```

### B. Add test images and run test suite

1. From root, run the following commands to add the proper image data to your files:
    ```
    mkdir data
    cd data
    mkdir original
    mkdir test_readable
    ```
2. Add any wanted test images to `/data/original/`.

3. Run the following program to convert any wanted images to .txt files that the C test suite can read directly.
    ```
    python3 PhotoHive_DSP/src/test/image_utils.py png2txt
    ```
4. Run the test suite to ensure that all works properly:
    ```
    ./PhotoHive_DSP/build/test_suite
    ```
    The output should resemble that given below:
    ```
    hello there!
    Enter an RGB image filename:
    30
    data/test_readable/30.txt
    There are 12 cores available to the C program.
    downsample rgb took 0.000000 seconds to execute
    rgb2hsv took 0.277513 seconds to execute
    rgb2pgm took 0.068783 seconds to execute
    rgb statistics took 0.246769 seconds to execute
    hsv average took 0.049124 seconds to execute
    color palette took 0.515915 seconds to execute
    sharpness took 0.000003 seconds to execute
    finding the FFT max value took 0.009982 seconds to execute normalizing FFT values took 0.045176 seconds to execute
    computing the magnitude fft took 0.184537 seconds to execute cartesian to polar conversion took 0.063377 seconds to execute accumulating pixels to blur_profile bins took 0.324324 seconds to execute calculating average of blur_profile bins took 0.000011 seconds to execute calculate_blur_profile took 0.324390 seconds to execute
    freeing FFT structures took 0.000001 seconds to execute
    blur profile took 0.572365 seconds to execute
    compile full report took 0.000000 seconds to execute
    free images took 0.001200 seconds to execute
    ```

### C. Build release version and run install

1. From the library root, remove old builds and create a new build:
    ```
    cd PhotoHive_DSP
    rm -rf build && mkdir build
    cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release
    ```
2. Run the install command on the library:
    ```
    pip install PhotoHive_DSP/
    ```
3. Import necessary library components to your program, including Pillow:
    ```python
    from PIL import Image
    from PhotoHive_DSP import Report, get_report, set_bounding_boxes
    ```

### D. Your program should now be ready to run!

## Accessing Features, and Viewing Report Data

### Code Interface Essentials

#### Using `get_report()`

`get_report()` is the basis of the entire PhotoHive_DSP library. With it, a user may request all data that the library can offer. `get_report()` cleans all data and abstracts the use of `ctypes` out. This ensures that all data can be easily recovered by the machine learning expert.

**Inputs:**
- `pil_image` (required): a Python Pillow.Image object of an RGB image.
- `salient_characters` (default None): a library-specific data-type, `Crop_Boundaries`. NoneType is also allowed, resulting in the library assuming that there are no crop boundaries to be used.
    - A `Crop_Boundaries` object can be obtained using the `set_bounding_boxes()` function, discussed below.
- `h_partitions` (Default: 18): Number of partitions for the hue component in the HSV color space. Increasing this number can result in more granular color detection but may also increase computational load.
- `s_partitions` (Default: 2): Number of partitions for the saturation component in the HSV color space. Higher values allow for more detailed saturation analysis.
- `v_partitions` (Default: 3): Number of partitions for the value (brightness) component in the HSV color space. Adjusting this affects the granularity of brightness analysis.
- `black_thresh` (Default: 0.1): Threshold for determining black color presence. Pixels with value below this threshold are considered black.
- `gray_thresh` (Default: 0.1): Threshold for identifying gray colors. This affects how the algorithm distinguishes between colored and grayscale areas.
- `coverage_thresh` (Default: 0.95): The minimum percentage of the image's pixels that parent colors should cover. Used in color quantization to decide major colors in the image.
- `linked_list_size` (Default: 1000): Size of the linked list used internally for data storage. This can impact the memory usage and speed of the operation.
- `downsample_rate` (Default: 1): Rate at which the image is downsampled. A higher rate decreases resolution and may speed up processing, but can reduce accuracy.
- `radius_partitions` (Default: 40): Number of partitions for radius in the blur profile analysis. This defines the granularity of radial blur detection.
- `angle_partitions` (Default: 72): Number of angular partitions in the blur profile analysis. More partitions mean finer angular resolution for blur direction detection.
- `quantity_weight` (Default: 0.1): Weight of the quantity of pixels in calculating the saliency of colors. Affects how much the number of pixels influences color importance.
- `saturation_value_weight` (Default: 0.9): Weight in the calculation of color saliency, emphasizing the product of saturation and value (brightness).
- `fft_streak_thresh` (Default: 1.20): Threshold used in the FFT process to detect streaks indicating motion blur. Higher values require stronger directional frequency components to consider as streaks.
- `magnitude_thresh` (Default: 0.3): Threshold for the magnitude of frequencies considered significant in blur analysis. This helps in identifying the extent of blurring.
- `blur_cutoff_ratio_denom` (Default: 2): Denominator of the ratio used to determine the cutoff point in blur profile analysis, affecting how quickly the algorithm identifies the transition from sharp to blurred regions.

**Outputs:**
- `rgb_stats`: This is a Python object directly mirroring the C structure, containing fields for brightness and contrast for each RGB channel (Br, Bg, Bb, Cr, Cg, Cb). The image dimensions are also included as height and width. Access these fields directly, for example, `report.rgb_stats.Br` for red brightness.
- `color_palette`: This is a Python list where each item is a tuple representing a color and its quantity in the image. Each tuple contains HSV values and the corresponding percentage, accessible via indexing. For example, `report.color_palette[0]` gives the first color in the palette as (H, S, V, Percentage).
- `blur_profile`: This contains a 2D list representing the blur intensity across different angles and radii. Access this data via two-dimensional indexing, like `report.blur_profile.bins[angle_index][radius_index]`. Each element in this 2D array represents the blur intensity at that specific angle and radius.
- `blur_vectors`: This is a list of objects, each with angle and magnitude attributes representing the directional blur vectors in the image. Access these vectors and their attributes like `report.blur_vectors[0].angle` and `report.blur_vectors[0].magnitude`.
- `average_saturation`: A simple float value representing the average saturation of the image. Access it directly, e.g., `report.average_saturation`.
- `sharpnesses`: A list of floats, each representing a sharpness value derived from different parts or features of the image. Access each sharpness value like `report.sharpnesses[0]`.

### Using `set_bounding_boxes()`

The `set_bounding_boxes()` function enables the definition of regions of interest within the image, allowing for targeted analysis of specified areas. Each bounding box gets its separate sharpness measured, as a way of estimating the camera’s accurate focus on the focus of a photo. The processor only allows for 10 bounding boxes per photo.

**Inputs:**
- `bounding_boxes`: This parameter should be a list where each item represents a bounding box. A bounding box can be defined using a tuple or dictionary with the keys 'top', 'bottom', 'left', and 'right'. These represent the coordinates of the bounding box in the image:
    - 'top': The top edge coordinate of the bounding box.
    - 'bottom': The bottom edge coordinate of the bounding box.
    - 'left': The left edge coordinate of the bounding box.
    - 'right': The right edge coordinate of the bounding box.
    - Each coordinate should be an integer representing the pixel position in the image.

**Outputs:**
- Returns a `Crop_Boundaries` object, which is an internal library-specific data-type used to represent the collection of bounding boxes for subsequent processing in the `get_report()` function. This object is not meant to be manipulated directly by the user but passed as an argument to `get_report()` where needed. The `Crop_Boundaries` object can also be “attached” to a report object, so that they are displayed in the output window when the `Report.display_all()` function is called.

**Example Usage:**
To use the `set_bounding_boxes()` function, you would typically define a list of bounding boxes, each specifying the area to be considered as salient or of interest. For example:

```python
bounding_boxes = [
    {'top': 50, 'bottom': 150, 'left': 30, 'right': 200},
    {'top': 200, 'bottom': 300, 'left': 50, 'right': 250}
]
crop_boundaries = set_bounding_boxes(bounding_boxes)
```

### Using the `Report()` Class

The `Report` class in the PhotoHive_DSP library serves as a container for the comprehensive data obtained from image analysis. This class provides several methods to interact with and visualize the data.

**Methods:**

- `generate_color_palette_image()`: Generates an image representation (a `Pillow.Image` object) of the color palette detected in the source image. The generated image can be viewed or saved using `display_color_palette_image()` or accessed directly via the `color_palette_image` attribute of the Report object.
- `display_color_palette_image()`: Displays the color palette image in a new window. The image should be generated beforehand using `generate_color_palette_image()`.
- `generate_blur_profile_image()`: Creates a visual image representation (a `Pillow.Image` object) of the blur profile. The resulting image can be displayed using `display_blur_profile()`.
- `display_blur_profile(path)`: Displays the blur profile image in a new window and saves the image to the specified path. The blur profile image is generated by calling `generate_blur_profile_image()` prior to this method.
- `display_all()`: Launches a comprehensive display window showing the analyzed image, color palette, blur vectors, and statistical data. Before calling this method, ensure that the following has been done:
    - `generate_color_palette_image()` is executed.
    - The `self.image` attribute is set to the source image.
    - The `self.bounding_boxes` attribute is set to the generated `Crop_Boundaries`, if any.
- `to_json()`: Converts the report data into a JSON format, suitable for use in data analysis or machine learning applications. This tool is essentially “a must” for ML engineers who wish to integrate with any machine learning libraries, such as Pandas, Pytorch, or TensorFlow.

### Viewing and Interpreting Report Data

#### Viewing Blur Profile

The blur profile is relatively simple. It is simply an FFT that has been sent through polar-coordinate simple-averaging downsampler. The results are only a 2D array of `[angle][radius]`, but the image rotates that array around a circle, looking like an un-shifted FFT. It can be deciphered as a polar form downsampled fft.

The blur profile viewer is very useful for visually calibrating the following hyperparameters:

- `radius_partitions` and `angle_partitions`
- `fft_streak_thresh`
- `magnitude_thresh`
- `blur_cutoff_ratio_denom`

#### Viewing Color Palette

The color palette is quite self-explanatory. It shows the colors extracted from the image, in order from most salient in the top-left corner to the least salient in the bottom-right corner of the image. Each color has a percentage listed, for the proportion of the image that the color accounts for.

The color palette viewer is very useful for visually calibrating the following hyperparameters:

- `h_partitions`, `s_partitions`, and `v_partitions`
- `black_thresh` and `gray_thresh`
- `coverage_thresh`
- `downsample_rate`
- `quantity_weight` and `saturation_value_weight`

#### Viewing `display_all()` Window

The `display_all` window contains all data that has been prepared to be sent to the machine learning model.

- On the left, you’ll see your original image.
- If you passed bounding boxes to the algorithm and properly prepared the `Report()` object, you should see those bounding boxes, with the sharpness measured within each bounding box. This can help you ensure that you are extracting and setting up the bounding boxes properly.
- If the algorithm discovered any motion blur in the image, you would see a red arrow at the center of the image, pointing in the direction of blur. It does NOT measure blur to the left, as it is the same as a blur to the right and gets measured as such.
- The length of the arrow indicates the strength of the blur effect. However, take caution in your visual analysis, as the length of the arrow and the amount of blur in the image is an inversely related ration; a shorter arrow means more blur in the arrow direction, and a longer arrow means less blur in the arrow direction.
- On the top right, you’ll see the rgb_statistics of your photo, as discussed in the "What Does PhotoHive Measure?" Section.
    - Red, Blue, and Green Brightness are equivalent to `Br`, `Bg`, `Bb`.
    - Red, Blue, and Green Contrast are equivalent to `Cr`, `Cg`, `Cb`.
    - Saturation is equivalent to S-bar, or the mean of saturation values across all pixels in the image.
- On the right side, you’ll see the color palette image, exactly as is seen in the color palette viewer. No further explanation necessary.
