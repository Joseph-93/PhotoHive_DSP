#include "image_processing.h"
#include "utilities.h"
#include "filtering.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SATURATION 0.999999
#define MAX_VALUE 0.999999


/******************************************************************************
 * create_rgb_image creates an image of given width and height, being a pointer
 *   to an Image_RGB data-type.
 *  -Function initializes all pixels to 0.
 *  -Function returns NULL if any allocations fail
 *  -Function returns NULL if height or width are equal to 0.
******************************************************************************/
Image_RGB* create_rgb_image(int width, int height) {
    if (!height) {  // If height is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_rgb_image received height==0.\n");
        return NULL;
    } else if (!width) {   // If width is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_rgb_image received width == 0.\n");
        return NULL;
    }
    Image_RGB* image = (Image_RGB*)calloc(1, sizeof(Image_RGB));
    if (!image) {   // Throw error and return NULL if malloc failed
        fprintf(stderr, "ERROR: Memory allocation for Image_RGB failed.\n");
        return NULL;
    }
    image->width = width;
    image->height = height;
    // allocate memory for rgb channels, initialized to 0.0
    image->r = (Pixel*)calloc(height * width, sizeof(Pixel));
    image->g = (Pixel*)calloc(height * width, sizeof(Pixel));
    image->b = (Pixel*)calloc(height * width, sizeof(Pixel));
    if (!image->r || !image->g || !image->b){   // Return NULL if callocs failed
        fprintf(stderr, "ERROR: Memory allocation for RGB channel pixels failed.\n");
        return NULL;
    }
    return image;
}


/******************************************************************************
 * create_hsv_image creates an image of given width and height, being a pointer
 *   to an Image_HSV data-type.
 *  -Function initializes all pixels to 0.
 *  -Function returns NULL if any allocations fail
 *  -Function returns NULL if height or width are equal to 0.
******************************************************************************/
Image_HSV* create_hsv_image(int width, int height) {
    if (!height) {  // If height is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_hsv_image received height==0.\n");
        return NULL;
    } else if (!width) {   // If width is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_hsv_image received width == 0.\n");
        return NULL;
    }
    Image_HSV* image = (Image_HSV*)calloc(1, sizeof(Image_HSV));
    if (!image) {   // Throw error and return NULL if malloc failed
        fprintf(stderr, "ERROR: Memory allocation for Image_RGB failed.\n");
        return NULL;
    }
    image->width = width;
    image->height = height;
    // allocate memory for rgb channels, initialized to 0.0
    image->pixels = (Pixel_HSV*)calloc(height * width, sizeof(Pixel_HSV));
    if (!image->pixels){   // Return NULL if callocs failed
        fprintf(stderr, "ERROR: Memory allocation for HSV pixel array failed.\n");
        return NULL;
    }
    return image;
}


/******************************************************************************
 * create_pgm_image creates an image of given width and height, being a pointer
 *   to an Image_PGM data-type.
 *  -Function initializes all pixels to 0.
 *  -Function returns NULL if any allocations fail
 *  -Function returns NULL if height or width are equal to 0.
******************************************************************************/
Image_PGM* create_pgm_image(int width, int height) {
    if (!height) {  // If height is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_pgm_image received height==0.\n");
        return NULL;
    } else if (!width) {    // If width is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_pgm_image received width == 0.\n");
        return NULL;
    }
    Image_PGM* image = (Image_PGM*)malloc(sizeof(Image_PGM));
    if (!image) {   // If malloc failed, throw error and return NULL
        fprintf(stderr, "ERROR: Memory allocation for Image_Gray failed");
        return NULL;
    }
    image->width = width;
    image->height = height;
    // allocate memory for image pixels, initialized to 0.0
    image->data = (Pixel*)calloc(height * width, sizeof(Pixel));
    if (!image->data){  // Throw error and return NULL if calloc failed
        fprintf(stderr, "ERROR: Memory allocation for image->pixels failed.\n");
        return NULL;
    }
    return image;
}


/******************************************************************************
 * read_image reads in a .txt file and returns an Image_RGB pointer.
 *  -All values stored to the Image_RGB file are doubles normalized to [0,1]
 *  -read_image assumes the image is in the format created by the supplemental
 *   python program imageDecompressor.py. If the file is incorrectly formatted,
 *   there is currently no check for it, so the program may experience failure,
 *   or may have undefined behavior.
 *  -File format from .txt file assumes that all values are between 0 and 255.
 *   If any other values are found, the program safely returns NULL.
 *  -If file fails to open, function returns NULL
 *  -If width or height read from file are invalid, function returns NULL.
******************************************************************************/
Image_RGB* read_image(const char* filepath) {
    FILE* file = fopen(filepath, "r");
    if (file == NULL) {
        fprintf(stderr, "Error occured opening file");
        return NULL;
    }
    int width, height;
    fscanf(file, "%d %d", &width, &height);
    if (width < 1) {
        fprintf(stderr, "ERROR: Width read from file was <1, which is not valid");
        return NULL;
    }
    else if (height < 1) {
        fprintf(stderr, "ERROR: Height read from file was <1, which is not valid");
        return NULL;
    }

    // Allocate memory for the image
    Image_RGB* image = (Image_RGB*)malloc(sizeof(Image_RGB));
    if (!image) {
        fprintf(stderr, "ERROR: Memory allocation for image failed.\n");
        return NULL;
    }
    image->height = height;
    image->width = width;
    image->r = (Pixel*)calloc(height * width, sizeof(Pixel));
    image->g = (Pixel*)calloc(height * width, sizeof(Pixel));
    image->b = (Pixel*)calloc(height * width, sizeof(Pixel));
    if (!image->r || !image->g || !image->b) {
        fprintf(stderr, "ERROR: Memory allocation for image color channels failed.\n");
        return NULL;
    }

    // Read in the RGB values for the image
    int i_tot = height * width;
    int r, g, b;
    for (int i=0; i<i_tot; i++) {
        // Scan in next r,g,b values
        fscanf(file, "%d %d %d", &r, &g, &b);
        // Check for incorrect values
        if (r<0 || r>255 || g<0 || g>255 || b<0 || b>255) {
            fprintf(stderr, "ERROR: file contained pixel values outside range [0,255].");
            return NULL;
        }
        // Normalize r,g,b values to the range [0,1]
        image->r[i] = (double) r / 255.0f;
        image->g[i] = (double) g / 255.0f;
        image->b[i] = (double) b / 255.0f;
    }
    fclose(file);
    return image;
}


/******************************************************************************
 * write_image_to_file converts an Image_RGB into a txt file at the given path.
 *  -Path passed must be of the format "*.txt".
 *  -All values will be scaled to [0,255] and converted to integers.
 *  -Given image must be scaled to range [0,1].
 *  -If file already exists, it WILL be rewritten!
 *  -Image may be viewed by running the python program imageSaver.py, inputting
 *   the name of this file given. imageSaver.py will create a valid .png file.
******************************************************************************/
void write_image_to_file(Image_RGB* image, const char* path) {
    // Open file for writing
    FILE* file = fopen(path, "w");
    if (file == NULL) {fprintf(stderr, "Error opening file.\n"); return;}
    fprintf(file, "%d %d\n", image->width, image->height);

    // Begin copying file
    int r, g, b;
    int i_tot = image->height * image->width;
    for (int i=0; i<i_tot; i++) {
        r = (int)(image->r[i] * 255.0f);
        g = (int)(image->g[i] * 255.0f);
        b = (int)(image->b[i] * 255.0f);
        fprintf(file, "%d %d %d\n", r, g, b);
    }
    fclose(file);
}


/******************************************************************************
 * crop_pgm returns an Image_PGM pointer, with only the indices encased by the
 *   boundaries (right, left, bottom, top).
 *  -Right and left are self-explanatory from the perspective of looking at the
 *   image itself. Left should always be less than right.
 *  -Top and bottom are self-explanatory from the perspective of looking at the
 *   image itself. Top should always be less than bottom.
 *  -No boundary indices may be negative. If any are, NULL will be returned.
******************************************************************************/
Image_PGM* crop_pgm(Image_PGM* image, int right, int left, int bottom, int top) {
    // Check that boundaries are valid
    if (right > image->width | left > image->width | bottom > image->height | top > image->height |
        left < 0 | right < 0 | top < 0 | bottom < 0) {
        fprintf(stderr, "Error: crop boundaries outside of image boundaries.");
        return NULL;
    }
    // Calculate cropped image width and height, then create the image
    int width = right - left;
    int height = bottom - top;
    Image_PGM* cropped_image = create_pgm_image(width, height);

    // Copy the portion of the original to be kept, onto the cropped image
    for (int x=0; x<width; x++) {
        for (int y=0; y<height; y++) {
            cropped_image->data[y*width + x] = image->data[(y+top)*image->width + x + left];
        }
    }
    return cropped_image;
}


/******************************************************************************
 * crop_image returns an Image_RGB pointer, with only the indices encased by
 *   the boundaries (right, left, bottom, top).
 *  -Right and left are self-explanatory from the perspective of looking at the
 *   image itself. Left should always be less than right.
 *  -Top and bottom are self-explanatory from the perspective of looking at the
 *   image itself. Top should always be less than bottom.
 *  -No boundary indices may be negative. If any are, NULL will be returned.
******************************************************************************/
Image_RGB* crop_image(Image_RGB* image, int right, int left, int bottom, int top) {
    // Check that boundaries are valid
    if (right > image->width | left > image->width | bottom > image->height | top > image->height |
        left < 0 | right < 0 | top < 0 | bottom < 0) {
        fprintf(stderr, "Error: crop boundaries outside of image boundaries.");
        return NULL;
    }

    // Calculate cropped image width and height, then create the image
    int width = right - left;
    int height = bottom - top;
    Image_RGB* cropped_image = create_rgb_image(width, height);

    // Copy the portion of the original to be kept, onto the cropped image
    for (int x=0; x<width; x++) {
        for (int y=0; y<height; y++) {
            cropped_image->r[y*width + x] = image->r[(y+top)*image->width + x + left];
            cropped_image->g[y*width + x] = image->g[(y+top)*image->width + x + left];
            cropped_image->b[y*width + x] = image->b[(y+top)*image->width + x + left];
        }
    }
    return cropped_image;
}


/******************************************************************************
 * free_crop_boundaries frees a Crop_Boundaries object given its pointer.
******************************************************************************/
void free_crop_boundaries(Crop_Boundaries* cb) {
    free(cb->top); cb->top = NULL;
    free(cb->bottom); cb->bottom = NULL;
    free(cb->left); cb->left = NULL;
    free(cb->right); cb->right = NULL;
    free(cb); cb = NULL;
}


/******************************************************************************
 * free_image_rgb frees all memory taken by an Image_RGB* and the image itself.
 *  -In the end it nullifies image as well, so there is no need for image=NULL.
******************************************************************************/
void free_image_rgb(Image_RGB* image) {
    free(image->r); image->r = NULL;
    free(image->g); image->g = NULL;
    free(image->b); image->b = NULL;
    free(image); image = NULL;
}


/******************************************************************************
 * free_image_hsv frees all memory taken by an Image_HSV* and the image itself.
 *  -In the end it nullifies image as well, so there is no need for hsv=NULL.
******************************************************************************/
void free_image_hsv(Image_HSV* hsv) {
    free(hsv->pixels); hsv->pixels = NULL;
    free(hsv); hsv = NULL;
}


/******************************************************************************
 * save_rgb writes an rgb image to user's desired image file (in .txt format 
 *   for python to save via imageSaver.py).
 *  -Function is designed for development, as it prompts users for command-line
 *   input.
******************************************************************************/
void save_rgb(Image_RGB* rgb) {
    const char* q = "\nEnter an RGB image filename: ";
    char* output_filename = create_path("images/output/", q, ".txt");
    write_image_to_file(rgb, output_filename);
    free(output_filename), output_filename = NULL;
}


/******************************************************************************
 * read_image_from_files reads in user's desired image from the list of
 *   readable .txt files
 *  -Function is designed for development, as it prompts users for command-line
 *   input.
******************************************************************************/
Image_RGB* read_image_from_files() {
    const char* q = "\nEnter an RGB image filename: ";
    char* input_filename = create_path("images/readable/", q, ".txt");
    Image_RGB* image = read_image(input_filename);
    free(input_filename), input_filename = NULL;
    return image;
}


/******************************************************************************
 * downsample_rgb returns a downsampled image from the image given.
 *  -Function uses the simple subsampling method
 *  -Function does NOT anti-alias the image prior to downsampling. If the image
 *   is meant to be displayed or operated on in the spacial domain or frequency
 *   domain, user must low-pass filter the image prior to downsampling.
 *  -Image_RGB* image is the rgb image to be downsampled
 *  -N is the interval at which pixels are sampled, in 2D. One pixel after
 *   downsampling attempts to describe an NxN square of pixels.
******************************************************************************/
Image_RGB* downsample_rgb(Image_RGB* image, short N) {
    // Create image for decimation
    int height = image->height; int width = image->width;
    int new_height = image->height/N; int new_width = image->width/N;
    Image_RGB* new_image = create_rgb_image(new_width, new_height);

    // Add every NxN pixel from image to new image
    int y_old_increment = (N - 1)*width - new_width*N;
    int i_old = 0;
    int i_new = 0;
    for(int y=0; y<new_height; y++) {
        for (int x=0; x<new_width; x++) {
            new_image->r[i_new] = image->r[i_old];
            new_image->g[i_new] = image->g[i_old];
            new_image->b[i_new] = image->b[i_old];
            i_old += N;
            i_new++;
        }
        i_old += y_old_increment;
    }

    return new_image;
}


/******************************************************************************
 * rgb2hsv creates a new hsv image from the pixels of an rgb image.
******************************************************************************/
Image_HSV* rgb2hsv(Image_RGB* rgb) {
    if (!rgb) {
        fprintf(stderr, "ERROR: rgb2hsv() recieved a null pointer to rgb.");
        return NULL;
    }
    // Set height and width to be used
    short height = rgb->height;
    short width = rgb->width;
    // Initialize 3-channel image as hsv, with all 0s
    Image_HSV* hsv = (Image_HSV*)create_hsv_image(width, height);
    // Iterate through every pixel
    int i_tot = height*width;
    for (int i=0; i<i_tot; i++) {
        // Set necessary values
        Pixel_HSV* cur = &hsv->pixels[i];
        Pixel r = rgb->r[i], g = rgb->g[i], b = rgb->b[i];
        Pixel max = fmax(fmax(r, g), b);
        Pixel min = fmin(fmin(r, g), b);
        Pixel delta = max - min;

        // Set H value
        Pixel h;
        if (delta == 0) {h = 0;}    // Textbook equation for rgb->h conversion
        else if (max == r) {h = 60 * ((g-b)/delta);}
        else if (max == g) {h = 60 * (2+(b-r)/delta);}
        else {h = 60 * (4+(r-g)/delta);}
        if (0<h && h<360) {;}   // if within range, continue forward
        else if (h<0) {     // if h is too low
            while (h<0) {h += 360;} // bring it up until within range
        }
        else if (h>360) {   // if h is too high
            while (h > 360) {h -= 360;} // bring it down until within range
        }
        cur->h = h;
        
        // Set V value
        if (max == 1) {cur->v = MAX_VALUE;}
        else {cur->v = max;}

        // Set S value
        if (max==0) {cur->s = 0;}
        else if (delta == max) {cur->s = MAX_SATURATION;}
        else {cur->s = delta/max;}
    }
    return hsv;
}


/******************************************************************************
 * hsv2rgb creates a new rgb image from the pixels of an hsv image.
******************************************************************************/
Image_RGB* hsv2rgb(Image_HSV* hsv) {
    if (!hsv) {
        fprintf(stderr, "ERROR: hsv2rgb() recieved a null pointer to rgb.");
        return NULL;
    }
    short height = hsv->height;
    short width = hsv->width;

    // Initialize 3-channel image as rgb, assuming create_rgb_image is similar to create_hsv_image
    Image_RGB* rgb = (Image_RGB*)create_rgb_image(width, height);
    if (!rgb) return NULL; // Check if memory allocation failed

    int i_tot = height * width;
    for (int i = 0; i < i_tot; i++) {
        Pixel_HSV cur = hsv->pixels[i];
        Pixel h = cur.h;
        Pixel s = cur.s;
        Pixel v = cur.v;

        Pixel c = v * s; // Chroma
        Pixel x = c * (1 - fabsf(fmodf(h / 60.0, 2) - 1));
        Pixel m = v - c;
        
        Pixel rs, gs, bs;
        if (h >= 0 && h < 60) {
            rs = c; gs = x; bs = 0;
        } else if (h >= 60 && h < 120) {
            rs = x; gs = c; bs = 0;
        } else if (h >= 120 && h < 180) {
            rs = 0; gs = c; bs = x;
        } else if (h >= 180 && h < 240) {
            rs = 0; gs = x; bs = c;
        } else if (h >= 240 && h < 300) {
            rs = x; gs = 0; bs = c;
        } else {
            rs = c; gs = 0; bs = x;
        }
        
        // Assign computed values to the RGB image
        rgb->r[i] = rs + m;
        rgb->g[i] = gs + m;
        rgb->b[i] = bs + m;
    }

    return rgb;
}


/******************************************************************************
 * Validates the given x and y coordinates against the dimensions of the given 
 *   image array.
 * -If the coordinates are out of bounds, an error message is printed, and the 
 *   error count is incremented.
 * 
 * @param array The image array against which the coordinates are validated.
 * @param x The x-coordinate to validate.
 * @param y The y-coordinate to validate.
 * @param errors Pointer to an integer that tracks the total number of errors.
******************************************************************************/
void validate_coordinates(Image_PGM* array, int x, int y, int* errors) {
    if (x>=array->width) {
        *errors +=1;
        fprintf(stderr, "ERROR %3d: x is above the array width: x=%3d, width=%3d.\n", *errors, x, array->width);
    }
    else if (x<0) {
        *errors +=1;
        fprintf(stderr, "ERROR %3d: x is below 0. x=%3d.\n", *errors, x);
    }
    if (y>=array->height) {
        *errors +=1;
        fprintf(stderr, "ERROR %3d: y is above the array height: y=%3d, height=%3d.\n", *errors, y, array->height);
    }
    else if (y<0) {
        *errors +=1;
        fprintf(stderr, "ERROR %3d: y is below zero. y=%3d.\n", *errors, y);
    }
}


/******************************************************************************
 * rgb2pgm takes in an rgb image and returns a grayscale image object.
******************************************************************************/
Image_PGM* rgb2pgm(Image_RGB* rgb) {
    Image_PGM* pgm = create_pgm_image(rgb->width, rgb->height);
    int length = rgb->height * rgb->width;
    for (int i=0; i<length; i++) {
        pgm->data[i] = 0.299*rgb->r[i] + 0.587*rgb->g[i] + 0.114*rgb->b[i];
    }
    return pgm;
}


Image_RGB* pgm2rgb(Image_PGM* pgm) {
    Image_RGB* rgb = create_rgb_image(pgm->width, pgm->height);
    int length = pgm->width*pgm->height;
    for (int i=0; i<length; i++) {
        rgb->r[i] = pgm->data[i];
        rgb->g[i] = pgm->data[i];
        rgb->b[i] = pgm->data[i];
    }
    return rgb;
}


void free_image_pgm(Image_PGM* image) {
    free(image->data); image->data = NULL;
    free(image); image = NULL;
}


Pixel get_hsv_average(Image_HSV* hsv) {
    double accumulator = 0.0;
    int length = hsv->height*hsv->width;
    for (int i=0; i<length; i++) {
        accumulator += (double)hsv->pixels[i].s;
    }
    return (Pixel)(accumulator/(double)length);
}


RGB_Statistics* get_rgb_statistics(Image_RGB* image) {
    RGB_Statistics* rgb_stats = (RGB_Statistics*)calloc(1, sizeof(RGB_Statistics));
    int length = image->height*image->width;
    rgb_stats->Br = get_average(image->r, length);
    rgb_stats->Bg = get_average(image->g, length);
    rgb_stats->Bb = get_average(image->b, length);
    rgb_stats->Cr = sqrt(get_variance(image->r, length, rgb_stats->Br));
    rgb_stats->Cg = sqrt(get_variance(image->g, length, rgb_stats->Bg));
    rgb_stats->Cb = sqrt(get_variance(image->b, length, rgb_stats->Bb));
    return rgb_stats;
}
