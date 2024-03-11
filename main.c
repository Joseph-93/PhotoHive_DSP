#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>

#define DEBUG
#define FFT_NORMALIZER_LOOKUP_SIZE 256
#define PI 3.14159265
#define RAD2DEG (180 / PI)

typedef double Pixel;   // Official data type for pixels.


/******************************************************************************
 * Lookup_1D is a simple lookup table structure, used to store values for any
 *   lookup table where a 1D input represented by Pixels maps to a 1D output
 *   represented by Pixels.
******************************************************************************/
typedef struct Lookup_1D {
    unsigned int length;
    Pixel* input;
    Pixel* output;
} Lookup_1D;


/******************************************************************************
 * Image_RGB holds a 2D image of RGB pixels.
 *  -Image_RGB is meant to hold r,g,b values as decimal numbers range [0, 1].
 *  -Effective as both an FFT representation and a literal image.
 *  -Cannot hold complex values.
 *  -The separation of R,G, and B channels is built for efficiency in filtering
 *   and FFT computation, understanding the slightly higher complexity in image
 *   representation.
 *  -The image is logically 2D, height by width, but is referenced as 1D, for
 *   increased efficiency in reference. Proper indexing: img.color[y*width + x]
******************************************************************************/
typedef struct Image_RGB {
    unsigned int height, width;
    Pixel* r;
    Pixel* g;
    Pixel* b;
} Image_RGB;


/******************************************************************************
 * Pixel_HSV contains a pixel for an HSV image, specifically designed for
 *   octree grouping in color quantization.
 *  -parent_id is the id of the parent octree node for which the pixel belongs.
******************************************************************************/
typedef struct Pixel_HSV {
    int parent_id;
    double h;
    double s;
    double v;
} Pixel_HSV;


/******************************************************************************
 * Image_HSV operates the same as Image_RGB, except that it is an HSV format,
 *   obviously!
 *  -The reason for separating the types is to facilitate static typing, as
 *   many functions are built to operate for HSV or RGB images exclusively.
 *  -The separation of H, S, and V channels is built for computational
 *   efficiency, understanding the slightly higher complexity in image
 *   representation.
 *  -The image is logically 2D, height by width, but is referenced as 1D, for
 *   increased efficiency in reference. Proper indexing: img.color[y*width + x]
******************************************************************************/
typedef struct Image_HSV {
    unsigned int height, width;
    Pixel_HSV* pixels;
} Image_HSV;


/******************************************************************************
 * Octree is necessary for keeping track of octree groups contained.
 *  -valid_parents: an array of valid parent_id values. valid parents are
 *   parents which are part of the color palette.
 *  -groups: pointer to an array of all octree groups
 *  -Lh: length of each normal (non gray or black) group in hue dimension.
 *  -Ls: length of each normal (non gray or black) group in saturation dimension.
 *  -Lv: length of each normal (non gray or black) group in value dimension.
 *  -num_h: number of normal (non gray or black) group partitions in hue
 *   dimension. In other words, it is the maximum h index available.
 *  -num_s: number of normal (non gray or black) group partitions in saturation
 *   dimension. In other words, it is the maximum s index available.
 *  -num_v: number of normal (non gray or black) group partitions in value
 *   dimension. In other words, it is the maximum v index available.
 *  -num_grays: the number of gray groups built into the octree. Gray groups are
 *   automatically placed after all normal groups.
 *  -total_length: the total length of the groups array, including gray and
 *   black. always set total_length = Lh*Ls*Lv + num_grays. 
 *  -black_thresh: the value which automatically puts a pixel into the black
 *   octree groups.
 *  -gray_thresh: the saturation which automatically puts a pixel into the gray
 *   octree groups, so long as said pixel's value > black_thresh.
 *  -NOTE: There is always a black node, which will always be located at the
 *   back of the array.
 *  -NOTE: Group indexing is a complex task. For reference, to reference the
 *   i-th index in normal groups, i = h_ind*(Ls*Lv) + s_ind*Lv + v_ind.
 *   To index gray groups, i=Lh*Ls*Lv+gray_i OR i=total_length-num_grays+gray_i
 *   To index the black node, i = Lh*Ls*Lv + num_grays OR i=total_length-1.
******************************************************************************/
typedef struct Octree {
    Octree_Group* groups;
    int* valid_parents;
    int len_valid_parents;
    double Lh, Ls, Lv;
    int num_h, num_s, num_v;
    int num_grays;
    int total_length;
    Pixel black_thresh, gray_thresh;
} Octree;


/******************************************************************************
 * Octree_Group is the struct that Octree.groups points to.
 *  -id is self-explanatory
 *  -quantity is the number of pixels that belong to the group.
 *  -h, s, v are designed to be the "average" color described by the nodes
 *   within the group. At first it is automatically assigned, but later must
 *   be set to the average H,S,V values of its children.
 *  -head is the pointer to the HSV_Linked_List head.
 *  -is_valid_parent should default to false, but is set to true if the group
 *   is decided to be a valid parent.
******************************************************************************/
typedef struct Octree_Group {
    int id;
    int quantity;
    double h, s, v;
    HSV_Linked_List* head;
    HSV_Linked_List* cur;
    bool is_valid_parent;
} Octree_Group;


/******************************************************************************
 * HSV_Linked_List is made to be pointed to by the Octree_Group. It allows the
 *   the octree group to contain an undefined number of pixels for relatively
 *   low computation and memory. 
 *  -When the end of the linked list is reached, the user must allocate a new
 *   HSV_Linked_List instance, point to it from this instance, and move the
 *   pointer to the next instance.
 *  -pixels is the array object
 *  -i is the current iteration of the array which has NOT yet been filled.
 *  -int num_pixels is the array size
 *  -next is a pointer to the next linked list node
******************************************************************************/
typedef struct HSV_Linked_List {
    Pixel_HSV* pixels;
    int num_pixels;
    int i;
    HSV_Linked_List* next;
} HSV_Linked_List;


/******************************************************************************
 * Image_PGM holds a 2D image of grayscale pixels.
 *  -Effective as both an FFT representation and a literal image
 *  -Cannot store complex values
 *  -The image is logically 2D, height x width, but is referenced as 1D, for
 *   increased efficiency in reference Proper indexing: img.data[y*width + x].
******************************************************************************/
typedef struct Image_PGM {
    unsigned int height, width;
    Pixel* data;
} Image_PGM;


/******************************************************************************
 * Polar_Coord contains a single set of polar coordinates
 *  -r_sq is equivalent to radius^2, for computational efficiency
 *  -phi is the angle arctan(y/x), limited to values between -90 and 90 degrees
 *   due to the reflective property of the FFT of a real-valued 2D image.
******************************************************************************/
typedef struct Polar_Coord {
    int r_sq;
    double phi;
} Polar_Coord;


/******************************************************************************
 * Cartesian_To_Polar an x,y coordinate of cartesian-polar conversions is
 *  -Essentially a replacement for a lookup table, because of the necessity for
 *   adherence to dynamic size requirements.
 *  -Table is logically 2D height by width, though it is literally 1D for
 *   efficient indexing.
******************************************************************************/
typedef struct Cartesian_To_Polar {
    unsigned int height, width;
    Polar_Coord* data;
} Cartesian_To_Polar;


// Official data type used for bin sums in Blur_Profiles.
typedef double Bin;


/******************************************************************************
 * Blur_Profile_RGB represents the strength of directional blur, per angle, by
 *   sum of all pixels within an FFT subsurface bounded by radii and angles.
 *  -num_angle_bins and num_radius_bins are parameters given by the user, to be
 *   tested by the machine learning model for optimal performance and accuracy.
 *  -angle_bin_size is found by total angular range divided by num_angle_bins.
 *  -radius_bin_size is determined by FFT width divided by num_radius_bins.
 *  -Bin* r, g, and b bins are 2D logically (angle by radius), but are
 *   literally 1D. Properly indexed as {color}_bins[angle*num_angle_bins + rad]
******************************************************************************/
typedef struct Blur_Profile_RGB {
    int num_angle_bins, num_radius_bins;
    int angle_bin_size, radius_bin_size;
    Bin** r_bins; // bins should be referenced as [angle][radius]
    Bin** g_bins; // bins should be referenced as [angle][radius]
    Bin** b_bins; // bins should be referenced as [angle][radius]
} Blur_Profile_RGB;     // Holds blur profile of bins for R,G,B


/******************************************************************************
 * GLOBAL VARIABLES
******************************************************************************/
int num_cores;


/******************************************************************************
 * cartesian_to_polar_conversion creates a mapping for cartesian to polar conversion.
 *  -Built to dynamically compute for a given FFT height and width
 *  -Returns Cartesian_To_Polar object pointer
 *  -CAUTION: Assumes that FFT image is NOT shifted, such that the DC component
 *   is found in the top and bottom-left corners of the image, rather than in
 *   the center-left side.
******************************************************************************/
Cartesian_To_Polar* cartesian_to_polar_conversion(unsigned int width, unsigned int height) {
    Cartesian_To_Polar* conversion = (Cartesian_To_Polar*)malloc(sizeof(Cartesian_To_Polar));
    if (!conversion) {
        fprintf(stderr, "ERROR: Memory allocation for Cartesion_To_Polar conversion table failed.\n");
        return NULL;
    }
    conversion->height = height, conversion->width = width;
    conversion->data = (Polar_Coord*)malloc(width * height * sizeof(Polar_Coord));
    if (!conversion->data) {
        fprintf(stderr, "ERROR: Memory allocation for data in conversion table failed.\n");
        return NULL;
    }
    int half_height = height/2;
    int y_times_width = 0;
    for (int y=0; y<height/2; y++) {    // height/2 for symmetry across x axis
        for (int x=0; x<width; x++) {
            int r_sq = x*x+y*y;
            double phi = atan2(y, x);
            // Not fft shifted, so for the top half, phi=-atan2(y, x)
            conversion->data[y_times_width + x].phi = -phi;
            conversion->data[y_times_width + x].r_sq = r_sq;
            // Not fft shifted, so for the bottom half phi=atan2(y, x)
            conversion->data[(height-1-y)*width + x].phi = phi;
            conversion->data[(height-1-y)*width + x].r_sq = r_sq;
        }
        y_times_width += width;
    }
    return conversion;
}


/******************************************************************************
 * newton_int_sqrt returns the nearest integer approximation to the square root
 *   of a given double.
 *  -The Newton square-root algorithm is used, but manually set for accuracy
 *   only good enough to reach the nearest integer and then stop and return the
 *   truncated integer value
******************************************************************************/
int newton_int_sqrt(double val) {
    if (val==0) {return 0;}
    double x = val;
    double sqrt = x;
    while(1) {
        sqrt = 0.5 * (x + (val/x));
        if (fabs(sqrt-x)<1) {return (int)sqrt;}
        x = sqrt;
    }
}


/******************************************************************************
 * calculate_blur_profile returns 2D array blur profile for an RGB image by
 *   mapping the image FFT data into polar coordinates and averaging the
 *   blurs in specified angular and radial bins.
 *  -Requires that the fft given is not fft-shifted and has no negative
 *   components, such that the DC component can be found in the top and bottom
 *   left corners.
 *  -Requires a Cartesian_To_Polar* pre-calculated. This can be done using the
 *   cartesian_to_polar_conversion() function. It is not enveloped into this
 *   function for the sake of efficiency if calculating several different FFTs
 *   of equal size.
 *  -Bin spacing and size is uniform in both radius and angle, determined by
 *   num_radius_bins and num_angle_bins.
 *  -The profile struct and bins themselves have been left as a 2D architecture
 *   because they will likely be small enough that the entire structure will be
 *   cached and quickly retrieved at compile time. It should not greatly affect
 *   performance.
 *  -Conversion and fft must have the same height and width. If they do not,
 *   the function will return NULL.
 *  -If the memory allocation for Blur_Profile_RGB or profile bins fail, the 
 *   function will return NULL.
******************************************************************************/
Blur_Profile_RGB* calculate_blur_profile(
                    const Cartesian_To_Polar* conversion, 
                    const Image_RGB* fft, 
                    int num_radius_bins, 
                    int num_angle_bins) {
    
    // Issue necessary warnings
    if (conversion->height != fft->height || conversion->width != fft->width) {
        printf("ERROR: conversion and fft dimensions do not match.\n");
        printf("\tconversion height: %4d\t conversion width: %4d\n", conversion->height, conversion->width);
        printf("\tfft height: %4d\t fft width: %4d\n", fft->height, fft->width);
        printf("Conversion and fft dimensions should always match for proper blur profile calculation\n\n");
        return NULL;
    }

    // Basic blur profile setup
    Blur_Profile_RGB* profile = (Blur_Profile_RGB*)malloc(sizeof(Blur_Profile_RGB));
    if (!profile) {
        fprintf(stderr, "ERROR: Memory allocation for blur profile failed.\n");
        return NULL;
    }
    profile->num_angle_bins = num_angle_bins;
    profile->num_radius_bins = num_radius_bins;
    profile->angle_bin_size = (double)(180 / num_angle_bins);
    double max_radius = sqrt(fft->width*fft->width + fft->height*fft->height/4);
    profile->radius_bin_size = (double)(max_radius / num_radius_bins);
    double radius_bin_size_sq = profile->radius_bin_size*profile->radius_bin_size;

    // Allocate space for empty bins
    profile->r_bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    profile->g_bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    profile->b_bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    if (!profile->r_bins || !profile->g_bins || !profile->b_bins) {
        fprintf(stderr, "ERROR: Memory allocation failed for r,g, or b bins.\n");
        return NULL;
    }
    for (int i=0; i<num_angle_bins; i++) {
        profile->r_bins[i] = (Bin*)calloc(num_radius_bins, sizeof(Bin));
        profile->g_bins[i] = (Bin*)calloc(num_radius_bins, sizeof(Bin));
        profile->b_bins[i] = (Bin*)calloc(num_radius_bins, sizeof(Bin));
    }

    // Sum up blurs in each bin space
    int conversion_height = conversion->height; // No need to perform member access every iteration
    int conversion_width = conversion->width;
    int i_tot = fft->height * fft->width;       // Store once to speed up by 1 multiply & 2 member accesses
    const Pixel* r = fft->r;    // locally store array pointers to assure high performance in dereferencing
    const Pixel* g = fft->g;
    const Pixel* b = fft->b;
    double max_phi = 0;

    // bin_quant counts number of pixels in each bin so the avg can be calculated
    int** bin_quant = (int**)malloc(num_angle_bins* sizeof(int*));
    for (int i=0; i<num_angle_bins; i++) {
        bin_quant[i] = (int*)calloc(num_radius_bins, sizeof(int));
    }

    for (int i=0; i<i_tot; i++) {
        // Store r_sq and phi for readability
        int r_sq = conversion->data[i].r_sq;
        double phi = conversion->data[i].phi;
        if (max_phi < phi) max_phi = phi;
        // phi_bin is found by int division of the phi value by the bin size
        // +num_angle_bins/2 to normalize the phi range of [-90deg, 90deg) to [0, num_angle_bins]
        int phi_bin = (int)((phi+PI*0.5f)/PI * (double)(profile->num_angle_bins-1));
        // r_bin is sqrt(r^2/r_bin_size^2), but a newton int approximation
        int r_bin = newton_int_sqrt((double)r_sq/radius_bin_size_sq);
        bin_quant[phi_bin][r_bin]++;    // Increment the counter for how many pixels have been assigned to this bin
        profile->r_bins[phi_bin][r_bin] += r[i];    // Profile not locally stored like fft members because
        profile->g_bins[phi_bin][r_bin] += g[i];    // due to small size, it will likely all be stored in
        profile->b_bins[phi_bin][r_bin] += b[i];    // CPU cache anyway. Access should be inexpensive
    }

    // Divide by the number of bins
    Bin bin_quantity;
    for (int phi_bin=0; phi_bin<num_angle_bins; phi_bin++) {
        for (int r_bin=0; r_bin<num_radius_bins; r_bin++) {
            bin_quantity = (Bin)bin_quant[phi_bin][r_bin];
            profile->r_bins[phi_bin][r_bin] /= bin_quantity;
            profile->g_bins[phi_bin][r_bin] /= bin_quantity;
            profile->b_bins[phi_bin][r_bin] /= bin_quantity;
        }
    }
    return profile;
}


/******************************************************************************
 * create_rgb_image creates an image of given width and height, being a pointer
 *   to an Image_RGB data-type.
 *  -Function initializes all pixels to 0.
 *  -Function returns NULL if any allocations fail
 *  -Function returns NULL if height or width are equal to 0.
******************************************************************************/
Image_RGB* create_rgb_image(unsigned int width, unsigned int height) {
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
Image_HSV* create_hsv_image(unsigned int width, unsigned int height) {
    if (!height) {  // If height is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_rgb_image received height==0.\n");
        return NULL;
    } else if (!width) {   // If width is 0, throw error and return NULL
        fprintf(stderr, "ERROR: create_rgb_image received width == 0.\n");
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
Image_PGM* create_pgm_image(unsigned int width, unsigned int height) {
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
    unsigned int width = right - left;
    unsigned int height = bottom - top;
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
 * rgb_fft calculates the FFT of an RGB image given.
 *  -If no fftw_plan* is given, then the function will automatically generate
 *   an fftw_plan* to be used, based on the inputted image height.
******************************************************************************/
Image_RGB* rgb_fft(Image_RGB* image) {
    // Allow threading
    fftw_init_threads();
    fftw_plan_with_nthreads(num_cores);

    // Create output and input
    fftw_complex* output = fftw_alloc_complex((image->width/2+1) * image->height);
    double* in = fftw_alloc_real(image->height * image->width);

    if (!output || !in) {
        fprintf(stderr, "ERROR: Memory allocation for fftw_malloc failed.\n");
        fftw_free(output);
        fftw_free(in);
        return NULL;
    }

    fftw_plan plan = fftw_plan_dft_r2c_2d(image->height, image->width, in, output, FFTW_ESTIMATE);
    if (!plan) {
        fprintf(stderr, "ERROR: fftw_plan generation failed.");
        fftw_destroy_plan(plan);
        return NULL;
    }

    int fft_width = image->width/2+1;
    Image_RGB* fft_image = create_rgb_image(fft_width, image->height);
    int i_tot = fft_image->width * fft_image->height;

    // Execute fftw for red, collect magnitudes^2 in fft_image
    memcpy(in, image->r, sizeof(double) * image->width * image->height);
    fftw_execute(plan);
    for (int i=0; i<i_tot; i++) {
        fft_image->r[i] = output[i][0]*output[i][0] + output[i][1]*output[i][1];
    }

    // Execute fftw for green, collect magnitudes^2 in fft_image
    memcpy(in, image->g, sizeof(double)*image->width*image->height);
    fftw_execute(plan);
    for (int i=0; i<i_tot; i++) {
        fft_image->g[i] = output[i][0]*output[i][0] + output[i][1]*output[i][1];
    }

    // Execute fftw for green, collect magnitudes^2 in fft_image
    memcpy(in, image->b, sizeof(double)*image->width*image->height);
    fftw_execute(plan);
    for (int i=0; i<i_tot; i++) {
        fft_image->b[i] = output[i][0]*output[i][0] + output[i][1]*output[i][1];
    }

    // Destroy plan, free memory, clean up threads and traces of FFTW
    fftw_destroy_plan(plan);
    fftw_free(in);
    in = NULL;
    fftw_free(output);
    output = NULL;
    fftw_cleanup_threads();
    fftw_cleanup();

    // Return final magnitude^2 fft image
    return fft_image;
}


/******************************************************************************
 * create_path returns a cleaned file path, retrieved from user input.
 *  -path must be the pre-filename path from the project root to the folder
 *   that is meant to hold the file
 *  -request_string is the message displayed to the user, requesting an input.
 *  -filetype is the file extension to be used.
 *  -No restraints are set on input variables or user input, as interface is
 *   only meant for developer testing and needs no safety constraints.
******************************************************************************/
char* create_path(const char* path, const char* request_string, const char* filetype) {
    // Request, get, and clean the filename
    char filename[40];
    printf("%s\n", request_string);
    fgets(filename, sizeof(filename), stdin);
    int len = strlen(filename);
    if (len > 0 && filename[len-1] == '\n') {
        filename[len-1] = '\0';
    }
    // Concatenate the filetype to the filename
    strcat(filename, filetype);

    // allocate memory for full path
    char* full_path = (char*)malloc(sizeof(char)*(strlen(path)+strlen(filename)+1));
    if (!full_path) {
        fprintf(stderr, "ERROR: Memory allocation for full_path failed.\n");
        return NULL;
    }

    // Assemble full path from filename and path
    strcpy(full_path, path);
    strcat(full_path, filename);
    printf("%s\n", full_path);
    return full_path;
}


/******************************************************************************
 * get_fft_normalizer_lookup returns a double array pointer to the values in
 *   the lookup table fft_normalizer.txt.
 *  -Returns NULL if fft_normalizer.txt is unable to be opened.
 *  -Returns NULL if fft_normalizer.txt is not a 1D lookup table.
 *  -Returns NULL if at any point the table does not follow proper formatting.
 *  -proper lookup_table formatting is defined as:
 *      {dimensionalidy}D\n
 *      {length}\n
 *      {input[0]} {output[0]}\n
 *      {input[0]} {output[0]}\n
 *      ...
 *      ... 
 *      {input[length-1]} {output[length-1]}\n
******************************************************************************/
Lookup_1D* get_fft_normalizer_lookup() {
    FILE* file = fopen("lookups/fft_normalizer.txt", "r");
    if (file == NULL) {
        fprintf(stderr, "Error occured opening lookups/fft_normalizer.txt.\n");
        free(file);
        return NULL;
    }

    char dimensionality[3];  // Buffer to hold the dimensionality string

    // Read the dimensionality string (expecting "1D\n")
    if (fscanf(file, "%2s\n", dimensionality) != 1) {
        fprintf(stderr, "Error reading dimensionality.\n");
        return NULL;
    }

    // Check if the dimensionality is "1D"
    if (dimensionality[0] != '1' || dimensionality[1] != 'D') {
        fprintf(stderr, "Invalid dimensionality: %s\n", dimensionality);
        fprintf(stderr, "get_fft_normalizer_lookup requires a 1Dx1D lookup table.");
        return NULL;
    }

    // Initialize Lookup_1D table to proper size
    Lookup_1D* fft_normalizer_lookup = malloc(sizeof(Lookup_1D));
    if (fscanf(file, "%u\n", &fft_normalizer_lookup->length) != 1) { // Read the length
        fprintf(stderr, "Error reading length.\n");
        return false;
    }
    fft_normalizer_lookup->input = (Pixel*)malloc(fft_normalizer_lookup->length * sizeof(Pixel));
    fft_normalizer_lookup->output = (Pixel*)malloc(fft_normalizer_lookup->length * sizeof(Pixel));
    if (!fft_normalizer_lookup->input || !fft_normalizer_lookup->output) {
        fprintf(stderr, "ERROR: Memory allocation for fft_normalizer conversion table failed.\n");
        return NULL;
    }

    for (int i=0; i<fft_normalizer_lookup->length; i++) {
        if (fscanf(file, "%lf %lf", &fft_normalizer_lookup->input[i], &fft_normalizer_lookup->output[i]) != 2) {
            fprintf(stderr, "Error reading value at index %d\n", i);
            free(fft_normalizer_lookup);
            fclose(file);
            return NULL;
        }
    }
    #ifdef DEBUG    // Print the values that were returned
        for (int i=0; i<fft_normalizer_lookup->length; i++) {
            fprintf(stderr, "input[%d]: %lf\toutput[%d]: %lf\n", 
            i, fft_normalizer_lookup->input[i], i, fft_normalizer_lookup->output[i]);
        }
    #endif

    fclose(file);
    return fft_normalizer_lookup;
}


// Frees a Lookup_1D table and nullifies its pointer
void free_1D_lookup_table(Lookup_1D* table) {
    free(table->input), table->input = NULL;
    free(table->output), table->input = NULL;
    free(table), table = NULL;
}


/******************************************************************************
 * rgb_normalize_fft takes in an fft with 3 channels and normalizes all values
 *   to within the range [0,1].
 *  -Function calculates one global max and applies it to all 3 channels.
 *  -Normalization uses a mathematical model found here:
 *   https://www.desmos.com/calculator/bm5nnnk7oo
 *      -h(x; G_s) is the function which defines what discrete value within [0, 1]
 *       is assigned given x and G_s.
 *      -x represents any value within the original FFT image
 *      -G_s represents the maximum value in the original FFT image
 *  -If DEBUG is set to 1, speed performance will suffer, as info messages for
 *   maximums and scaled values will be shown.
******************************************************************************/
void rgb_normalize_fft(Image_RGB* fft, Lookup_1D* fft_normalizer_lookup) {
    double max = fft->r[0];     // Set Maximum as any valid value
    int i_tot = fft->height * fft->width;   // Get stopping point
    int max_index = 0;

    // Find maximum input
    for (int i=0; i<i_tot; i++) {
        double r_pixel = fft->r[i]; // Set temporary values for efficiency
        double g_pixel = fft->g[i];
        double b_pixel = fft->b[i];
        if (max < r_pixel) {max = r_pixel; max_index = i;}
        if (max < g_pixel) {max = g_pixel;}
        if (max < b_pixel) {max = b_pixel;}
    }

    #ifdef DEBUG
        printf("Maximum value found in the raw FFT: %f\n", max);
    #endif
    
    // Set G_s
    Pixel G_s = 1/(2*log(sqrt(max) + 1));
    
    // calculate normalized values
    for (int i=0; i<i_tot; i++) {
        if (fft->r[i] < 1) fft->r[i] = 0;
        else fft->r[i] = log(fft->r[i]) * G_s;
        if (fft->g[i] < 1) fft->g[i] = 0;
        else fft->g[i] = log(fft->g[i]) * G_s;
        if (fft->b[i] < 1) fft->b[i] = 0;
        else fft->b[i] = log(fft->b[i]) * G_s;
    }

    #ifdef DEBUG
        // Set Maximum as any valid value
        max = fft->r[0];

        // Set maximum values
        for (int i=0; i<i_tot; i++) {
            double r_pixel = fft->r[i]; // Set temporary values for efficiency
            double g_pixel = fft->g[i];
            double b_pixel = fft->b[i];
            if (max < r_pixel) {max = r_pixel;}
            if (max < g_pixel) {max = g_pixel;}
            if (max < b_pixel) {max = b_pixel;}
        }
        printf("Maximum value found in the normalized FFT: %f\n", max);
    #endif
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
void validate_coordinates(Image_RGB* array, int x, int y, int* errors) {
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
 * fft_shift performs a quadrant swap and rotation on an FFT image.
 * 
 * This function creates a new FFT image with dimensions (2*width-1, height)
 * and then performs a quadrant swap and 180-degree rotation to shift the FFT.
 * The swap and rotation are applied to each of the R, G, and B channels.
 * 
 * Parameters:
 *  - fft: Pointer to the original FFT image.
 * 
 * Returns:
 *  - Pointer to the new FFT image with quadrant swap and rotation applied.
 * 
 * The function iterates through one quarter of the required image area (half
 * of the original image dimensions). It handles even and odd-sized images
 * differently, using compensators for height and width. The function swaps
 * quadrants Q1 with Q4 and Q2 with Q3, and also rotates Q2 and Q3 by 180 degrees.
 * Coordinate validation is performed to ensure that array access is within bounds.
 ******************************************************************************/
Image_RGB* fft_shift(Image_RGB* fft) {
    Image_RGB* fft_image = create_rgb_image(fft->width*2 - 1, fft->height);
    // Iterate through one quarter of needed image, via one half of the image
    printf("new_image size: (%3d, %3d)\t old_image size: (%3d, %3d)\n", fft_image->width, fft_image->height, fft->width, fft->height);
    
    // If even image, height_compensator is 0. if odd image, compensator is 1
    int height_compensator = (int)(fft_image->height%2==1);
    int width_compensator = (int)(fft_image->width%2==1); // width comp is always 1, but I want to make it dynamic.
    printf("height_equalizer: %d\t width_equalizer: %d\n", height_compensator, width_compensator);

    int errors = 0;
    Pixel r_temp;
    Pixel g_temp;
    Pixel b_temp;
    int half_fft_height = fft->height/2;
    int y_val, x_val;
    int half_height_plus_compensator = half_fft_height+height_compensator;
    int fft_width = fft->width;
    for (int y=0; y<half_height_plus_compensator; y++) {
        for (int x=0; x<fft_width; x++) {
            /*** Swap Q1 and Q4 ***/
            y_val = y+half_fft_height, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // New Q1 gets old Q4
            r_temp = fft->r[y_val*fft_width + x_val];
            g_temp = fft->g[y_val*fft_width + x_val];
            b_temp = fft->b[y_val*fft_width + x_val];
            if (y-height_compensator != -1) { // If height is odd, middle row goes in a theoretical pixels[-1]
                y_val = y-height_compensator, x_val = x+fft->width-width_compensator, validate_coordinates(fft_image, x_val, y_val, &errors);
                fft_image->r[y_val*fft_width + x_val] = r_temp;
                fft_image->g[y_val*fft_width + x_val] = g_temp;
                fft_image->b[y_val*fft_width + x_val] = b_temp;
            }
            y_val = y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // New Q4 gets old Q1
            r_temp = fft->r[y_val * fft_width + x_val];
            g_temp = fft->g[y_val * fft_width + x_val];
            b_temp = fft->b[y_val * fft_width + x_val];
            y_val = y+half_fft_height, x_val = x+fft->width-width_compensator, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->r[y_val * fft_width + x_val] = r_temp;
            fft_image->g[y_val * fft_width + x_val] = g_temp;
            fft_image->b[y_val * fft_width + x_val] = b_temp;

            /*** Rotate to get Q2 and Q3 ***/
            y_val = y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // Q2 is 180 deg rotation of Q4
            r_temp = fft->r[y_val * fft_width + x_val];
            g_temp = fft->g[y_val * fft_width + x_val];
            b_temp = fft->b[y_val * fft_width + x_val];
            y_val = half_fft_height-y, x_val = fft->width-1-x, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->r[y_val * fft_width + x_val] = r_temp; //width-1 because width-1=last_index
            fft_image->g[y_val * fft_width + x_val] = g_temp;
            fft_image->b[y_val * fft_width + x_val] = b_temp;
            y_val = half_fft_height+y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // Q3 is 180 deg rotation of Q1
            r_temp = fft->r[y_val * fft_width + x_val];
            g_temp = fft->g[y_val * fft_width + x_val];
            b_temp = fft->b[y_val * fft_width + x_val];
            y_val = fft->height-1-y, x_val = fft->width-1-x, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->r[y_val * fft_width + x_val] = r_temp; // height-1 because height-1=last_index
            fft_image->g[y_val * fft_width + x_val] = g_temp;
            fft_image->b[y_val * fft_width + x_val] = b_temp;
        }
    }
    return fft_image;
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
 * get_blur_profile_visual creates an Image_RGB* that represents the contents
 *   of a measured blur_profile visually.
 *  -height should be the fft height
 *  -width should be the fft widht
 *  -The user may choose to change the height and width, but should make sure
 *   that their aspect ratio stays the same
 *  -The outputted image should look like a radially pixelated version of the
 *   actual fft, and as num_radius_bins and num_angle_bins increases it should
 *   come to approximate the fft itself.
******************************************************************************/
Image_RGB* get_blur_profile_visual(Blur_Profile_RGB* blur_profile, 
                                   Cartesian_To_Polar* conversion, 
                                   int height, int width) {
    // Allocate memory for an RGB image, with height and width
    Image_RGB* output_image = create_rgb_image(width, height);
    if (!output_image) {
        fprintf(stderr, "Error creating output image.\n");
        return NULL;
    }
    int y_times_width = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double deltaX = x;
            double deltaY;
            if (y<height/2) deltaY = -y;
            else deltaY = height - y;

            double r = sqrt(deltaX * deltaX + deltaY * deltaY);
            double phi = atan2(deltaY, deltaX);
            // if (phi < 0) phi += 360; // Normalize phi to [0, 360) range

            int r_bin = (int)(r / blur_profile->radius_bin_size);
            if (r_bin >= blur_profile->num_radius_bins) r_bin = blur_profile->num_radius_bins - 1;

            // +num_angle_bins/2 to normalize possible angles from [-90deg, 90deg) to [0, num_angle_bins]
            int phi_bin = (int)((phi+PI*0.5f)/PI * (double)(blur_profile->num_angle_bins-1));

            if (phi_bin >= blur_profile->num_angle_bins) {
                phi_bin = blur_profile->num_angle_bins - 1;
            }
            if (phi_bin < 0){
                phi_bin = 0;
            }

            // Direct assignment without inappropriate mirroring
            output_image->r[y_times_width + x] = blur_profile->r_bins[phi_bin][r_bin];
            output_image->g[y_times_width + x] = blur_profile->g_bins[phi_bin][r_bin];
            output_image->b[y_times_width + x] = blur_profile->b_bins[phi_bin][r_bin];
        }
        y_times_width += width;
    }
    return output_image;
}


/******************************************************************************
 * view_2d_array_contents simply prints out the values of a 2D array after
 *   normalizing its values to the range [0,255].
 *  -Assumes that the array is in the range [0,1].
 *  -Height and width must be accurate (duh).
******************************************************************************/
void view_2d_array_contents(Bin** array, int height, int width) {
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            printf("%d, ", (int)(array[i][j]*255.0f));
        }
        printf("\n");
    }
    printf("\n");
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
    unsigned short height = rgb->height;
    unsigned short width = rgb->width;
    // Initialize 3-channel image as hsv, with all 0s
    Image_HSV* hsv = (Image_HSV*)create_hsv_image(width, height);
    // Iterate through every pixel
    unsigned int i_tot = height*width;
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
        cur->v = max;

        // Set S value
        if (max==0) {cur->s = 0;}
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
    unsigned short height = hsv->height;
    unsigned short width = hsv->width;

    // Initialize 3-channel image as rgb, assuming create_rgb_image is similar to create_hsv_image
    Image_RGB* rgb = (Image_RGB*)create_rgb_image(width, height);
    if (!rgb) return NULL; // Check if memory allocation failed

    unsigned int i_tot = height * width;
    for (unsigned int i = 0; i < i_tot; i++) {
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
 * free_cartesian_to_polar frees the memory of Cartesian_To_Polar objects.
******************************************************************************/
void free_cartesian_to_polar(Cartesian_To_Polar* c2p) {
    if (c2p) {
        free(c2p->data);
        c2p->data = NULL;
        free(c2p);
    }
}


/******************************************************************************
 * blur_profile_tests takes an fft and saves the blur_profile visualization.
 *  -Function is used for development, as it prompts the user for command-line
 *   inputs.
******************************************************************************/
void blur_profile_tests(Image_RGB* fft) {
    // Calculate bin-approximated FFT
    Cartesian_To_Polar* conversion;   // Create cartesian to polar conversion
    conversion = cartesian_to_polar_conversion(fft->width, fft->height);
    int num_radius_bins = 4;    // Set up pseudo-hyperparameters
    int num_angle_bins = 18;
    Blur_Profile_RGB* blur_profile; // Calculate the blur profile
    blur_profile = calculate_blur_profile(conversion, fft, num_radius_bins, num_angle_bins);

    // Create image representation of bin-approximated FFT
    Image_RGB* blur_profile_visualization = get_blur_profile_visual(blur_profile, conversion, fft->height, fft->width);

    // Save image representation of bin-approximated FFT
    const char* path = "images/visualizations/";
    const char* q = "\nEnter visualization filename: ";
    char* visualization_filename = create_path(path, q, ".txt");
    write_image_to_file(blur_profile_visualization, visualization_filename);

    // Free all function data
    free_cartesian_to_polar(conversion);
    free(visualization_filename), visualization_filename = NULL;
    free_image_rgb(blur_profile_visualization);
}


/******************************************************************************
 * compute_magnitude_fft returns an Image_RGB* to a normalized fft, given an
 *   rgb image.
******************************************************************************/
Image_RGB* compute_magnitude_fft(Image_RGB* image) {
    Image_RGB* fft = rgb_fft(image);
    rgb_normalize_fft(fft, NULL);
    return fft;
}

/******************************************************************************
 * save_fft writes an fft image to user's desired image file (in .txt format 
 *   for python to save via imageSaver.py).
 *  -Function is designed for development, as it prompts users for command-line
 *   input.
******************************************************************************/
void save_fft(Image_RGB* fft) {
    const char* q = "Enter the name of the .txt file you wish to write to: ";
    char* output_filename = create_path("images/output/", q, ".txt");
    write_image_to_file(fft, output_filename);
    free(output_filename), output_filename = NULL;
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


// Initialization to work with CPU cores
int threading_setup() {
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
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
Image_RGB* downsample_rgb(Image_RGB* image, unsigned short N) {
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
 * initialize_octree sets up an octree according to hyperparamters and returns
 *   a pointer to the octree.
 *  -h_parts, s_parts, v_parts, and num_grays may not be 0.
 *  -Function does NOT initialize valid_parents array, as it is not known upon
 *   initialization. valid_parents will be NULL.
******************************************************************************/
Octree* initialize_octree(int h_parts, 
                          int s_parts,
                          int v_parts, 
                          int num_grays, 
                          double black_thresh,
                          double gray_thresh) {
    if (h_parts==0 || s_parts==0 || v_parts==0 || num_grays == 0) {
        fprintf(stderr, "ERROR: initialize_octree() expects nonzero values for h_parts, s_parts, v_parts, and num_grays.\n");
        return NULL;
    }

    // Create octree, set all measurement values
    Octree* octree = (Octree*)malloc(sizeof(Octree));
    if (!octree) {
        fprintf(stderr, "ERROR: Malloc() of Octree failed.\n");
        return NULL;
    }
    octree->total_length = h_parts*s_parts*v_parts+num_grays+1;
    octree->num_h = h_parts;
    octree->Lh = 360/h_parts;
    octree->num_s = s_parts;
    octree->Ls = (1-gray_thresh)/s_parts; // length(s) = (max_s - gray_thresh)/(number_s_parts)
    octree->num_v = v_parts;
    octree->Lv = (1-black_thresh)/v_parts; // length(v) = (max-black_thresh)/(number_v_parts)
    octree->num_grays = num_grays;
    octree->black_thresh = black_thresh;
    octree->gray_thresh = gray_thresh;

    // Create Octree Groups
    octree->groups = (Octree_Group*)malloc(octree->total_length*sizeof(Octree_Group));
    if (!octree->groups) {
        fprintf(stderr, "ERROR: Malloc() of Octree_Group failed.\n");
        return NULL;
    }

    // Initialize hue groups from h, s, and v partitions.
    double half_h = octree->Lh/2;
    double s_offs = octree->Ls/2 + gray_thresh;     // s_offs = len(s)/2 + gray_thresh
    double v_offs = octree->Lv/2 + black_thresh;    // v_offs = len(v)/2 + black_thresh
    int i;
    for (int h=0; h<h_parts; h++) {
        for (int s=0; s<s_parts; s++) {
            for (int v=0; v<v_parts; v++) {
                i = h*s_parts*v_parts + s*v_parts + v;
                octree->groups[i].h = h*octree->Lh + half_h;
                octree->groups[i].s = s*octree->Ls + s_offs;
                octree->groups[i].v = v*octree->Lv + v_offs;
                octree->groups[i].id = i;
                octree->groups[i].quantity = 0;
                octree->groups[i].head = NULL;
            }
        }
    }

    // Initialize gray groups
    double L_gray = (1.0f-black_thresh)/(double)num_grays;
    int total_length = octree->total_length;
    for (int i=total_length-num_grays; i<total_length-1; i++) {
        octree->groups[i].h = 0;
        octree->groups[i].s = 0;
        octree->groups[i].v = L_gray*i + black_thresh;
        octree->groups[i].id = i;
        octree->groups[i].quantity = 0;
        octree->groups[i].head = NULL;
    }
    
    // Initialize black group
    octree->groups[total_length-1].h = 0;
    octree->groups[total_length-1].s = 0;
    octree->groups[total_length-1].v = 0;
    octree->groups[total_length-1].id = total_length-1;
    octree->groups[total_length-1].quantity = 0;
    octree->groups[total_length-1].head = NULL;

    return octree;
}


/******************************************************************************
 * get_hsv_linked_list_node creates a new HSV_Linked_List with an array of 
 *   arr_size. if it is unable to allocate memory, it prints an error to stderr
 *   and returns NULL.
 *  -arr_size > 0 is required. if arr_size < 0, Function returns an stderr and
 *   a NULL pointer.
******************************************************************************/
HSV_Linked_List* get_hsv_linked_list_node(int arr_size) {
    if (arr_size == 0) {
        fprintf(stderr, "ERROR: get_hsv_linked_list_node() expected a positive arr_size");
        return NULL;
    }
    HSV_Linked_List* node = (HSV_Linked_List*)calloc(1, sizeof(HSV_Linked_List));
    if (!node) {
        fprintf(stderr, "ERROR: failed to calloc node in get_hsv_linked_list_node()");
        return NULL;
    }
    node->num_pixels = arr_size;
    node->pixels = (Pixel_HSV*)calloc(arr_size, sizeof(Pixel_HSV));
    if (!node->pixels) {
        fprintf(stderr, "ERROR: failed to calloc pixels array in get_hsv_linked_list_node()");
        return NULL;
    }
    return node;
}


/*****************************************************************************
 * arm_octree assigns all pixels of an HSV image into the proper groups of the 
 *   octree given.
******************************************************************************/
void arm_octree(Image_HSV* hsv, Octree* octree, int HSV_Linked_List_Size) {
    // Create array of current octree group linked list nodes for use in the initialization.
    HSV_Linked_List** list_curs = (HSV_Linked_List**)malloc(octree->total_length * sizeof(HSV_Linked_List*));
    if (!list_curs) {
        fprintf(stderr, "ERROR: arm_octree() failed to malloc list_curs.");
        return NULL;
    }
    for (int i=0; i<octree->total_length; i++) {
        // Set array values all to the heads of their respective linked lists.
        list_curs[i] = octree->groups[i].head;
    }

    // Loop through all pixels
    int i_tot = hsv->height * hsv->width;
    for (int i=0; i<i_tot; i++) {
        // Vi, Si, and Hi are index i trackers for H,S,V values.
        int Vi, Si, Hi, g;
        // If value is under black threshold, set g to the black grouping.
        if (hsv->pixels[i].v < octree->black_thresh){
            g = octree->total_length;
        }
        // if current pixel falls under gray groupings, set g to its proper gray grouping.
        else if (hsv->pixels[i].v < octree->gray_thresh) {
            Vi=(int)(hsv->pixels[i].v-octree->black_thresh)*octree->num_grays/(1-octree->black_thresh);
            g = octree->total_length-octree->num_grays+Vi;
        }
        // if current pixel is not gray or black, set g to the proper color grouping.
        else {
            Vi=(int)((hsv->pixels[i].v-octree->black_thresh)/octree->Lv);
            Si=(int)((hsv->pixels[i].s-octree->gray_thresh)/octree->Ls);
            Hi=(int)(hsv->pixels[i].h/octree->Lh);
            g = (Hi*octree->num_s+Si)*octree->num_v + Vi;
        }

        // If the current linked list is NULL, then create one for its space.
        if (!list_curs[g]) {
            list_curs[g] = get_hsv_linked_list_node(HSV_Linked_List_Size);
        }
        // If the current linked list is filled, allocate a new list and move cur to it.
        else if (list_curs[g]->i == list_curs[g]->num_pixels) {
            list_curs[g]->next = get_hsv_linked_list_node(HSV_Linked_List_Size);
            list_curs[g] = list_curs[g]->next;
        }
        // The first unfilled pixel in the g-th list_cur gets the current hsv pixel
        // The i inside of list_curs[g] is NOT the same i as the for loop uses!!
        list_curs[g]->pixels[list_curs[g]->i++] = hsv->pixels[i];
        octree->groups[g].quantity++;
    }
}


/******************************************************************************
 * compare_quantity is used to compare the quantities of two octree groups.
******************************************************************************/
int compare_quantities(const int* x, const int* y, const Octree_Group* groups) {
    return groups[*y].quantity - groups[*x].quantity; // For descending order
}


/******************************************************************************
 * find_valid_parents calculates the quantities of pixels in each octree group
 *   and assigns valid_parent=true to all octree groups that are required to 
 *   reach the hyperparameter threshold.
 *  -Function also adds group id to octree->valid_parents array.
 *  -total_pixels is the total number of pixels in the image.
 *  -coverage_threshold is the hyperparameter for the threshold for how many of
 *   the pixels should be a part of the main colors in the color palette.
 *   should be a number in the range [0, 1]
******************************************************************************/
void find_valid_octree_parents(Octree* octree, int total_pixels, double coverage_threshold) {
    // sort an array of ids in order from least to greatest
    int* sorted_ids = malloc(octree->total_length * sizeof(int));
    for (int i=0; i<octree->total_length; i++) {
        sorted_ids[i] = i;
    }
    qsort_r(sorted_ids, octree->total_length, sizeof(int), compare_quantities, &octree->groups);

    // Fill out values of valid_parents array until coverage_threshold is hit
    int goal_num_pixels = (int)((double)total_pixels*coverage_threshold);
    int i;
    for (i=0; i<octree->total_length; i++) {
        goal_num_pixels -= octree->groups[sorted_ids[i]].quantity;
        // if goal_num_pixels is reached, set up valid_parents as 0 to i, and return.
        if (goal_num_pixels <= 0) {
            int* valid_parents = (int*)malloc(i*sizeof(int));
            for(int j=0; j<i; j++) {
                valid_parents[j] = sorted_ids[j];
            }
            octree->valid_parents = valid_parents;
            octree->len_valid_parents = i;
            return;
        }
    }
    fprintf(stderr, "ERROR: find_valid_octree_parents should not reach the end of its valid_parents loop");
}


int main() {
    // Get image from files
    Image_RGB* image = read_image_from_files();

    // Downsample image haphazardly
    // Image_RGB* downsampled = downsample_rgb(image, 10);

    // Save new rgb
    // save_rgb(downsampled);

    // Compute the FFT for RGB image
    // Image_RGB* fft = compute_magnitude_fft(image);

    // Save fft image to .txt file
    // save_fft(fft);
    
    // Run blur_profile code
    // blur_profile_tests(fft);

    // Free data
    // free_image_rgb(fft), fft = NULL;

    // free_image_rgb(downsampled);
    free_image_rgb(image);
    return 0;
}
