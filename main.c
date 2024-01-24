#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>

#define DEBUG
#define FFT_NORMALIZER_LOOKUP_SIZE 256


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
 *  -Effective as both an FFT representation and a literal image.
 *  -Cannot hold complex values
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
int max_index;


/******************************************************************************
 * cartesian_to_polar_conversion creates a mapping for cartesian to polar conversion.
 *  -Built to dynamically compute for a given FFT height and width
 *  -Returns Cartesian_To_Polar object pointer
 *  -CAUTION: Assumes that FFT image is NOT shifted, such that the DC component
 *   is found in the top and bottom-left corners of the image, rather than in
 *   the center-left side.
******************************************************************************/
Cartesian_To_Polar* cartesian_to_polar_conversion(unsigned int width, unsigned int height) {
    Polar_Coord coord;
    Cartesian_To_Polar* conversion = (Cartesian_To_Polar*)malloc(sizeof(Cartesian_To_Polar*));
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
    for (int y=0; y<height/2; y++) {    // height/2 for symmetry across x axis
        for (int x=0; x<width; x++) {
            int r_sq = x*x+y*y;
            double phi = atan(y/x);
            // Not fft shifted, so for the top half, phi=-atan(y/x)
            conversion->data[y + x].phi = -phi;
            conversion->data[y + x].r_sq = r_sq;
            // Not fft shifted, so for the bottom half phi=atan(y/x)
            conversion->data[y + x].phi = phi;
            conversion->data[y-half_height + x].r_sq = r_sq;
        }
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
    profile->radius_bin_size = (double)(fft->width / num_radius_bins);
    double radius_bin_size_sq = profile->radius_bin_size*profile->radius_bin_size;

    // Allocate space for empty bins
    profile->r_bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    profile->g_bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    profile->b_bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    if (!profile->r_bins || !profile->g_bins || !profile->b_bins) {
        fprintf(stderr, "ERROR: Memory allocation failed for r,g, or b bins.\n");
        return NULL;
    }
    for (int i=0; i<num_radius_bins; i++) {
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
    for (int i=0; i<i_tot; i++) {
        // Store r_sq and phi for readability
        int r_sq = conversion->data[i].r_sq;
        double phi = conversion->data[i].phi;
        // phi_bin is found by int division of the phi value by the bin size
        int phi_bin = (int)phi/profile->angle_bin_size;
        // r_bin is sqrt(r^2/r_bin_size^2), but a newton int approximation
        int r_bin = newton_int_sqrt((double)r_sq/radius_bin_size_sq);
        profile->r_bins[phi_bin][r_bin] += r[i];    // Profile not locally stored like fft members because
        profile->g_bins[phi_bin][r_bin] += g[i];    // due to small size, it will likely all be stored in
        profile->b_bins[phi_bin][r_bin] += b[i];    // CPU cache anyway. Access should be inexpensive
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
    Image_RGB* image = (Image_RGB*)malloc(sizeof(Image_RGB));
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
        printf("Maximum value found in the normalized FFT: %f\n", max);
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


int main() {
    num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    fprintf(stderr, "num_cores: %d\n", num_cores);
    Image_RGB* image;
    {   // read in user's desired image from the list of readable .txt files
        char* input_filename = create_path("images/readable/", "Enter the name of the .txt file you wish to use: ", ".txt");
        image = read_image(input_filename);
        free(input_filename), input_filename = NULL;
        if (image == NULL) {
            return -1;
        }
    }

    Image_RGB* fft;
    {   // Compute the FFT for RGB image
        fft = rgb_fft(image);
        rgb_normalize_fft(fft, NULL);
        free_image_rgb(image), image=NULL;
    }

    {   // write fft_image to user's desired image file (in .txt format for python to save via imageSaver.py)
        char* output_filename = create_path("images/output/", "Enter the name of the .txt file you wish to write to: ", ".txt");
        write_image_to_file(fft, output_filename);
        free(output_filename), output_filename = NULL;
    }

    // Free original image last since it's used periodically during runtime.
    free_image_rgb(fft), fft = NULL;

    return 0;
}
