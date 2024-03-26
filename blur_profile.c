#include "image_processing.h"
#include "blur_profile.h"
#include "utilities.h"
#include "fft_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

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
    double radius_bin_size_sq = (double)(profile->radius_bin_size*profile->radius_bin_size);

    radius_bin_size_sq = (double)((fft->width*fft->width + fft->height*fft->height/4) / (num_radius_bins*num_radius_bins));

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

    START_TIMING(accumulate_bins_time);
    for (int i=0; i<i_tot; i++) {
        // Store r_sq and phi for readability
        int r_sq = conversion->data[i].r_sq;
        double phi = conversion->data[i].phi;
        if (max_phi < phi) max_phi = phi;
        // phi_bin is found by int division of the phi value by the bin size
        // +num_angle_bins/2 to normalize the phi range of [-90deg, 90deg) to [0, num_angle_bins]
        int phi_bin = (int)((phi+PI*0.5f)/PI * (double)(profile->num_angle_bins-1));
        // r_bin is sqrt(r^2/r_bin_size^2), but a newton int approximation
        int r_bin = newton_int_sqrt(((double)r_sq)/radius_bin_size_sq);
        bin_quant[phi_bin][r_bin]++;    // Increment the counter for how many pixels have been assigned to this bin
        profile->r_bins[phi_bin][r_bin] += r[i];    // Profile not locally stored like fft members because
        profile->g_bins[phi_bin][r_bin] += g[i];    // due to small size, it will likely all be stored in
        profile->b_bins[phi_bin][r_bin] += b[i];    // CPU cache anyway. Access should be inexpensive
    }
    END_TIMING(accumulate_bins_time, "accumulating pixels to blur_profile bins");

    // Divide by the number of bins
    START_TIMING(average_bin_time);
    Bin bin_quantity;
    for (int phi_bin=0; phi_bin<num_angle_bins; phi_bin++) {
        for (int r_bin=0; r_bin<num_radius_bins; r_bin++) {
            bin_quantity = (Bin)bin_quant[phi_bin][r_bin];
            if (bin_quantity != 0) {
                profile->r_bins[phi_bin][r_bin] /= bin_quantity;
                profile->g_bins[phi_bin][r_bin] /= bin_quantity;
                profile->b_bins[phi_bin][r_bin] /= bin_quantity;
            }
            else {
                profile->r_bins[phi_bin][r_bin] = 0;
                profile->g_bins[phi_bin][r_bin] = 0;
                profile->b_bins[phi_bin][r_bin] = 0;
            }
        }
    }
    END_TIMING(average_bin_time, "calculating average of blur_profile bins");
    return profile;
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

    #ifdef DEBUG
    // Create image representation of bin-approximated FFT
    Image_RGB* blur_profile_visualization = get_blur_profile_visual(blur_profile, fft->height, fft->width);

    // Save image representation of bin-approximated FFT
    const char* path = "images/visualizations/";
    const char* q = "\nEnter visualization filename: ";
    char* visualization_filename = create_path(path, q, ".txt");
    write_image_to_file(blur_profile_visualization, visualization_filename);
    free(visualization_filename), visualization_filename = NULL;
    free_image_rgb(blur_profile_visualization);
    #endif

    // Free all function data
    free_cartesian_to_polar(conversion);
}


/******************************************************************************
 * get_blur_profile is the top-level function for production implementation for
 *   the blur_profile.c file. It takes in image and returns the blur profile
 *   object.
 *  -image is an RGB image of which the blur profile is desired.
 *  -num_radius_bins and num_angle_bins are hyperparameters, together deciding
 *   the size of each polar-coordinate FFT bin, and therefore the granularity
 *   of the FFT measurement.
******************************************************************************/
Blur_Profile_RGB* get_blur_profile(Image_RGB* image, int num_radius_bins, int num_angle_bins) {
    // Calculate FFT
    START_TIMING(compute_mag_fft_time);
    Image_RGB* fft = compute_magnitude_fft(image);
    END_TIMING(compute_mag_fft_time, "computing the magnitude fft");

    // Create cartesian to polar conversion
    START_TIMING(c2p_time);
    Cartesian_To_Polar* conversion;
    conversion = cartesian_to_polar_conversion(fft->width, fft->height);
    END_TIMING(c2p_time, "cartesian to polar conversion");

    // Calculate the blur profile
    START_TIMING(calc_blur_profile_time);
    Blur_Profile_RGB* blur_profile;
    blur_profile = calculate_blur_profile(conversion, fft, num_radius_bins, num_angle_bins);
    END_TIMING(calc_blur_profile_time, "calculate_blur_profile");
    
    #ifdef DEBUG
    // Create image representation of bin-approximated FFT
    Image_RGB* blur_profile_visualization = get_blur_profile_visual(blur_profile, fft->height, fft->width);

    // Save image representation of bin-approximated FFT
    const char* path = "images/visualizations/";
    const char* q = "\nEnter visualization filename: ";
    char* visualization_filename = create_path(path, q, ".txt");
    write_image_to_file(blur_profile_visualization, visualization_filename);
    free(visualization_filename), visualization_filename = NULL;
    free_image_rgb(blur_profile_visualization);
    #endif

    // Clean up
    START_TIMING(free_fft_time);
    free_image_rgb(fft);
    free_cartesian_to_polar(conversion);
    END_TIMING(free_fft_time, "freeing FFT structures");
    return blur_profile;
}


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
    int height_bound;
    if (height %2 == 1) height_bound = half_height+1;
    else height_bound = half_height;
    int y_times_width = 0;
    for (int y=0; y<height_bound; y++) {    // height/2 for symmetry across x axis
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
 * free_cartesian_to_polar frees the memory of Cartesian_To_Polar objects.
******************************************************************************/
void free_cartesian_to_polar(Cartesian_To_Polar* c2p) {
    if (c2p) {
        free(c2p->data);
        c2p->data = NULL;
        free(c2p);
    }
}


void free_blur_profile_rgb(Blur_Profile_RGB* bp) {
    free(bp->b_bins);
    free(bp->g_bins);
    free(bp->r_bins);
    free(bp);
    bp = NULL;
}