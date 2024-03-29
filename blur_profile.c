#include "image_processing.h"
#include "blur_profile.h"
#include "utilities.h"
#include "fft_processing.h"
#include "filtering.h"
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
Blur_Profile* calculate_blur_profile(
                    const Cartesian_To_Polar* conversion, 
                    const Image_PGM* fft, 
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
    Blur_Profile* profile = (Blur_Profile*)malloc(sizeof(Blur_Profile));
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
    profile->bins = (Bin**)malloc(num_angle_bins * sizeof(Bin*));
    if (!profile->bins) {
        fprintf(stderr, "ERROR: Memory allocation failed for r,g, or b bins.\n");
        return NULL;
    }
    for (int i=0; i<num_angle_bins; i++) {
        profile->bins[i] = (Bin*)calloc(num_radius_bins, sizeof(Bin));
    }

    // Sum up blurs in each bin space
    int conversion_height = conversion->height; // No need to perform member access every iteration
    int conversion_width = conversion->width;
    int i_tot = fft->height * fft->width;       // Store once to speed up by 1 multiply & 2 member accesses
    const Pixel* data = fft->data;    // locally store array pointers to assure high performance in dereferencing
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
        if (r_bin == num_radius_bins) r_bin--;
        bin_quant[phi_bin][r_bin]++;    // Increment the counter for how many pixels have been assigned to this bin
        profile->bins[phi_bin][r_bin] += data[i];    // Profile not locally stored like fft members because
    }
    END_TIMING(accumulate_bins_time, "accumulating pixels to blur_profile bins");

    // Divide by the number of bins
    START_TIMING(average_bin_time);
    Bin bin_quantity;
    for (int phi_bin=0; phi_bin<num_angle_bins; phi_bin++) {
        for (int r_bin=0; r_bin<num_radius_bins; r_bin++) {
            bin_quantity = (Bin)bin_quant[phi_bin][r_bin];
            if (bin_quantity != 0) {
                profile->bins[phi_bin][r_bin] /= bin_quantity;
            }
            else {
                profile->bins[phi_bin][r_bin] = 0;
            }
        }
    }    
    END_TIMING(average_bin_time, "calculating average of blur_profile bins");

    // Clean up
    for (int i=0; i < num_angle_bins; i++) {
        free(bin_quant[i]); bin_quant[i] = NULL;
    }
    free(bin_quant); bin_quant = NULL;

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
Image_PGM* get_blur_profile_visual(Blur_Profile* blur_profile, 
                                   int height, int width) {
    // Allocate memory for an RGB image, with height and width
    Image_PGM* output_image = create_pgm_image(width, height);

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

            // Direct assignment without unnecessary mirroring
            output_image->data[y_times_width + x] = blur_profile->bins[phi_bin][r_bin];
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
void blur_profile_tests(Image_PGM* fft) {
    // Calculate bin-approximated FFT
    Cartesian_To_Polar* conversion;   // Create cartesian to polar conversion
    conversion = cartesian_to_polar_conversion(fft->width, fft->height);
    int num_radius_bins = 4;    // Set up pseudo-hyperparameters
    int num_angle_bins = 18;
    Blur_Profile* blur_profile; // Calculate the blur profile
    blur_profile = calculate_blur_profile(conversion, fft, num_radius_bins, num_angle_bins);

    #ifdef DEBUG
    // Create image representation of bin-approximated FFT
    Image_PGM* blur_profile_visualization = get_blur_profile_visual(blur_profile, fft->height, fft->width);
    Image_RGB* vis = pgm2rgb(blur_profile_visualization);
    // Save image representation of bin-approximated FFT
    const char* path = "images/visualizations/";
    const char* q = "\nEnter visualization filename: ";
    char* visualization_filename = create_path(path, q, ".txt");
    write_image_to_file(vis, visualization_filename);
    free(visualization_filename), visualization_filename = NULL;
    free_image_rgb(vis); vis = NULL;
    free_image_pgm(blur_profile_visualization); blur_profile_visualization = NULL;
    #endif

    // Free all function data
    free_cartesian_to_polar(conversion); conversion = NULL;
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
Blur_Profile* get_blur_profile(Image_PGM* pgm, int num_radius_bins, int num_angle_bins) {
    // Calculate FFT
    START_TIMING(compute_mag_fft_time);
    Image_PGM* fft = compute_magnitude_fft(pgm);
    END_TIMING(compute_mag_fft_time, "computing the magnitude fft");

    // Create cartesian to polar conversion
    START_TIMING(c2p_time);
    Cartesian_To_Polar* conversion;
    conversion = cartesian_to_polar_conversion(fft->width, fft->height);
    END_TIMING(c2p_time, "cartesian to polar conversion");

    // Calculate the blur profile
    START_TIMING(calc_blur_profile_time);
    Blur_Profile* blur_profile;
    blur_profile = calculate_blur_profile(conversion, fft, num_radius_bins, num_angle_bins);
    END_TIMING(calc_blur_profile_time, "calculate_blur_profile");
    
    #ifdef DEBUG
    // Create image representation of bin-approximated FFT
    Image_PGM* blur_profile_visualization = get_blur_profile_visual(blur_profile, fft->height, fft->width);
    Image_RGB* vis = pgm2rgb(blur_profile_visualization);

    // Save image representation of bin-approximated FFT
    const char* path = "images/visualizations/";
    const char* q = "\nEnter visualization filename: ";
    char* visualization_filename = create_path(path, q, ".txt");
    write_image_to_file(vis, visualization_filename);
    free(visualization_filename), visualization_filename = NULL;
    free_image_rgb(vis); vis = NULL;
    free_image_pgm(blur_profile_visualization); blur_profile_visualization = NULL;
    #endif

    // Clean up
    START_TIMING(free_fft_time);
    free_image_pgm(fft); fft = NULL;
    free_cartesian_to_polar(conversion); conversion = NULL;
    END_TIMING(free_fft_time, "freeing FFT structures");
    return blur_profile;
}


// Create an instance of Blur_Vector_RGB
Blur_Vector_Group* initialize_blur_vector_group(int len_vectors) {
    Blur_Vector_Group* blur_vector_group = (Blur_Vector_Group*)calloc(1, sizeof(Blur_Vector_Group));
    blur_vector_group->blur_vectors = (Blur_Vector*)calloc(len_vectors, sizeof(Blur_Vector));
    blur_vector_group->len_vectors = len_vectors;
    return blur_vector_group;
}


// Free Blur_Vector_RGB object
void free_blur_vectors_rgb(Blur_Vector_Group* bv) {
    free(bv->blur_vectors);
    bv->blur_vectors = NULL;
    free(bv);
}


/******************************************************************************
 * vectorize_blur_profile returns an array Blur_Vector structs, representing
 *   the angle and magnitude of blur in an image.
 *  -blur_profile contains the polar-form quantized FFT of an image
 *  -error_thresh is the value that represents how much higher an FFT streak
 *   must be in order to count as an effect of reasonable blurring
 *  -mag_thresh is the magnitude threshold (in the range [0,1]) in the
 *   Blur_Profile spectrum that marks what is considered the magnitude of blur.
 *  -cutoff_ration_denom represents the proportion (1/cutoff) of the total
 *   radius bins used to calculate the average that to vets possible maxima.
******************************************************************************/
Blur_Vector_Group* vectorize_blur_profile(Blur_Profile* blur_profile,
                                        Pixel error_thresh,
                                        Pixel mag_thresh,
                                        int cutoff_ratio_denom) {
    Blur_Vector_Group* bv = initialize_blur_vector_group(10);
    Blur_Vector* blur_vectors = bv->blur_vectors;

    // Store common values for easy use
    int num_angle_bins = blur_profile->num_angle_bins;
    int num_radius_bins = blur_profile->num_radius_bins;

    Bin** bins = blur_profile->bins;

    // Assign the current color of bins to use
    Pixel* tot = (Pixel*)calloc(num_angle_bins, sizeof(Pixel));

    // Get average of bins, and total of first half of the radii, per angle
    Pixel avg = 0;
    int radius_cutoff = num_radius_bins/cutoff_ratio_denom;
    for (int i=0; i<num_angle_bins; i++) {
        for (int j=0; j<radius_cutoff; j++) {
            tot[i] += bins[i][j];
        }
        avg += tot[i];
    }
    avg /= num_angle_bins;

    // Smooth the data to reduce noise
    int smooth_size = 5;
    Pixel* smoother = initialize_1d_smooth_filter(smooth_size);
    Pixel* smooth = convolve_1d(tot, smoother, num_angle_bins, smooth_size);
    
    // Find maxima in the smoothed data
    Blur_Vector maxima[10];
    int maxima_idx = 0;
    if (smooth[0] > smooth[num_angle_bins-1] && smooth[0] > smooth[1]) {
        if (smooth[0] > avg * error_thresh && maxima_idx < 10) {
            maxima[maxima_idx].angle = 0;
            maxima[maxima_idx++].magnitude = tot[0];
        }
    }
    for (int i=1; i<num_angle_bins-1; i++) {
        if (smooth[i] > smooth[i-1] && smooth[i] > smooth[i+1]) {
            if (smooth[i] > avg * error_thresh && maxima_idx < 10) {
                maxima[maxima_idx].angle = i;
                maxima[maxima_idx++].magnitude = tot[i]/radius_cutoff;
            }
        }
    }
    if (smooth[num_angle_bins-1] > smooth[num_angle_bins-2] && smooth[num_angle_bins-1] > smooth[0]) {
        if (smooth[num_angle_bins-1] > avg * error_thresh && maxima_idx < 10) {
            maxima[maxima_idx].angle = num_angle_bins-1;
            maxima[maxima_idx++].magnitude = tot[num_angle_bins-1]/radius_cutoff;
        }
    }
    free(tot);

    // for max_angle in maxima:
    for (int i=0; i<maxima_idx; i++) {
        // Adjust angle_index to be within the valid range
        int angle_idx = (maxima[i].angle + num_angle_bins/2) % num_angle_bins;
        // Get the signal for the current angle
        Pixel* cur_sig = bins[angle_idx];

        // Find the index where the signal falls below the threshold
        int cur_max_radius = num_radius_bins;
        for (int j=0; j<num_radius_bins; j++) {
            if (cur_sig[j] < mag_thresh) {
                cur_max_radius = j;
                break;
            }
        }

        // Calculate and store relative magnitude and angle
        blur_vectors[i].magnitude = ((float)cur_max_radius / (float)num_radius_bins);
        blur_vectors[i].angle = (int)(180 * ((float)angle_idx / (float)num_angle_bins) - 90);
    }
    return bv;
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


void free_blur_profile(Blur_Profile* bp) {
    free_2d_array(bp->bins, bp->num_angle_bins);
    bp->bins = NULL;
    free(bp);
    bp = NULL;
}


void free_blur_vector_group(Blur_Vector_Group* bv) {
    free(bv->blur_vectors); bv->blur_vectors = NULL;
    free(bv);
}