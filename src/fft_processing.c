// Contains: FFT related functions and structures like rgb_fft, compute_magnitude_fft, fft_shift, and the inclusion of the FFTW3 library.

#include <fftw3.h>
#include "fft_processing.h"
#include "image_processing.h"
#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

/******************************************************************************
 * rgb_fft calculates the FFT of an RGB image given.
 *  -If no fftw_plan* is given, then the function will automatically generate
 *   an fftw_plan* to be used, based on the inputted image height.
******************************************************************************/
Image_PGM* pgm_fft(Image_PGM* pgm) {
    // Allow threading
    fftw_init_threads();
    fftw_plan_with_nthreads(num_cores);

    // Create output and input
    fftw_complex* output = fftw_alloc_complex((pgm->width/2+1) * pgm->height);
    double* in = fftw_alloc_real(pgm->height * pgm->width);

    if (!output || !in) {
        fprintf(stderr, "ERROR: Memory allocation for fftw_malloc failed.\n");
        fftw_free(output);
        fftw_free(in);
        return NULL;
    }

    fftw_plan plan = fftw_plan_dft_r2c_2d(pgm->height, pgm->width, in, output, FFTW_ESTIMATE);
    if (!plan) {
        fprintf(stderr, "ERROR: fftw_plan generation failed.");
        fftw_destroy_plan(plan);
        return NULL;
    }

    int fft_width = pgm->width/2+1;
    Image_PGM* fft_image = create_pgm_image(fft_width, pgm->height);
    int i_tot = fft_image->width * fft_image->height;

    // Execute fftw for red, collect magnitudes^2 in fft_image
    memcpy(in, pgm->data, sizeof(double) * pgm->width * pgm->height);
    fftw_execute(plan);
    for (int i=0; i<i_tot; i++) {
        fft_image->data[i] = output[i][0]*output[i][0] + output[i][1]*output[i][1];
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
 * compute_magnitude_fft returns an Image_RGB* to a normalized fft, given an
 *   rgb image.
******************************************************************************/
Image_PGM* compute_magnitude_fft(Image_PGM* pgm) {
    Image_PGM* fft = pgm_fft(pgm);
    pgm_normalize_fft(fft, NULL);
    return fft;
}

/******************************************************************************
 * save_fft writes an fft image to user's desired image file (in .txt format 
 *   for python to save via imageSaver.py).
 *  -Function is designed for development, as it prompts users for command-line
 *   input.
******************************************************************************/
void save_fft(Image_PGM* fft) {
    Image_RGB* rgb = pgm2rgb(fft);
    const char* q = "Enter the name of the .txt file you wish to write to: ";
    char* output_filename = create_path("images/output/", q, ".txt");
    write_image_to_file(rgb, output_filename);
    free(output_filename), output_filename = NULL;
    free_image_rgb(rgb);
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
Image_PGM* fft_shift(Image_PGM* fft) {
    Image_PGM* fft_image = create_pgm_image(fft->width*2 - 1, fft->height);
    // Iterate through one quarter of needed image, via one half of the image
    printf("new_image size: (%3d, %3d)\t old_image size: (%3d, %3d)\n", fft_image->width, fft_image->height, fft->width, fft->height);
    
    // If even image, height_compensator is 0. if odd image, compensator is 1
    int height_compensator = (int)(fft_image->height%2==1);
    int width_compensator = (int)(fft_image->width%2==1); // width comp is always 1, but I want to make it dynamic.
    printf("height_equalizer: %d\t width_equalizer: %d\n", height_compensator, width_compensator);

    int errors = 0;
    Pixel temp;
    int half_fft_height = fft->height/2;
    int y_val, x_val;
    int half_height_plus_compensator = half_fft_height+height_compensator;
    int fft_width = fft->width;
    for (int y=0; y<half_height_plus_compensator; y++) {
        for (int x=0; x<fft_width; x++) {
            /*** Swap Q1 and Q4 ***/
            y_val = y+half_fft_height, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // New Q1 gets old Q4
            temp = fft->data[y_val*fft_width + x_val];
            if (y-height_compensator != -1) { // If height is odd, middle row goes in a theoretical pixels[-1]
                y_val = y-height_compensator, x_val = x+fft->width-width_compensator, validate_coordinates(fft_image, x_val, y_val, &errors);
                fft_image->data[y_val*fft_width + x_val] = temp;
            }
            y_val = y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // New Q4 gets old Q1
            temp = fft->data[y_val * fft_width + x_val];
            y_val = y+half_fft_height, x_val = x+fft->width-width_compensator, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->data[y_val * fft_width + x_val] = temp;

            /*** Rotate to get Q2 and Q3 ***/
            y_val = y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // Q2 is 180 deg rotation of Q4
            temp = fft->data[y_val * fft_width + x_val];
            y_val = half_fft_height-y, x_val = fft->width-1-x, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->data[y_val * fft_width + x_val] = temp; //width-1 because width-1=last_index
            y_val = half_fft_height+y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            // Q3 is 180 deg rotation of Q1
            temp = fft->data[y_val * fft_width + x_val];
            y_val = fft->height-1-y, x_val = fft->width-1-x, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->data[y_val * fft_width + x_val] = temp; // height-1 because height-1=last_index
        }
    }
    return fft_image;
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
void pgm_normalize_fft(Image_PGM* fft, Lookup_1D* fft_normalizer_lookup) {
    int guess_max_location = fft->height*fft->width/2;
    double max = fft->data[guess_max_location];     // Set Maximum as any valid value
    int i_tot = fft->height * fft->width;   // Get stopping point
    int max_index = 0;

    // Find maximum input
    START_TIMING(find_max_time);
    for (int i=0; i<i_tot; i++) {
        double p = fft->data[i]; // Set temporary values for efficiency
        if (max < p) {max = p; max_index = i;}
    }
    END_TIMING(find_max_time, "finding the FFT max value");

    #ifdef VERBOSE
        printf("Maximum value found in the raw FFT: %f\n", max);
    #endif
    
    // Set G_s
    Pixel G_s = 1/(2*log(sqrt(max) + 1));
    
    START_TIMING(normalize_fft_time);
    // calculate normalized values
    for (int i=0; i<i_tot; i++) {
        if (fft->data[i] < 1) fft->data[i] = 0;
        else fft->data[i] = log(fft->data[i]) * G_s;
    }
    END_TIMING(normalize_fft_time, "normalizing FFT values");

    #ifdef VERBOSE
        // Set Maximum as any valid value
        max = fft->data[0];

        // Set maximum values
        for (int i=0; i<i_tot; i++) {
            double p = fft->data[i]; // Set temporary values for efficiency
            if (max < p) {max = p;}
        }
        printf("Maximum value found in the normalized FFT: %f\n", max);
    #endif
}

