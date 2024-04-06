#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/time.h>

#include "image_processing.h"
#include "fft_processing.h"
#include "color_quantization.h"
#include "blur_profile.h"
#include "debug.h"
#include "filtering.h"
#include "types.h"
#include "utilities.h"


Full_Report_Data* get_full_report_data(Image_RGB* image, Crop_Boundaries* crop_boundaries,
                  int h_partitions, int s_partitions, int v_partitions,
                  double black_thresh, double gray_thresh,
                  double coverage_thresh, int linked_list_size,
                  int downsample_rate, int radius_partitions, int angle_partitions,
                  float quantity_weight, float saturation_value_weight,
                  double fft_streak_thresh, double magnitude_thresh, int blur_cutoff_ratio_denom
                  ) {
    if (pre_compute_error_checks(image)) return NULL;

    #ifdef DEBUG
    setbuf(stdout, NULL);
    #endif

    threading_setup();
    printf("\n There are %d cores available to the C program.\n\n", num_cores);

    // Downsample RGB image
    START_TIMING(downsample_time);
    Image_RGB* downsampled;
    if(downsample_rate > 1) downsampled = downsample_rgb(image, downsample_rate);
    else downsampled = image;
    END_TIMING(downsample_time, "downsample rgb");

    // Convert to HSV and PGM
    START_TIMING(hsv_time);
    Image_HSV* hsv = rgb2hsv(downsampled);
    END_TIMING(hsv_time, "rgb2hsv");

    START_TIMING(pgm_time);
    Image_PGM* pgm = rgb2pgm(image);
    END_TIMING(pgm_time, "rgb2pgm");

    // Get brightness and contrast measurements
    START_TIMING(rgb_stats_time);
    RGB_Statistics* rgb_stats = get_rgb_statistics(image);
    END_TIMING(rgb_stats_time, "rgb statistics");

    // Get average saturation from hsv
    START_TIMING(s_bar_time);
    Pixel S_bar = get_hsv_average(hsv);
    END_TIMING(s_bar_time, "hsv average");

    // Get color palette from hsv
    START_TIMING(color_palette_time);
    Color_Palette* cp = get_color_palette(hsv, linked_list_size, coverage_thresh,
                                          h_partitions, s_partitions, v_partitions,
                                          black_thresh, gray_thresh,
                                          quantity_weight, saturation_value_weight);
    END_TIMING(color_palette_time, "color palette");

    // get the variance sharpness measurement
    START_TIMING(sharpness_time);
    Sharpnesses* variance_sharpness = get_variance_sharpness(pgm, crop_boundaries);
    END_TIMING(sharpness_time, "sharpness");

    // Get blur profile
    START_TIMING(blur_profile_time);
    double avg = (rgb_stats->Br + rgb_stats->Bg + rgb_stats->Bb)/3.0;
    Blur_Profile* bp = get_blur_profile(pgm, radius_partitions, angle_partitions, avg);
    Blur_Vector_Group* bv = vectorize_blur_profile(bp, fft_streak_thresh, magnitude_thresh, blur_cutoff_ratio_denom);
    END_TIMING(blur_profile_time, "blur profile");

    START_TIMING(compile_report_time);
    Full_Report_Data* full_report = compile_full_report(rgb_stats, cp, bp, bv, S_bar, variance_sharpness);
    END_TIMING(compile_report_time, "compile full report");

    // Free data
    START_TIMING(free_image_time);
    free_image_pgm(pgm); pgm = NULL;
    free_image_hsv(hsv); hsv = NULL;
    if (downsample_rate > 1) free_image_rgb(downsampled); downsampled = NULL;
    END_TIMING(free_image_time, "free images");
    return full_report;
}


void free_full_report(Full_Report_Data** report) {
    free_color_palette((*report)->color_palette); (*report)->color_palette = NULL;
    if ((*report)->blur_profile != NULL) {
        free_blur_profile((*report)->blur_profile); (*report)->blur_profile = NULL;
    }
    free_blur_vector_group((*report)->blur_vectors); (*report)->blur_vectors = NULL;
    if ((*report)->sharpness != NULL) {
        free((*report)->sharpness->sharpness); (*report)->sharpness->sharpness = NULL;
        free((*report)->sharpness); (*report)->sharpness = NULL;
    }
    free((*report)->rgb_stats); (*report)->rgb_stats = NULL;
    
    free((*report)); *report = NULL;
    report = NULL;
}
