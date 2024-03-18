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

#define LINKED_LIST_SIZE 1000
#define COVERAGE_THRESHOLD 0.90
#define DOWNSAMPLE_RATE 5
#define H_PARTS 6
#define S_PARTS 7
#define V_PARTS 7
#define BLACK_THRESH 0.1
#define GRAY_THRESH 0.1
#define NUM_RADIUS_BINS 4
#define NUM_ANGLE_BINS 18


Full_Report_Data* get_full_report_data(Image_RGB* image) {
    setbuf(stdout, NULL);
    // Convert to HSV and PGM
    Image_HSV* hsv = rgb2hsv(image);
    printf("Got past hsv conversion.\n");
    Image_PGM* pgm = rgb2pgm(image);
    printf("got past the pgm conversion.\n");

    // Get brightness and contrast measurements
    RGB_Statistics* rgb_stats = get_rgb_statistics(image);
    printf("got past rgb_statistics.\n");

    // Get average saturation from hsv
    Pixel S_bar = get_hsv_average(hsv);
    printf("got past S_bar\n");

    // Get color palette from hsv
    Color_Palette* cp = get_color_palette(hsv, LINKED_LIST_SIZE, COVERAGE_THRESHOLD,
                                          H_PARTS, S_PARTS, V_PARTS,
                                          BLACK_THRESH, GRAY_THRESH);
    printf("got past color_palette.\n");    

    // get the variance sharpness measurement
    Pixel variance_sharpness = get_variance_sharpness(pgm->data, pgm->height, pgm->width);
    printf("got past variance_sharpness.\n");

    // Get blur profile
    Blur_Profile_RGB* bp = get_blur_profile(image, NUM_RADIUS_BINS, NUM_ANGLE_BINS);
    printf("got past Blur Profile.\n");

    Full_Report_Data* full_report = compile_full_report(rgb_stats, cp, bp, S_bar, variance_sharpness);
    printf("Got past full_data_report*.\n");

    // Free data
    free_image_pgm(pgm); pgm = NULL;
    printf("got past freeing the pgm.\n");
    free_image_hsv(hsv); hsv = NULL;
    printf("got past freeing hsv.\n");
    return full_report;
}


void free_full_report(Full_Report_Data** report) {
    free_color_palette((*report)->color_palette); (*report)->color_palette = NULL;
    free_blur_profile_rgb((*report)->blur_profile); (*report)->blur_profile = NULL;
    free((*report)); *report = NULL;
    report = NULL;
}


int main() {

    #ifdef DEBUG
    setbuf(stdout, NULL);
    #endif

    // Get image from files
    Image_RGB* image = read_image_from_files();

    // Arm full report
    Full_Report_Data* full_report_data = get_full_report_data(image);

    print_full_report(full_report_data);

    return 0;
}
