#include "../image_processing.h"
#include "../fft_processing.h"
#include "../color_quantization.h"
#include "../blur_profile.h"
#include "../debug.h"
#include "../filtering.h"
#include "../types.h"
#include "../utilities.h"
#include "../interface.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

// Test configuration parameters
#define RUN_FAILING_TESTS
#define RUN_PASSING_TESTS
#define MASSIVE_WIDTH 120000
#define MASSIVE_HEIGHT 10000
#define WIDE_WIDTH 2001
#define WIDE_HEIGHT 400
#define TALL_WIDTH 400
#define TALL_HEIGHT 2001
#define SMALL_WIDTH 349
#define SMALL_HEIGHT 350

// Default parameters for get_full_report_data
#define DEFAULT_H_PARTITIONS 18
#define DEFAULT_S_PARTITIONS 2
#define DEFAULT_V_PARTITIONS 3
#define DEFAULT_BLACK_THRESH 0.1
#define DEFAULT_GRAY_THRESH 0.1
#define DEFAULT_COVERAGE_THRESH 0.95
#define DEFAULT_LINKED_LIST_SIZE 1000
#define DEFAULT_DOWNSAMPLE_RATE 1
#define DEFAULT_RADIUS_PARTITIONS 40
#define DEFAULT_ANGLE_PARTITIONS 72
#define DEFAULT_QUANTITY_WEIGHT 0.1
#define DEFAULT_SATURATION_VALUE_WEIGHT 0.9
#define DEFAULT_FFT_STREAK_THRESH 1.20
#define DEFAULT_MAGNITUDE_THRESH 0.3
#define DEFAULT_BLUR_CUTOFF_RATIO_DENOM 2


#define OUTPUT_BUFFER_SIZE 10000

char output_buffer[OUTPUT_BUFFER_SIZE];

void append_to_output(const char* format, ...) {
    char buffer[256];
    va_list args;

    // Initialize the va_list variable with the last known fixed parameter
    va_start(args, format);
    // Use vsnprintf for safe formatting with variable arguments
    vsnprintf(buffer, sizeof(buffer) - 1, format, args);
    va_end(args); // Clean up the va_list variable

    // Append newline character to the buffer
    size_t buf_len = strlen(buffer);
    if (buf_len < sizeof(buffer) - 1) {
        buffer[buf_len] = '\n'; // Add newline character
        buffer[buf_len + 1] = '\0'; // Ensure string is null-terminated
    }

    // Concatenate the formatted string to the main output buffer
    if (strlen(output_buffer) + strlen(buffer) < OUTPUT_BUFFER_SIZE) {
        strcat(output_buffer, buffer);
    } else {
        printf("Output buffer full. Truncating or handling overflow...\n");
    }
}


Full_Report_Data* run_full_report_data(Image_RGB* image, Crop_Boundaries* crop_boundaries) {
    return get_full_report_data(image, crop_boundaries,
                                DEFAULT_H_PARTITIONS, DEFAULT_S_PARTITIONS, DEFAULT_V_PARTITIONS, 
                                DEFAULT_BLACK_THRESH, DEFAULT_GRAY_THRESH, DEFAULT_COVERAGE_THRESH, 
                                DEFAULT_LINKED_LIST_SIZE, DEFAULT_DOWNSAMPLE_RATE, 
                                DEFAULT_RADIUS_PARTITIONS, DEFAULT_ANGLE_PARTITIONS, 
                                DEFAULT_QUANTITY_WEIGHT, DEFAULT_SATURATION_VALUE_WEIGHT,
                                DEFAULT_FFT_STREAK_THRESH, DEFAULT_MAGNITUDE_THRESH, DEFAULT_BLUR_CUTOFF_RATIO_DENOM);
}


bool test_large_file_size() {
    Image_RGB* image = create_rgb_image(MASSIVE_WIDTH, MASSIVE_HEIGHT);
    Crop_Boundaries* crop_boundaries = NULL;  // Or set appropriately if needed

    Full_Report_Data* report = run_full_report_data(image, crop_boundaries);
    if (report == NULL) {
        append_to_output("Test Passed: System handled large image.");
    } else {
        append_to_output("Test Failed: System could not handle large image.");
    }
    if (image) free_image_rgb(image);
    if (report != NULL) free_full_report(&report);
}


bool test_unusual_aspect_ratios() {
    Image_RGB* wide_image = create_rgb_image(WIDE_WIDTH, WIDE_HEIGHT);
    Image_RGB* tall_image = create_rgb_image(TALL_WIDTH, TALL_HEIGHT);
    Crop_Boundaries* crop_boundaries = NULL;  // Or set appropriately if needed

    Full_Report_Data* wide_report = run_full_report_data(wide_image, crop_boundaries);
    Full_Report_Data* tall_report = run_full_report_data(tall_image, crop_boundaries);

    if (wide_report == NULL && tall_report == NULL) {
        append_to_output("Test Passed: System handled images with unusual aspect ratios.");
    } else {
        append_to_output("Test Failed: System could not handle images with unusual aspect ratios.");
    }

    if (wide_image) (wide_image);
    if (tall_image) free_image_rgb(tall_image);
    if (wide_report != NULL) free_full_report(&wide_report);
    if (tall_report != NULL) free_full_report(&tall_report);
}


bool test_minimum_size_constraint() {
    Image_RGB* small_image = create_rgb_image(SMALL_WIDTH, SMALL_HEIGHT);
    Crop_Boundaries* crop_boundaries = NULL;  // Or set appropriately if needed

    Full_Report_Data* report = run_full_report_data(small_image, crop_boundaries);
    if (report == NULL) {
        append_to_output("Test Passed: System rejected image below minimum size constraint.");
    } else {
        append_to_output("Test Failed: System did not reject image below minimum size constraint.");
    }
    if (small_image) free_image_rgb(small_image);
    if (report) free_full_report(&report);
}


bool run_time_test() {
    // Get image from files
    Image_RGB* image = read_image_from_files();

    Crop_Boundaries* crop_boundaries = (Crop_Boundaries*)malloc(sizeof(Crop_Boundaries));
    crop_boundaries->N = 0;

    struct timeval start = start_timing();
    Full_Report_Data* full_report_data = run_full_report_data(image, crop_boundaries);
    double elapsed_time = end_timing(start);

    if (full_report_data == NULL) {
        append_to_output("Test Failed: system could not finish report on a proper image.");
    }
    else if (elapsed_time > 0.5) {
        append_to_output("Test Failed: full report data took %lf seconds to compute.", elapsed_time);
    }
    else {
        append_to_output("Test Passed: full report data took under 0.5 seconds (took %lf seconds).", elapsed_time);
        free_full_report(&full_report_data);
    }

    // Clean up
    if (image) free_image_rgb(image);
    free_crop_boundaries(crop_boundaries);

    return 0;
}


int main() {
    memset(output_buffer, 0, OUTPUT_BUFFER_SIZE);
    printf("Hello there!\n");
    printf("\nGeneral Kenobi....\n");
    printf("Running test suite...\n");
    
    #ifdef RUN_FAILING_TESTS
    test_large_file_size();
    test_unusual_aspect_ratios();
    test_minimum_size_constraint();
    #endif

    #ifdef RUN_PASSING_TESTS
    run_time_test();
    #endif

    printf("\nTest suite completed. Test Results:\n%s\n", output_buffer);

    return 0;
}
