#include "../image_processing.h"
#include "../fft_processing.h"
#include "../color_quantization.h"
#include "../blur_profile.h"
#include "../debug.h"
#include "../filtering.h"
#include "../types.h"
#include "../utilities.h"
#include "../interface.h"


int main() {
    printf("hello there!");

    // Get image from files
    Image_RGB* image = read_image_from_files();

    float coverage_thresh=.4;
    int downsample_rate=1;
    int radius_partitions=40;
    int angle_partitions=72;
    float quantity_weight=0.9;
    float saturation_value_weight=0.1;
    int h_partitions=18;
    int s_partitions=2;
    int v_partitions=3;
    float black_thresh=0.1;
    float gray_thresh=0.1;
    int linked_list_size=1000;
    double fft_streak_thresh=1.1;
    double magnitude_thresh=0.3;
    int blur_cutoff_ratio_denom=2;

    // Arm full report
    Full_Report_Data* full_report_data = get_full_report_data(image,
                                         h_partitions, s_partitions, v_partitions, 
                                         black_thresh, gray_thresh, coverage_thresh, 
                                         linked_list_size, downsample_rate, 
                                         radius_partitions, angle_partitions, 
                                         quantity_weight, saturation_value_weight,
                                         fft_streak_thresh, magnitude_thresh, blur_cutoff_ratio_denom);


    // Clean up
    free_image_rgb(image);
    free_full_report(&full_report_data);

    return 0;
}

