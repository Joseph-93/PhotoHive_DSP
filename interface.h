#include "image_processing.h"
#include "utilities.h"

#define LINKED_LIST_SIZE 5000
#define COVERAGE_THRESHOLD 0.90
#define DOWNSAMPLE_RATE 5
#define H_PARTS 10
#define S_PARTS 5
#define V_PARTS 5
#define BLACK_THRESH 0.15
#define GRAY_THRESH 0.1
#define NUM_RADIUS_BINS 4
#define NUM_ANGLE_BINS 18


Full_Report_Data* get_full_report_data(Image_RGB* image,
                                       int h_partitions, int s_partitions, int v_partitions,
                                       double black_thresh, double gray_thresh,
                                       double coverage_thresh, int linked_list_size,
                                       int downsample_rate, int radius_partitions, int angle_partitions,
                                       float quantity_weight, float saturation_value_weight
                                       );


void free_full_report(Full_Report_Data** report);