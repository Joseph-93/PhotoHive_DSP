#include "types.h"
#include "color_quantization.h"
#include "stdio.h"
#include "blur_profile.h"

#ifndef UTILITIES_H
#define UTILITIES_H

extern int num_cores;

typedef struct Full_Report_Data {
    RGB_Statistics* rgb_stats;
    Color_Palette* color_palette;
    Blur_Profile_RGB* blur_profile;
    Pixel average_saturation;
    Pixel sharpness;
} Full_Report_Data;

int newton_int_sqrt(double val);

char* create_path(const char* path, const char* request_string, const char* filetype);

void threading_setup();

void custom_sort(void *base, size_t nmemb, size_t size,
                 int (*compar)(const void *, const void *, void *),
                 void *arg);

bool ask_yes_no_question(const char *question);

void normalize_array(Pixel* array, int length);

Pixel get_variance_sharpness(Pixel* input, int height, int width);

Full_Report_Data* compile_full_report(RGB_Statistics* rgb_stats,
                                  Color_Palette* color_palette,
                                  Blur_Profile_RGB* blur_profile,
                                  Pixel average_saturation,
                                  Pixel sharpness);

void print_full_report(Full_Report_Data* data);

#endif