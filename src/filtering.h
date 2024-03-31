#include "image_processing.h"
#include "utilities.h"

#ifndef FILTERING_H
#define FILTERING_H

typedef struct Filter {
    int height;
    int width;
    double* coefs;
} Filter;


Pixel* convolve_1d(Pixel* x, Pixel* h, int x_len, int y_len);

Pixel* initialize_1d_smooth_filter(int size);

Filter* initialize_3x3_laplacian();

Pixel get_variance(Pixel* input, int length, Pixel average);

Pixel get_average(Pixel* input, int length);

Pixel* filter_image(Filter* filt, Pixel* input, int height, int width);

double sharpness_avg(double* input, int length);

Image_RGB* create_filtered_RGB(Image_RGB* input, Filter* filt);

Sharpnesses* get_variance_sharpness(Image_PGM* image, Crop_Boundaries* bounds);

Pixel get_average_sharpness(Pixel* input, int height, int width);

#endif
