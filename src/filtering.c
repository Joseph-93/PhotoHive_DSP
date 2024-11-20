#include "filtering.h"
#include <stdlib.h>
#include <math.h>


#define THRESHOLD 0.2


/******************************************************************************
 * convolve_1d is exactly what it sounds like.
******************************************************************************/
Pixel* convolve_1d(Pixel* x, Pixel* h, int x_len, int h_len) {
    Pixel* result = (Pixel*)calloc(x_len, sizeof(Pixel));  // Result array initialized with zeros
    for (int i = 0; i < x_len; i++) {
        for (int j = 0; j < h_len; j++) {
            int wrapped_index = (i - j + x_len) % x_len;  // Wrap-around index calculation
            result[i] += x[wrapped_index] * h[j];
        }
    }
    for (int i=0; i<x_len; i++) {
        result[i] /= h_len;
    }
    return result;
}


// Initialize a filter used for smoothing the signal
Pixel* initialize_1d_smooth_filter(int size) {
    Pixel* filter = (Pixel*)malloc(size * sizeof(Pixel));
    for(int i=0; i<size; i++) {
        filter[i] = 1;
    }
    return filter;
}


/******************************************************************************
 * initialize_3x3_laplacian simply builds a laplacian filter of size 3x3.
******************************************************************************/
Filter* initialize_3x3_laplacian() {
    Filter* filt = (Filter*)malloc(sizeof(Filter));
    filt->height = 3;
    filt->width = 3;

    filt->coefs = (double*)malloc(9 * sizeof(double));
    filt->coefs[0] = -1.0; filt->coefs[1] = -1.0; filt->coefs[2] = -1.0;
    filt->coefs[3] = -1.0; filt->coefs[4] = 8; filt->coefs[5] = -1.0;
    filt->coefs[6] = -1.0; filt->coefs[7] = -1.0; filt->coefs[8] = -1.0;
    return filt;
}


/******************************************************************************
 * sharpness_avg calculates a heuristic for sharpness, by averaging the values
 *   of a laplacian-filtered image that are above a given THRESHOLD. Function
 *   attempts to replace the laplacian-variance sharpness measure.
******************************************************************************/
Pixel sharpness_avg(Pixel* input, int length) {
    Pixel p;
    Pixel p_total = 0.0;
    int num_p = 0;
    for (int i=0; i<length; i++) {
        p = input[i];
        if (p > THRESHOLD) {
            p_total += p;
            num_p++;
        }
    }
    Pixel avg = p_total/((Pixel)(num_p));
    return avg;
}


/******************************************************************************
 * filter_image runs a Pixel array through a filter and returns a pointer to
 *   its output array.
 *  -filt is the MxN filter object. Can be any size.
 *  -input is the input array.
 *  -height and width are the image's height and width.
******************************************************************************/
Pixel* filter_image(Filter* filt, Pixel* input, int height, int width) {
    Pixel* output = (Pixel*)malloc(height * width * sizeof(Pixel));
    int filt_h = filt->height;
    int filt_w = filt->width;
    int yoffs = filt_h / 2;
    int xoffs = filt_w / 2;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double dotp = 0.0;
            for (int fy = 0; fy < filt_h; fy++) {
                for (int fx = 0; fx < filt_w; fx++) {
                    int iy = y + fy - yoffs;
                    int ix = x + fx - xoffs;

                    if (iy >= 0 && iy < height && ix >= 0 && ix < width) {
                        int img_idx = iy * width + ix;
                        int filt_idx = fy * filt_w + fx;
                        dotp += input[img_idx] * filt->coefs[filt_idx];
                    }
                }
            }
            output[y * width + x] = dotp;
        }
    }
    return output;
}


Image_RGB* create_filtered_RGB(Image_RGB* input, Filter* filt) {
    Image_RGB* output = create_rgb_image(input->width, input->height);
    output->r = filter_image(filt, input->r, output->height, output->width);
    output->g = filter_image(filt, input->g, output->height, output->width);
    output->b = filter_image(filt, input->b, output->height, output->width);
    return output;
}


/******************************************************************************
 * get_average gets the average value from an array of inputs and its length.
 *  -Function converts from Pixel to double to prevent round-to-nearest errors
 *   in accumulation between the massive accumulator and comparatively tiny
 *   input values.
******************************************************************************/
Pixel get_average(Pixel* input, int length) {
    double accumulator = 0;
    for (int i=0; i<length; i++) {
        accumulator += (double)input[i];
    }
    return (Pixel)(accumulator/(double)length);
}


/******************************************************************************
 * get_variance gets the variance of an array, given its length and average.
 *  -Function converts from Pixel to double to prevent round-to-nearest errors
 *   in accumulation between the massive accumulator and comparatively tiny
 *   input values.
******************************************************************************/
Pixel get_variance(Pixel* input, int length, Pixel average) {
    double accumulator = 0;
    for (int i=0; i<length; i++) {
        double diff = (double)(input[i]-average);
        accumulator += diff*diff;
    }
    accumulator = accumulator/(double)length;
    return accumulator;
}


Sharpnesses* get_variance_sharpness(Image_PGM* image, Crop_Boundaries* bounds) {
    if (bounds == NULL) {
        return NULL;
    }
    // run laplacian kernel filter through image
    Filter* filt = initialize_3x3_laplacian();

    Sharpnesses* sharpnesses = (Sharpnesses*)malloc(sizeof(Sharpnesses));
    sharpnesses->N = bounds->N;
    sharpnesses->sharpness = calloc(sharpnesses->N, sizeof(Pixel));

    for (int i=0; i<bounds->N; i++) {
        Image_PGM* cropped = crop_pgm(image, 
                                      bounds->right[i], bounds->left[i], 
                                      bounds->bottom[i], bounds->top[i]);
        
        Pixel* filtered = filter_image(filt, cropped->data, cropped->height, cropped->width);

        // get variance
        int length = cropped->height*cropped->width;
        Pixel avg = get_average(filtered, length);

        // Get the sharpness. NOTE: devide by avg to get a relatively scale-invariant measure
        sharpnesses->sharpness[i] = get_variance(filtered, length, avg) / avg;

        free(filtered);
        free_image_pgm(cropped);
    }
    // Clean up memory
    free(filt->coefs);
    free(filt);
    return sharpnesses;
}


Pixel get_average_sharpness(Pixel* input, int height, int width) {
    // Run laplacian kernel filter through image
    Filter* filt = initialize_3x3_laplacian();
    Pixel* filtered = filter_image(filt, input, height, width);

    // Get sharpness average
    int length = height*width;
    Pixel sharpness_average = sharpness_avg(filtered, length);

    // Clean up memory
    free(filt->coefs);
    free(filt);
    free(filtered);
    return sharpness_average;
}
