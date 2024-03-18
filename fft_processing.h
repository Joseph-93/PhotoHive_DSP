
#include "types.h"
#include "image_processing.h"

#ifndef FFT_PROCESSING_H
#define FFT_PROCESSING_H

/******************************************************************************
 * Lookup_1D is a simple lookup table structure, used to store values for any
 *   lookup table where a 1D input represented by Pixels maps to a 1D output
 *   represented by Pixels.
******************************************************************************/
typedef struct Lookup_1D {
    unsigned int length;
    Pixel* input;
    Pixel* output;
} Lookup_1D;


Image_RGB* rgb_fft(Image_RGB* image);

Image_RGB* compute_magnitude_fft(Image_RGB* image);

void save_fft(Image_RGB* fft);

Image_RGB* fft_shift(Image_RGB* fft);

void rgb_normalize_fft(Image_RGB* fft, Lookup_1D* fft_normalizer_lookup);

#endif