#include "types.h"

#ifndef IMAGE_PROCESSING_H
#define IMAGE_PROCESSING_H


/******************************************************************************
 * Pixel_HSV contains a pixel for an HSV image, specifically designed for
 *   octree grouping in color quantization.
 *  -parent_id is the id of the parent octree node for which the pixel belongs.
******************************************************************************/
typedef struct Pixel_HSV {
    int parent_id;
    double h;
    double s;
    double v;
} Pixel_HSV;


/******************************************************************************
 * Image_RGB holds a 2D image of RGB pixels.
 *  -Image_RGB is meant to hold r,g,b values as decimal numbers range [0, 1].
 *  -Effective as both an FFT representation and a literal image.
 *  -Cannot hold complex values.
 *  -The separation of R,G, and B channels is built for efficiency in filtering
 *   and FFT computation, understanding the slightly higher complexity in image
 *   representation.
 *  -The image is logically 2D, height by width, but is referenced as 1D, for
 *   increased efficiency in reference. Proper indexing: img.color[y*width + x]
******************************************************************************/
typedef struct Image_RGB {
    int height, width;
    Pixel* r;
    Pixel* g;
    Pixel* b;
} Image_RGB;


/******************************************************************************
 * Image_HSV operates the same as Image_RGB, except that it is an HSV format,
 *   obviously!
 *  -The reason for separating the types is to facilitate static typing, as
 *   many functions are built to operate for HSV or RGB images exclusively.
 *  -The separation of H, S, and V channels is built for computational
 *   efficiency, understanding the slightly higher complexity in image
 *   representation.
 *  -The image is logically 2D, height by width, but is referenced as 1D, for
 *   increased efficiency in reference. Proper indexing: img.color[y*width + x]
******************************************************************************/
typedef struct Image_HSV {
    int height, width;
    Pixel_HSV* pixels;
} Image_HSV;


/******************************************************************************
 * Image_PGM holds a 2D image of grayscale pixels.
 *  -Effective as both an FFT representation and a literal image
 *  -Cannot store complex values
 *  -The image is logically 2D, height x width, but is referenced as 1D, for
 *   increased efficiency in reference Proper indexing: img.data[y*width + x].
******************************************************************************/
typedef struct Image_PGM {
    int height, width;
    Pixel* data;
} Image_PGM;


/******************************************************************************
 * RGB_Statistics is simply a structure that holds values for the basic RGB
 *   measurements- brightness and contrast of RGB channels.
******************************************************************************/
typedef struct RGB_Statistics {
    Pixel Br;
    Pixel Bg;
    Pixel Bb;
    Pixel Cr;
    Pixel Cg;
    Pixel Cb;
} RGB_Statistics;


/******************************************************************************
 * Crop_Boundaries simply holds the integer values for the top, bottom, left,
 *   and right boundaries for image crops.
 *  -N is the number of crops
 *  -top is a pointer to all N top boundaries
 *  -bottom is a pointer to all N bottom boundaries
 *  -left is a pointer to all N left boundaries
 *  -right is a pointer to all N right boundaries
******************************************************************************/
typedef struct Crop_Boundaries {
    int N;
    int* top;
    int* bottom;
    int* left;
    int* right;
} Crop_Boundaries;


Image_RGB* create_rgb_image(int width, int height);

Image_HSV* create_hsv_image(int width, int height);

Image_PGM* create_pgm_image(int width, int height);

Image_RGB* read_image(const char* filepath);

void write_image_to_file(Image_RGB* image, const char* path);

Image_PGM* crop_pgm(Image_PGM* image, int right, int left, int bottom, int top);

Image_RGB* crop_image(Image_RGB* image, int right, int left, int bottom, int top);

Image_HSV* rgb2hsv(Image_RGB* rgb);

Image_RGB* hsv2rgb(Image_HSV* hsv);

void free_image_rgb(Image_RGB* image);

void free_image_hsv(Image_HSV* hsv);

Image_RGB* read_image_from_files();

Image_RGB* downsample_rgb(Image_RGB* image, short N);

void validate_coordinates(Image_PGM* array, int x, int y, int* errors);

Image_PGM* rgb2pgm(Image_RGB* rgb);

Image_RGB* pgm2rgb(Image_PGM* pgm);

void free_crop_boundaries(Crop_Boundaries* cb);

void free_image_pgm(Image_PGM* image);

Pixel get_hsv_average(Image_HSV* hsv);

RGB_Statistics* get_rgb_statistics(Image_RGB* image);

#endif