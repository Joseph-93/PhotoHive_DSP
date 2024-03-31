#include <stdio.h>
#include "color_quantization.h"
#include "utilities.h"

#ifndef DEBUG_H
#define DEBUG_H

// #define VERBOSE

void display_octree_settings(Octree* octree);

void compare_rgb_to_hsv(Image_HSV* hsv, Image_RGB* rgb);

Image_RGB* create_test_rgb(int height, int width);

void verify_arm_octree(Octree* octree, Image_HSV* hsv);

void validate_octree_parents(Octree* octree);

void generate_color_palette_image(Octree* octree, Image_HSV* hsv);

void report_color_palette(Color_Palette* cp);

#endif