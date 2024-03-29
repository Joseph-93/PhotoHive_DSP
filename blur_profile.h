
#include "types.h"

#ifndef BLUR_PROFILE_H
#define BLUR_PROFILE_H

// Official data type used for bin sums in Blur_Profiles.
typedef double Bin;

/******************************************************************************
 * Blur_Profile_RGB represents the strength of directional blur, per angle, by
 *   sum of all pixels within an FFT subsurface bounded by radii and angles.
 *  -num_angle_bins and num_radius_bins are parameters given by the user, to be
 *   tested by the machine learning model for optimal performance and accuracy.
 *  -angle_bin_size is found by total angular range divided by num_angle_bins.
 *  -radius_bin_size is determined by FFT width divided by num_radius_bins.
 *  -Bin* r, g, and b bins are 2D logically (angle by radius), but are
 *   literally 1D. Properly indexed as {color}_bins[angle*num_angle_bins + rad]
******************************************************************************/
typedef struct Blur_Profile_RGB {
    int num_angle_bins, num_radius_bins;
    int angle_bin_size, radius_bin_size;
    Bin** r_bins; // bins should be referenced as [angle][radius]
    Bin** g_bins; // bins should be referenced as [angle][radius]
    Bin** b_bins; // bins should be referenced as [angle][radius]
} Blur_Profile_RGB;     // Holds blur profile of bins for R,G,B


/******************************************************************************
 * Polar_Coord contains a single set of polar coordinates
 *  -r_sq is equivalent to radius^2, for computational efficiency
 *  -phi is the angle arctan(y/x), limited to values between -90 and 90 degrees
 *   due to the reflective property of the FFT of a real-valued 2D image.
******************************************************************************/
typedef struct Polar_Coord {
    int r_sq;
    double phi;
} Polar_Coord;


/******************************************************************************
 * Cartesian_To_Polar an x,y coordinate of cartesian-polar conversions is
 *  -Essentially a replacement for a lookup table, because of the necessity for
 *   adherence to dynamic size requirements.
 *  -Table is logically 2D height by width, though it is literally 1D for
 *   efficient indexing.
******************************************************************************/
typedef struct Cartesian_To_Polar {
    unsigned int height, width;
    Polar_Coord* data;
} Cartesian_To_Polar;


typedef struct Blur_Vector {
    int angle;
    float magnitude;
} Blur_Vector;


typedef struct Blur_Vector_RGB {
    int len_vectors;
    Blur_Vector** blur_vectors_rgb;
} Blur_Vector_RGB;


Blur_Vector_RGB* vectorize_blur_profile(Blur_Profile_RGB* blur_profile,
                                     Pixel error_thresh,
                                     Pixel mag_thresh,
                                     int cutoff_ratio_denom);

Blur_Profile_RGB* calculate_blur_profile(
                    const Cartesian_To_Polar* conversion, 
                    const Image_RGB* fft, 
                    int num_radius_bins, 
                    int num_angle_bins);

Image_RGB* get_blur_profile_visual(Blur_Profile_RGB* blur_profile, 
                                   int height, int width);

void view_2d_array_contents(Bin** array, int height, int width);

void blur_profile_tests(Image_RGB* fft);

Cartesian_To_Polar* cartesian_to_polar_conversion(unsigned int width, unsigned int height);

void free_cartesian_to_polar(Cartesian_To_Polar* c2p);

Blur_Profile_RGB* get_blur_profile(Image_RGB* image, int num_radius_bins, int num_angle_bins);

void free_blur_profile_rgb(Blur_Profile_RGB* bp);

void free_blur_vectors_rgb(Blur_Vector_RGB* bv);

#endif