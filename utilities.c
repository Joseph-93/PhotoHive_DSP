#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "utilities.h"
#include "image_processing.h"
#include "color_quantization.h"


int num_cores;


/******************************************************************************
 * newton_int_sqrt returns the nearest integer approximation to the square root
 *   of a given double.
 *  -The Newton square-root algorithm is used, but manually set for accuracy
 *   only good enough to reach the nearest integer and then stop and return the
 *   truncated integer value
******************************************************************************/
int newton_int_sqrt(double val) {
    if (val==0) {return 0;}
    double x = val;
    double sqrt = x;
    while(1) {
        sqrt = 0.5 * (x + (val/x));
        if (fabs(sqrt-x)<1) {return (int)sqrt;}
        x = sqrt;
    }
}


/******************************************************************************
 * create_path returns a cleaned file path, retrieved from user input.
 *  -path must be the pre-filename path from the project root to the folder
 *   that is meant to hold the file
 *  -request_string is the message displayed to the user, requesting an input.
 *  -filetype is the file extension to be used.
 *  -No restraints are set on input variables or user input, as interface is
 *   only meant for developer testing and needs no safety constraints.
******************************************************************************/
char* create_path(const char* path, const char* request_string, const char* filetype) {
    // Request, get, and clean the filename
    char filename[40];
    printf("%s\n", request_string);
    fgets(filename, sizeof(filename), stdin);
    int len = strlen(filename);
    if (len > 0 && filename[len-1] == '\n') {
        filename[len-1] = '\0';
    }
    // Concatenate the filetype to the filename
    strcat(filename, filetype);

    // allocate memory for full path
    char* full_path = (char*)malloc(sizeof(char)*(strlen(path)+strlen(filename)+1));
    if (!full_path) {
        fprintf(stderr, "ERROR: Memory allocation for full_path failed.\n");
        return NULL;
    }

    // Assemble full path from filename and path
    strcpy(full_path, path);
    strcat(full_path, filename);
    printf("%s\n", full_path);
    return full_path;
}


// Initialization to work with CPU cores
void threading_setup() {
    num_cores = sysconf(_SC_NPROCESSORS_ONLN);
}


void custom_sort(void *base, size_t nmemb, size_t size,
                 int (*compar)(const void *, const void *, void *),
                 void *arg) {
    for (size_t i = 1; i < nmemb; i++) {
        for (size_t j = i; j > 0; j--) {
            char *elem_j = (char *)base + j * size;
            char *elem_j1 = (char *)base + (j - 1) * size;

            // Using the comparator function with the additional context argument
            if (compar(elem_j, elem_j1, arg) < 0) {
                // Swap elements
                for (size_t k = 0; k < size; k++) {
                    char temp = elem_j[k];
                    elem_j[k] = elem_j1[k];
                    elem_j1[k] = temp;
                }
            } else {
                break;
            }
        }
    }
}


// Return True or False for user-inputted yes/no question
bool ask_yes_no_question(const char *question) {
    char response[10]; // Buffer to store user input

    while (true) {
        printf("%s (yes/no): ", question);
        scanf("%9s", response); // Read the input from user and discard whitespace

        // Clear input buffer to remove any leftover characters, such as newline
        int c;
        while ((c = getchar()) != '\n' && c != EOF);

        // Check the response
        if (strncmp(response, "yes", 3) == 0) {
            return true; // User responded yes
        } else if (strncmp(response, "no", 2) == 0) {
            return false; // User responded no
        } else {
            // User gave an unrecognized response
            printf("Unrecognized response. Please type 'yes' or 'no'.\n");
        }
    }
}




void normalize_array(Pixel* array, int length) {
    Pixel min = array[0];
    Pixel max = array[0];
    for (int i=0; i<length; i++) {
        // Update min and max
        if (min > array[i]) {
            min = array[i];
        }
        if (max < array[i]) {
            max = array[i];
        }
    }

    // Scale Laplacian
    double range_multiplier = 1/(max - min);
    for (int i=0; i<length; i++) {
        array[i] = (array[i] - min)*range_multiplier;
    }
}


Full_Report_Data* compile_full_report(RGB_Statistics* rgb_stats,
                                 Color_Palette* color_palette,
                                 Blur_Profile_RGB* blur_profile,
                                 Pixel average_saturation,
                                 Pixel sharpness) {
    Full_Report_Data* full_report_data = (Full_Report_Data*)malloc(sizeof(Full_Report_Data));

    full_report_data->rgb_stats = rgb_stats;
    full_report_data->color_palette = color_palette;
    full_report_data->blur_profile = blur_profile;
    full_report_data->average_saturation = average_saturation;
    full_report_data->sharpness = sharpness;

    return full_report_data;
}


void print_full_report(Full_Report_Data* data) {
    char* path = create_path("reports/", "Enter a report filename: ", ".txt");
    FILE* f = fopen(path, "w");
    free(path);
    fprintf(f, "FULL REPORT:\n");
    fprintf(f, "Variance-based Sharpness: %lf\n", data->sharpness);
    fprintf(f, "Average Saturation: %lf\n", data->average_saturation);
    fprintf(f, "Brightness of RGB: (%lf,%lf,%lf)\n", data->rgb_stats->Br, data->rgb_stats->Bg, data->rgb_stats->Bb);
    fprintf(f, "Contrast of RGB; (%lf,%lf,%lf)\n", data->rgb_stats->Cr, data->rgb_stats->Cg, data->rgb_stats->Cb);
    
    fprintf(f, "\nColor Palette Contents:\n");
    for (int i=0; i<data->color_palette->N; i++) {
        Pixel_HSV p = data->color_palette->averages[i];
        fprintf(f, "%d\tHSV: (%3d,%3d,%3d), Portion of image accounted for: %lf\n", 
                i+1, (int)p.h, (int)(p.s*100), (int)(p.v*100), data->color_palette->percentages[i]);
    }

    fprintf(f, "\nBlur Profile:\n");
    Blur_Profile_RGB* bp = data->blur_profile;
    for(int i=0; i<bp->num_angle_bins; i++) {
        for (int j=0; j<bp->num_radius_bins; j++) {
            float norm_freq = ((float)j)/((float)bp->num_radius_bins);
            int angle = bp->angle_bin_size * i;
            fprintf(f, "angle: %3d, frequency: %.3f\t\t RGB Bins: (%lf,%lf,%lf)\n", 
                    angle, norm_freq, bp->r_bins[i][j], bp->g_bins[i][j], bp->b_bins[i][j]);
        }
    }
    fprintf(f, "\n\nEND OF REPORT.\n");
}