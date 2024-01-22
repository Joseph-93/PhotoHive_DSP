#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fftw3.h>
#include <math.h>

typedef struct {
    double r;
    double g;
    double b;
} Pixel_RGB;    // Holds the RGB values of a pixel

typedef struct {
    int height, width;
    Pixel_RGB** pixels;
} Image_RGB;    // Holds a 2D image of RGB values

typedef struct {
    double v;
} Pixel_Gray;   // Holds a single value for a grayscale pixel

typedef struct {
    int height, width;
    Pixel_Gray** pixels;
} Image_Gray;   // Holds a 2D image of grayscale pixels

typedef struct {
    int r_sq;
    double phi;
} Polar_Coord; // Struct contains a polar representation of cartesian coords

typedef struct {
    int height, width;
    Polar_Coord** data;
} Cartesian_To_Polar;   // An x,y coordinate of cartesian-polar conversions

typedef struct {
    double sum;
} Bin;

typedef struct {
    int angle_bin_size, radius_bin_size;
    int num_angle_bins, num_radius_bins;
    Bin** r_bins; // bins should be referenced as [angle][radius]
    Bin** g_bins; // bins should be referenced as [angle][radius]
    Bin** b_bins; // bins should be referenced as [angle][radius]
} Blur_Profile_RGB;     // Holds blur profile of bins for R,G,B


// create a mapping for cartesian to polar conversions, built for a given fft
// height and width. CAUTION: is built to assume fft is notshifted
Cartesian_To_Polar* cartesian_to_polar_conversion(int width, int height) {
    Polar_Coord coord;
    Cartesian_To_Polar* conversion = (Cartesian_To_Polar*)malloc(sizeof(Cartesian_To_Polar*));
    conversion->height = height, conversion->width = width;
    conversion->data = (Polar_Coord**)malloc(height * sizeof(Polar_Coord*));
    for (int y=0; y<height/2; y++) {    // height/2 for symmetry across x axis
        conversion->data[y] = (Polar_Coord*)malloc(width * sizeof(Polar_Coord));
        for (int x=0; x<width; x++) {
            int r_sq = x*x+y*y;
            double phi = atan(y/x);
            // Not fft shifted, so for the top half, phi=-atan(y/x)
            conversion->data[y][x].phi = -phi;
            conversion->data[y][x].r_sq = r_sq;
            // Not fft shifted, so for the bottom half phi=atan(y/x)
            conversion->data[y][x].phi = phi;
            conversion->data[y-height/2][x].r_sq = r_sq;
        }
    }
    return conversion;
}


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


// return calculated blur profile, given the conversion table and the fft
Blur_Profile_RGB* calculate_blur_profile(
                    const Cartesian_To_Polar* conversion, 
                    const Image_RGB* fft, 
                    const int* num_radius_bins, 
                    const int* num_angle_bins) {
    // Basic blur profile setup
    Blur_Profile_RGB* profile = (Blur_Profile_RGB*)malloc(sizeof(Blur_Profile_RGB*));
    profile->num_angle_bins = *num_angle_bins;
    profile->num_radius_bins = *num_radius_bins;
    profile->angle_bin_size = (double)(180 / *num_angle_bins);
    profile->radius_bin_size = (double)(fft->width / *num_radius_bins);
    double radius_bin_size_sq = profile->radius_bin_size*profile->radius_bin_size;

    // Allocate space for empty bins
    profile->r_bins = (Bin**)malloc(*num_angle_bins * sizeof(Bin*));
    profile->g_bins = (Bin**)malloc(*num_angle_bins * sizeof(Bin*));
    profile->b_bins = (Bin**)malloc(*num_angle_bins * sizeof(Bin*));
    for (int i=0; i<*num_angle_bins; i++) {
        profile->r_bins[i] = (Bin*)malloc(*num_radius_bins * sizeof(Bin*));
        profile->g_bins[i] = (Bin*)malloc(*num_radius_bins * sizeof(Bin*));
        profile->b_bins[i] = (Bin*)malloc(*num_radius_bins * sizeof(Bin*));
    }

    // Sum up blurs in each bin space
    for (int y=0; y<fft->height; y++) {     //TODO: assure that fft_shift is accomodated properly
        for (int x=0; x<fft->width; x++) {
            // Store r_sq and phi for readability
            int r_sq = conversion->data[y][x].r_sq;
            double phi = conversion->data[y][x].phi;
            // phi_bin is found by int division of the phi value by the bin size
            int phi_bin = (int)phi/profile->angle_bin_size;
            // r_bin is sqrt(r^2/r_bin_size^2), but a newton int approximation
            int r_bin = newton_int_sqrt((double)r_sq/radius_bin_size_sq);
            profile->r_bins[phi_bin][r_bin].sum += fft->pixels[y][x].r;
            profile->g_bins[phi_bin][r_bin].sum += fft->pixels[y][x].g;
            profile->b_bins[phi_bin][r_bin].sum += fft->pixels[y][x].b;
        }
    }
}


Image_RGB* create_rgb_image(int width, int height) {
    Image_RGB* image = (Image_RGB*)malloc(sizeof(Image_RGB*));
    image->width = width;
    image->height = height;
    image->pixels = (Pixel_RGB**)malloc(height * sizeof(Pixel_RGB*));
    for (int i=0; i<height; i++) {
        image->pixels[i] = (Pixel_RGB*)malloc(width * sizeof(Pixel_RGB));
    }
    return image;
}


Image_Gray* create_gray_image(int width, int height) {
    Image_Gray* image = (Image_Gray*)malloc(sizeof(Image_Gray));
    image->width = width;
    image->height = height;
    image->pixels = (Pixel_Gray**)malloc(height * sizeof(Pixel_Gray*));
    for (int i=0; i<height; i++) {
        image->pixels[i] = (Pixel_Gray*)malloc(width * sizeof(Pixel_Gray));
    }
    return image;
}


Image_RGB* read_image(const char* filepath) {
    FILE* file = fopen(filepath, "r");
    if (file == NULL) {
        perror("Error occured opening file");
        return NULL;
    }
    int width, height;
    fscanf(file, "%d %d", &width, &height);

    // Allocate memory for the image array
    Pixel_RGB** pixel_array = (Pixel_RGB**)malloc(height * sizeof(Pixel_RGB*));
    for (int i=0; i<height; i++) {
        pixel_array[i] = (Pixel_RGB*)malloc(width * sizeof(Pixel_RGB));
    }

    // Read in the RGB values for the image
    for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
            int r, g, b;
            fscanf(file, "%d %d %d", &r, &g, &b);
            pixel_array[y][x].r = (double) r / 255.0f;
            pixel_array[y][x].g = (double) g / 255.0f;
            pixel_array[y][x].b = (double) b / 255.0f;
        }
    }
    Image_RGB* image = (Image_RGB*)malloc(sizeof(Image_RGB*));
    image->pixels = pixel_array;
    image->height = height;
    image->width = width;
    fclose(file);
    return image;
}


void write_image_to_file(Image_RGB* image, const char* path) {
    // Open file for writing
    FILE* file = fopen(path, "w");
    if (file == NULL) {perror("Error opening file."); return;}
    fprintf(file, "%d %d\n", image->width, image->height);

    // Begin copying file
    int r, g, b;
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<image->width; x++) {
            r = (int)(image->pixels[y][x].r * 255);
            g = (int)(image->pixels[y][x].g * 255);
            b = (int)(image->pixels[y][x].b * 255);
            fprintf(file, "%d %d %d\n", r, g, b);
        }
    }
    fclose(file);
}


Image_RGB* crop_image(Image_RGB* image, int right, int left, int bottom, int top) {
    if (right > image->width | left > image->width | bottom > image->height | top > image->height |
        left < 0 | right < 0 | top < 0 | bottom < 0) {
        perror("Error: crop boundaries outside of image boundaries.");
        return NULL;
    }

    // Calculate cropped image width and height, then create the image
    int width = right - left;
    int height = bottom - top;
    Image_RGB* cropped_image = create_rgb_image(width, height);

    // Copy the portion of the original to be kept, onto the cropped image
    for (int x=0; x<width; x++) {
        for (int y=0; y<height; y++) {
            cropped_image->pixels[y][x] = image->pixels[y + top][x + left];
        }
    }
    return cropped_image;
}


void separate_rgb_channels(Image_RGB* image, double** r_channel, double** g_channel, double** b_channel) {
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<image->width; x++) {
            r_channel[y][x] = image->pixels[y][x].r;
            g_channel[y][x] = image->pixels[y][x].g;
            b_channel[y][x] = image->pixels[y][x].b;
        }
    }
}


double* convert_to_1d(double** array, int height, int width) {
    double* array_1d = (double*)malloc(height * width * sizeof(double));
    for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
            array_1d[y*width + x] = array[y][x];
        }
    }
    return array_1d;
}


Image_RGB* rgb_fft(Image_RGB* image, double** r_channel, double** g_channel, double** b_channel) {
    // Separate rgb channels from the image
    separate_rgb_channels(image, r_channel, g_channel, b_channel);

    double* r_1d = convert_to_1d(r_channel, image->height, image->width);
    double* g_1d = convert_to_1d(g_channel, image->height, image->width);
    double* b_1d = convert_to_1d(b_channel, image->height, image->width);

    fftw_complex* output = fftw_malloc(sizeof(fftw_complex) * image->width * image->height);
    double* in = fftw_malloc(sizeof(double) * image->height * image->width);

    // Create a plan
    fftw_plan plan = fftw_plan_dft_r2c_2d(image->height, image->width, in, output, FFTW_ESTIMATE);

    int fft_width = image->width/2+1;
    Image_RGB* fft_image = create_rgb_image(fft_width, image->height);

    // Execute fftw for red, collect magnitudes in fft_image
    memcpy(in, r_1d, sizeof(double)*image->width*image->height);
    fftw_execute(plan);
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<fft_width; x++) {
            fftw_complex fft_point = {output[y*fft_width + x][0], output[y*fft_width + x][1]};
            double mag_sq = fft_point[0]*fft_point[0] + fft_point[1]*fft_point[1];
            fft_image->pixels[y][x].r = mag_sq;
        }
    }

    // Execute fftw for green, collect magnitudes in fft_image
    memcpy(in, g_1d, sizeof(double)*image->width*image->height);
    fftw_execute(plan);
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<fft_width; x++) {
            fftw_complex fft_point = {output[y*fft_width + x][0], output[y*fft_width + x][1]};
            double mag_sq = fft_point[0]*fft_point[0] + fft_point[1]*fft_point[1];
            fft_image->pixels[y][x].g = mag_sq;
        }
    }

    // Execute fftw for green, collect magnitudes in fft_image
    memcpy(in, b_1d, sizeof(double)*image->width*image->height);
    fftw_execute(plan);
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<fft_width; x++) {
            fftw_complex fft_point = {output[y*fft_width + x][0], output[y*fft_width + x][1]};
            double mag_sq = fft_point[0]*fft_point[0] + fft_point[1]*fft_point[1];
            fft_image->pixels[y][x].b = mag_sq;
        }
    }

    return fft_image;
}


double** allocate_2d_array(int width, int height) {
    double** array = malloc(height * sizeof(double*));
    for (int y=0; y<height; y++) {
        array[y] = malloc(width * sizeof(double));
    }
    return array;
}


void free_2d_array(double** array, int height) {
    for (int y=0; y<height; y++) {
        free(array[y]);
    }
    free(array);
}


char* create_path(const char* path, const char* request_string, const char* filetype) {
    // Request, get, and clean the filename
    char filename[40];
    printf(request_string);
    fgets(filename, sizeof(filename), stdin);
    int len = strlen(filename);
    if (len > 0 && filename[len-1] == '\n') {
        filename[len-1] = '\0';
    }
    // Concatenate the filetype to the filename
    strcat(filename, filetype);

    // full_path = path + filename
    char* full_path = (char*)malloc(sizeof(char)*(strlen(path)+strlen(filename)+1));
    if (full_path == NULL) {perror("Error allocating memory"); return NULL;}
    strcpy(full_path, path);
    strcat(full_path, filename);
    printf("%s\n", full_path);
    return full_path;
}


// Function normalizes an RGB image, optimizing for color balance preservation
void rgb_normalize_photo(Image_RGB* image) {
    // Set Minimums and Maximums as any valid value
    double max = image->pixels[0][0].r;
    double min = max;

    // Set maximum and minimum values
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<image->width; x++) {
            Pixel_RGB pixel = image->pixels[y][x];            
            if (max < pixel.r) {max = pixel.r;}
            if (min > pixel.r) {min = pixel.r;}
            if (max < pixel.g) {max = pixel.g;}
            if (min > pixel.g) {min = pixel.g;}
            if (max < pixel.b) {max = pixel.b;}
            if (min > pixel.b) {min = pixel.b;}
        }
    }

    printf("Maximum value found in the NOT normalized FFT: %f\n", max);
    printf("Minimum value found in the NOT normalized FFT: %f\n", min);

    // Pull the logarithm of the square root of each value. The +1 is to avoid dividing by zero.
    double mag;
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<image->width; x++) {
            // Replace image red pixel with its own calculated magnitude
            mag = sqrt(image->pixels[y][x].r);
            mag = log(mag+1);
            image->pixels[y][x].r = mag;

            // Replace image red pixel with its own calculated magnitude
            mag = sqrt(image->pixels[y][x].g);
            mag = log(mag+1);
            image->pixels[y][x].g = mag;

            // Replace image red pixel with its own calculated magnitude
            mag = sqrt(image->pixels[y][x].b);
            mag = log(mag+1);
            image->pixels[y][x].b = mag;
        }
    }

    // Set Minimums and Maximums as any valid value
    max = image->pixels[0][0].r;
    min = max;

    // Set maximum and minimum values
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<image->width; x++) {
            Pixel_RGB pixel = image->pixels[y][x];            
            if (max < pixel.r) {max = pixel.r;}
            if (min > pixel.r) {min = pixel.r;}
            if (max < pixel.g) {max = pixel.g;}
            if (min > pixel.g) {min = pixel.g;}
            if (max < pixel.b) {max = pixel.b;}
            if (min > pixel.b) {min = pixel.b;}
        }
    }

    printf("Maximum value found in the normalized FFT: %f\n", max);
    printf("Minimum value found in the normalized FFT: %f\n", min);

    // Normalize values
    double ratio = 1 / (max - min);
    for (int y=0; y<image->height; y++) {
        for (int x=0; x<image->width; x++) {
            // Bring up values to zero, then normalize with the ratio
            image->pixels[y][x].r = (image->pixels[y][x].r - min) * ratio;
            image->pixels[y][x].g = (image->pixels[y][x].g - min) * ratio;
            image->pixels[y][x].b = (image->pixels[y][x].b - min) * ratio;
        }
    }
}


void validate_coordinates(Image_RGB* array, int x, int y, int* errors) {
    if (x>=array->width) {
        *errors +=1;
        printf("ERROR %3d: x is above the array width: x=%3d, width=%3d.\n", *errors, x, array->width);
    }
    else if (x<0) {
        *errors +=1;
        printf("ERROR %3d: x is below 0. x=%3d.\n", *errors, x);
    }
    if (y>=array->height) {
        *errors +=1;
        printf("ERROR %3d: y is above the array height: y=%3d, height=%3d.\n", *errors, y, array->height);
    }
    else if (y<0) {
        *errors +=1;
        printf("ERROR %3d: y is below zero. y=%3d.\n", *errors, y);
    }
}


Image_RGB* fft_shift(Image_RGB* fft) {
    Image_RGB* fft_image = create_rgb_image(fft->width*2 - 1, fft->height);
    // Iterate through one quarter of needed image, via one half of the image
    printf("new_image size: (%3d, %3d)\t old_image size: (%3d, %3d)\n", fft_image->width, fft_image->height, fft->width, fft->height);
    
    // If even image, height_compensator is 0. if odd image, compensator is 1
    int height_compensator = (int)(fft_image->height%2==1);
    int width_compensator = (int)(fft_image->width%2==1); // width comp is always 1, but I want to make it dynamic.
    printf("height_equalizer: %d\t width_equalizer: %d\n", height_compensator, width_compensator);

    int errors = 0;
    Pixel_RGB temp;
    int half_fft_height = fft->height/2;
    int y_val, x_val;
    for (int y=0; y<half_fft_height+height_compensator; y++) {
        for (int x=0; x<fft->width; x++) {
            if (x==161 && y==161) {
                printf("(%d, %d)\n", x, y);
            }
            // printf("(%d, %d)\n", x, y);

            // Swap Q1 and Q4
            y_val = y+half_fft_height, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            temp = fft->pixels[y_val][x_val]; // New Q1 gets old Q4
            if (y-height_compensator != -1) { // If height is odd, middle row goes in a theoretical pixels[-1]
                y_val = y-height_compensator, x_val = x+fft->width-width_compensator, validate_coordinates(fft_image, x_val, y_val, &errors);
                fft_image->pixels[y_val][x_val] = temp;
            }
            y_val = y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            temp = fft->pixels[y_val][x_val]; // New Q4 gets old Q1
            y_val = y+half_fft_height, x_val = x+fft->width-width_compensator, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->pixels[y_val][x_val] = temp;
            
            // Rotate to get Q2 and Q3
            y_val = y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            temp = fft->pixels[y_val][x_val]; // Q2 is 180 deg rotation of Q4
            y_val = half_fft_height-y, x_val = fft->width-1-x, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->pixels[y_val][x_val] = temp; //width-1 because width-1=last_index
            y_val = half_fft_height+y, x_val = x, validate_coordinates(fft, x_val, y_val, &errors);
            temp = fft->pixels[y_val][x_val]; // Q3 is 180 deg rotation of Q1
            y_val = fft->height-1-y, x_val = fft->width-1-x, validate_coordinates(fft_image, x_val, y_val, &errors);
            fft_image->pixels[y_val][x_val] = temp; // height-1 because height-1=last_index
        }
    }
    return fft_image;
}





int main() {
    // read in user's desired image from the list of readable .txt files
    char* input_filename = create_path("images/readable/", "Enter the name of the .txt file you wish to use: ", ".txt");
    Image_RGB* image = read_image(input_filename);
    if (image == NULL) {
        return -1;
    }

    // Allocate arrays for each channel
    double** r_channel = allocate_2d_array(image->width, image->height);
    double** g_channel = allocate_2d_array(image->width, image->height);
    double** b_channel = allocate_2d_array(image->width, image->height);

    // Compute the FFT for RGB image
    Image_RGB* fft = rgb_fft(image, r_channel, g_channel, b_channel);
    rgb_normalize_photo(fft);
    Image_RGB* fft_image = fft_shift(fft);

    // write fft_image to user's desired image file (in .txt format for python to save via imageSaver.py)
    char* output_filename = create_path("images/output/", "Enter the name of the .txt file you wish to write to: ", ".txt");
    write_image_to_file(fft_image, output_filename);

    return 0;
}
