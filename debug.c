#include <stdlib.h>
#include "types.h"

#include "debug.h"

// #define VERBOSE

void display_octree_settings(Octree* octree) {
    printf("<=================OCTREE SETTINGS====================>\n\n");

    // print out octree settings
    printf("num_h: %d\t", octree->num_h);
    printf("num_s: %d\t", octree->num_s);
    printf("num_v: %d\t", octree->num_v);
    printf("num_grays: %d\t", octree->num_grays);
    printf("Total length: %d\n", octree->total_length);

    printf("black_thresh: %lf\n", octree->black_thresh);
    printf("gray_thresh: %lf\n", octree->gray_thresh);

    printf("Length of h: %lf\n", octree->Lh);
    printf("Length of s: %lf\n", octree->Ls);
    printf("Length of v: %lf\n\n", octree->Lv);

    // Print out group information
    printf("GROUPS:\n");
    int i;
    printf("Colored groups:\n");
    for (i=0; i<octree->total_length; i++) {
        if (i == octree->total_length-(octree->num_grays+1)) printf("\nGrays:\n");
        if (i == octree->total_length-1) printf("\nBlacks\n");
        printf("\tnode id: %d\n", octree->groups[i].id);
        printf("\t\tHSV values: (%lf, %lf, %lf)\n", octree->groups[i].h, octree->groups[i].s, octree->groups[i].v);
        printf("\t\tquantity: %d\n", octree->groups[i].quantity);
        printf("\t\tis_valid_parent: %d\n", octree->groups[i].is_valid_parent);
        printf("\t\tGroup head pointer: %p\n", octree->groups[i].head);
    }
    return;
}


void compare_rgb_to_hsv(Image_HSV* hsv, Image_RGB* rgb) {
    for (int i=0; i<rgb->height*rgb->width; i++) {
        double r=rgb->r[i]*255, g=rgb->g[i]*255, b=rgb->b[i]*255;
        double h=hsv->pixels[i].h;
        double s=hsv->pixels[i].s*100;
        double v=hsv->pixels[i].v*100;
        printf("i: %5d\trgb: (%lf, %lf, %lf)\thsv: (%lf, %lf, %lf)\n", i, r, g, b, h, s, v);
    }
}


Image_RGB* create_test_rgb(int height, int width) {
    Image_RGB* image = create_rgb_image(width, height);
    for (int i=0; i<height*width; i++) {
        image->r[i] = 1/(double)(i+1);
        image->g[i] = 1/(double)(i*2+1);
        image->b[i] = 1/(double)(i*4+1);
    }
    return image;
}


void verify_arm_octree(Octree* octree, Image_HSV* hsv) {
    // Verify that octree and hsv have the same number of pixels
    printf("<=============================Verify_Arm_Octree========================>\n\n");
    bool error_found = false;
    int num_octree_pixels = 0;
    for (int i=0; i<octree->total_length; i++) {
        HSV_Linked_List* list = octree->groups[i].head;
        while (list) {
            if (list->num_pixels > list->array_size) {
                error_found = true;
                printf("-------------list->num_pixels > list->array_size.-------------\n");
            }
            num_octree_pixels += list->num_pixels;
            list = list->next;
        }
    }
    printf("height: %d\t width: %d\n", hsv->height, hsv->width);
    int num_hsv_pixels = hsv->height * hsv->width;
    printf("Number of octree pixels: %d\t Number of hsv pixels: %d\n", num_octree_pixels, num_hsv_pixels);

    for (int i=0; i<octree->total_length-(octree->num_grays+1); i++) {
        HSV_Linked_List* list = octree->groups[i].head;
        while (list) {
            int* coords = get_octree_hsv_coords(octree, i);
            for (int j=0; j<list->num_pixels; j++) {
                Pixel_HSV p = list->pixels[j];
                Pixel h_min = octree->Lh*(Pixel)coords[0];
                Pixel h_max = h_min + octree->Lh;
                Pixel s_min = octree->gray_thresh + octree->Ls*(Pixel)coords[1];
                Pixel s_max = s_min + octree->Ls;
                Pixel v_min = octree->black_thresh + octree->Lv*(Pixel)coords[2];
                Pixel v_max = v_min + octree->Lv;
                if (p.h < h_min) {
                    error_found = true;
                    printf("Group i=%d, pixel j=%d found error: pixel.h < h_min.\n", i, j);
                    printf("\tpixel.h: %lf\t h_min: %lf\n\n", p.h, h_min);
                }
                if (p.h > h_max) {
                    error_found = true;
                    printf("Group i=%d, pixel j=%d found error: pixel.h > h_max.\n", i, j);
                    printf("\tpixel.h: %lf\t h_max: %lf\n\n", p.h, h_max);
                }
                if (p.s < s_min) {
                    error_found = true;
                    printf("Group i=%d, pixel j=%d found error: pixel.s < s_min.\n", i, j);
                    printf("\tpixel.h: %lf\t s_min: %lf\n\n", p.s, s_min);
                }
                if (p.s > s_max) {
                    error_found = true;
                    printf("Group i=%d, pixel j=%d found error: pixel.s > s_max.\n", i, j);
                    printf("\tpixel.s: %lf\t s_min: %lf\n\n", p.s, s_min);
                }
                if (p.v < v_min) {
                    error_found = true;
                    printf("Group i=%d, pixel j=%d found error: pixel.v < v_min.\n", i, j);
                    printf("\tpixel.v: %lf\t v_min: %lf\n\n", p.v, v_min);
                }
                if (p.v > v_max) {
                    error_found = true;
                    printf("Group i=%d, pixel j=%d found error: pixel.v > v_min.\n", i, j);
                    printf("\tpixel.v: %lf\t v_max: %lf\n\n", p.v, v_min);
                }
            }
            list = list->next;
        }
    }
    if (!error_found) {printf("\n======>No errors were found in arm_octree().\n\n");}
}


void validate_octree_parents(Octree* octree) {
    printf("<====================Validate_Octree_Parents=======================>\n");
    #ifdef VERBOSE
    for (int i=0; i<octree->len_valid_parents; i++) {
        printf("i:%3d,\tgroup:%4d,\tquantity:%4d\n", i, octree->valid_parents[i], octree->groups[octree->valid_parents[i]].quantity);
    }
    #endif

    bool error_found = false;
    int old_quantity = octree->groups[octree->valid_parents[0]].quantity;
    printf("length of valid parents: %d\n", octree->len_valid_parents);
    for (int i=1; i<octree->len_valid_parents; i++) {
        int new_quantity = octree->groups[octree->valid_parents[i]].quantity;
        if (old_quantity < new_quantity) {
            error_found = true;
            printf("i: %3d,\tGroup[%3d].quantity: %4d,\tGroup[%3d].quantity: %4d\n", 
                i, octree->valid_parents[i], new_quantity, octree->valid_parents[i-1], old_quantity);
        }
        old_quantity = new_quantity;
    }
    if (!error_found) {
        printf("=====>No errors were found inside of validate_octree_sort().\n");
    }
}


/******************************************************************************
 * Report the results of group_irregular_pixels().
 *  -Function reports all valid_parents, with their HSV values and quantities.
******************************************************************************/
void generate_color_palette_image(Octree* octree, Image_HSV* hsv) {
    #ifdef VERBOSE
    printf("<====================group_irregular_pixels REPORT================>\n\n");
    int total_pixels = 0;
    printf("Valid parents:\n");
    for (int i=0; i<octree->len_valid_parents; i++) {
        Octree_Group parent = octree->groups[octree->valid_parents[i]];
        total_pixels += parent.quantity;
        printf("i=%d\t\tid:%4d,\tHSV: (%f,%f,%f),\tQuantity:%5d\n", i, parent.id, parent.h, parent.s, parent.v, parent.quantity);
    }
    printf("\nTOTALS:\n");
    printf("hsv pixels:%5d\toctree pixels:%5d\n\n", hsv->height*hsv->width, total_pixels);
    #endif

    // Set up dimensions. Each color is represented as a 200x200 block.
    int num_colors = octree->len_valid_parents;
    int block_size = 50; // Each color will be a block of 200x200 pixels
    int height = newton_int_sqrt(num_colors) * block_size;
    int width = ((num_colors + height/block_size - 1) / (height/block_size)) * block_size; // Ensure full rows

    // Create image
    Image_HSV* hsv_out = create_hsv_image(width, height);
    if (hsv_out == NULL) {
        // handle error
    }

    for (int color_index = 0; color_index < num_colors; color_index++) {
        Octree_Group parent = octree->groups[octree->valid_parents[color_index]];
        
        // Calculate starting position for current block
        int start_row = (color_index / (width / block_size)) * block_size;
        int start_col = (color_index % (width / block_size)) * block_size;

        for (int row = start_row; row < start_row + block_size; row++) {
            for (int col = start_col; col < start_col + block_size; col++) {
                int pixel_index = row * width + col;
                hsv_out->pixels[pixel_index].h = parent.h;
                hsv_out->pixels[pixel_index].s = parent.s;
                hsv_out->pixels[pixel_index].v = parent.v;
            }
        }
    }

    Image_RGB* rgb = hsv2rgb(hsv_out);
    // Create output image
    const char* q = "\nEnter a color palette image filename: ";
    char* output_filename = create_path("images/output/", q, ".txt");
    write_image_to_file(rgb, output_filename);
    // Clean up
    free(output_filename); output_filename = NULL;
    free_image_hsv(hsv_out);
    free_image_rgb(rgb);
}


void report_color_palette(Color_Palette* cp) {
    printf("\n<======================Color Palette Report===========================>\n\n");
    bool errors_found = false;
    double p_tot = 0;
    for (int i=0; i<cp->N; i++) {
        int h = (int)cp->averages[i].h;
        int s = (int)cp->averages[i].s*100.0;
        int v = (int)cp->averages[i].v*100.0;
        double p = cp->percentages[i]*100.0;
        if (h > 360 || h < 0) {
            errors_found = true;
            printf("ERROR: h was out of bounds. i=%3d,\th=%3d\n", i, h);
        }
        if (s > 100 || v < 0) {
            errors_found = true;
            printf("ERROR: s was out of bounds. i=%3d,\ts=%3d\n", i, s);
        }
        if (v > 100 || v < 0) {
            errors_found = true;
            printf("ERROR: v was out of bounds. i=%3d,\tv=%3d\n", i, v);
        }
        p_tot += p;
    }
    if ((int)p_tot > 100) {
        errors_found = true;
        printf("ERROR: total of percentages was not 100. percentage total: %5.2lf", p_tot);
    }
    if (errors_found) return; // Do not create full report if errors found

    for (int i=0; i<cp->N; i++) {
        int h = (int)cp->averages[i].h;
        int s = (int)(cp->averages[i].s*100.0);
        int v = (int)(cp->averages[i].v*100.0);
        double p = cp->percentages[i]*100.0;
        printf("N=%2d,\tHSV: (%3d,%3d,%3d)\tPercentage of photo: %5.2lf\n", i, h, s, v, p);
    }
}

