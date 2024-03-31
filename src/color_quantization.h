#include "types.h"
#include "stdbool.h"
#include "image_processing.h"
#ifndef COLOR_QUANTIZATION_H
#define COLOR_QUANTIZATION_H


/******************************************************************************
 * Color_Palette is the deliverable object for the measured hsv color palette.
******************************************************************************/
typedef struct Color_Palette {
    int N;
    Pixel_HSV* averages;
    double* percentages;
} Color_Palette;


/******************************************************************************
 * HSV_Linked_List is made to be pointed to by the Octree_Group. It allows the
 *   the octree group to contain an undefined number of pixels for relatively
 *   low computation and memory. 
 *  -When the end of the linked list is reached, the user must allocate a new
 *   HSV_Linked_List instance, point to it from this instance, and move the
 *   pointer to the next instance.
 *  -pixels is the array object
 *  -i is the current iteration of the array which has NOT yet been filled.
 *  -next is a pointer to the next linked list node
 *  -int num_pixels is the number of pixels the array is to hold
 *  -int array_size is the total amount of memory in the array.
 *  -The difference between num_pixels and array_size is that num_pixels is
 *   used to keep track of pixels in an image, and if the array is not all used
 *   num_pixels keeps track of the functional end of the array, whereas
 *   array_size keeps track of the total size of the array, so that all of its
 *   memory can be properly freed.
******************************************************************************/
typedef struct HSV_Linked_List HSV_Linked_List;

struct HSV_Linked_List {
    Pixel_HSV* pixels;
    int num_pixels;
    int array_size;
    int i;
    HSV_Linked_List* next;
};


/******************************************************************************
 * Octree_Group is the struct that Octree.groups points to.
 *  -id is self-explanatory
 *  -quantity is the number of pixels that belong to the group.
 *  -h, s, v are designed to be the "average" color described by the nodes
 *   within the group. At first it is automatically assigned, but later must
 *   be set to the average H,S,V values of its children.
 *  -head is the pointer to the HSV_Linked_List head.
 *  -is_valid_parent should default to false, but is set to true if the group
 *   is decided to be a valid parent.
 *  -crossed_zero is used to track what colors allowed consolidation across the
 *   H=0 boundary, so that averaging them can be done properly.
******************************************************************************/
typedef struct Octree_Group {
    int id;
    int quantity;
    double h, s, v;
    HSV_Linked_List* head;
    HSV_Linked_List* cur;
    bool is_valid_parent;
    bool crossed_zero;
} Octree_Group;


/******************************************************************************
 * Octree is necessary for keeping track of octree groups contained.
 *  -valid_parents: an array of valid parent_id values. valid parents are
 *   parents which are part of the color palette.
 *  -groups: pointer to an array of all octree groups
 *  -Lh: length of each normal (non gray or black) group in hue dimension.
 *  -Ls: length of each normal (non gray or black) group in saturation dimension.
 *  -Lv: length of each normal (non gray or black) group in value dimension.
 *  -num_h: number of normal (non gray or black) group partitions in hue
 *   dimension. In other words, it is the maximum h index available.
 *  -num_s: number of normal (non gray or black) group partitions in saturation
 *   dimension. In other words, it is the maximum s index available.
 *  -num_v: number of normal (non gray or black) group partitions in value
 *   dimension. In other words, it is the maximum v index available.
 *  -total_length: the total length of the groups array, including gray and
 *   black. always set total_length = Lh*Ls*Lv + num_grays.
 *  -black_thresh: the value which automatically puts a pixel into the black
 *   octree groups.
 *  -gray_thresh: the saturation which automatically puts a pixel into the gray
 *   octree groups, so long as said pixel's value > black_thresh.
 *  -NOTE: There is always a black node, which will always be located at the
 *   back of the array.
 *  -NOTE: Group indexing is a complex task. For reference, to reference the
 *   i-th index in normal groups, i = h_ind*(Ls*Lv) + s_ind*Lv + v_ind.
 *   To index gray groups, i=Lh*Ls*Lv+gray_i OR i=total_length-num_grays+gray_i
 *   To index the black node, i = Lh*Ls*Lv + num_grays OR i=total_length-1.
******************************************************************************/
typedef struct Octree {
    Octree_Group* groups;
    int* valid_parents;
    int len_valid_parents;
    double Lh, Ls, Lv;
    int num_h, num_s, num_v;
    int num_grays;
    int total_length;
    Pixel black_thresh, gray_thresh;
} Octree;


Octree* initialize_octree(int h_parts, 
                          int s_parts,
                          int v_parts, 
                          double black_thresh,
                          double gray_thresh);


void arm_octree(Image_HSV* hsv, Octree* octree, int HSV_Linked_List_Size);

void find_valid_octree_parents(Octree* octree, int total_pixels, double coverage_threshold);

void consolidate_valid_octree_parents(Octree* octree);

int* get_octree_hsv_coords(Octree* octree, int id);

Pixel get_node_distance_heuristic(Octree* octree, int groupid, int parentid);

double get_distance_pixel_to_parent(Pixel_HSV* pixel, Octree* octree, int parentid);

void move_all_pixels_to_new_group(Octree* octree,
                                  HSV_Linked_List** cur_groups, 
                                  HSV_Linked_List* from_list, 
                                  int closest_parent, int group_id);

void group_irregular_pixels(Octree* octree);

HSV_Linked_List* get_hsv_linked_list_node(int arr_size);

Color_Palette* calculate_avg_hsv(Octree* octree, Image_HSV* hsv);

int compare_quantities(const void *a, const void *b, void *arg);

void free_color_palette(Color_Palette* cp);

void free_octree(Octree* octree);

float saliency(const Octree_Group* group);

Color_Palette* get_color_palette(Image_HSV* hsv,
                                 const int linked_list_size,
                                 const double coverage_threshold,
                                 const int h_parts, const int s_parts, const int v_parts,
                                 const double black_thresh, const double gray_thresh,
                                 float quantity_weight,
                                 float saturation_value_weight);

#endif