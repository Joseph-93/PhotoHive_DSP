#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "color_quantization.h"
#include "image_processing.h"
#include "utilities.h"
#include "debug.h"

#define HUE_NORMALIZER (1.0)/(360.0)

float QUANTITY_WEIGHT = 0.1;
float SATURATION_VALUE_WEIGHT = 0.9;

/******************************************************************************
 * initialize_octree sets up an octree according to hyperparamters and returns
 *   a pointer to the octree.
 *  -h_parts, s_parts, v_parts, and num_grays may not be 0.
 *  -Function does NOT initialize valid_parents array, as it is not known upon
 *   initialization. valid_parents will be NULL.
******************************************************************************/
Octree* initialize_octree(int h_parts, 
                          int s_parts,
                          int v_parts, 
                          double black_thresh,
                          double gray_thresh) {
    int num_grays = v_parts;
    if (h_parts==0 || s_parts==0 || v_parts==0 || num_grays == 0) {
        fprintf(stderr, "ERROR: initialize_octree() expects nonzero values for h_parts, s_parts, v_parts, and num_grays.\n");
        return NULL;
    }

    // Create octree, set all measurement values
    Octree* octree = (Octree*)malloc(sizeof(Octree));
    if (!octree) {
        fprintf(stderr, "ERROR: Malloc() of Octree failed.\n");
        return NULL;
    }
    octree->total_length = h_parts*s_parts*v_parts+num_grays+1;
    octree->num_h = h_parts;
    octree->Lh = 360/h_parts;
    octree->num_s = s_parts;
    octree->Ls = (1-gray_thresh)/s_parts; // length(s) = (max_s - gray_thresh)/(number_s_parts)
    octree->num_v = v_parts;
    octree->Lv = (1-black_thresh)/v_parts; // length(v) = (max-black_thresh)/(number_v_parts)
    octree->num_grays = num_grays;
    octree->black_thresh = black_thresh;
    octree->gray_thresh = gray_thresh;

    // Create Octree Groups
    octree->groups = (Octree_Group*)malloc(octree->total_length*sizeof(Octree_Group));
    if (!octree->groups) {
        fprintf(stderr, "ERROR: Malloc() of Octree_Group failed.\n");
        return NULL;
    }

    // Initialize hue groups from h, s, and v partitions.
    double half_h = octree->Lh/2;
    double s_offs = octree->Ls/2 + gray_thresh;     // s_offs = len(s)/2 + gray_thresh
    double v_offs = octree->Lv/2 + black_thresh;    // v_offs = len(v)/2 + black_thresh
    int i;
    for (int h=0; h<h_parts; h++) {
        for (int s=0; s<s_parts; s++) {
            for (int v=0; v<v_parts; v++) {
                i = h*s_parts*v_parts + s*v_parts + v;
                octree->groups[i].h = h*octree->Lh + half_h;
                octree->groups[i].s = s*octree->Ls + s_offs;
                octree->groups[i].v = v*octree->Lv + v_offs;
                octree->groups[i].id = i;
                octree->groups[i].quantity = 0;
                octree->groups[i].head = NULL;
                octree->groups[i].crossed_zero = false;
            }
        }
    }

    // Initialize gray groups
    double L_gray = (1.0f-black_thresh)/(double)num_grays;
    for (int j=0; j<num_grays; j++) {
        i++;
        octree->groups[i].h = 0;
        octree->groups[i].s = 0;
        octree->groups[i].v = L_gray*j + v_offs;
        octree->groups[i].id = i;
        octree->groups[i].quantity = 0;
        octree->groups[i].head = NULL;
        octree->groups[i].crossed_zero = false;
    }
    
    // Initialize black group
    i++;
    octree->groups[i].h = 0;
    octree->groups[i].s = 0;
    octree->groups[i].v = 0;
    octree->groups[i].id = i;
    octree->groups[i].quantity = 0;
    octree->groups[i].head = NULL;
    octree->groups[i].crossed_zero = false;

    return octree;
}


/*****************************************************************************
 * arm_octree assigns all pixels of an HSV image into the proper groups of the 
 *   octree given.
******************************************************************************/
void arm_octree(Image_HSV* hsv, Octree* octree, int HSV_Linked_List_Size) {
    // Create array of current octree group linked list nodes for use in the initialization.
    HSV_Linked_List** list_curs = (HSV_Linked_List**)calloc(octree->total_length, sizeof(HSV_Linked_List*));
    if (!list_curs) {
        fprintf(stderr, "ERROR: arm_octree() failed to malloc list_curs.");
        return;
    }
    for (int i=0; i<octree->total_length; i++) {
        // Set array values all to the heads of their respective linked lists.
        octree->groups[i].head = get_hsv_linked_list_node(HSV_Linked_List_Size);
        list_curs[i] = octree->groups[i].head;
        // if (list_curs[i] == 0x55c37fef01c0) {
        //     printf("")
        // }
    }

    // Loop through all pixels
    int i_tot = hsv->height * hsv->width;
    int old_g = 0;
    for (int i=0; i<i_tot; i++) {
        // Vi, Si, and Hi are index i trackers for H,S,V values.
        int Vi, Si, Hi, g;
        // If value is under black threshold, set g to the black grouping.
        if (hsv->pixels[i].v < octree->black_thresh){
            g = octree->total_length-1;
        }
        // if current pixel falls under gray groupings, set g to its proper gray grouping.
        else if (hsv->pixels[i].s < octree->gray_thresh) {
            Vi=(int)(hsv->pixels[i].v-octree->black_thresh)*octree->num_grays/(1-octree->black_thresh);
            g = octree->total_length-(octree->num_grays+1)+Vi;
        }
        // if current pixel is not gray or black, set g to the proper color grouping.
        else {
            Vi=((int)((hsv->pixels[i].v-octree->black_thresh)/octree->Lv));
            Si=(int)((hsv->pixels[i].s-octree->gray_thresh)/octree->Ls);
            Hi=(int)(hsv->pixels[i].h/octree->Lh);
            g = (Hi*octree->num_s+Si)*octree->num_v + Vi;
        }
        if (!list_curs[g]) {
            list_curs[g] = get_hsv_linked_list_node(HSV_Linked_List_Size);
        }
        // If the current linked list is filled, allocate a new list and move cur to it.
        else if (list_curs[g]->i == list_curs[g]->array_size) {
            list_curs[g]->next = get_hsv_linked_list_node(HSV_Linked_List_Size);
            list_curs[g] = list_curs[g]->next;
        }
        // The first unfilled pixel in the g-th list_cur gets the current hsv pixel
        // The i inside of list_curs[g] is NOT the same i as the for loop uses!!
        list_curs[g]->pixels[list_curs[g]->i++] = hsv->pixels[i];
        list_curs[g]->num_pixels++;
        octree->groups[g].quantity++;
    }
    free(list_curs); list_curs = NULL;
}


/******************************************************************************
 * find_valid_parents calculates the quantities of pixels in each octree group
 *   and assigns valid_parent=true to all octree groups that are required to 
 *   reach the hyperparameter threshold.
 *  -Function also adds group id to octree->valid_parents array.
 *  -total_pixels is the total number of pixels in the image.
 *  -coverage_threshold is the hyperparameter for the threshold for how many of
 *   the pixels should be a part of the main colors in the color palette.
 *   should be a number in the range [0, 1]
******************************************************************************/
void find_valid_octree_parents(Octree* octree, int total_pixels, double coverage_threshold) {
    // sort an array of ids in order from least to greatest
    int* sorted_ids = malloc(octree->total_length * sizeof(int));
    for (int i=0; i<octree->total_length; i++) {
        sorted_ids[i] = i;
    }

    custom_sort(sorted_ids, octree->total_length, sizeof(int), compare_quantities, octree->groups);

    // Fill out values of valid_parents array until coverage_threshold is hit
    int goal_num_pixels = (int)((double)total_pixels*coverage_threshold);
    int i;
    for (i=0; i<octree->total_length; i++) {
        goal_num_pixels -= octree->groups[sorted_ids[i]].quantity;
        // if goal_num_pixels is reached, set up valid_parents as 0 to i, and return.
        if (goal_num_pixels <= 0) {
            int* valid_parents = (int*)malloc(i*sizeof(int));
            for(int j=0; j<i; j++) {
                valid_parents[j] = sorted_ids[j];
            }
            octree->valid_parents = valid_parents;
            octree->len_valid_parents = i;
            free(sorted_ids);
            return;
        }
    }

    fprintf(stderr, "ERROR: find_valid_octree_parents should not reach the end of its valid_parents loop");
}


/******************************************************************************
 * consolidate_valid_octree_parents brings parents that are very similar into
 *   one. If two parents appear to be the same value, it forces their pixels
 *   into one. It only uses indices and manhatten distance to decide what is
 *   close enough.
******************************************************************************/
void consolidate_valid_octree_parents(Octree* octree) {
    return;
}


/******************************************************************************
 * 
******************************************************************************/
int* get_octree_hsv_coords(Octree* octree, int id) {
    int* coords = (int*)malloc(3*sizeof(int));
    if (id < (octree->total_length - octree->num_grays)) {
        coords[0] = id / (octree->num_s * octree->num_v);
        coords[1] = (id % (octree->num_s * octree->num_v)) / octree->num_v;
        coords[2] = id % octree->num_v;
    }
    else if (id < octree->total_length - 1) { // if id shows node to be gray
        coords[0] = 0;  // Gray values are rounded to have no hue
        coords[1] = 0;  // Gray values are rounded to have no saturation
        coords[2] = id % octree->num_grays + 1; // +1 because black takes up the zero space
    }
    else {  // Otherwise, node must be black
        coords[0] = 0;  // black hue is 0
        coords[0] = 0;  // black saturation is 0
        coords[0] = 0;  // black value is 0
    }
    return coords;
}


/******************************************************************************
 * get_node_distance_heuristic calculates the distance from a group node to a
 *   parent node, with the caveats of the octree grouping system, gray nodes
 *   and the black node.
 *  -Function returns integer distance, similar to the square of the euclidean
 *   distance.
 *  -octree is the octree used.
 *  -groupid is the id attribute of the group in question
 *  -parentid is the id attribute of the parent group in question
 *  -Function converts hue values to the range [0,1] and acts as if the hSV
 *   space was linear, for efficiency.
******************************************************************************/
Pixel get_node_distance_heuristic(Octree* octree, int groupid, int parentid) {
    // Gray_start is the index that starts the gray area
    int gray_start = octree->total_length - (octree->num_grays+1);
    int black_id = octree->total_length-1;

    // Get hsv coordinates
    Octree_Group* group = &octree->groups[groupid];
    Octree_Group* parent =&octree->groups[parentid];

    Pixel dist;
    // if group is color and parent is color
    if (groupid < gray_start && parentid < gray_start) {
        // dist = h^2 + s^2 + v^2
        Pixel h_diff = fabs(group->h - parent->h);
        if (h_diff > 180) h_diff = 360 - h_diff;
        h_diff *= HUE_NORMALIZER;
        Pixel s_diff = group->s - parent->s;
        Pixel v_diff = group->v - parent->v;
        dist = h_diff*h_diff + s_diff*s_diff + v_diff*v_diff;
    }
    // if group is gray and parent is color OR group is color and parent is gray
    else if ((gray_start<=groupid && groupid<black_id && parentid<gray_start) || 
             (gray_start<=parentid && parentid<black_id && groupid<gray_start)) {        
        // dist = s^2 + v^2
        Pixel s_diff = group->s - parent->s;
        Pixel v_diff = group->v - parent->v;
        dist = s_diff*s_diff + v_diff*v_diff;
    }
    // else
    else {
        // dist = v^2
        Pixel v_diff = group->v - parent->v;
        dist = v_diff*v_diff;
    }
    return dist;
}


/******************************************************************************
 * get_distance_pixel_to_parent gets the distance from a pixel to an octree 
 *   parent.
 *  -pixel is the pixel to compare against the octree parent coords.
 *  -parentid is the id of the parent, so it can be located in the octree.
 *  -octree is the structure where the parent resides.
 *  -Function assumes that the hsv space can be approximated linearly, because
 *   pixel and possible parents are assumed to be close, which creates a nearly
 *   linear relationship.
 *  -Function normalizes hue so that no distance above 180 degrees is possible,
 *   because in polar form this is true.
******************************************************************************/
double get_distance_pixel_to_parent(Pixel_HSV* pixel, Octree* octree, int parentid) {
    Octree_Group group = octree->groups[parentid];
    double h = fabs(pixel->h - group.h);
    if (h > 180) h = 360-h;
    h *= HUE_NORMALIZER;
    double s = pixel->s - group.s;
    double v = pixel->v - group.v;
    double dist = h*h + s*s + v*v;
}


/******************************************************************************
 * move_all_pixels_to_new_group moves the entire HSV_Linked_List list from one
 *   octree group to another, updating quantities, cleaning as necessary, and
 *   assuring that all pixels are moved from the from_list to the new list
 *   properly.
 *  -octree is the main octree object pointer
 *  -cur_groups is the linked list that points to the farthest grouping of
 *   HSV_Linked_Lists for each group in octree->groups.
 *  -from_list is the pointer to the head of the HSV_Linked_List that needs
 *   reassignation.
 *  -closest_parent is the id of the octree group that from_list needs to be
 *   appended to. it's the id of the destination.
 *  -group_id is the id of the group that contains from_list.
******************************************************************************/
void move_all_pixels_to_new_group(Octree* octree, HSV_Linked_List** cur_groups, HSV_Linked_List* from_list, int closest_parent, int group_id);


// Check if the difference between group.h and h crosses zero and update answer
void update_crosses_zero(Octree_Group* group, Pixel h) {
    if (abs(group->h - h) > 180) {
        group->crossed_zero = true;
    }
}


/******************************************************************************
 * group_irregular_pixels assigns pixels that are not part of a valid_parent.
******************************************************************************/
void group_irregular_pixels(Octree* octree) {
    // initialize cur_groups to the latest HSV_Linked_List in each group
    HSV_Linked_List** cur_groups = (HSV_Linked_List**)calloc(octree->total_length, sizeof(HSV_Linked_List*));
    for (int i=0; i<octree->total_length; i++) {
        cur_groups[i] = octree->groups[i].head;
        while(cur_groups[i]->next) {
            cur_groups[i] = cur_groups[i]->next;
        }
    }
    // for node in octree:
    for (int i=0; i<octree->total_length; i++) {
        // If the current value is already a valid_parent or is empty
        if (octree->groups[i].quantity == 0) {
            continue;
        }
        bool already_parent = false;
        for (int j=0; j<octree->len_valid_parents; j++) {
            if (octree->valid_parents[j] == i) {
                already_parent = true;
                break;
            }
        }
        if (already_parent) continue;

        // CHECK FOR NEAREST NODES:
        // initialize cur_min_dist to high number
        Pixel cur_min_dist = (Pixel)octree->total_length*octree->total_length; // a real min_dist cannot be higher
        // initialize number_minimums to 0
        int num_mins=0;
        // initialize a parent_distances array
        Pixel* parent_distances = (Pixel*)malloc(octree->len_valid_parents * sizeof(Pixel));
        // for parent in valid_parents:
        for (int j=0; j<octree->len_valid_parents; j++) {
            // calculate squared distance heuristic
            Pixel distance = get_node_distance_heuristic(octree, i, octree->valid_parents[j]);
            // if distance < cur_min_dist:
            if (distance < cur_min_dist) {
                // update cur_min_dist
                cur_min_dist = distance;
                // num_minimums = 1
                num_mins = 1;
            }
            // if distance == cur_min_dist:
            else if (distance == cur_min_dist) {
                // num_minimums++
                num_mins++;
            }
            // update parent_distances array
            parent_distances[j] = distance;
        }
        // Malloc a closest[] array, num_minimums
        int num_closest_parents = num_mins;
        int* closest_parents = (int*)malloc(num_closest_parents*sizeof(int));
        {
            int k = 0;
            // for parent_distance in parent_distances:
            for (int j=0; j<octree->len_valid_parents; j++) {
                // if parent_distance == cur_min_dist
                if (parent_distances[j] == cur_min_dist) {
                    // append valid_parents[j] to closest_parents
                    closest_parents[k++] = octree->valid_parents[j];
                }
            }
        }
        free(parent_distances);

        HSV_Linked_List* from_list = octree->groups[i].head;
        // ASSIGN TO NEAREST
        // if several equally nearest nodes:
        if (num_closest_parents > 1) {
            // USE EUCLIDEAN^2 DISTANCE FOR EACH PIXEL
            // for pixel in node:
            while(from_list) {
                for(int j=0; j<from_list->num_pixels; j++) {
                    // initialize cur_min_dist to something HIGH or something...
                    double cur_min_dist = (double)octree->total_length;
                    // initialize cur_min_parent_id to 0
                    int cur_min_parent_id = 0;
                    // for parent in valid_parents
                    for(int k=0; k<num_closest_parents; k++) {
                        // euc_dist = (parentx - pixelx)^2 + (parenty-pixely)^2
                        double parent_id = closest_parents[k];
                        Pixel_HSV* pixel = &from_list->pixels[j];
                        double distance = get_distance_pixel_to_parent(pixel, octree, parent_id);
                        // if euc_dist < cur_min_dist
                        if (distance < cur_min_dist) {
                            // update cur_min_dist
                            cur_min_dist = distance;
                            // update cur_min_parent_id
                            cur_min_parent_id = parent_id;
                        }
                    }
                    // assign pixel to the octree->groups[cur_min_parent_id]
                    HSV_Linked_List* cur_group = cur_groups[cur_min_parent_id];
                    // If the current linked list is filled, allocate a new list and move cur to it.
                    if (cur_group->i == cur_group->array_size) {
                        cur_group->next = get_hsv_linked_list_node(cur_group->array_size);
                        cur_group = cur_group->next;
                    }
                    // Add pixel to proper cur_group, increment counters
                    cur_group->pixels[cur_group->i++] = from_list->pixels[j];
                    cur_group->num_pixels++;
                    octree->groups[cur_min_parent_id].quantity++;
                    // Check if the pixel update causes the crosses zero, unless it already has
                    update_crosses_zero(&octree->groups[cur_min_parent_id], from_list->pixels[j].h);
                }
                from_list = from_list->next;
            }
        }

        // if there is only one nearest node
        else {
            int closest_parent = closest_parents[0];
            update_crosses_zero(&octree->groups[closest_parent], octree->groups[i].h);
            // every pixel goes to the nearest node
            HSV_Linked_List* cur_group = cur_groups[closest_parent];
            cur_group->num_pixels = cur_group->i;
            if (from_list->num_pixels != 0) {
                cur_group->next = from_list;
            }
            // valid_parent.quantity += node.quantity
            octree->groups[closest_parent].quantity += octree->groups[i].quantity;
            // node.quantity = 0
            octree->groups[i].quantity = 0;
            // get cur_group to the very back of cur_group+from_list
            while(cur_groups[closest_parent]->next) {
                cur_groups[closest_parent] = cur_groups[closest_parent]->next;
            }
            // empty octree->groups[i] because it's been transferred to cur_group
            from_list = NULL;
            octree->groups[i].head = NULL;
            // cur_groups[i] = NULL;
        }
        free(closest_parents);
    }
    free(cur_groups);
}


/******************************************************************************
 * get_hsv_linked_list_node creates a new HSV_Linked_List with an array of 
 *   arr_size. if it is unable to allocate memory, it prints an error to stderr
 *   and returns NULL.
 *  -arr_size > 0 is required. if arr_size < 0, Function returns an stderr and
 *   a NULL pointer.
******************************************************************************/
HSV_Linked_List* get_hsv_linked_list_node(int arr_size) {
    if (arr_size == 0) {
        fprintf(stderr, "ERROR: get_hsv_linked_list_node() expected a positive arr_size");
        return NULL;
    }
    HSV_Linked_List* node = (HSV_Linked_List*)calloc(1, sizeof(HSV_Linked_List));
    if (!node) {
        fprintf(stderr, "ERROR: failed to calloc node in get_hsv_linked_list_node()");
        return NULL;
    }
    node->num_pixels = 0;
    node->array_size = arr_size;
    node->pixels = (Pixel_HSV*)calloc(arr_size, sizeof(Pixel_HSV));
    if (!node->pixels) {
        fprintf(stderr, "ERROR: failed to calloc pixels array in get_hsv_linked_list_node()");
        return NULL;
    }
    return node;
}


Color_Palette* calculate_avg_hsv(Octree* octree, Image_HSV* hsv) {
    // Initialize Color Palette
    Color_Palette* cp = (Color_Palette*)malloc(sizeof(Color_Palette));
    cp->N = octree->len_valid_parents;
    cp->averages = (Pixel_HSV*)malloc(cp->N*sizeof(Pixel_HSV));
    cp->percentages = (Pixel*)malloc(cp->N*sizeof(Pixel));

    // Get 1/(number of pixels) for fast multiplication
    Pixel inverse_of_quantity = 1.0/(hsv->height * hsv->width);

    // Iterate through valid parents of octree
    for (int i=0; i<octree->len_valid_parents; i++) {
        int id = octree->valid_parents[i];
        bool crossed_zero = octree->groups[id].crossed_zero;
        HSV_Linked_List* list = octree->groups[id].head;
        Pixel h_total = 0, s_total=0, v_total=0;
        int total_pixels = 0;
        // Iterate through all linked lists
        while(list) {
            // accumulate list's number of pixels to the total pixels
            total_pixels += list->num_pixels;
            // Iterate through pixels of linked lists
            for(int j=0; j<list->num_pixels; j++) {
                // Get totals of the h,s,v values.
                // NOTE: Use large data types, else precision will be lost here.
                if (crossed_zero && list->pixels[j].h > 180) h_total -= 360;
                h_total += list->pixels[j].h;
                s_total += list->pixels[j].s;
                v_total += list->pixels[j].v;
            }
            list = list->next;
        }
        // divide h,s,v totals by the number of pixels (calculate averages)
        // Store values in color palette structure
        double list_inverse_of_quantity = 1.0/(double)total_pixels;
        cp->averages[i].h = h_total * list_inverse_of_quantity;
        if (crossed_zero && cp->averages[i].h < 0) {cp->averages[i].h += 360;}
        cp->averages[i].s = s_total * list_inverse_of_quantity;
        cp->averages[i].v = v_total * list_inverse_of_quantity;
        cp->percentages[i] = (double)total_pixels * inverse_of_quantity;
    }
    // return the color palette
    return cp;
}


/******************************************************************************
 * saliency measures a heuristic of how "salient" a color in an image is. 
 *  -The heuristic is based on Saturation, Value, and quantity of pixels.
 *  -The equation is as follows:
 *      h(x,y)=If(x<THRESH, Y_WEIGHT*y, Y_WEIGHT*y + S_V_WEIGHT*(x-THRESH))
 *      where g is saliency, x is S*V, and y is quantity,
 *      Y_WEIGHT is QUANTITY WEIGHT, S_V_WEIGHT is SATURATION_VALUE_WEIGHT
 *      Visual of the model: https://www.geogebra.org/calculator/vuabsjrb
******************************************************************************/
float saliency(const Octree_Group* group) {
    float saliency;
    float s_v = group->s * group->v;

    saliency = (float)group->quantity*(QUANTITY_WEIGHT+SATURATION_VALUE_WEIGHT*s_v);

    return saliency*1000; // simply multiply by 1000 to reduce issues in truncation
}


/******************************************************************************
 * compare_quantity is used to compare the quantities of two octree groups.
******************************************************************************/
int compare_quantities(const void *a, const void *b, void *arg) {
    const int index_a = *(const int *)a;
    const int index_b = *(const int *)b;
    const Octree_Group *groups = (const Octree_Group *)arg;

    float saliency_a = saliency(&groups[index_a]);
    float saliency_b = saliency(&groups[index_b]);

    return (int)(saliency_b - saliency_a);
    // return groups[index_b].quantity - groups[index_a].quantity;
}


void free_color_palette(Color_Palette* cp) {
    free(cp->averages); cp->averages = NULL;
    free(cp->percentages); cp->percentages = NULL;
    free(cp); cp = NULL;
}


void free_octree(Octree* octree) {
    if (octree == NULL) return;

    // Free valid_parents
    free(octree->valid_parents);

    // Free each group and its linked lists
    for (int i = 0; i < octree->total_length; i++) {
        HSV_Linked_List* current_list = octree->groups[i].head;
        int j = -1;
        while (current_list != NULL) {
            j++;
            HSV_Linked_List* next_list = current_list->next;

            printf("------->(%d,%d)", i, j);
            printf("%p,\t%d\t", current_list->pixels, current_list->num_pixels);

            printf("next: %p", current_list->next);
            printf("<--------\n");

            // Free the pixels array of the current linked list node
            if (i!=69 || false) {
                printf("frick at least I did it...\n");
                free(current_list->pixels);
                free(current_list);
            }

            // Free the current linked list node itself

            current_list = next_list;
        }
    }

    // Free the groups array itself
    free(octree->groups);

    // Finally, free the octree structure
    free(octree);
}


Color_Palette* get_color_palette(Image_HSV* hsv,
                                 const int linked_list_size,
                                 const double coverage_threshold,
                                 const int h_parts, const int s_parts, const int v_parts,
                                 const double black_thresh, const double gray_thresh,
                                 float quantity_weight,
                                 float saturation_value_weight) {

    // Set up constants
    QUANTITY_WEIGHT = quantity_weight;
    SATURATION_VALUE_WEIGHT = saturation_value_weight;

    // Initialize Octree
    Octree* octree = initialize_octree(h_parts, s_parts, v_parts, black_thresh, gray_thresh);

    // Put hsv pixels into octree structure
    arm_octree(hsv, octree, linked_list_size);

    // Find the octree parents required to make up the COVERAGE_THRESHOLD
    find_valid_octree_parents(octree, hsv->height*hsv->width, coverage_threshold);

    // Group un-parented pixels to their nearest parents
    group_irregular_pixels(octree);

    // Create a color palette
    Color_Palette* cp = calculate_avg_hsv(octree, hsv);

    // Clean up
    free_octree(octree);

    // Return color palette
    return cp;
}