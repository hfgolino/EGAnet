#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>

// Structure for `modularity_values`
struct modularity_result {
    double* positive_modularity_values;
    double* negative_modularity_values;
    double positive_sum_flag;
    double negative_sum_flag;
    double positive_total_sum;
    double negative_total_sum;
};

// Function to compute modularity values
struct modularity_result modularity_values(double* network, int cols) {

    // Initialize iterators
    int i, j, network_offset;
    double edge;
    int count = 0;

    // Initialize sums
    double* positive_column_sums = (double*)calloc(cols, sizeof(double));
    double* negative_column_sums = (double*)calloc(cols, sizeof(double));
    double positive_sum = 0.0, negative_sum = 0.0;

    // Loop over to get sums
    for(i = 0; i < cols; i++) {

        // Compute network offset
        network_offset = i * cols;

        for(j = i; j < cols; j++) {

            // Get edge
            edge = network[network_offset + j];

            // Compute based on sign
            if(edge > 0) {

                // Add to sums
                positive_column_sums[i] += edge;

                // Check for off-diagonal edge
                if(i != j) {
                    positive_column_sums[j] += edge;
                }

            }else if(edge < 0) {

                // Add to sums
                negative_column_sums[i] += edge;

                // Check for off-diagonal edge
                if(i != j) {
                    negative_column_sums[j] += edge;
                }

            }

        }
    }

    // Compute sums
    for(i = 0; i < cols; i++) {
        positive_sum += positive_column_sums[i];
        negative_sum += negative_column_sums[i];
    }

    // Compute total sum
    double total_sum = positive_sum + negative_sum;

    // Get indices for modularity values
    int value_indices = (((cols * (cols - 1)) / 2) + cols);

    // Initialize modularity values
    double* positive_modularity_values = (double*)calloc(value_indices, sizeof(double));
    double* negative_modularity_values = (double*)calloc(value_indices, sizeof(double));

    // Check if positive_sum and negative_sum are not zero
    int positive_sum_flag = positive_sum != 0;
    int negative_sum_flag = negative_sum != 0;

    // Loop over to compute modularity
    for (i = 0; i < cols; i++) {

        // Compute network offset
        network_offset = i * cols;

        for (j = i; j < cols; j++) {

            // Obtain edges
            edge = network[network_offset + j];

            // Check for positive sum
            if(positive_sum_flag) {

                // Update positive modularity
                positive_modularity_values[count] += (
                    ((edge > 0) ? edge : 0) - // positive edge
                    positive_column_sums[i] *
                    positive_column_sums[j] /
                    positive_sum
                ) / positive_sum;

            }

            // Check for negative sum
            if(negative_sum_flag) {

                // Update negative modularity
                negative_modularity_values[count] += (
                    ((edge < 0) ? edge : 0) - // negative edge
                    negative_column_sums[i] *
                    negative_column_sums[j] /
                    negative_sum
                ) / negative_sum;

            }

            // Increase count
            count++;
        }
    }

    // Free memory
    free(positive_column_sums);
    free(negative_column_sums);

    // Set pointers to NULL
    positive_column_sums = NULL;
    negative_column_sums = NULL;

    // Pre-compute division of sums
    double positive_total_sum = positive_sum / total_sum;
    double negative_total_sum = negative_sum / total_sum;

    // Set up result
    struct modularity_result result = {
        positive_modularity_values,
        negative_modularity_values,
        positive_sum_flag,
        negative_sum_flag,
        positive_total_sum,
        negative_total_sum
    };

    // Return result
    return(result);

}

// Signed modularity function
double signed_modularity(struct modularity_result Q_values, int* membership, int cols) {

    // Initialize iterators
    int i, j;
    int count = 0;

    // Initialize positive and negative modularity
    double Q_positive = 0.0, Q_negative = 0.0;

    // Initialize positive and negative value
    double positive_value, negative_value;

    // Loop over to compute modularity
    for (i = 0; i < cols; i++) {

        for (j = i; j < cols; j++) {

            // Check for memberships
            if (membership[i] == membership[j]) {

                // Values
                positive_value = ((Q_values.positive_sum_flag) ? Q_values.positive_modularity_values[count] : 0);
                negative_value = ((Q_values.negative_sum_flag) ? Q_values.negative_modularity_values[count] : 0);

                // Update modularity
                Q_positive += positive_value * 2;
                Q_negative += negative_value * 2;

                // Check for diagonal edge
                if (i == j) {
                    Q_positive -= positive_value;
                    Q_negative -= negative_value;
                }

            }

            // Increase count
            count++;
        }

    }

    // Return modularity
    return Q_values.positive_total_sum * Q_positive - Q_values.negative_total_sum * Q_negative;

}

// Function to compute modularity gain
double modularity_gain(
    struct modularity_result Q_values,
    int target_node,
    int target_membership, int neighbor_membership,
    int* membership, int* new_memberships,
    double update_modularity, int cols
) {

    // Initialize iterators
    int i, index;

    // Initialize positive and negative value
    double positive_value, negative_value;

    // Initialize gain in modularity
    double gain_positive = 0.0, gain_negative = 0.0;

    // Initialize row start
    int row_start = 0;

    // Compute row start
    for(i = 0; i < target_node; i++) {
        row_start += (cols - i);
    }

    // Compute row end
    int row_end = row_start + (cols - target_node);

    // Initialize count
    int count = 0;

    // Loop over to get values
    for(i = row_start; i < row_end; i++) {

        // Get index
        index = target_node + count;

        // Check for memberships
        if (target_membership == membership[index]) {

            // Values
            positive_value = ((Q_values.positive_sum_flag) ? Q_values.positive_modularity_values[i] : 0);
            negative_value = ((Q_values.negative_sum_flag) ? Q_values.negative_modularity_values[i] : 0);

            // Update modularity
            gain_positive -= positive_value * 2;
            gain_negative -= negative_value * 2;

            // Check for diagonal edge
            if (i == row_start) {
                gain_positive += positive_value;
                gain_negative += negative_value;
            }

        }

        if(neighbor_membership == new_memberships[index]) {

            // Values
            positive_value = ((Q_values.positive_sum_flag) ? Q_values.positive_modularity_values[i] : 0);
            negative_value = ((Q_values.negative_sum_flag) ? Q_values.negative_modularity_values[i] : 0);

            // Update modularity
            gain_positive += positive_value * 2;
            gain_negative += negative_value * 2;

            // Check for diagonal edge
            if (i == row_start) {
                gain_positive -= positive_value;
                gain_negative -= negative_value;
            }

        }

        // Increase count
        count++;

    }

    // Check for whether target node equals zero
    if(target_node != 0) {

        // Initialize target index
        int column_index = target_node;

        // Reset count
        count = target_node;

        // Loop over to get columns
        for(i = 0; i < target_node; i++) {

            // Get index
            index = target_node - count;

            // Check for memberships
            if (target_membership == membership[index]) {

                // Values
                positive_value = ((Q_values.positive_sum_flag) ? Q_values.positive_modularity_values[column_index] : 0);
                negative_value = ((Q_values.negative_sum_flag) ? Q_values.negative_modularity_values[column_index] : 0);

                // Update modularity
                gain_positive -= positive_value * 2;
                gain_negative -= negative_value * 2;

            }

            if(neighbor_membership == new_memberships[index]) {

                // Values
                positive_value = ((Q_values.positive_sum_flag) ? Q_values.positive_modularity_values[column_index] : 0);
                negative_value = ((Q_values.negative_sum_flag) ? Q_values.negative_modularity_values[column_index] : 0);

                // Update modularity
                gain_positive += positive_value * 2;
                gain_negative += negative_value * 2;

            }

            // Update target index
            column_index += (cols - i - 1);

            // Decrease count
            count--;

        }

    }

    // Return modularity
    return Q_values.positive_total_sum * gain_positive - Q_values.negative_total_sum * gain_negative;

}

/* Functions for Main Louvain */

// Structure for `make_higher_order`
struct higher_order_result {
    double* higher_order;
    int number;
};

// Function to make higher-order network
struct higher_order_result make_higher_order(double* network, int* membership, int cols) {

    // Initialize iterators
    int i, j, network_offset;
    double edge;
    int first_position, second_position;

    // Initialize number of communities
    int number = 0;

    // Loop over to find highest value
    for(i = 0; i < cols; i++) {

        // Check for new highest value
        if(membership[i] > number) {
            number = membership[i];
        }

    }

    // Increase number by 1
    number += 1;

    // Initialize higher-order network
    double* higher_order = (double*)calloc(number * number, sizeof(double));

    // Populate network
    for(i = 0; i < cols; i++) {

        // Compute network offset
        network_offset = i * cols;

        for(j = i; j < cols; j++) {

            // Obtain edge
            edge = network[network_offset + j];

            // Obtain positions
            first_position = membership[i] * number + membership[j];
            second_position = membership[j] * number + membership[i];

            // Add edge to current position of memberships
            higher_order[first_position] += edge;
            higher_order[second_position] += edge;

        }
    }

    // Set up result
    struct higher_order_result result = {
        higher_order,
        number
    };

    // Return result
    return(result);

}

// Function to re-index membership
int* reindex_membership(int* membership, int cols) {

    // Initialize iterators
    int i;
    int count = 0;

    // Initialize index array
    int* index_array = (int*)calloc(cols, sizeof(int));

    // Loop over to determine index correspondence
    for(i = 0; i < cols; i++) {

        // Check if index already exists
        if(index_array[membership[i]] == 0) {

            // Increase count
            count++;

            // Set index array
            index_array[membership[i]] = count;

        }

    }

    // Update membership
    for(i = 0; i < cols; i++) {

        // Check value in index array
        membership[i] = index_array[membership[i]] - 1;

    }

    // Free memory
    free(index_array);

    // Return membership
    return(membership);

}

// Get clock time in nanoseconds
uint64_t get_time_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

// Fisher-Yates (or Knuth) Shuffle
void shuffle_nodes(int *arr, int cols) {

    // Initialize iterators
    int i, j, temp;

    // Get clock time
    uint64_t current_time_ns = get_time_ns();

    // Seed the random number generator
    srand(current_time_ns);

    // Iterate through the array from the last element to the first
    for (i = cols - 1; i > 0; i--) {

        // Generate a random index between 0 and i (inclusive)
        j = rand() % (i + 1);

        // Swap the current element with the randomly selected element
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;

    }

}

// Main Louvain function
void main_louvain(
    double* network,
    struct modularity_result Q_values,
    int* membership_copy, double previous_modularity,
    int cols, int original_cols
) {

    // Initialize iterators
    int i, j, order, network_offset;

    // Initialize membership
    int* membership = (int*)malloc(cols * sizeof(int));
    int* index = (int*)malloc(cols * sizeof(int));

    // Populate membership
    for(i = 0; i < cols; i++) {
        membership[i] = i;
        index[i] = i;
    }

    // Initialize update to modularity
    double update_modularity = signed_modularity(Q_values, membership, cols);

    // Initialize commonly used values
    int target_membership, neighbor_membership, new_membership;
    double best_increase;

    // Initialize higher order structure
    struct higher_order_result higher_order;

    // Initialize new memberships
    int* new_memberships = (int*)malloc(cols * sizeof(int));

    // `while` until there is no gain
    while(1) {

        // Permutate index
        shuffle_nodes(index, cols);

        // Set gain to zero
        int gain = 0;

        // Loop over nodes
        for(i = 0; i < cols; i++) {

            // Initialize gain vector
            double* gain_vector = (double*)calloc(cols, sizeof(double));

            // Reset best increase
            best_increase = 0.0;

            // Get index order
            order = index[i];

            // Obtain membership of target node
            target_membership = membership[order];

            // Initialize new membership
            new_membership = target_membership;

            // Set 1D network offset
            network_offset = order * cols;

            // Loop over nodes
            for(j = 0; j < cols; j++){

                // Check if node is a neighbor
                if(network[network_offset + j] != 0) {

                    // Get neighbor's membership
                    neighbor_membership = membership[j];

                    // Check whether gain is zero
                    if(gain_vector[neighbor_membership] == 0 && target_membership != neighbor_membership) {

                        // Make copy of membership
                        memcpy(new_memberships, membership, cols * sizeof(int));

                        // Update with neighbor membership
                        new_memberships[order] = neighbor_membership;

                        // Compute gain
                        gain_vector[neighbor_membership] = modularity_gain(
                            Q_values, order, target_membership, neighbor_membership,
                            membership, new_memberships, update_modularity, cols
                        );

                        // Check for increase
                        if(gain_vector[neighbor_membership] > best_increase) {

                            // Update best increase
                            best_increase = gain_vector[neighbor_membership];

                            // Make copy of membership
                            new_membership = neighbor_membership;

                        }

                    }

                }

            }

            // Check if there was an increase
            if(best_increase > 0) {

                // Update membership
                membership[order] = new_membership;

                // Update modularity
                update_modularity += best_increase;

                // Reset gain
                gain = 1;

            }

            // Free memory
            free(gain_vector);

        }

        // Check for higher-order network
        if(previous_modularity != 0) {

            // Check if modularity has changed
            if(previous_modularity == update_modularity) {

                // Set gain to zero (break)
                gain = 0;

            } else {

                // Reset previous modularity
                previous_modularity = update_modularity;

                // Re-index membership
                membership = reindex_membership(membership, cols);

                // Update previous membership
                for(i = 0; i < cols; i++) {
                    for(j = 0; j < original_cols; j++) {

                        // Update each membership
                        if(membership_copy[j] == i) {
                            membership_copy[j] = membership[i];
                        }

                    }
                }

                // Obtain new higher-order
                higher_order = make_higher_order(
                    network, membership, cols
                );

                // Update network and number of membership
                memcpy(network, higher_order.higher_order, higher_order.number * higher_order.number * sizeof(double));
                cols = higher_order.number;

                // Free memory
                free(higher_order.higher_order);
                higher_order.higher_order = NULL;

                // Reallocate memory for membership
                membership = (int*)realloc(membership, cols * sizeof(int));
                index = (int*)realloc(index, cols * sizeof(int));

                // Populate membership
                for(i = 0; i < cols; i++) {
                    membership[i] = i;
                    index[i] = i;
                }

                // Update modularity values
                Q_values = modularity_values(network, cols);

                // Update modularity
                update_modularity = signed_modularity(Q_values, membership, cols);

            }

        }

        // Check gain for break condition
        if(gain == 0){
            break;
        }

    }

    // Check for level
    if(previous_modularity == 0) {


        // Check for unidimensionality
        if(cols == 1) {

            // Index membership
            for(i = 0; i < original_cols; i++){
                membership_copy[i] = 0;
            }

        } else {

            // Re-index membership
            membership = reindex_membership(membership, original_cols);

            // Copy to membership copy
            memcpy(membership_copy, membership, original_cols * sizeof(int));

        }

    } else {

        // Check for unidimensionality
        if(cols == 1) {

            // Index membership
            for(i = 0; i < original_cols; i++){
                membership_copy[i] = 0;
            }

        } else {

            // Re-index membership
            membership_copy = reindex_membership(membership_copy, original_cols);

        }

        // Free memory
        free(Q_values.positive_modularity_values);
    	free(Q_values.negative_modularity_values);
    	Q_values.positive_modularity_values = NULL;
        Q_values.negative_modularity_values = NULL;

    }

    // Free memory
    free(index);
    free(membership);
    free(new_memberships);
    index = NULL;
    membership = NULL;
    new_memberships = NULL;

}

/* Functions for Signed Louvain */

// Structure for `signed_louvain`
struct signed_louvain_result {
    int** memberships;
    double* modularities;
    int level;
};

// Signed Louvain function
struct signed_louvain_result signed_louvain(double* original_network, int original_cols) {

    // Initialize iterators
    int level = 0;

    // Initialize membership and modularity results
    int** memberships = (int**)malloc(original_cols * sizeof(int*));
    int* membership_copy = (int*)calloc(original_cols, sizeof(int));
    double* modularities = (double*)malloc(original_cols * sizeof(double));

    // Compute original modularity values
    struct modularity_result original_Q_values = modularity_values(original_network, original_cols);

    // Compute first level
    main_louvain(
        original_network, original_Q_values,
        membership_copy, 0.0,
        original_cols, original_cols
    );

    // Compute modularity
    modularities[level] = signed_modularity(
        original_Q_values, membership_copy, original_cols
    );

    // Populate membership
    memberships[level] = (int*)malloc(original_cols * sizeof(int));
    memcpy(memberships[level], membership_copy, original_cols * sizeof(int));

    // Initialize higher order and modularity values
    struct higher_order_result higher_order;
    struct modularity_result Q_values;

    // Compute next levels
    while(1) {

        // Increase level
        level++;

        // Obtain higher-order network
        higher_order = make_higher_order(
            original_network, membership_copy, original_cols
        );

        // Compute modularity values
        Q_values = modularity_values(
            higher_order.higher_order, higher_order.number
        );

        // Compute next level
        main_louvain(
            higher_order.higher_order, Q_values,
            membership_copy, modularities[level - 1],
            higher_order.number, original_cols
        );

        // Compute modularity
        modularities[level] = signed_modularity(
            original_Q_values, membership_copy, original_cols
        );

        // Populate membership
        memberships[level] = (int*)malloc(original_cols * sizeof(int));
        memcpy(memberships[level], membership_copy, original_cols * sizeof(int));

        // Free memory
        free(Q_values.positive_modularity_values);
        free(Q_values.negative_modularity_values);
        free(higher_order.higher_order);

        // Set pointers to NULL
        Q_values.positive_modularity_values = NULL;
        Q_values.negative_modularity_values = NULL;
        higher_order.higher_order = NULL;

        // Check for any change
        if(modularities[level] <= modularities[level - 1]) {
            break;
        }

    }

    // Free memory
    free(original_Q_values.positive_modularity_values);
    free(original_Q_values.negative_modularity_values);
    free(membership_copy);

    // Set pointers to NULL
    original_Q_values.positive_modularity_values = NULL;
    original_Q_values.negative_modularity_values = NULL;
    membership_copy = NULL;

    // Set up result
    struct signed_louvain_result result = {
        memberships,
        modularities,
        level
    };

    // Return result
    return(result);

}

SEXP r_signed_louvain(SEXP r_input_network) {

    // Initialize iterators
    int i, j;

    // Obtain columns
    int cols = ncols(r_input_network);

    // Call the C function
    struct signed_louvain_result c_result = signed_louvain(REAL(r_input_network), cols);

    // Create R output list
    SEXP r_output = PROTECT(allocVector(VECSXP, 2));

    // Convert the memberships result to an R matrix
    SEXP r_memberships = PROTECT(allocMatrix(INTSXP, c_result.level, cols));
    for (i = 0; i < c_result.level; i++) {
        for (j = 0; j < cols; j++) {
            INTEGER(r_memberships)[j * c_result.level + i] = c_result.memberships[i][j] + 1;
        }
    }

    // Convert the modularities result to an R numeric vector
    SEXP r_modularities = PROTECT(allocVector(REALSXP, c_result.level));
    for (i = 0; i < c_result.level; i++) {
        REAL(r_modularities)[i] = c_result.modularities[i];
        free(c_result.memberships[i]);
    }

    // Set output list elements
    SET_VECTOR_ELT(r_output, 0, r_memberships);
    SET_VECTOR_ELT(r_output, 1, r_modularities);

    // Assign names to the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("memberships"));
    SET_STRING_ELT(names, 1, mkChar("modularities"));
    setAttrib(r_output, R_NamesSymbol, names);

    // Free memory
    free(c_result.memberships);
    free(c_result.modularities);

    // Release protected SEXP objects
    UNPROTECT(4);

    // Return the result
    return r_output;

}
