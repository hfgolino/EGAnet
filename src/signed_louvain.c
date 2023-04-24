#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

/* Function for signed modularity */

// Signed modularity
double signed_modularity(double* network, int *membership, int cols) {

    // Initialize iterators and index
    int i, j;
    double edge; // Initialize edge value

    // Initialize return values
    double *positive_column_sum = (double *)calloc(cols, sizeof(double));
    double *negative_column_sum = (double *)calloc(cols, sizeof(double));

    // Initialize matrix sums
    double positive_sum = 0, negative_sum = 0;

    // Loop over network to obtain matrices
    for (i = 0; i < cols; i++) {
        for (j = 0; j < cols; j++) {

            // Obtain edge
            edge = network[i * cols + j];

            // Check for positive value
            if (edge > 0) {

                // Add to positive sums
                positive_column_sum[j] += edge;
                positive_sum += edge;


            } else if (edge < 0) { // Check for negative value

                // Add to negative sums
                negative_column_sum[j] += edge;
                negative_sum += edge;

            }

        }
    }

    // Compute total sum
    double total_sum = positive_sum + negative_sum;

    // Initialize positive and negative modularity
    double Q_positive = 0;
    double Q_negative = 0;

    // Loop over matrices
    for (i = 0; i < cols; i++) {
        for (j = 0; j < cols; j++) {

            // Obtain edge
            edge = network[i * cols + j];

            // Check for positive values in network
            if (positive_sum != 0) {

                // Check if memberships match, if yes, then update positive modularity
                if (membership[i] == membership[j]) {

                    Q_positive += (
                        ((edge > 0) ? edge : 0) -
                        positive_column_sum[i] *
                        positive_column_sum[j] /
                        positive_sum
                    ) / positive_sum;

                }

            }

            // Check for negative values in network
            if (negative_sum != 0) {

                // Check if memberships match, if yes, then update negative modularity
                if (membership[i] == membership[j]) {

                    Q_negative += (
                        ((edge < 0) ? edge : 0) -
                        negative_column_sum[i] *
                        negative_column_sum[j] /
                        negative_sum
                    ) / negative_sum;

                }

            }

        }
    }

    // Compute modularity
    double Q = positive_sum * Q_positive / total_sum -
               negative_sum * Q_negative / total_sum;

    // Free memory
    free(positive_column_sum);
    free(negative_column_sum);

    // Return modularity
    return Q;

}

/* Functions for main louvain */

// Function to initialize communities
int* initialize_communities(int cols){

    // Allocate memory for communities
    int* communities = (int*)malloc(cols * sizeof(int));

    // Loop over to populate communities
    for(int i = 0; i < cols; i++){
        communities[i] = i + 1;
    }

    // Return communities
    return(communities);

}

// Function to find neighboring nodes
void find_neighbors(
    double *network, int index, int cols,
    int **target_neighbors, int *num_target_neighbors
) {

    // Count the number of neighbors
    int count = 0;

    // Find neighbors
    for (int j = 0; j < cols; j++) {
        if (network[index * cols + j] != 0) {
            count++;
        }
    }

    // Allocate memory for target_neighbors array
    *target_neighbors = (int *)malloc(count * sizeof(int));

    // Fill the target_neighbors array
    int neighbor = 0;
    for (int j = 0; j < cols; j++) {
        if (network[index * cols + j] != 0) {
            (*target_neighbors)[neighbor] = j;
            neighbor++;
        }
    }

    // Set the number of neighbors
    *num_target_neighbors = count;

}

// Re-index communities function
void reindex_comm(int *communities, int size) {
    int *comm_freq = (int *)calloc(size, sizeof(int));
    int i;

    // Get frequencies
    for (i = 0; i < size; i++) {
        comm_freq[communities[i] - 1]++;
    }

    // Order frequencies in descending order
    int *freq_order = (int *)malloc(size * sizeof(int));
    for (i = 0; i < size; i++) {
        freq_order[i] = i;
    }

    int j, temp;
    for (i = 0; i < size; i++) {
        for (j = i + 1; j < size; j++) {
            if (comm_freq[freq_order[i]] < comm_freq[freq_order[j]]) {
                temp = freq_order[i];
                freq_order[i] = freq_order[j];
                freq_order[j] = temp;
            }
        }
    }

    // Re-index them by matching
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (communities[i] == freq_order[j] + 1) {
                communities[i] = j + 1;
                break;
            }
        }
    }

    free(comm_freq);
    free(freq_order);
}

// Structure for `make_higher_order`
struct higher_order_result {
    double* new_network;
    int new_cols;
};

// Function to make higher-order Louvain results
struct higher_order_result make_higher_order(double* network, int* communities, int cols){

    // Initialize iterators
    int i, j;

    // Re-index communities
    reindex_comm(communities, cols);

    // Obtain number of communities
    int num_communities = 0;

    // Loop over to get maximum communities
    for(i = 0; i < cols; i++){

        // Check for greater value
        if(num_communities < communities[i]){
            num_communities = communities[i];
        }

    }

    // Initialize unique communities (automatically sorts)
    int* unique_communities = (int *) calloc(num_communities, sizeof(int)); // freed later

    // Populate unique communities
    for(i = 0; i < num_communities; i++){
        unique_communities[i] = i + 1;
    }

    // Initialize new network
    double* new_network = (double *) calloc(num_communities * num_communities, sizeof (double));

    // Populate new network with sums of edges in each community
    for(i = 0; i < cols; i++){
        for(j = 0; j < cols; j++){

            // Check for i less than j
            if(i <= j){

                // Obtain target communities
                int target_community_i = communities[i] - 1;
                int target_community_j = communities[j] - 1;

                // Obtain target edge
                double target_edge = network[i * cols + j];

                // Populate new network
                new_network[
                    target_community_i * num_communities + target_community_j
                ] += target_edge;

                // Fill other side
                new_network[
                    target_community_j * num_communities + target_community_i
                ] += target_edge;

            }

        }
    }


    // Set up return structure
    struct higher_order_result results = {
        new_network,
        num_communities
    };

    // Free memory
    free(unique_communities);

    // Return results
    return(results);

}

// Define the main louvain function
struct louvain_result {
    int *communities;
    double modularity;
};

// Main Louvain function
struct louvain_result main_louvain(
    double *network, double *original_network,
    int *previous_communities, double previous_modularity,
    int cols, int original_cols
) {

    // Initialize iterators
    int i, j;

    // Initialize communities and modularity
    int* update_communities = (int*)malloc(original_cols * sizeof(int)); // not yet freed
    double update_modularity = 0.0;

    // Update communities and modularity
    if(previous_communities != NULL){

        // Copy values
        memcpy(update_communities, previous_communities, original_cols * sizeof(int));
        update_modularity = previous_modularity; // not yet freed

    }

    // Initialize communities
    int* communities = initialize_communities(cols); // freed by return

    // Initialize modularity
    double Q = 0.0;

    // Initialize gain and iterations
    int gain = 1;
    int iter = 0;

    // Set up `while` loop for modularity improvement
    while(gain == 1){

        // Initialize gain
        gain = 0;

        // Loop over nodes
        for(i = 0; i < cols; i++){

            // Initialize best increase
            double best_increase = 0.0;

            // Initialize gain vector
            double* gain_vector = (double *) calloc(cols + 1, sizeof(double)); // freed later

            // Set target community
            int target_community = communities[i];

            // Initialize new community
            int new_community = target_community;
            int new_community_copy = new_community;

            // Find neighbors of node
            int *target_neighbors = NULL;
            int num_target_neighbors;
            find_neighbors(
                network, i, cols,
                &target_neighbors, // freed later
                &num_target_neighbors
            );

            // Check that there are neighbors
            if (num_target_neighbors != 0) {

                // Loop over neighbors
                for (j = 0; j < num_target_neighbors; j++) {

                    // Determine neighbor's community
                    int neighbor_community = communities[target_neighbors[j]];

                    // Check for gain with neighbor's community
                    if (gain_vector[neighbor_community] == 0) {

                        // Check for whether current community equals neighbor's community
                        if (communities[i] != neighbor_community) {

                            // Allocate memory for new communities
                            int *new_communities = (int *) malloc(cols * sizeof(int)); // freed later

                            // Copy original communities to new communities
                            memcpy(new_communities, communities, cols * sizeof(int));

                            // Update current node's community with neighboring community
                            new_communities[i] = neighbor_community;

                            // Compute difference
                            double gain_difference = signed_modularity(network, new_communities, cols) -
                                                     signed_modularity(network, communities, cols);

                            // Compute difference and add to gains
                            gain_vector[neighbor_community] = gain_difference;

                            // Check if gain is better than previous best increase
                            if (gain_difference > best_increase) {
                                best_increase = gain_difference; // Update best increase
                                new_community_copy = neighbor_community; // Update community membership
                            }

                            // Free memory
                            free(new_communities);

                        }

                    }
                }
            }

            // Free memory
            free(target_neighbors);
            free(gain_vector);

            // Check for whether best increase is greater than 0
            // That is, there is an improvement in modularity
            if(best_increase > 0){

                // Update new community
                new_community = new_community_copy;

            }

            // Check if there are previous communities
            if(previous_communities == NULL){

                // Update current communities
                communities[i] = new_community;

            }else{

                // Allocate memory for improve communities
                int* improve_communities = (int*)calloc(original_cols, sizeof(int)); // freed later

                // Make copy of previous communities
                memcpy(improve_communities, previous_communities, original_cols * sizeof(int));

                // Update communities
                for(j = 0; j < original_cols; j++){

                    // Check for target community
                    if(improve_communities[j] == target_community){
                        // Update communities
                        improve_communities[j] = new_community;
                    }

                }

                // Compute modularity
                double improve_modularity = signed_modularity(
                    original_network, improve_communities, original_cols
                );

                // Check for update
                if(improve_modularity > previous_modularity){

                    // Update communities
                    communities[i] = new_community;

                    // Copy into update communities
                    for(j = 0; j < original_cols; j++){
                        update_communities[j] = improve_communities[j];
                    }

                    // Update modularity
                    update_modularity = improve_modularity;

                }

                // Free memory
                free(improve_communities);

            }

            // Reset gain
            if(new_community != target_community){
                gain = 1; // there was an improvement, keep going
            }

        }

        // Check for whether there was an original network
        // This part works on higher levels after first pass
        if(original_network == NULL){

            // Update modularity
            Q = signed_modularity(network, communities, cols);

        }else{

            // Re-index communities
            reindex_comm(previous_communities, original_cols);
            reindex_comm(update_communities, original_cols);

            // Check if there is no change in modularity
            if(previous_modularity == update_modularity){

                // Set gain to zero
                gain = 0;

                // Re-initialize communities
                communities = update_communities;

                // Update modularity
                Q = update_modularity;

            }else{

                // Obtain modularity
                Q = signed_modularity(original_network, update_communities, original_cols);

                // Make network higher order
                struct higher_order_result higher_order = make_higher_order(
                    original_network,
                    update_communities,
                    original_cols
                );

                // Update columns
                cols = higher_order.new_cols;

                // Re-initialize communities
                communities = initialize_communities(cols); // freed by return

                // Free network
                free(network);

                // Re-initialize network
                // double* // apparently, don't need to re-intialize
                network = higher_order.new_network;

                // Free new network
                // free(higher_order.new_network); // apparently, don't need to free?

                // Re-initialize previous communities
                previous_communities = update_communities;

                // Update modularity
                previous_modularity = update_modularity;

            }

        }

        // Increase iterations
        iter++;

    }

    // Re-index communities for final return
    reindex_comm(communities, original_cols);

    // Free communities
    free(update_communities);

    // Set up results
    struct louvain_result results = {
        communities, Q
    };

    // Return results
    return(results);

}

// Define the signed louvain function
struct signed_louvain_result {
    int *communities;
    double *modularities;
    int count;
};

// Signed Louvain function
struct signed_louvain_result signed_louvain(double* network, int cols) {

    // Initialize iterators
    int i, j;

    // Allocate memory for communities
    int* communities = (int*)malloc(cols * cols * sizeof(int)); // freed later

    // Allocate memory for modularities
    double* modularities = (double*)malloc(cols * sizeof(double)); // freed later

    // Initialize count
    int count = 0;

    // Obtain lowest order Louvain results
    struct louvain_result louvain_initial = main_louvain(
        network, NULL,  // double* network, original_network
        NULL,           // int *previous_communities
        0,              // double previous_modularity
        cols, cols      // int cols, original_cols
    );

    // Store results
    for (i = 0; i < cols; i++) {
        communities[count * cols + i] = louvain_initial.communities[i];
    }
    modularities[count] = louvain_initial.modularity;

    // Free memory
    free(louvain_initial.communities); // frees first run communities

    // Declare re-used structures
    struct louvain_result louvain_next_run;
    struct higher_order_result higher_order;

    // Enter `while` loop for next levels
    while (1) {

        // Increase count
        count++;

        // Obtain higher-order network
        higher_order = make_higher_order(
            network, &communities[(count - 1) * cols], cols
        );

        // Obtain first-level Louvain results
        louvain_next_run = main_louvain(
            higher_order.new_network,       // double *network
            network,                        // double *original_network
            &communities[(count - 1) * cols], // int *previous_communities
            modularities[count - 1],        // double previous_modularity
            higher_order.new_cols,          // int cols
            cols                            // int original_cols
        );

        // Check for break condition
        if (louvain_next_run.modularity == modularities[count - 1]) {

            // Break out of loop
            break;

        } else {

            // Store results
            for (i = 0; i < cols; i++) {
                communities[count * cols + i] = louvain_next_run.communities[i];
            }
            modularities[count] = louvain_next_run.modularity;

        }

    }

    // Allocate memory for final communities
    int* final_communities = (int*)malloc(count * cols * sizeof(int)); // freed with return

    // Shrink matrix
    for (i = 0; i < count; i++) {
        for (j = 0; j < cols; j++) {
            final_communities[i * cols + j] = communities[i * cols + j];
        }
    }

    // Free communities
    free(communities);

    // Allocate memory for final modularities
    double* final_modularities = (double*)malloc(count * sizeof(double)); // freed with return

    // Shrink vector
    for (i = 0; i < count; i++) {
        final_modularities[i] = modularities[i];
    }

    // Free modularities
    free(modularities);

    // Set up result
    struct signed_louvain_result result = {
        final_communities,
        final_modularities,
        count
    };

    // Return result
    return(result);

}

SEXP r_signed_louvain(SEXP r_input_network) {

    // Initialize iterators
    int i;

    // Obtain columns
    int cols = ncols(r_input_network);

    // Call the C function
    struct signed_louvain_result c_result = signed_louvain(REAL(r_input_network), cols);

    // Create R output list
    SEXP r_output = PROTECT(allocVector(VECSXP, 2));

    // Convert the communities result to an R matrix
    SEXP r_communities = PROTECT(allocMatrix(INTSXP, cols, c_result.count));
    for (i = 0; i < c_result.count * cols; i++) {
        INTEGER(r_communities)[i] = c_result.communities[i];
    }

    // Convert the modularities result to an R numeric vector
    SEXP r_modularities = PROTECT(allocVector(REALSXP, c_result.count));
    for (i = 0; i < c_result.count; i++) {
        REAL(r_modularities)[i] = c_result.modularities[i];
    }

    // Set output list elements
    SET_VECTOR_ELT(r_output, 0, r_communities);
    SET_VECTOR_ELT(r_output, 1, r_modularities);

    // Assign names to the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("communities"));
    SET_STRING_ELT(names, 1, mkChar("modularities"));
    setAttrib(r_output, R_NamesSymbol, names);

    // Free memory
    free(c_result.communities);
    free(c_result.modularities);

    // Release protected SEXP objects
    UNPROTECT(4);

    // Return the result
    return r_output;

}

///* Final output */
//int main() {
//
//    // Initialize iterators
//    int i;
//
//    // Initialize columns
//    int cols = 10;
//
//    // Set up network
//    double network_data[] = {
//        0.00, 0.21, 0.03, 0.21, 0.10, 0.00, 0.00, 0.00, 0.00, 0.06,
//        0.21, 0.00, 0.36, 0.02, 0.23, 0.06, 0.00, 0.03, 0.05, 0.09,
//        0.03, 0.36, 0.00, 0.08, 0.00, 0.15, 0.01, 0.02, 0.00, 0.06,
//        0.21, 0.02, 0.08, 0.00, 0.10, 0.07, 0.08, 0.08, 0.02, 0.12,
//        0.10, 0.23, 0.00, 0.10, 0.00, 0.01, 0.06, 0.01, 0.00, 0.03,
//        0.00, 0.06, 0.15, 0.07, 0.01, 0.00, 0.16, 0.11, 0.10, 0.13,
//        0.00, 0.00, 0.01, 0.08, 0.06, 0.16, 0.00, 0.18, 0.00, 0.06,
//        0.00, 0.03, 0.02, 0.08, 0.01, 0.11, 0.18, 0.00, 0.14, 0.00,
//        0.00, 0.05, 0.00, 0.02, 0.00, 0.10, 0.00, 0.14, 0.00, 0.20,
//        0.06, 0.09, 0.06, 0.12, 0.03, 0.13, 0.06, 0.00, 0.20, 0.00
//    };
//
//    // Initialize network
//    double *network = (double *)malloc(cols * cols * sizeof(double));
//
//    // Copy network_data to network
//    for (i = 0; i < cols * cols; i++) {
//        network[i] = network_data[i];
//    }
//
//    // Obtain signed Louvain result
//    struct signed_louvain_result results = signed_louvain(network, cols);
//
//    // Calculate the number of levels
//    int levels = 0;
//    while (results.modularities[levels] != 0) {
//        levels++;
//    }
//
//    // Print communities
//    for (int i = 0; i < levels; i++) {
//        for (int j = 0; j < cols; j++) {
//            printf("%d ", results.communities[i * cols + j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//
//    // Print modularity
//    for (int i = 0; i < levels; i++) {
//        printf("%f ", results.modularities[i]);
//    }
//    printf("\n");
//
//    // Free memory
//    free(network);
//
//    // Free output
//    free(results.communities);
//    free(results.modularities);
//
//    // Return success
//    return 0;
//
//}

