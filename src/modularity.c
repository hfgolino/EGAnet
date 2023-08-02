#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "modularity.h"

// Function to compute modularity values
struct modularity_result modularity_values(double* network, int cols, double resolution) {

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
                    (resolution * positive_column_sums[i] * positive_column_sums[j] / positive_sum)
                ) / positive_sum;

            }

            // Check for negative sum
            if(negative_sum_flag) {

                // Update negative modularity
                negative_modularity_values[count] += (
                    ((edge < 0) ? edge : 0) - // negative edge
                    (resolution * negative_column_sums[i] * negative_column_sums[j] / negative_sum)
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
