#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include "xoshiro.h"
#include "nanotime.h"

/* Unsigned (standard) modularity Louvain

   A leaner companion to `signed_louvain.c`: rather than splitting
   modularity into separate positive/negative components (Gomez et al.,
   2009), this operates on a single non-negative weight matrix and
   optimizes standard Newman modularity. Callers are expected to pass
   in non-negative weights (e.g., `abs()`'d correlations), matching how
   {EGAnet}'s `community.detection()` already treats every other
   "louvain"-style algorithm.

   Uses the standard incremental Louvain bookkeeping (Blondel et al.,
   2008) rather than a naive dense reference: a sparse adjacency list
   is built once per aggregation level, and each candidate move's
   modularity gain is evaluated in O(1) from three things -- a node's
   own weighted degree, its edge weight into each candidate community
   (found via a single O(degree) neighbor scan), and each community's
   running total degree (updated in O(1) as nodes move). This keeps
   per-pass cost at roughly O(cols + edges) instead of O(cols^2). Two
   smaller, unconditionally value-preserving choices ride along:

     - No hypothetical-membership copy: since a candidate move only
       ever differs from the current membership at the one node being
       moved, there's no need for an O(cols) copy per candidate --
       that one index is just treated specially in O(1).
     - No per-node `calloc`/`free` of a size-`cols` scratch vector: an
       "epoch stamp" (`touched_epoch`) marks which communities have
       been seen for the *current* node in O(1), avoiding an O(cols)
       reset (or allocation) on every node visited.

   `modularity_values()`/`compute_modularity()` do the dense O(cols^2)
   pairwise computation and are used only for the *reported*
   modularity at each (external) level and for the internal
   convergence bookkeeping's starting point -- not a hot path, since
   they're called once per level rather than once per node per pass.

   Move acceptance requires the incremental gain to exceed
   `GAIN_EPSILON`, not just be positive. Structurally symmetric
   networks (e.g., two nodes with identical connectivity to every
   other node) produce *exact* modularity ties between candidate
   moves; the dense pairwise computation happens to round these to
   precisely 0 (the same precomputed terms are summed for both
   directions), but the incremental formula reaches the same value via
   different floating-point operations and can round to a tiny
   nonzero epsilon (observed: ~5e-19) instead of exact 0. A strict
   `> 0` check then accepts that "gain" every time it's evaluated --
   including the reverse move on the next pass, which rounds the same
   way -- so a node can swap back and forth between two communities
   forever without a full pass ever reporting zero moves.
*/

// Minimum modularity gain for a candidate move to be accepted -- filters
// out floating-point noise around exact ties (see above) while being far
// below any gain that would matter for realistically-scaled networks
#define GAIN_EPSILON 1e-10

// Structure for `modularity_values`
struct modularity_result {
    double* values;
    double total_sum;
    int nonzero_flag;
};

// Function to compute modularity values
static struct modularity_result modularity_values(double* network, int cols, double resolution) {

    // Initialize iterators
    int i, j, network_offset;
    double edge;
    int count = 0;

    // Initialize sums
    double* column_sums = (double*)calloc(cols, sizeof(double));
    double total_sum = 0.0;

    // Loop over to get sums
    for(i = 0; i < cols; i++) {

        // Compute network offset
        network_offset = i * cols;

        for(j = i; j < cols; j++) {

            // Get edge
            edge = network[network_offset + j];

            if(edge != 0) {

                // Add to sums
                column_sums[i] += edge;

                // Check for off-diagonal edge
                if(i != j) {
                    column_sums[j] += edge;
                }

            }

        }
    }

    // Compute total sum
    for(i = 0; i < cols; i++) {
        total_sum += column_sums[i];
    }

    // Get indices for modularity values
    int value_indices = (((cols * (cols - 1)) / 2) + cols);

    // Initialize modularity values
    double* values = (double*)calloc(value_indices, sizeof(double));

    // Check if total_sum is not zero (guards against divide-by-zero for
    // empty/disconnected networks)
    int nonzero_flag = total_sum != 0;

    // Loop over to compute modularity
    for(i = 0; i < cols; i++) {

        // Compute network offset
        network_offset = i * cols;

        for(j = i; j < cols; j++) {

            if(nonzero_flag) {

                // Obtain edge
                edge = network[network_offset + j];

                // Update modularity
                values[count] = (
                    edge - (resolution * column_sums[i] * column_sums[j] / total_sum)
                ) / total_sum;

            }

            // Increase count
            count++;

        }
    }

    // Free memory
    free(column_sums);
    column_sums = NULL;

    // Set up result
    struct modularity_result result = {
        values, total_sum, nonzero_flag
    };

    // Return result
    return(result);

}

// Standard (unsigned) modularity function
static double compute_modularity(struct modularity_result Q_values, int* membership, int cols) {

    // Initialize iterators
    int i, j;
    int count = 0;

    // Initialize modularity
    double Q = 0.0;
    double value;

    // Loop over to compute modularity
    for(i = 0; i < cols; i++) {

        for(j = i; j < cols; j++) {

            // Check for memberships
            if(membership[i] == membership[j]) {

                value = (Q_values.nonzero_flag) ? Q_values.values[count] : 0;

                // Update modularity
                Q += value * 2;

                // Check for diagonal edge
                if(i == j) {
                    Q -= value;
                }

            }

            // Increase count
            count++;

        }

    }

    // Return modularity
    return Q;

}

/* Functions for Main Louvain */

// Get maximum membership
static int maximum_membership(int* memberships, int cols) {

    // Initialize max membership
    int max_membership = 0;

    // Loop over to get maximum membership
    for(int i = 0; i < cols; i++) {

        // Check for greater membership
        if(memberships[i] > max_membership) {
            max_membership = memberships[i];
        }

    }

    // Return maximum membership
    return max_membership;

}

// Structure for `make_higher_order`
struct higher_order_result {
    double* higher_order;
    int number;
};

// Function to make higher-order network
static struct higher_order_result make_higher_order(double* network, int* membership, int cols) {

    // Initialize iterators
    int i, j, network_offset;
    double edge;
    int first_position, second_position;

    // Get maximum membership (number of communities)
    int number = maximum_membership(membership, cols);

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
static int* reindex_membership(int* membership, int cols) {

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

// Fisher-Yates (or Knuth) Shuffle
static void shuffle_nodes(int *arr, int cols, xoshiro256_state* rng_state) {

    // Initialize iterators
    int i, j, temp;

    // Iterate through the array from the last element to the first
    for(i = cols - 1; i > 0; i--) {

        // Generate a random index between 0 and i (inclusive)
        j = next(rng_state) % (i + 1);

        // Swap the current element with the randomly selected element
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;

    }

}

/* Sparse adjacency list (CSR), built once per aggregation level and
   reused across every local-moving pass at that level -- replaces a
   repeated O(cols) dense-row neighbor scan */
struct adjacency_list {
    int* neighbors;    // flattened neighbor indices, size = nnz
    double* weights;   // corresponding edge weights, size = nnz
    int* row_start;    // size cols + 1; node i's neighbors are
                       // neighbors[row_start[i] .. row_start[i+1])
    double* degree;    // size cols; k_i (row sum), matches `column_sums`
};

static struct adjacency_list build_adjacency_list(double* network, int cols) {

    int i, j, offset, nnz = 0;

    int* counts = (int*)calloc(cols, sizeof(int));
    double* degree = (double*)calloc(cols, sizeof(double));

    // First pass: node degrees and neighbor counts (for CSR offsets).
    // `degree[i]` includes the diagonal (self-loop) term, matching
    // `column_sums` in `modularity_values` -- but self-loops are
    // deliberately *excluded* from the neighbor list built below: a
    // node's own self-loop contributes identically to Q regardless of
    // which community that node is in, so it cancels out of every
    // move's gain exactly. Treating it as a normal neighbor edge here
    // would double-count it as "evidence" for the node's *current*
    // community -- this matters because aggregated (level 2+)
    // networks carry self-loops representing internal community
    // weight, even when the original input network has none
    for(i = 0; i < cols; i++) {

        offset = i * cols;

        for(j = 0; j < cols; j++) {

            double edge = network[offset + j];
            degree[i] += edge;

            if(edge != 0 && j != i) {
                counts[i]++;
                nnz++;
            }

        }

    }

    int* row_start = (int*)malloc((cols + 1) * sizeof(int));
    row_start[0] = 0;

    for(i = 0; i < cols; i++) {
        row_start[i + 1] = row_start[i] + counts[i];
    }

    int* neighbors = (int*)malloc((nnz > 0 ? nnz : 1) * sizeof(int));
    double* weights = (double*)malloc((nnz > 0 ? nnz : 1) * sizeof(double));
    int* fill_pos = (int*)malloc(cols * sizeof(int));
    memcpy(fill_pos, row_start, cols * sizeof(int));

    // Second pass: fill neighbor/weight arrays in ascending-index order,
    // excluding the diagonal (see above)
    for(i = 0; i < cols; i++) {

        offset = i * cols;

        for(j = 0; j < cols; j++) {

            double edge = network[offset + j];

            if(edge != 0 && j != i) {
                int p = fill_pos[i]++;
                neighbors[p] = j;
                weights[p] = edge;
            }

        }

    }

    free(counts);
    free(fill_pos);

    struct adjacency_list result = {neighbors, weights, row_start, degree};
    return(result);

}

static void free_adjacency_list(struct adjacency_list* adj) {
    free(adj->neighbors); adj->neighbors = NULL;
    free(adj->weights); adj->weights = NULL;
    free(adj->row_start); adj->row_start = NULL;
    free(adj->degree); adj->degree = NULL;
}

// Main Louvain function
static void main_louvain(
    double* network,
    struct modularity_result Q_values,
    int* membership_copy, double previous_modularity,
    int cols, int original_cols,
    double resolution, xoshiro256_state* rng_state,
    int is_first_level
) {

    // Initialize iterators
    int i, j;

    // Tracks whether `Q_values` currently points to memory this function
    // allocated itself (via `modularity_values`) rather than memory owned
    // by the caller -- determines what gets freed below
    int Q_values_owned = 0;

    // Initialize membership
    int* membership = (int*)malloc(cols * sizeof(int));
    int* index = (int*)malloc(cols * sizeof(int));

    // Populate membership
    for(i = 0; i < cols; i++) {
        membership[i] = i;
        index[i] = i;
    }

    // Initialize update to modularity
    double update_modularity = compute_modularity(Q_values, membership, cols);

    // Initialize commonly used values
    int target_membership, neighbor_membership, new_membership;
    double best_increase;

    struct higher_order_result higher_order;

    // Adjacency list + running community-degree totals for this level --
    // `sigma_tot[c]` is Sigma_tot(c), the running total weighted degree of
    // every node currently assigned to community `c`
    struct adjacency_list adj = build_adjacency_list(network, cols);
    double* sigma_tot = (double*)malloc(cols * sizeof(double));
    for(i = 0; i < cols; i++) {
        sigma_tot[i] = adj.degree[i];
    }
    double total_sum = Q_values.total_sum;

    // Epoch-cached per-node scratch: `k_in[c]` accumulates the edge weight
    // from the node currently being visited into community `c`;
    // `touched_epoch[c] == epoch` marks `c` as already seen for *this*
    // node, standing in for an O(cols) reset/allocation every node visit
    double* k_in = (double*)calloc(cols, sizeof(double));
    int* touched_epoch = (int*)malloc(cols * sizeof(int));
    for(i = 0; i < cols; i++) {
        touched_epoch[i] = -1;
    }
    int* touched = (int*)malloc(cols * sizeof(int));
    int epoch = 0;

    // `while` until there is no gain
    while(1) {

        // Allow user to interrupt long-running optimization
        R_CheckUserInterrupt();

        // Permutate index
        shuffle_nodes(index, cols, rng_state);

        // Set gain to zero
        int gain = 0;

        // Loop over nodes
        for(i = 0; i < cols; i++) {

            // Reset best increase
            best_increase = 0.0;

            // Get index order
            int order = index[i];

            // Obtain membership of target node
            target_membership = membership[order];

            // Initialize new membership
            new_membership = target_membership;

            double k_v = adj.degree[order];

            // Scan `order`'s neighbors once, accumulating the edge weight
            // into each distinct community encountered (`k_in`), and
            // recording the order communities are first seen (`touched`)
            int n_touched = 0;
            epoch++;

            for(int e = adj.row_start[order]; e < adj.row_start[order + 1]; e++) {

                int neighbor = adj.neighbors[e];
                int neighbor_comm = membership[neighbor];

                if(touched_epoch[neighbor_comm] != epoch) {
                    touched_epoch[neighbor_comm] = epoch;
                    k_in[neighbor_comm] = 0.0;
                    touched[n_touched++] = neighbor_comm;
                }

                k_in[neighbor_comm] += adj.weights[e];

            }

            // `order`'s edge weight into its own (current) community --
            // 0 if none of its neighbors are currently co-members
            double k_in_target = (touched_epoch[target_membership] == epoch) ?
                k_in[target_membership] : 0.0;
            double sigma_target = sigma_tot[target_membership];

            // Evaluate each distinct candidate community exactly once, in
            // the order first encountered -- standard incremental
            // modularity-gain formula (Blondel et al., 2008)
            for(int t = 0; t < n_touched; t++) {

                neighbor_membership = touched[t];

                // Moving to the current community is a no-op
                if(neighbor_membership == target_membership) {
                    continue;
                }

                // Note the `+ k_v`: `sigma_target` (Sigma_tot(C)) still
                // includes `order` itself at this point (it's only
                // removed from the running total once a move is actually
                // accepted, below), so the null-model term needs Sigma_tot(C)
                // *after* removal, i.e. `sigma_target - k_v`; folded into
                // the difference this becomes `+ k_v` on the other side
                double gain_value = (
                    2.0 * (k_in[neighbor_membership] - k_in_target) -
                    (2.0 * resolution * k_v / total_sum) *
                        (sigma_tot[neighbor_membership] - sigma_target + k_v)
                ) / total_sum;

                if(gain_value > best_increase) {
                    best_increase = gain_value;
                    new_membership = neighbor_membership;
                }

            }

            // Check if there was an increase (see `GAIN_EPSILON` above --
            // guards against floating-point noise around exact ties)
            if(best_increase > GAIN_EPSILON) {

                // Update running community degree totals
                sigma_tot[target_membership] -= k_v;
                sigma_tot[new_membership] += k_v;

                // Update membership
                membership[order] = new_membership;

                // Update modularity
                update_modularity += best_increase;

                // Reset gain
                gain = 1;

            }

        }

        // Check for higher-order network
        if(!is_first_level) {

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

                // Free previously (internally) allocated modularity values
                // before overwriting -- avoids leaking every intermediate
                // aggregation level's `Q_values` arrays
                if(Q_values_owned) {
                    free(Q_values.values);
                }

                // Update modularity values
                Q_values = modularity_values(network, cols, resolution);
                Q_values_owned = 1;

                // Update modularity
                update_modularity = compute_modularity(Q_values, membership, cols);

                // Rebuild adjacency list / community bookkeeping for the
                // new (aggregated) network
                free_adjacency_list(&adj);
                adj = build_adjacency_list(network, cols);

                sigma_tot = (double*)realloc(sigma_tot, cols * sizeof(double));
                for(i = 0; i < cols; i++) {
                    sigma_tot[i] = adj.degree[i];
                }
                total_sum = Q_values.total_sum;

                k_in = (double*)realloc(k_in, cols * sizeof(double));
                touched_epoch = (int*)realloc(touched_epoch, cols * sizeof(int));
                for(i = 0; i < cols; i++) {
                    touched_epoch[i] = -1;
                }
                touched = (int*)realloc(touched, cols * sizeof(int));
                epoch = 0;

            }

        }

        // Check gain for break condition
        if(gain == 0) {
            break;
        }

    }

    // Check for level
    if(is_first_level) {

        // Check for unidimensionality
        if(cols == 1) {

            // Index membership
            for(i = 0; i < original_cols; i++) {
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
            for(i = 0; i < original_cols; i++) {
                membership_copy[i] = 0;
            }

        } else {

            // Re-index membership
            membership_copy = reindex_membership(membership_copy, original_cols);

        }

        // Free memory (only if this function allocated it -- otherwise
        // these are the caller's arrays and the caller owns freeing them)
        if(Q_values_owned) {
            free(Q_values.values);
            Q_values.values = NULL;
        }

    }

    // Free memory
    free_adjacency_list(&adj);
    free(sigma_tot);
    free(k_in);
    free(touched_epoch);
    free(touched);

    free(index);
    free(membership);
    index = NULL;
    membership = NULL;

}

/* Functions for Louvain */

// Structure for `louvain`
struct louvain_result {
    int** memberships;
    double* modularities;
    int level; // index of the last *allocated* level (may be a rejected level)
    int valid_levels; // number of levels, starting at 0, that are safe to return
};

// Louvain function
static struct louvain_result louvain(
    double* original_network, int original_cols,
    double resolution, double seed, int lower_order
) {

    // Seed the random number generator (xoshiro256++, consistent with the
    // rest of {EGAnet}) -- a `seed` of 0 falls back to clock time, matching
    // the convention used by `r_xoshiro_*` in xoshiro.c
    uint64_t seed_value = (uint64_t) seed;
    if(seed_value == 0) {
        seed_value = get_time_ns();
    }
    xoshiro256_state rng_state;
    seed_xoshiro256(&rng_state, seed_value);

    // Initialize iterators
    int i, level = 0;

    // Level 0 (the finest partition, computed below) is always valid to
    // return -- `valid_levels` only grows further when a coarser level
    // both improves modularity and isn't a full collapse to one community
    int valid_levels = 1;

    // Initialize membership and modularity results
    int** memberships = (int**)malloc(original_cols * sizeof(int*));
    int* membership_copy = (int*)calloc(original_cols, sizeof(int));
    double* modularities = (double*)malloc(original_cols * sizeof(double));

    // Compute original modularity values
    struct modularity_result original_Q_values = modularity_values(
        original_network, original_cols, resolution
    );

    // Compute first level
    main_louvain(
        original_network, original_Q_values,
        membership_copy, 0.0,
        original_cols, original_cols,
        resolution, &rng_state, 1 /* is_first_level */
    );

    // Compute modularity
    modularities[level] = compute_modularity(
        original_Q_values, membership_copy, original_cols
    );

    // Populate membership
    memberships[level] = (int*)malloc(original_cols * sizeof(int));
    memcpy(memberships[level], membership_copy, original_cols * sizeof(int));

    // Initialize higher order and modularity values
    struct higher_order_result higher_order;
    struct modularity_result Q_values;

    // Get maximum membership
    int max_communities = maximum_membership(memberships[level], original_cols);

    // `order = "lower"` stops here: level 0 (the first, un-aggregated
    // local-moving pass over the original network) is the entire result,
    // with no coarsening/merging of communities into higher-order nodes.
    // `order = "higher"` (default) continues to repeatedly aggregate and
    // re-optimize until modularity stops improving
    while(!lower_order && max_communities != (original_cols - 1)) {

        // Allow user to interrupt long-running optimization
        R_CheckUserInterrupt();

        // Increase level
        level++;

        // Obtain higher-order network
        higher_order = make_higher_order(
            original_network, membership_copy, original_cols
        );

        // Compute modularity values
        Q_values = modularity_values(
            higher_order.higher_order, higher_order.number, resolution
        );

        // Compute next level
        main_louvain(
            higher_order.higher_order, Q_values,
            membership_copy, modularities[level - 1],
            higher_order.number, original_cols,
            resolution, &rng_state, 0 /* is_first_level */
        );

        // Compute modularity
        modularities[level] = compute_modularity(
            original_Q_values, membership_copy, original_cols
        );

        // Populate membership
        memberships[level] = (int*)malloc(original_cols * sizeof(int));
        memcpy(memberships[level], membership_copy, original_cols * sizeof(int));

        // Free memory
        free(Q_values.values);
        free(higher_order.higher_order);

        // Set pointers to NULL
        Q_values.values = NULL;
        higher_order.higher_order = NULL;

        // Check for any change
        if(modularities[level] <= modularities[level - 1]) {
            break;
        }

        // Initialize all zero count (unidimensionality)
        int all_zero = 0;

        // Count number of zeros
        for(i = 0; i < original_cols; i++) {
            all_zero += memberships[level][i] == 0;
        }

        // If number of zeros is the same as length of memberships, then break
        if(all_zero == original_cols) {
            break;
        }

        // This level improved modularity and isn't a full collapse -- safe
        // to return it (and any level up to it) to the caller
        valid_levels = level + 1;

        // Update communities
        max_communities = maximum_membership(memberships[level], original_cols);

    }

    // Any allocated level at index `valid_levels` or beyond was rejected
    // (didn't improve modularity, or collapsed to one community) and is
    // never read by the caller (`r_louvain` only returns levels
    // `0` .. `valid_levels - 1`) -- free it here or it leaks. At most one
    // such level can exist, since the loop above breaks immediately upon
    // rejecting a level
    if(valid_levels <= level) {
        free(memberships[level]);
        memberships[level] = NULL;
    }

    // Free memory
    free(original_Q_values.values);
    free(membership_copy);

    // Set pointers to NULL
    original_Q_values.values = NULL;
    membership_copy = NULL;

    // Set up result
    struct louvain_result result = {
        memberships,
        modularities,
        level,
        valid_levels
    };

    // Return result
    return(result);

}

// `r_lower_order`: 0 = order = "higher" (default; full multilevel
// aggregation), 1 = order = "lower" (first local-moving pass only,
// no community merging)
SEXP r_louvain(SEXP r_input_network, SEXP r_resolution, SEXP r_seed, SEXP r_lower_order) {

    // Initialize iterators
    int i, j;

    // Obtain columns
    int cols = ncols(r_input_network);

    // Call the C function
    struct louvain_result c_result = louvain(
        REAL(r_input_network), cols, REAL(r_resolution)[0], REAL(r_seed)[0],
        INTEGER(r_lower_order)[0]
    );

    // Create R output list
    SEXP r_output = PROTECT(allocVector(VECSXP, 2));

    // Convert the memberships result to an R matrix
    SEXP r_memberships = PROTECT(allocMatrix(INTSXP, c_result.valid_levels, cols));
    for(i = 0; i < cols; i++) {
        for(j = 0; j < c_result.valid_levels; j++) {
            INTEGER(r_memberships)[i * c_result.valid_levels + j] = c_result.memberships[j][i] + 1;
        }
    }

    // Convert the modularities result to an R numeric vector
    SEXP r_modularities = PROTECT(allocVector(REALSXP, c_result.valid_levels));
    for(i = 0; i < c_result.valid_levels; i++) {
        REAL(r_modularities)[i] = c_result.modularities[i];
        free(c_result.memberships[i]);
    }

    // Set output list elements
    SET_VECTOR_ELT(r_output, 0, r_memberships);
    SET_VECTOR_ELT(r_output, 1, r_modularities);

    // Assign names to the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("memberships"));
    SET_STRING_ELT(names, 1, mkChar("modularity"));
    setAttrib(r_output, R_NamesSymbol, names);

    // Release protected SEXP objects
    UNPROTECT(4);

    // Free memory
    free(c_result.memberships);
    free(c_result.modularities);

    // Return the result
    return r_output;

}
