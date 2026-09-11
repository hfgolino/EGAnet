#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

// Tolerance for the branch-and-bound prune comparison. Sparse partial-
// correlation networks (e.g., glasso output) routinely have exact-zero
// edges, which can make many distinct sign assignments exactly (or
// near-exactly) tied for the optimum. Two numerically-equivalent ways of
// summing the same edge weights -- e.g., one running total across all
// pairs vs. per-node partial sums added together -- do not generally
// agree to the last bit (floating-point addition isn't associative), so
// a tied branch's `remaining_bound` can land a few ULPs above exact
// zero instead of at it, letting the strict `<=` comparison let it
// through as a spurious "improvement" over an already-optimal score.
// With enough near-tied nodes those spurious improvements compound
// combinatorially -- this was observed directly (a real, non-adversarial
// network whose R prototype finished in 99 visited nodes took millions
// in a first C translation using plain `<=`). A small absolute epsilon
// -- several orders of magnitude above worst-case floating-point
// accumulation error for any community size this is used on, and
// utterly negligible next to any `sign_objective` difference a user
// could care about -- prunes true ties (and near-ties well below
// measurement noise) without affecting which distinct-valued solution
// is reported: the true optimum is still found, since anything more
// than `PRUNE_EPSILON` better than the current best is never pruned
#define PRUNE_EPSILON 1e-9

/* Exact sign search (branch and bound) ----

   A direct C translation of the R prototype `exact_signs()` in
   `net.loads.R` -- same node-visit order and pruning rule, differing
   only in the `PRUNE_EPSILON` tolerance described above (validated
   against the R prototype in the package's development scripts, not
   shipped as a test since it's an internal, deterministic equivalence
   check rather than user-facing behavior: on every tested network, the
   two reach the same `sign_objective` value -- the same certified
   optimum -- even on the rare input where they land on different
   tied-optimal sign vectors that achieve it). Moving the search itself
   to C removes R's per-call function-dispatch and allocation overhead
   from what can be millions of recursive calls; node reordering (by
   descending absolute degree) stays in R, which already does it in
   vectorized O(nodes^2) time.

   `network` is the *already reordered* nodes x nodes matrix (symmetric,
   zero diagonal), passed in column-major (R's native layout). Node 0
   is fixed to sign +1 (a global sign flip doesn't change the
   objective, `sum_{i<j} A[i,j] * s[i] * s[j]`), and the search assigns
   signs to nodes 1..n-1 depth-first, pruning any branch whose best
   possible completion (`current_score + remaining_bound`) can't beat
   the best complete assignment found so far. `remaining_bound` is the
   sum of |A[i,j]| over all pairs with at least one endpoint still
   unassigned -- a valid upper bound since every such edge can
   contribute at most its absolute weight regardless of final signs. */

typedef struct {
    const double *row_major;       // reordered network, row-major, n x n
    int n;
    const double *bound_decrement; // bound_decrement[i] = sum_{j<i} |A[i,j]|
    signed char *signs;            // current branch's sign assignment
    double *signs_d;                // same, pre-cast to double (see below)
    signed char *best_signs;       // best complete assignment found so far
    double best_score;
    double max_visited;
    double visited;
} bnb_state;

static void search(bnb_state *state, int index, double current_score, double remaining_bound) {
    state->visited += 1.0;

    // Prune: even the best case from here can't beat the current best
    // (by more than a negligible tolerance; see `PRUNE_EPSILON` above)
    if (current_score + remaining_bound <= state->best_score + PRUNE_EPSILON || state->visited > state->max_visited) {
        return;
    }

    // Complete assignment
    if (index == state->n) {
        state->best_score = current_score;
        memcpy(state->best_signs, state->signs, (size_t) state->n * sizeof(signed char));
        return;
    }

    // Edges to already-assigned nodes resolve exactly; edges to nodes
    // not yet assigned remain bounded by their absolute weight. This dot
    // product is, by a wide margin, the hot spot of the whole search --
    // profiling (callgrind) attributes 96% of all instructions executed
    // by the search to this one loop. It's summed as 4 independent
    // running totals rather than 1 (using `signs_d`, a `double`-typed
    // copy of `signs` kept alongside it purely so this loop isn't also
    // re-casting `signed char` to `double` on every element): a single
    // running total forces each addition to wait on the previous one,
    // while 4 independent ones can execute (and, with `-O2`, do get
    // auto-vectorized as) largely in parallel. This reorders the
    // summation relative to a naive left-to-right loop, which is exactly
    // the class of floating-point non-determinism `PRUNE_EPSILON` above
    // already makes the search robust to -- measured at ~2-3x fewer
    // instructions and ~2-3x faster wall-clock time (n=50 to n=200)
    // with no change in the objective value reached
    const double *row = state->row_major + (size_t) index * state->n;
    const double *sd = state->signs_d;
    double c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;
    int j = 0;
    for (; j + 4 <= index; j += 4) {
        c0 += row[j] * sd[j];
        c1 += row[j + 1] * sd[j + 1];
        c2 += row[j + 2] * sd[j + 2];
        c3 += row[j + 3] * sd[j + 3];
    }
    double contribution = c0 + c1 + c2 + c3;
    for (; j < index; j++) {
        contribution += row[j] * sd[j];
    }

    double new_bound = remaining_bound - state->bound_decrement[index];

    // Try the locally-preferred sign first (tighter bound, sooner)
    int preferred = (contribution >= 0.0) ? 1 : -1;

    state->signs[index] = (signed char) preferred;
    state->signs_d[index] = (double) preferred;
    search(state, index + 1, current_score + preferred * contribution, new_bound);

    state->signs[index] = (signed char) -preferred;
    state->signs_d[index] = (double) -preferred;
    search(state, index + 1, current_score - preferred * contribution, new_bound);
}

SEXP r_exact_signs(SEXP r_ordered_network, SEXP r_max_visited) {
    int n = ncols(r_ordered_network);
    const double *network = REAL(r_ordered_network);
    double max_visited = REAL(r_max_visited)[0];

    // Row-major copy for contiguous per-node access in the hot loop
    // (the source matrix is symmetric, but is stored column-major),
    // plus each node's bound decrement, computed once up front since
    // both depend only on node index, not on any particular branch
    double *row_major = (double *) malloc((size_t) n * (size_t) n * sizeof(double));
    double *bound_decrement = (double *) malloc((size_t) n * sizeof(double));

    for (int i = 0; i < n; i++) {
        double abs_sum = 0.0;
        for (int j = 0; j < n; j++) {
            double value = network[i + (size_t) j * n]; // column-major source
            row_major[(size_t) i * n + j] = value;
            if (j < i) {
                abs_sum += fabs(value);
            }
        }
        bound_decrement[i] = abs_sum;
    }

    double total_bound = 0.0;
    for (int i = 1; i < n; i++) {
        total_bound += bound_decrement[i];
    }

    signed char *signs = (signed char *) malloc((size_t) n * sizeof(signed char));
    double *signs_d = (double *) malloc((size_t) n * sizeof(double));
    signed char *best_signs = (signed char *) malloc((size_t) n * sizeof(signed char));
    for (int i = 0; i < n; i++) {
        signs[i] = 1;
        signs_d[i] = 1.0;
        best_signs[i] = 1;
    }

    bnb_state state;
    state.row_major = row_major;
    state.n = n;
    state.bound_decrement = bound_decrement;
    state.signs = signs;
    state.signs_d = signs_d;
    state.best_signs = best_signs;
    state.best_score = R_NegInf;
    state.max_visited = max_visited;
    state.visited = 0.0;

    // Node 0 fixed at +1; search assigns nodes 1..n-1
    if (n > 1) {
        search(&state, 1, 0.0, total_bound);
    } else {
        state.best_score = 0.0;
    }

    // Build return value: the best complete sign assignment found, whether
    // the search exhausted itself (a certified global optimum) or was cut
    // off by `max_visited` first (the best found by then -- there is no
    // fallback method to hand an unfinished search off to, so this is
    // always what `exact_signs()` in R returns)
    SEXP r_signs = PROTECT(allocVector(INTSXP, n));
    int *out_signs = INTEGER(r_signs);
    for (int i = 0; i < n; i++) {
        out_signs[i] = (int) state.best_signs[i];
    }

    UNPROTECT(1);

    free(row_major);
    free(bound_decrement);
    free(signs);
    free(signs_d);
    free(best_signs);

    return r_signs;
}
