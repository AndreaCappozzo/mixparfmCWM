# ICM classification step (iterated conditional modes) of the ICM-CEM
# algorithm.
#
# Updates one patient label at a time while retaining the joint
# hospital-by-component shared-frailty contribution. Every accepted move
# strictly increases the classification log-likelihood
# (.mixparfm_classification_loglik), and the frailty term is recomputed
# for both the source and the destination (group, component) cells, so
# no per-patient marginal shortcut is used.
#
# Visiting order is the row order of the data; ties retain the current
# label. The endpoint (not the monotonicity) may depend on the visiting
# order, which is why multiple starts are recommended.
#
# Robustness measures against degenerate component parameters:
# - the (group, component) cell sums are rebuilt exactly at the start of
#   every sweep, so incremental drift cannot accumulate across sweeps;
# - if a score evaluates to something non-finite (catastrophic
#   cancellation in the incrementally maintained sums), the affected
#   group's cell sums are rebuilt exactly and the scores recomputed;
# - a label change must strictly improve the score by a relative
#   tolerance (move_tolerance), so rounding-level "improvements"
#   cannot cause spurious moves or cycling.

# state is a list with at least:
#   G, class, frailty, group_index (integer vector of cell indices),
#   terms (list per component, as from .mixparfm_component_terms),
#   additive (n x G matrix from .mixparfm_additive_log_matrix).
.mixparfm_conditional_cstep <- function(
  state,
  max_sweeps = 100L,
  move_tolerance = 1e-9,
  minimum_component_size = 1L,
  minimum_component_events = 0L,
  objective_tolerance = 1e-7
) {
  class <- as.integer(state$class)
  n <- length(class)
  G <- state$G
  group_index <- state$group_index
  events <- state$terms[[1L]]$event
  counts <- tabulate(class, nbins = G)
  # events per component; a patient cannot leave a component if that
  # would leave it with fewer than minimum_component_events events
  # (degenerate zero-event or near-zero-event components are the
  # source of runaway parameter estimates)
  counts_events <- tabulate(class[events == 1L], nbins = G)

  # score[g] is the exact change in the classification log-likelihood
  # of adding patient i to cell (j, g), relative to the common
  # configuration with patient i removed; the frailty terms of both
  # affected cells are recomputed. D and S are the cell sums of events
  # and cumulative hazards in the patient-removed configuration.
  compute_scores <- function(i, j, old) {
    score <- rep(-Inf, G)
    for (g in seq_len(G)) {
      if (
        g != old &&
          (counts[old] <= minimum_component_size ||
            counts_events[old] - events[i] < minimum_component_events)
      ) {
        next
      }
      score[g] <- state$additive[i, g] +
        .mixparfm_frailty_logLT(
          state$frailty,
          D[j, g] + events[i],
          S[j, g] + state$terms[[g]]$risk[i],
          state$terms[[g]]$theta
        ) -
        .mixparfm_frailty_logLT(
          state$frailty,
          D[j, g],
          S[j, g],
          state$terms[[g]]$theta
        )
    }
    score
  }

  objective_before <- .mixparfm_classification_loglik(state)
  sweep_records <- vector("list", max_sweeps)

  for (sweep in seq_len(max_sweeps)) {
    moves <- 0L

    # (group, component) cell sums, rebuilt from scratch at the start of
    # every sweep (O(N), negligible) and maintained incrementally only
    # within the sweep.
    n_cells <- max(group_index)
    D <- matrix(0L, nrow = n_cells, ncol = G)
    S <- matrix(0, nrow = n_cells, ncol = G)
    for (i in seq_len(n)) {
      j <- group_index[i]
      g <- class[i]
      D[j, g] <- D[j, g] + events[i]
      S[j, g] <- S[j, g] + state$terms[[g]]$risk[i]
    }

    for (i in seq_len(n)) {
      j <- group_index[i]
      old <- class[i]
      D[j, old] <- D[j, old] - events[i]
      S[j, old] <- S[j, old] - state$terms[[old]]$risk[i]
      if (is.finite(S[j, old]) && abs(S[j, old]) < 1e-12) {
        S[j, old] <- 0
      }

      score <- compute_scores(i, j, old)

      if (!is.finite(score[old])) {
        # The incrementally maintained cell sums have lost precision
        # (catastrophic cancellation when a component's parameters run
        # away), which would produce NaN scores. Rebuild this group's
        # cell sums exactly from the current labels and recompute.
        others <- which(group_index == j)
        others <- others[others != i]
        for (g in seq_len(G)) {
          in_g <- class[others] == g
          D[j, g] <- sum(events[others][in_g])
          S[j, g] <- sum(state$terms[[g]]$risk[others][in_g])
        }
        score <- compute_scores(i, j, old)
      }
      # any remaining non-finite score is treated as -Inf: the label is
      # never chosen on the basis of a NaN comparison
      score[is.nan(score)] <- -Inf

      best <- which.max(score)
      if (score[best] == -Inf) {
        # no admissible label yields a finite score; keep the incumbent
        best <- old
      } else {
        # strict improvement, relative to the magnitude of the scores
        tolerance <- move_tolerance *
          (if (is.finite(score[old])) 1 + abs(score[old]) else 1)
        if (score[old] >= score[best] - tolerance) best <- old
      }
      if (best != old) {
        moves <- moves + 1L
        counts[old] <- counts[old] - 1L
        counts[best] <- counts[best] + 1L
        counts_events[old] <- counts_events[old] - events[i]
        counts_events[best] <- counts_events[best] + events[i]
        class[i] <- best
      }
      D[j, best] <- D[j, best] + events[i]
      S[j, best] <- S[j, best] + state$terms[[best]]$risk[i]
    }

    sweep_records[[sweep]] <- data.frame(sweep = sweep, moves = moves)
    if (moves == 0L) break
  }

  state_after <- state
  state_after$class <- class
  objective_after <- .mixparfm_classification_loglik(state_after)
  # Single-implementation sanity check: any violation here can only come
  # from floating-point drift, never from a mismatched likelihood.
  if (objective_after < objective_before - objective_tolerance) {
    stop(
      "Internal error: the ICM classification step decreased the ",
      "classification log-likelihood by ",
      sprintf("%.12g", objective_before - objective_after),
      call. = FALSE
    )
  }

  sweeps <- do.call(
    rbind,
    sweep_records[!vapply(sweep_records, is.null, logical(1))]
  )

  list(
    class = class,
    sweeps = sweeps,
    total_moves = sum(sweeps$moves),
    c_step_converged = tail(sweeps$moves, 1L) == 0L,
    objective_before = objective_before,
    objective_after = objective_after
  )
}
