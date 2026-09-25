# Regression tests for the numerical robustness of the conditional
# classification step (audit 07, findings N1 and N5).

test_that("the C-step survives NaN-level cancellation in the cell sums", {
  set.seed(71)
  n <- 48
  d <- data.frame(
    time = exp(stats::rnorm(n, 0, 1)),
    status = stats::rbinom(n, 1, 0.7),
    x = stats::rnorm(n)
  )
  d$family <- factor(rep(seq_len(12), each = 4))
  design <- stats::model.matrix(~x, d)[, -1L, drop = FALSE]

  G <- 2
  class <- rep(1:2, each = n / 2)
  pars <- list(c(log(691), 0, 0), c(log(1), 0, 0))
  terms <- lapply(seq_len(G), function(g) {
    mixparfmCWM:::.mixparfm_component_terms(
      pars[[g]],
      design,
      d$time,
      d$status,
      "weibull",
      "gamma"
    )
  })

  # group 1, cell (1, 1): three small risks and one enormous one. The
  # cell sum rounds to the largest summand, so when the big patient is
  # removed the small ones are lost to cancellation; subtracting a small
  # patient afterwards used to leave a NEGATIVE cell sum, producing NaN
  # scores in the incrementally maintained implementation.
  terms[[1]]$risk[1:4] <- c(8.6e53, 8.6e53, 1e40, 1e70)
  d$status[1:4] <- 1
  terms[[1]]$event[1:4] <- 1L
  class[1:4] <- 1

  additive <- mixparfmCWM:::.mixparfm_additive_log_matrix(
    tau = c(0.5, 0.5),
    event_log_hazard = vapply(terms, function(t) t$event_log_hazard, numeric(n))
  )
  # lure the big patient away so its cell is left with a stale sum
  additive[4, 2] <- additive[4, 2] + 100

  state <- list(
    G = G,
    class = class,
    frailty = "gamma",
    group_index = as.integer(d$family),
    terms = terms,
    additive = additive
  )

  # the adversarial state intentionally produces NaNs inside the frailty
  # term; the guard detects and repairs them, so warnings are expected
  res <- suppressWarnings(mixparfmCWM:::.mixparfm_conditional_cstep(
    state,
    max_sweeps = 100L,
    move_tolerance = 1e-9,
    minimum_component_size = 1L,
    minimum_component_events = 0L,
    objective_tolerance = 1e-7
  ))

  expect_true(is.finite(res$objective_after))
  expect_gte(res$objective_after, res$objective_before - 1e-7)
  expect_true(isTRUE(res$c_step_converged) || res$total_moves > 0)
})
