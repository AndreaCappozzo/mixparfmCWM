# Baseline hazards, frailty Laplace-transform terms and the classification
# log-likelihood of the shared-frailty cluster-weighted model.
#
# The baseline hazard and frailty routines are adapted from the parfm R
# package (version 2.7.6): Federico Rotolo, Marco Munda and Antonio
# Legrand (2012), "Parfm: parametric frailty models in R", Journal of
# Statistical Software 51(11), 1-30. parfm is licensed under GPL-2, which
# is why this package is distributed under GPL-2 as well.
#
# There is a single likelihood implementation in this package: the ICM
# classification step, the component-wise M-step optimizer and
# the objective value reported in the trace all evaluate the same
# functions below.

################################################################################
# Baseline hazard distributions
################################################################################

# Supported baselines and their number of parameters (on the natural scale).
.mixparfm_baseline_npar <- function(baseline) {
  switch(
    baseline,
    exponential = 1L,
    weibull = 2L,
    inweibull = 2L,
    frechet = 2L,
    gompertz = 2L,
    lognormal = 2L,
    loglogistic = 2L,
    stop("Unsupported baseline: ", baseline, call. = FALSE)
  )
}

# Cumulative baseline hazard H0(t) and log baseline hazard log h0(t).
# pars is on the natural scale; the parameter conventions follow parfm.
.mixparfm_baseline_H <- function(baseline, pars, t) {
  switch(
    baseline,
    exponential = pars[1] * t,
    weibull = pars[2] * t^pars[1],
    inweibull = -log(1 - exp(-pars[2] * t^(-pars[1]))),
    frechet = -log(1 - exp(-pars[2] * t^(-pars[1]))),
    gompertz = pars[2] / pars[1] * (exp(pars[1] * t) - 1),
    lognormal = -stats::plnorm(
      t,
      meanlog = pars[1],
      sdlog = pars[2],
      lower.tail = FALSE,
      log.p = TRUE
    ),
    loglogistic = log(1 + exp(pars[1]) * t^pars[2]),
    stop("Unsupported baseline: ", baseline, call. = FALSE)
  )
}

.mixparfm_baseline_lh <- function(baseline, pars, t) {
  switch(
    baseline,
    exponential = log(pars[1]),
    weibull = log(pars[1]) + log(pars[2]) + (pars[1] - 1) * log(t),
    inweibull = log(pars[1]) +
      log(pars[2]) -
      (pars[1] + 1) * log(t) -
      log(exp(pars[2] * t^(-pars[1])) - 1),
    frechet = log(pars[1]) +
      log(pars[2]) -
      (pars[1] + 1) * log(t) -
      log(exp(pars[2] * t^(-pars[1])) - 1),
    gompertz = log(pars[2]) + pars[1] * t,
    lognormal = stats::dlnorm(
      t,
      meanlog = pars[1],
      sdlog = pars[2],
      log = TRUE
    ) -
      stats::plnorm(
        t,
        meanlog = pars[1],
        sdlog = pars[2],
        lower.tail = FALSE,
        log.p = TRUE
      ),
    loglogistic = pars[1] +
      log(pars[2]) +
      (pars[2] - 1) * log(t) -
      log(1 + exp(pars[1]) * t^pars[2]),
    stop("Unsupported baseline: ", baseline, call. = FALSE)
  )
}

# Natural-scale baseline parameters from the optimizer (transformed) scale.
.mixparfm_baseline_from_transformed <- function(baseline, raw) {
  switch(
    baseline,
    exponential = exp(raw),
    weibull = exp(raw),
    inweibull = exp(raw),
    frechet = exp(raw),
    gompertz = exp(raw),
    lognormal = c(raw[1], exp(raw[2])),
    loglogistic = c(raw[1], exp(raw[2])),
    stop("Unsupported baseline: ", baseline, call. = FALSE)
  )
}

################################################################################
# Frailty distributions
################################################################################

# Supported frailties and their number of parameters.
.mixparfm_frailty_npar <- function(frailty) {
  switch(
    frailty,
    none = 0L,
    gamma = 1L,
    ingau = 1L,
    stop("Unsupported frailty: ", frailty, call. = FALSE)
  )
}

# Natural-scale frailty parameter from the optimizer scale.
# Gamma and inverse Gaussian frailties have mean 1 and variance theta.
.mixparfm_frailty_from_transformed <- function(frailty, raw) {
  switch(
    frailty,
    none = NA_real_,
    gamma = exp(raw),
    ingau = exp(raw),
    stop("Unsupported frailty: ", frailty, call. = FALSE)
  )
}

# log[(-1)^k L^(k)(s)] where L is the Laplace transform of the frailty
# distribution and k a non-negative integer. This is the log of the
# group-level survival likelihood contribution
# int_0^Inf m^k exp(-s m) f_M(m) dm
# with k = number of events and s = sum of cumulative hazards in the
# (group, component) cell. Formulas follow parfm's fr.gamma and fr.ingau.
.mixparfm_frailty_logLT <- function(frailty, k, s, theta) {
  if (frailty == "none") {
    return(-s)
  }
  if (frailty == "gamma") {
    if (k == 0) {
      return(-1 / theta * log(1 + theta * s))
    }
    return(
      -(k + 1 / theta) *
        log(1 + theta * s) +
        sum(log(1 + (seq_len(k) - 1) * theta))
    )
  }
  if (frailty == "ingau") {
    # (1 - sqrt(1 + 2*theta*s)) / theta in cancellation-free form,
    # accurate also for very small theta*s
    base <- -2 * s / (1 + sqrt(1 + 2 * theta * s))
    if (k == 0) {
      return(base)
    }
    # z grows like 1/sqrt(theta); the exponentially scaled Bessel function
    # prevents the underflow of besselK(z, k - 0.5) for small theta, which
    # used to create a -Inf cliff that trapped the M-step optimizer
    z <- theta^(-0.5) * sqrt(2 * s + theta^(-1))
    return(
      -k /
        2 *
        log(2 * theta * s + 1) +
        log(besselK(z, k - 0.5, expon.scaled = TRUE)) -
        log(pi / (2 * z)) / 2 +
        base
    )
  }
  stop("Unsupported frailty: ", frailty, call. = FALSE)
}

################################################################################
# Component terms
################################################################################

# Per-observation quantities implied by a component's parameter vector.
# p is on the optimizer scale: [frailty parameter (if any),
# transformed baseline parameters, regression coefficients].
# design is the model matrix WITHOUT the intercept (column order is fixed
# once for the whole data set, so full-data and component-level
# evaluations can never disagree on the design).
.mixparfm_component_terms <- function(
  p,
  design,
  time,
  event,
  baseline,
  frailty
) {
  n_f <- .mixparfm_frailty_npar(frailty)
  n_b <- .mixparfm_baseline_npar(baseline)

  theta <- if (n_f > 0L) {
    .mixparfm_frailty_from_transformed(frailty, p[1])
  } else {
    NA_real_
  }
  baseline_pars <- .mixparfm_baseline_from_transformed(
    baseline,
    p[n_f + seq_len(n_b)]
  )
  beta <- if (n_f + n_b < length(p)) {
    p[(n_f + n_b + 1):length(p)]
  } else {
    numeric(0)
  }

  linear_predictor <- if (length(beta)) {
    as.numeric(design %*% beta)
  } else {
    rep(0, length(time))
  }

  H0 <- .mixparfm_baseline_H(baseline, baseline_pars, time)
  lh0 <- .mixparfm_baseline_lh(baseline, baseline_pars, time)

  list(
    theta = theta,
    baseline_pars = baseline_pars,
    beta = beta,
    risk = as.numeric(H0 * exp(linear_predictor)),
    event_log_hazard = as.numeric(event * (lh0 + linear_predictor)),
    event = as.integer(event)
  )
}

# Additive per-observation, per-component part of the classification
# log-likelihood: mixing weight + covariate densities + event log-hazard.
# event_log_hazard is a n x G matrix, covariate_logdens likewise (or NULL).
.mixparfm_additive_log_matrix <- function(
  tau,
  event_log_hazard,
  covariate_logdens = NULL
) {
  out <- sweep(event_log_hazard, MARGIN = 2L, STATS = log(tau), FUN = "+")
  if (!is.null(covariate_logdens)) {
    out <- out + covariate_logdens
  }
  out
}

################################################################################
# Survival blocks given a partition
################################################################################

# Survival contribution of component g to the classification
# log-likelihood, given the rows currently assigned to g:
# sum_i event_log_hazard_i + sum over (group, component) cells of
# log[(-1)^d L^(d)(s)]. This is the same expression that the ICM
# classification step uses, and the negative of it is the objective of the
# component-wise M-step optimizer.
.mixparfm_component_survival_loglik <- function(terms, group_index, frailty) {
  if (!length(terms$event_log_hazard)) {
    return(0)
  }
  additive <- sum(terms$event_log_hazard)
  frailty_part <- 0
  if (frailty == "none") {
    frailty_part <- .mixparfm_frailty_logLT(
      "none",
      0,
      sum(terms$risk),
      NA_real_
    )
  } else {
    groups <- unique(group_index)
    for (j in groups) {
      index <- group_index == j
      d <- sum(terms$event[index])
      s <- sum(terms$risk[index])
      frailty_part <- frailty_part +
        .mixparfm_frailty_logLT(
          frailty,
          d,
          s,
          terms$theta
        )
    }
  }
  additive + frailty_part
}

################################################################################
# Empirical-Bayes frailty prediction
################################################################################

# Posterior expectation and variance of the shared frailty of a
# (group, component) cell with d events and cumulative hazard sum s:
# E[M]    = exp(logLT(k = d + 1, s) - logLT(k = d, s))
# E[M^2]  = exp(logLT(k = d + 2, s) - logLT(k = d, s))
# This mirrors parfm::predict.parfm() in closed form for all supported
# frailties, since both are ratios of derivatives of the same Laplace
# transform.
.mixparfm_frailty_posterior <- function(frailty, d, s, theta) {
  logLT <- .mixparfm_frailty_logLT
  if (frailty == "none") {
    stop("Frailty predictions require a frailty distribution.", call. = FALSE)
  }
  base <- logLT(frailty, d, s, theta)
  expected <- exp(logLT(frailty, d + 1L, s, theta) - base)
  second_moment <- exp(logLT(frailty, d + 2L, s, theta) - base)
  list(mean = expected, var = second_moment - expected^2)
}

################################################################################
# Classification log-likelihood (the objective of the algorithm)
################################################################################

# Cell sums of events (D) and cumulative hazards (S) for every
# (group, component) pair at a given partition.
.mixparfm_cell_sums <- function(terms, group_index, class, G) {
  n_cells <- max(group_index)
  D <- matrix(0L, nrow = n_cells, ncol = G)
  S <- matrix(0, nrow = n_cells, ncol = G)
  for (i in seq_along(group_index)) {
    j <- group_index[i]
    g <- class[i]
    D[j, g] <- D[j, g] + terms[[g]]$event[i]
    S[j, g] <- S[j, g] + terms[[g]]$risk[i]
  }
  list(D = D, S = S)
}

# The classification log-likelihood of Equation (5) of the manuscript,
# evaluated at the partition state$class and the component parameters in
# state$terms. This is the only objective in the package: the ICM
# classification step, the trace, the BIC and the multi-start selection
# all call this function (or the component-wise blocks it is built from).
.mixparfm_classification_loglik <- function(state) {
  n <- nrow(state$additive)
  cl <- as.integer(state$class)
  value <- sum(state$additive[cbind(seq_len(n), cl)])
  G <- state$G

  if (state$frailty == "none") {
    for (g in seq_len(G)) {
      value <- value +
        .mixparfm_frailty_logLT(
          "none",
          0,
          sum(state$terms[[g]]$risk[cl == g]),
          NA_real_
        )
    }
    return(value)
  }

  cells <- .mixparfm_cell_sums(state$terms, state$group_index, cl, G)
  for (g in seq_len(G)) {
    theta <- state$terms[[g]]$theta
    for (j in seq_len(nrow(cells$D))) {
      if (cells$D[j, g] == 0L && cells$S[j, g] == 0) {
        next
      }
      value <- value +
        .mixparfm_frailty_logLT(
          state$frailty,
          cells$D[j, g],
          cells$S[j, g],
          theta
        )
    }
  }
  value
}

# Human-readable, natural-scale parameter vector of one component.
.mixparfm_transform_parameters <- function(p, baseline, frailty, beta_names) {
  n_f <- .mixparfm_frailty_npar(frailty)
  n_b <- .mixparfm_baseline_npar(baseline)
  baseline_pars <- .mixparfm_baseline_from_transformed(
    baseline,
    p[n_f + seq_len(n_b)]
  )
  baseline_names <- switch(
    baseline,
    exponential = "lambda",
    weibull = c("rho", "lambda"),
    inweibull = c("rho", "lambda"),
    frechet = c("rho", "lambda"),
    gompertz = c("gamma", "lambda"),
    lognormal = c("mu", "sigma"),
    loglogistic = c("alpha", "kappa"),
    stop("Unsupported baseline: ", baseline, call. = FALSE)
  )
  out <- numeric(0)
  if (n_f > 0L) {
    out <- c(
      out,
      setNames(
        .mixparfm_frailty_from_transformed(frailty, p[1]),
        "theta"
      )
    )
  }
  out <- c(out, setNames(baseline_pars, baseline_names))
  if (length(beta_names)) {
    beta <- p[(n_f + n_b + 1):length(p)]
    out <- c(out, setNames(beta, paste0("beta.", beta_names)))
  }
  out
}
