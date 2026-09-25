# Internal fitting of one component's shared-frailty survival model.
#
# Replaces the previous reliance on parfm::parfm(). The optimizer works on
# the same parameter scale and the same likelihood as the classification
# objective (see R/likelihood.R), so the M-step maximizes exactly the
# terms that drive the allocation step.
#
# Robustness measures added after the algorithm audit:
# - warm starts from the current component estimates;
# - a floor on the frailty variance (a warning is emitted when hit);
# - non-finite objective evaluations are penalized instead of crashing;
# - optimizer failures are reported, not thrown, so the caller can keep
#   the previous estimates.

# Deterministic default starting values on the optimizer scale.
# inip (optional) holds [transformed baseline parameters, regression
# coefficients]; iniFpar (optional) is the natural-scale frailty
# parameter. This mirrors the semantics of the old control_parfm().
.mixparfm_default_init <- function(
  time,
  event,
  design,
  baseline,
  frailty,
  inip = NULL,
  iniFpar = NULL
) {
  n_b <- .mixparfm_baseline_npar(baseline)
  n_r <- ncol(design)
  rate <- max(sum(event), 0.5) / max(sum(time), .Machine$double.eps)

  if (is.null(inip)) {
    baseline_raw <- switch(
      baseline,
      exponential = log(rate),
      weibull = c(0, log(rate)),
      inweibull = c(0, log(rate)),
      frechet = c(0, log(rate)),
      gompertz = c(log(0.1), log(rate)),
      lognormal = c(mean(log(time)), log(stats::sd(log(time)) + 0.1)),
      loglogistic = c(log(rate), 0),
      stop("Unsupported baseline: ", baseline, call. = FALSE)
    )
    inip <- c(baseline_raw, rep(0, n_r))
  } else {
    if (length(inip) != n_b + n_r) {
      stop(
        sprintf(
          "'inip' must have length %d (%d baseline + %d regression parameters)",
          n_b + n_r,
          n_b,
          n_r
        ),
        call. = FALSE
      )
    }
  }

  n_f <- .mixparfm_frailty_npar(frailty)
  frailty_raw <- if (n_f > 0L) {
    if (is.null(iniFpar)) {
      0 # theta = 1
    } else {
      log(iniFpar)
    }
  } else {
    numeric(0)
  }

  c(frailty_raw, inip)
}

# Fit one component given the rows assigned to it.
# Returns a list with elements:
#   ok          - TRUE if the optimizer returned usable estimates
#   estim_par   - parameter vector on the optimizer scale
#   loglik      - maximized component survival log-likelihood
#   convergence - optimizer convergence code (0 = success)
#   message     - optimizer message
#   theta_floored - TRUE if the frailty variance hit the floor
.mixparfm_fit_component <- function(
  design,
  time,
  event,
  group_index,
  baseline,
  frailty,
  start,
  method = c("nlminb", "Nelder-Mead", "BFGS"),
  maxit = 500,
  theta_floor = 1e-4
) {
  method <- match.arg(method)
  n_f <- .mixparfm_frailty_npar(frailty)

  # Negative component survival log-likelihood. The frailty parameter is
  # projected onto [log(theta_floor), Inf) so that every method enforces
  # the floor; the projection also guards the degenerate corner where the
  # frailty variance collapses to zero.
  nll <- function(p) {
    if (n_f > 0L) {
      p[1] <- max(p[1], log(theta_floor))
    }
    value <- tryCatch(
      {
        terms <- .mixparfm_component_terms(
          p,
          design,
          time,
          event,
          baseline,
          frailty
        )
        -.mixparfm_component_survival_loglik(terms, group_index, frailty)
      },
      error = function(e) NA_real_
    )
    if (!is.finite(value)) {
      return(1e10)
    }
    value
  }

  run <- tryCatch(
    {
      if (method == "nlminb") {
        lower <- c(
          if (n_f > 0L) log(theta_floor) else numeric(0),
          rep(-Inf, length(start) - n_f)
        )
        stats::nlminb(
          start = start,
          objective = nll,
          lower = lower,
          control = list(iter.max = maxit, eval.max = 10 * maxit)
        )
      } else {
        stats::optim(
          par = start,
          fn = nll,
          method = method,
          control = list(maxit = maxit)
        )
      }
    },
    error = function(e) {
      list(
        par = start,
        value = NA_real_,
        convergence = 100L,
        message = conditionMessage(e)
      )
    }
  )

  if (length(run$par) != length(start) || !all(is.finite(run$par))) {
    return(list(
      ok = FALSE,
      estim_par = start,
      loglik = NA_real_,
      convergence = 100L,
      message = "optimizer returned non-finite estimates",
      theta_floored = FALSE
    ))
  }
  if (n_f > 0L) {
    run$par[1] <- max(run$par[1], log(theta_floor))
  }

  loglik <- -nll(run$par)
  theta_floored <- n_f > 0L &&
    .mixparfm_frailty_from_transformed(frailty, run$par[1]) <=
      theta_floor * (1 + 1e-6)

  list(
    ok = is.finite(loglik),
    estim_par = run$par,
    loglik = loglik,
    convergence = run$convergence,
    message = run$message,
    theta_floored = theta_floored
  )
}
