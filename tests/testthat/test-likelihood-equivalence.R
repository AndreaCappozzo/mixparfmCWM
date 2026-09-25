# The internal likelihood must reproduce the quantities provided by the
# (modified) parfm package: per-patient marginal terms via
# parfm::loglikelihood_i(), and the group-level marginal likelihood via
# parfm::Mloglikelihood().

test_that("per-patient terms match parfm::loglikelihood_i to 1e-8", {
  skip_if_not_installed("parfm")
  has_lli <- exists(
    "loglikelihood_i",
    where = asNamespace("parfm"),
    inherits = FALSE
  )
  skip_if(!has_lli, "installed parfm lacks loglikelihood_i() (fork-only)")

  sim <- simulate_mixparfm_data(seed = 11)
  d <- sim$data
  formula <- survival::Surv(time, status) ~ x
  design <- stats::model.matrix(formula, d)[, -1L, drop = FALSE]
  n <- nrow(d)

  combos <- expand.grid(
    baseline = c(
      "exponential",
      "weibull",
      "inweibull",
      "gompertz",
      "lognormal",
      "loglogistic"
    ),
    frailty = c("none", "gamma", "ingau"),
    stringsAsFactors = FALSE
  )
  baseline_starts <- list(
    exponential = c(-1.2),
    weibull = c(0.3, -1.6),
    inweibull = c(0.4, -0.9),
    gompertz = c(-1.6, -2.3),
    lognormal = c(0.4, -0.2),
    loglogistic = c(-0.2, 0.3)
  )

  for (draw in seq_len(3)) {
    for (i in seq_len(nrow(combos))) {
      baseline <- combos$baseline[i]
      frailty <- combos$frailty[i]
      p <- c(
        log(0.3 + 0.4 * draw), # frailty variance (ignored when none)
        baseline_starts[[baseline]] + 0.1 * draw,
        0.4 + 0.2 * draw # beta
      )
      if (frailty == "none") {
        p <- p[-1L]
      }

      obs <- parfm::obsdata_creator(
        formula = formula,
        data = d,
        cluster = "family",
        strata = NULL,
        frailty = frailty,
        dist = baseline
      )
      expected <- suppressWarnings(parfm::loglikelihood_i(
        p = p,
        obs = obs,
        dist = baseline,
        frailty = frailty,
        correct = 0,
        transform = TRUE
      ))
      expected <- as.numeric(expected)

      terms <- mixparfmCWM:::.mixparfm_component_terms(
        p,
        design,
        d$time,
        d$status,
        baseline,
        frailty
      )
      theta <- terms$theta
      computed <- terms$event_log_hazard +
        vapply(
          seq_len(n),
          function(k) {
            mixparfmCWM:::.mixparfm_frailty_logLT(
              frailty,
              d$status[k],
              terms$risk[k],
              theta
            )
          },
          0
        )

      expect_true(
        max(abs(computed - expected)) < 1e-8,
        info = paste(baseline, frailty, "draw", draw)
      )
    }
  }
})

test_that("the G = 1 classification loglik matches -Mloglikelihood", {
  skip_if_not_installed("parfm")
  sim <- simulate_mixparfm_data(seed = 12)
  d <- sim$data
  formula <- survival::Surv(time, status) ~ x
  design <- stats::model.matrix(formula, d)[, -1L, drop = FALSE]
  n <- nrow(d)

  for (frailty in c("none", "gamma", "ingau")) {
    p_base <- c(-1.7, 0.3, 0.4) # log rho, log lambda, beta
    p <- if (frailty == "none") p_base else c(0.2, p_base)
    obs <- parfm::obsdata_creator(
      formula = formula,
      data = d,
      cluster = "family",
      strata = NULL,
      frailty = frailty,
      dist = "weibull"
    )
    expected <- -suppressWarnings(as.numeric(parfm::Mloglikelihood(
      p = p,
      obs = obs,
      dist = "weibull",
      frailty = frailty,
      correct = 0
    )))

    terms <- mixparfmCWM:::.mixparfm_component_terms(
      p,
      design,
      d$time,
      d$status,
      "weibull",
      frailty
    )
    group_index <- as.integer(factor(d$family))
    computed <- mixparfmCWM:::.mixparfm_component_survival_loglik(
      terms,
      group_index,
      frailty
    )

    expect_true(
      abs(computed - expected) < 1e-8 * (1 + abs(expected)),
      info = frailty
    )
  }
})
