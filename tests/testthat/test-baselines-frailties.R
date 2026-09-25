# Every supported baseline x frailty combination must fit without errors
# and produce a non-decreasing objective trace.

test_that("all baseline x frailty combinations fit and stay monotone", {
  sim <- simulate_mixparfm_data(seed = 61, groups_per_component = 3)
  baselines <- c(
    "weibull",
    "exponential",
    "inweibull",
    "frechet",
    "gompertz",
    "lognormal",
    "loglogistic"
  )
  frailties <- c("none", "gamma", "ingau")

  for (baseline in baselines) {
    for (frailty in frailties) {
      set.seed(61)
      fit <- suppressWarnings(mixparfmCWM::fit_mixparfm(
        formula = survival::Surv(time, status) ~ x,
        G = 2,
        class_init = sim$true_class,
        grouping_variable = if (frailty == "none") "family" else "family",
        data = sim$data,
        baseline = baseline,
        frailty = frailty,
        control_EM_algorithm = mixparfmCWM::control_EM(
          itermax = 3,
          tol = 1e-4
        ),
        control_parfm_algorithm = mixparfmCWM::control_parfm(maxit = 500)
      ))
      expect_true(
        is.finite(fit$classification_loglik),
        info = paste(baseline, frailty)
      )
      expect_true(
        all(fit$cem_trace$total_delta >= -1e-8),
        info = paste(baseline, frailty)
      )
    }
  }
})
