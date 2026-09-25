# Shared simulation helper for the test suite.
# Two well-separated components, Weibull baseline, gamma shared frailty,
# one survival covariate, two Gaussian and one multinomial cluster-weighted
# covariates.
simulate_mixparfm_data <- function(
  n_per_group = 4,
  groups_per_component = 6,
  sep = 2.5,
  theta = 0.5,
  seed = 1
) {
  set.seed(seed)
  dat_list <- list()
  true_class <- integer()
  for (g in 1:2) {
    for (j in seq_len(groups_per_component)) {
      m_j <- stats::rgamma(n_per_group, shape = 1 / theta, scale = theta)
      lam <- c(0.4, 0.12)[g]
      rho <- c(1.3, 1.7)[g]
      beta <- c(0.3, 0.5)[g]
      x <- stats::rnorm(n_per_group)
      t_ev <- ((-log(stats::runif(n_per_group))) /
        (m_j * lam * exp(beta * x)))^(1 / rho)
      cens <- stats::runif(n_per_group, 2, 9)
      dat_list[[length(dat_list) + 1]] <- data.frame(
        time = pmin(t_ev, cens),
        status = as.integer(t_ev <= cens),
        x = x,
        Z1 = stats::rnorm(n_per_group, c(-sep / 2, sep / 2)[g]),
        Z2 = stats::rnorm(n_per_group, c(1, -1)[g]),
        V = factor(
          stats::rbinom(n_per_group, 1, c(0.3, 0.6)[g])
        ),
        family = paste0("f", g, j)
      )
      true_class <- c(true_class, rep(g, n_per_group))
    }
  }
  d <- do.call(rbind, dat_list)
  d$family <- factor(d$family)
  list(data = d, true_class = true_class)
}

fit_default <- function(data, class_init, G = 2, ...) {
  mixparfmCWM::fit_mixparfm(
    formula = survival::Surv(time, status) ~ x,
    G = G,
    class_init = class_init,
    grouping_variable = "family",
    X_gaussian_variables = c("Z1", "Z2"),
    X_multinomial_variables = "V",
    data = data,
    ...
  )
}
