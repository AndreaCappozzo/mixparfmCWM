#' Fit Finite Mixtures of Parametric Frailty Models
#'
#' @description
#' Fits finite mixtures of parametric frailty models with possibly random (cluster-weighted)
#' covariates using an EM-type algorithm. The model accommodates hierarchical survival data
#' with a grouping structure (e.g., patients within hospitals) and allows for heterogeneity
#' across latent subpopulations.
#'
#' @param formula A formula specifying the survival model, with a \code{Surv} object on the
#'   left-hand side and covariates on the right-hand side, e.g.,
#'   \code{Surv(time, status) ~ x1 + x2}.
#' @param G Integer. Number of mixture components (clusters) to fit.
#' @param class_init Integer vector of length \code{nrow(data)} specifying initial cluster
#'   assignments. Each element should be an integer between 1 and \code{G}.
#' @param grouping_variable Character string specifying the name of the grouping variable
#'   in \code{data} that defines the hierarchical structure (e.g., hospital ID, family ID).
#'   This variable identifies the groups for which frailty effects are estimated. If NULL,
#'   no frailty structure is assumed.
#' @param strata Character string specifying the stratification variable name, or NULL
#'   for no stratification. Default is NULL.
#' @param X_gaussian_variables Character vector specifying names of continuous covariates
#'   to be modeled as cluster-weighted (random) via Gaussian distributions. If NULL
#'   (default), all covariates are treated as fixed.
#' @param X_multinomial_variables Character vector specifying names of categorical
#'   covariates to be modeled as cluster-weighted (random) via multinomial distributions.
#'   If NULL (default), categorical covariates are treated as fixed.
#' @param data A data frame containing the variables specified in \code{formula},
#'   \code{grouping_variable}, and the cluster-weighted covariate specifications.
#' @param baseline Character string specifying the baseline hazard distribution. Options are:
#'   \describe{
#'     \item{\code{"weibull"}}{Weibull distribution (default)}
#'     \item{\code{"exponential"}}{Exponential distribution}
#'     \item{\code{"inweibull"}}{Inverse Weibull distribution}
#'     \item{\code{"frechet"}}{Fréchet distribution}
#'     \item{\code{"gompertz"}}{Gompertz distribution}
#'     \item{\code{"loglogistic"}}{Log-logistic distribution}
#'     \item{\code{"lognormal"}}{Log-normal distribution}
#'     \item{\code{"logskewnormal"}}{Log-skew-normal distribution}
#'   }
#' @param frailty Character string specifying the frailty distribution. Options are:
#'   \describe{
#'     \item{\code{"none"}}{No frailty (independence model)}
#'     \item{\code{"gamma"}}{Gamma frailty (default)}
#'     \item{\code{"ingau"}}{Inverse Gaussian frailty}
#'     \item{\code{"possta"}}{Positive stable frailty}
#'     \item{\code{"lognormal"}}{Log-normal frailty}
#'     \item{\code{"loglogistic"}}{Log-logistic frailty}
#'   }
#' @param control_EM_algorithm A list of control parameters for the EM algorithm,
#'   as returned by \code{\link{control_EM}}. Default is \code{control_EM()}.
#' @param control_parfm_algorithm A list of control parameters for fitting the
#'   parametric frailty models, as returned by \code{\link{control_parfm}}.
#'   Default is \code{control_parfm()}.
#'
#' @return An object of class \code{mixparfm} (a list) with components:
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{parameters}{A list containing:
#'     \describe{
#'       \item{\code{tau}}{Vector of mixing proportions (length \code{G}).}
#'       \item{\code{AFT_parameters}}{Matrix of AFT model parameters (baseline and
#'         regression coefficients) with one column per mixture component.}
#'       \item{\code{X_gaussian_parameters}}{List with \code{mu} (means) and
#'         \code{sigma} (covariances) for Gaussian cluster-weighted covariates
#'         (only present if \code{X_gaussian_variables} is not NULL).}
#'       \item{\code{X_multinomial_parameters}}{List of probability matrices for
#'         categorical cluster-weighted covariates (only present if
#'         \code{X_multinomial_variables} is not NULL).}
#'     }
#'   }
#'   \item{frailty_effect}{Data frame with estimated frailty effects for each group
#'     and mixture component (NULL if \code{frailty = "none"}).}
#'   \item{frailty_var_effect}{Data frame with variances of frailty effects
#'     (NULL if \code{frailty = "none"}).}
#'   \item{z}{Matrix of posterior probabilities of cluster membership (N x G).}
#'   \item{class}{Vector of MAP cluster assignments (length N).}
#'   \item{bic}{Bayesian Information Criterion value.}
#'   \item{baseline}{The baseline hazard distribution used.}
#'   \item{frailty}{The frailty distribution used.}
#'   \item{fit_parfm}{List of length \code{G} containing the fitted \code{parfm}
#'     objects for each mixture component.}
#'   \item{loglik_vec}{Vector of log-likelihood values at each EM iteration.}
#'
#' @details
#' This function implements a finite mixture of parametric frailty models suitable for
#' clustered survival data. The model extends standard mixture models by:
#' \enumerate{
#'   \item Incorporating a frailty (random effect) structure within each mixture component
#'     to account for within-group correlation.
#'   \item Allowing covariates to be cluster-specific (cluster-weighted), meaning their
#'     distributions can vary across mixture components.
#' }
#'
#' The EM algorithm alternates between:
#' \itemize{
#'   \item \strong{E-step}: Computing posterior probabilities of cluster membership
#'   \item \strong{M-step}: Updating mixture component parameters by fitting separate
#'     parametric frailty models via the \code{parfm} package
#' }
#'
#' The function uses either Classification EM (CEM) with hard assignments or
#' Stochastic EM (SEM) with probabilistic assignments, as specified in
#' \code{control_EM_algorithm}.
#'
#' @note
#' \itemize{
#'   \item This package requires the modified \code{parfm} package available at
#'     \url{https://github.com/AndreaCappozzo/parfm}.
#'   \item Character variables in \code{data} are automatically converted to factors.
#'   \item The function handles cases where some factor levels may be absent in certain
#'     mixture components.
#' }
#'
#' @references
#' Caldera, Cappozzo, Masci, Forlani, Antonelli, Leoni, Paganoni, Ieva (2025+).
#' Cluster-weighted modeling of lifetime hierarchical data for profiling COVID-19
#' heart failure patients. \url{https://arxiv.org/abs/2507.12230}
#'
#' @seealso \code{\link{control_EM}}, \code{\link{control_parfm}}, \code{\link[parfm]{parfm}},
#'   \code{\link[survival]{Surv}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data (using frailtySurv package)
#' library(frailtySurv)
#' set.seed(46)
#'
#' # Generate two subpopulations with different characteristics
#' df1 <- genfrail(N = 10, K = 40, beta = c(2, 3),
#'                 frailty = "gamma", theta = 6,
#'                 covar.matrix = matrix(rnorm(400*2, c(2,3), c(1,2)), ncol=2),
#'                 Lambda_0 = function(t, lambda=0.2, rho=3) lambda*t^rho)
#'
#' df2 <- genfrail(N = 10, K = 50, beta = c(-2, -3),
#'                 frailty = "gamma", theta = 0.1,
#'                 covar.matrix = matrix(rnorm(500*2, c(4,5), c(2,1)), ncol=2),
#'                 Lambda_0 = function(t, lambda=0.2, rho=15) lambda*t^rho)
#'
#' df <- rbind(df1, df2)
#' true_class <- rep(1:2, c(nrow(df1), nrow(df2)))
#'
#' # Initialize cluster assignments with some noise
#' class_init <- true_class
#' class_init[1:450] <- sample(1:2, 450, replace = TRUE, prob = c(0.8, 0.2))
#' class_init[451:900] <- sample(1:2, 450, replace = TRUE, prob = c(0.2, 0.8))
#'
#' # Fit mixture model with fixed covariates
#' fit_fixed <- fit_mixparfm(
#'   formula = Surv(time, status) ~ Z1 + Z2,
#'   G = 2,
#'   data = df,
#'   class_init = class_init,
#'   grouping_variable = "family",
#'   baseline = "weibull",
#'   frailty = "gamma",
#'   control_EM_algorithm = control_EM(itermax = 50),
#'   control_parfm_algorithm = control_parfm(maxit = 1000)
#' )
#'
#' # Check convergence
#' plot(fit_fixed$loglik_vec, type = "l",
#'      xlab = "Iteration", ylab = "Log-likelihood")
#'
#' # Fit mixture model with cluster-weighted covariates
#' fit_cw <- fit_mixparfm(
#'   formula = Surv(time, status) ~ Z1 + Z2,
#'   G = 2,
#'   data = df,
#'   class_init = class_init,
#'   grouping_variable = "family",
#'   X_gaussian_variables = c("Z1", "Z2"),
#'   baseline = "weibull",
#'   frailty = "gamma",
#'   control_EM_algorithm = control_EM(itermax = 50),
#'   control_parfm_algorithm = control_parfm()
#' )
#'
#' # Compare classifications
#' table(Estimated = fit_cw$class, True = true_class)
#'
#' # View parameter estimates
#' fit_cw$parameters
#'
#' # View BIC
#' fit_cw$bic
#' }
fit_mixparfm <-
  function(
    formula,
    G,
    class_init,
    grouping_variable = NULL,
    strata = NULL,
    X_gaussian_variables = NULL,
    X_multinomial_variables = NULL,
    # modelXnorm = "VVV", # FIXME think if we want parsimonious covariance structures for continuous covariates
    data,
    baseline = c(
      "weibull",
      "inweibull",
      "frechet",
      "exponential",
      "gompertz",
      "loglogistic",
      "lognormal",
      "logskewnormal"
    ),
    frailty = c("none", "gamma", "ingau", "possta", "lognormal", "loglogistic"),
    control_EM_algorithm = control_EM(),
    control_parfm_algorithm = control_parfm()
  ) {
    if (missing(data)) {
      data <- eval(parse(
        text = paste(
          "data.frame(",
          paste(all.vars(formula), collapse = ", "),
          ")"
        )
      ))
    }
    data <- as.data.frame(unclass(data), stringsAsFactors = TRUE) # converts all character variables to factors
    # (this automatically handles situations in which for some clusters not all the levels of the qualitative variable are observed)

    N <- nrow(data)
    # Retrieve marginal variables
    # Continous covariates
    is_X_gaussian <- !is.null(X_gaussian_variables)
    is_X_multinomial <- !is.null(X_multinomial_variables)
    if (is_X_gaussian) {
      X_gaussian <- data[, X_gaussian_variables]
      p_gaussian <- length(X_gaussian_variables)
      mclust_model_name <- ifelse(p_gaussian == 1, "V", "VVV")
    } else {
      # if no continuous random covariates are modeled, the associated loglik contribution is set to 0 and will stay like this
      # throughout the algorithm
      X_gaussian_log_density <- 0
      log_density_X_gaussian_g <- 0
      n_par_X_gaussian <- 0
    }

    if (is_X_multinomial) {
      X_multinomial <- data[, X_multinomial_variables, drop = FALSE]

      p_multinomial <- length(X_multinomial_variables) # the number of variables modeled independently via multinomial dist (also denoted with M hereafter)
      C_M_vec <- apply(X_multinomial, 2, function(X_mult_col) {
        length(unique(X_mult_col))
      })
      X_multinomial_params <- vector(mode = "list", length = p_multinomial) # container for the C_m x G matrices of probabilities of variable m
      # C_m is the # of categories for variable m, m=1,...M aka p_multinomial (see notation in Section 3 of Variable selection for latent class analysis with application to low back pain diagnosis Fop et al)
      # I need to store them in a list since C_m (the total # of categories) may differ for the p_multinomial variables
      for (m in 1:p_multinomial) {
        X_multinomial_params[[m]] <- matrix(nrow = C_M_vec[m], ncol = G)
        rownames(X_multinomial_params[[m]]) <- names(table(X_multinomial[, m]))
      }
    } else {
      # if no multinomial random covariates are modeled, the associated loglik contribution is set to 0 and will stay like this
      # throughout the algorithm
      X_multinomial_log_density <- 0
      log_density_X_multinomial_g <- 0
      n_par_X_multinomial <- 0
    }

    # EM controls
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err
    E_step_update <- control_EM_algorithm$E_step_update

    # parfm controls
    inip = control_parfm_algorithm$inip
    iniFpar = control_parfm_algorithm$iniFpar
    method = control_parfm_algorithm$method
    maxit = control_parfm_algorithm$maxit
    Fparscale = control_parfm_algorithm$Fparscale
    showtime = control_parfm_algorithm$showtime
    correct = control_parfm_algorithm$correct

    # While loop controls
    iter <- 0
    loglik <- loglik_prev <- -.Machine$integer.max / 2
    loglik_vec <- NULL
    crit <- TRUE

    # Parameters initialization according to class_init
    z_init <- mclust::unmap(class_init)
    tau <- colMeans(z_init)
    parfm_params <- vector("list", length = G)
    fit_parfm_g_list <- vector("list", length = G) # container for the G estimated models via parfm (I need them in the end for the frailty calc)

    # obsdata object needed for log likelihood computation from the parfm package:
    # it will automatically account for different
    # dist and frailties
    obsdata <-
      parfm::obsdata_creator(
        formula = formula,
        data = data,
        cluster = grouping_variable,
        strata = strata,
        frailty = frailty,
        dist = baseline
      )

    for (g in 1:G) {
      data_g <- data[class_init == g, ]

      fit_parfm_g_list[[g]] <-
        parfm::parfm(
          formula = formula,
          cluster = grouping_variable,
          strata = strata,
          inip = inip,
          iniFpar = iniFpar,
          maxit = maxit,
          Fparscale = Fparscale,
          showtime = showtime,
          correct = correct,
          data = data_g,
          dist = baseline,
          frailty = frailty,
          method = method
        )
      parfm_params[[g]] <- attributes(fit_parfm_g_list[[g]])$estim_par

      if (is_X_multinomial) {
        X_multinomial_g <- X_multinomial[class_init == g, , drop = FALSE]
        for (m in 1:p_multinomial) {
          tab_X_multinomial_gm <- table(factor(
            X_multinomial_g[, m],
            levels = rownames(X_multinomial_params[[m]])
          ))
          # this handles the case in which for some levels
          # of the qualitative variable there are 0 entries in cluster g
          X_multinomial_params[[m]][, g] <- tab_X_multinomial_gm /
            sum(tab_X_multinomial_gm)
        }
      }
    }

    # 3 ways to compute the same quantities
    # X_gaussian_params  <- flexCWM:::.PX_norm(colXn=ncol(X_gaussian),Xnorm=X_gaussian,modelXnorm="VVV",z=z_init,k=G,n=nrow(data),eps=10e-6)
    # X_gaussian_params <- mclust::covw(X=X_gaussian,Z=z_init, normalize = FALSE)
    if (is_X_gaussian) {
      X_gaussian_params <-
        mclust::mstep(
          data = X_gaussian,
          z = z_init,
          modelName = mclust_model_name
        )
    }

    while (crit) {
      # E step ------------------------------------------------------------------

      time_log_density <-
        sapply(1:G, function(g) {
          parfm::loglikelihood_i(
            p = parfm_params[[g]],
            obs = obsdata,
            dist = baseline,
            frailty = frailty,
            correct = correct,
            transform = TRUE
          )
        })
      if (is_X_gaussian) {
        X_gaussian_log_density <- tryCatch(
          do.call(
            mclust::cdens,
            c(
              list(data = X_gaussian, logarithm = TRUE),
              X_gaussian_params
            )
          ),
          error = function(e) {
            -Inf
          }
        )
      }

      if (is_X_multinomial) {
        X_multinomial_density_list <- lapply(1:p_multinomial, function(m) {
          sapply(1:G, function(g) {
            X_multinomial_params[[m]][, g][as.character(X_multinomial[, m])]
          })
        })
        X_multinomial_log_density <-
          Reduce(f = "+", x = lapply(X_multinomial_density_list, log))
      }

      log_comp_mixture <-
        sweep(
          (time_log_density +
            X_gaussian_log_density +
            X_multinomial_log_density),
          MARGIN = 2,
          STATS = log(tau),
          FUN = "+"
        )

      z_max <- apply(log_comp_mixture, 1, max)
      log_density <-
        z_max + log(rowSums(exp(log_comp_mixture - z_max)))
      log_z <- log_comp_mixture - log_density
      z <- exp(log_z)

      # Classification (CEM) or stochastic (SEM) step
      if (E_step_update == "classification") {
        map_z <- mclust::map(z)
      } else if (E_step_update == "stochastic") {
        map_z <- mclust::map(t(apply(z, 1, function(prob) {
          stats::rmultinom(1, 1, prob)
        })))
        # FIXME check the "average" strategy of weibullRMM_SEM in mixtools \cite{Bordes2016} to see whether it may improve the estimates
      }

      z <- mclust::unmap(map_z)

      data_g_list <- split(data, map_z)

      # M step ------------------------------------------------------------------

      tau <- colMeans(z)

      for (g in 1:G) {
        fit_parfm_g_list[[g]] <-
          parfm::parfm(
            formula = formula,
            cluster = grouping_variable,
            strata = strata,
            inip = inip,
            iniFpar = iniFpar,
            maxit = maxit,
            Fparscale = Fparscale,
            showtime = showtime,
            correct = correct,
            data = data_g_list[[g]],
            dist = baseline,
            frailty = frailty,
            method = method
          )
        parfm_params[[g]] <- attributes(fit_parfm_g_list[[g]])$estim_par

        if (is_X_multinomial) {
          # MLE for categorical covariates
          X_multinomial_g_list <- split(X_multinomial, map_z)
          for (m in 1:p_multinomial) {
            tab_X_multinomial_gm <- table(factor(
              X_multinomial_g_list[[g]][, m],
              levels = rownames(X_multinomial_params[[m]])
            ))
            X_multinomial_params[[m]][, g] <- tab_X_multinomial_gm /
              sum(tab_X_multinomial_gm)
          }
        }
      }

      if (is_X_gaussian) {
        # MLE for continuous covariates
        X_gaussian_params <-
          mclust::mstep(data = X_gaussian, z = z, modelName = mclust_model_name)
      }

      # Log likelihood ----------------------------------------------

      log_density_parfm_g <- numeric(length = G)
      for (g in 1:G) {
        obsdata_g <-
          parfm::obsdata_creator(
            formula = formula,
            data = data_g_list[[g]],
            cluster = grouping_variable,
            strata = strata,
            frailty = frailty,
            dist = baseline
          )

        log_density_parfm_g[g] <-
          -parfm::Mloglikelihood(
            p = parfm_params[[g]],
            obs = obsdata_g,
            dist = baseline,
            frailty = frailty,
            correct = correct
          )
      }
      if (is_X_gaussian) {
        log_density_X_gaussian_g <-
          colSums(
            do.call(
              mclust::cdens,
              c(
                list(data = X_gaussian, logarithm = TRUE),
                X_gaussian_params
              )
            ) *
              z
          )
      }

      if (is_X_multinomial) {
        X_multinomial_density_list <- lapply(1:p_multinomial, function(m) {
          sapply(1:G, function(g) {
            X_multinomial_params[[m]][, g][as.character(X_multinomial[, m])]
          })
        })
        X_multinomial_log_dens <-
          (Reduce(f = "+", x = lapply(X_multinomial_density_list, log))) * z
        X_multinomial_log_dens_no_Nan <- ifelse(
          is.nan(X_multinomial_log_dens),
          no = X_multinomial_log_dens,
          yes = 0
        )
        # this handles the case in which for some levels
        # of the qualitative variable there are 0 entries in cluster g
        log_density_X_multinomial_g <- colSums(X_multinomial_log_dens_no_Nan)
      }

      n_g <- colSums(z)

      loglik <-
        sum(
          log_density_parfm_g +
            log_density_X_gaussian_g +
            log_density_X_multinomial_g +
            log(tau) * n_g
        )

      # Check convergence
      err <- abs(loglik - loglik_prev) / (1 + abs(loglik))
      loglik_prev <- loglik
      loglik_vec <- c(loglik_vec, loglik)
      iter <- iter + 1

      crit <- err > tol & iter < itermax
    }

    parfm_params_transformed <-
      lapply(1:G, function(g) {
        parfm::paramaters_transformator(
          estim_par = parfm_params[[g]],
          frailty = frailty,
          dist = baseline,
          obsdata = obsdata
        )
      })
    if (G == 1) {
      AFT_parameters <- as.matrix(Reduce(x = parfm_params_transformed, "rbind"))
    } else {
      AFT_parameters <- t(Reduce(x = parfm_params_transformed, "rbind"))
    }

    colnames(AFT_parameters) <- 1:G
    parameters <-
      list(
        tau = tau,
        AFT_parameters = AFT_parameters
      )

    if (is_X_gaussian) {
      parameters$X_gaussian_parameters = list(
        mu = X_gaussian_params$parameters$mean,
        sigma = X_gaussian_params$parameters$variance$sigma
      )
    }

    if (is_X_multinomial) {
      parameters$X_multinomial_parameters = X_multinomial_params
    }

    if (frailty == "none") {
      # in case no frailty is considered simply return NULL for frailty_effect (this is necessary to prevent errors in the call to predict.parfm)
      frailty_effect <- NULL
      frailty_var_effect <- NULL
    } else {
      frailty_effect_df_list <- vector(mode = "list", length = G)
      frailty_var_effect_df_list <- vector(mode = "list", length = G)

      for (g in 1:G) {
        # all this mess is done as we may have some clusters for which no obs fall within a given group, and I want to keep them all
        frailty_g <- parfm::predict.parfm(fit_parfm_g_list[[g]])[[1]]
        frailty_var <- parfm::predict.parfm(fit_parfm_g_list[[g]])[[2]]

        frailty_effect_df_list[[g]] <-
          data.frame(
            group = attributes(frailty_g)$names,
            frailty = c(frailty_g)
          )
        colnames(frailty_effect_df_list[[g]])[2] <- g

        frailty_var_effect_df_list[[g]] <-
          data.frame(
            group = attributes(frailty_var)$names,
            frailty_var = c(frailty_var)
          )
        colnames(frailty_var_effect_df_list[[g]])[2] <- g
      }

      frailty_effect = plyr::join_all(
        frailty_effect_df_list,
        by = "group",
        type = "full"
      )
      frailty_var_effect = plyr::join_all(
        frailty_var_effect_df_list,
        by = "group",
        type = "full"
      )

      # if no frailty_effect is present for a given group in a cluster (i.e., no obs from that group belong to the g-th cluster) it returns NA
      colnames(frailty_effect)[1] <-
        grouping_variable # I manually specify the grouping variable name

      colnames(frailty_var_effect)[1] <-
        grouping_variable
    }

    # Compute bic

    n_par_tau <- G - 1
    n_par_parfm <- length(parameters$AFT_parameters)
    if (is_X_gaussian) {
      n_par_X_gaussian <-
        mclust::nMclustParams(
          modelName = mclust_model_name,
          d = p_gaussian,
          G = G
        ) -
        n_par_tau # nMclustParams counts also the G-1 mixture weights, I subtract them otherwise I would count them twice
    }
    if (is_X_multinomial) {
      n_par_X_multinomial <- sum(sapply(X_multinomial_params, length) - G)
    }
    bic_final <- 2 *
      loglik -
      (n_par_tau + n_par_parfm + n_par_X_gaussian + n_par_X_multinomial) *
        log(N)

    res <-
      list(
        loglik = loglik,
        parameters = parameters,
        frailty_effect = frailty_effect,
        frailty_var_effect = frailty_var_effect,
        z = z,
        class = mclust::map(z),
        bic = bic_final,
        baseline = baseline,
        frailty = frailty,
        fit_parfm = fit_parfm_g_list, # #FIXME add droplevels + formula_g a list of dimension G containing the estimated parametric frailty models in the G clusters. Useful for assessing significance
        loglik_vec = loglik_vec
      )
  }
