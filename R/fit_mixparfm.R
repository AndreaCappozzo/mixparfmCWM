#' Fit Finite Mixtures of Parametric Frailty Models
#'
#' @description
#' Fits finite mixtures of parametric frailty models with possibly random (cluster-weighted)
#' covariates using an ICM-CEM algorithm: a classification EM (CEM) algorithm whose
#' classification step is carried out by iterated conditional modes (ICM). The model accommodates
#' hierarchical survival data with a grouping structure (e.g., patients within hospitals)
#' and allows for heterogeneity across latent subpopulations.
#'
#' @param formula A formula specifying the survival model, with a \code{Surv} object on the
#'   left-hand side and covariates on the right-hand side, e.g.,
#'   \code{Surv(time, status) ~ x1 + x2}. Only right-censored responses are
#'   supported; left-truncated (counting-process) \code{Surv} objects raise an error.
#' @param G Integer. Number of mixture components (clusters) to fit.
#' @param class_init Integer vector of length \code{nrow(data)} specifying initial cluster
#'   assignments. Each element should be an integer between 1 and \code{G}. Used for the
#'   first start; further starts (see \code{n_start} in \code{\link{control_EM}}) use
#'   random partitions.
#' @param grouping_variable Character string specifying the name of the grouping variable
#'   in \code{data} that defines the hierarchical structure (e.g., hospital ID, family ID).
#'   This variable identifies the groups for which frailty effects are estimated. Required
#'   (not NULL) unless \code{frailty = "none"}.
#' @param strata Character string specifying the stratification variable name, or NULL
#'   for no stratification. Default is NULL. Stratified models are currently not supported
#'   and a non-NULL value raises an error.
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
#'     \item{\code{"frechet"}}{Fréchet distribution (same as \code{"inweibull"})}
#'     \item{\code{"gompertz"}}{Gompertz distribution}
#'     \item{\code{"loglogistic"}}{Log-logistic distribution}
#'     \item{\code{"lognormal"}}{Log-normal distribution}
#'   }
#' @param frailty Character string specifying the frailty distribution. Options are:
#'   \describe{
#'     \item{\code{"none"}}{No frailty (independence model)}
#'     \item{\code{"gamma"}}{Gamma frailty with mean 1 and variance theta (default)}
#'     \item{\code{"ingau"}}{Inverse Gaussian frailty with mean 1 and variance theta}
#'   }
#' @param control_EM_algorithm A list of control parameters for the ICM-CEM
#'   algorithm, as returned by \code{\link{control_EM}}. Default is \code{control_EM()}.
#' @param control_parfm_algorithm A list of control parameters for the component-wise
#'   parametric frailty fits, as returned by \code{\link{control_parfm}}.
#'   Default is \code{control_parfm()}.
#'
#' @return A list with components:
#'   \item{loglik}{Final classification log-likelihood.}
#'   \item{classification_loglik}{Final classification log-likelihood (same as
#'     \code{loglik}; the name is kept to be explicit about the optimized criterion).}
#'   \item{parameters}{A list containing:
#'     \describe{
#'       \item{\code{tau}}{Vector of mixing proportions (length \code{G}).}
#'       \item{\code{AFT_parameters}}{Matrix of natural-scale parameters (frailty
#'         variance, baseline and regression coefficients) with one column per
#'         mixture component.}
#'       \item{\code{X_gaussian_parameters}}{List with \code{mu} (means) and
#'         \code{sigma} (covariances) for Gaussian cluster-weighted covariates
#'         (only present if \code{X_gaussian_variables} is not NULL).}
#'       \item{\code{X_multinomial_parameters}}{List of probability matrices for
#'         categorical cluster-weighted covariates
#'         (only present if \code{X_multinomial_variables} is not NULL).}
#'     }
#'   }
#'   \item{frailty_effect}{Data frame with posterior means of the frailty for each group
#'     and mixture component (NULL if \code{frailty = "none"}); groups with no patient
#'     assigned to a component have NA entries.}
#'   \item{frailty_var_effect}{Data frame with posterior variances of the frailty effects
#'     (NULL if \code{frailty = "none"}).}
#'   \item{z}{Matrix (N x G) of hard cluster-membership indicators.}
#'   \item{class}{Vector of hard cluster assignments (length N).}
#'   \item{bic}{Classification BIC, see Details.}
#'   \item{classification_bic}{Same as \code{bic}.}
#'   \item{baseline}{The baseline hazard distribution used.}
#'   \item{frailty}{The frailty distribution used.}
#'   \item{fit_parfm}{List of length \code{G} with one internal fit summary per component
#'     (estimates on the optimizer scale, component log-likelihood, convergence code,
#'     optimizer message and frailty-floor flag). These are plain lists, not
#'     \code{parfm} objects.}
#'   \item{loglik_vec}{Classification log-likelihood trace, one value per iteration.}
#'   \item{classification_loglik_vec}{Same as \code{loglik_vec}.}
#'   \item{cem_trace}{Data frame with the objective before the ICM
#'     classification step, after it, and after the M-step, at every iteration.}
#'   \item{objective_name}{Name of the optimized objective.}
#'   \item{algorithm}{Name of the fitted algorithm.}
#'   \item{converged}{Logical convergence indicator based on the configured tolerance.}
#'   \item{stopping_reason}{Character description of why the algorithm stopped.}
#'   \item{n_iter}{Number of completed iterations of the selected start.}
#'   \item{reached_itermax}{Logical indicator for stopping at the iteration limit.}
#'   \item{final_relative_change}{Final relative change of the classification
#'     log-likelihood.}
#'   \item{tol_zero_var}{Gaussian covariance eigenvalue tolerance used by the fit.}
#'   \item{n_start}{Number of starts requested.}
#'   \item{selected_start}{Index of the start attaining the highest classification
#'     log-likelihood.}
#'   \item{start_summary}{Data frame with the final objective, number of iterations,
#'     convergence flag, error message, degeneracy flag and degeneracy reason
#'     of every start. Degenerate starts are excluded from the selection.}
#'   \item{n_degenerate_starts}{Number of starts that ended in degenerate
#'     solutions and were excluded from the selection.}
#'   \item{n_m_step_rejections}{Number of M-step candidates rejected because they
#'     decreased the component survival log-likelihood.}
#'   \item{theta_floor}{The frailty-variance floor used in the component fits.}
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
#' The estimation criterion is the classification (complete-data) log-likelihood.
#' Because the integrated shared-frailty term couples the labels of patients in the
#' same group-by-component cell, the classification log-likelihood does not separate
#' across patients, and the usual C-step of CEM—assigning each unit by its marginal
#' posterior probability—is not available. The model is therefore fitted with an
#' \strong{ICM-CEM} algorithm: a classification EM (CEM) algorithm (Celeux & Govaert,
#' 1992) whose classification step is carried out by iterated conditional modes
#' (ICM; Besag, 1986). Each label is set in turn to the mode of its full conditional
#' distribution given the current labels of all other patients and the current
#' parameter estimates, and the new label is used immediately. The ICM-CEM algorithm
#' alternates between:
#' \itemize{
#'   \item \strong{ICM classification step}: Updating one hard assignment at a
#'     time. Each update recomputes the joint group-by-component shared-frailty
#'     contribution of both the source and the destination cells, so every accepted move
#'     increases the classification log-likelihood.
#'   \item \strong{M-step}: Updating mixture component parameters. Mixing proportions,
#'     Gaussian and multinomial covariate parameters have closed-form updates given the
#'     hard partition; the survival parameters of each component are refitted by
#'     numerical maximization of the same component survival log-likelihood used by the
#'     classification step. Refits are warm-started from the current estimates, and a
#'     candidate that does not improve a component's contribution is rejected (the
#'     previous estimates are retained).
#' }
#'
#' Every ICM update and every block of the M-step leaves the classification
#' log-likelihood non-decreasing, so ICM-CEM is a monotone block-coordinate ascent
#' algorithm that converges to a partial optimum: a partition that no single
#' reassignment improves at the final parameters, together with parameters that locally
#' maximise the classification log-likelihood given that partition. As with any
#' CEM-type algorithm the solution is local, so multiple random starts are used.
#'
#' The \code{bic} element is a \emph{classification} BIC: \code{2 * loglik - d * log(N)}
#' with the classification (not the observed-data marginal) log-likelihood and
#' \code{d} the total number of parameters. As such it behaves like an ICL-type criterion
#' and tends to favour fewer, well-separated components; comparisons across \code{G},
#' baselines or frailties are internally consistent, but it is not the Schwarz BIC of the
#' marginal mixture likelihood.
#'
#' The ICM classification step is monotone but converges to a coordinate-wise
#' local optimum, and its endpoint may depend on the visiting order; multiple random
#' starts are therefore recommended (\code{n_start} in \code{\link{control_EM}}, default 5),
#' and the reported fit is the one with the highest classification log-likelihood among
#' the non-degenerate starts. The classification log-likelihood has degenerate boundary
#' solutions (components with too few events, or with runaway survival parameters);
#' starts ending in such solutions are flagged and excluded from the selection, and the
#' flags are reported in \code{start_summary}. The event floor and the runaway-parameter
#' threshold are controlled by \code{minimum_component_events} and
#' \code{runaway_parameter_limit} in \code{\link{control_EM}}.
#'
#' The package fits the component-specific frailty models with its own optimizer
#' (see \code{\link{control_parfm}}); it does not depend on the \code{parfm} package at
#' runtime. The likelihood routines are adapted from parfm (GPL-2), which is why this
#' package is distributed under GPL-2.
#'
#' @note
#' Character variables in \code{data} are automatically converted to factors. The
#' function handles cases where some factor levels of the cluster-weighted categorical
#' covariates are absent in certain mixture components.
#'
#' @references
#' Besag, J. (1986). On the statistical analysis of dirty pictures.
#' \emph{Journal of the Royal Statistical Society: Series B}, 48(3), 259-302.
#'
#' Caldera, Cappozzo, Masci, Forlani, Antonelli, Leoni, Paganoni, Ieva (2025+).
#' Cluster-weighted modeling of lifetime hierarchical data for profiling COVID-19
#' heart failure patients. \url{https://arxiv.org/abs/2507.12230}
#'
#' Celeux, G. & Govaert, G. (1992). A classification EM algorithm for clustering
#' and two stochastic versions. \emph{Computational Statistics & Data Analysis},
#' 14(3), 315-332.
#'
#' Munda, M., Rotolo, F., Legrand, C. (2012). Parfm: parametric frailty models in R.
#' \emph{Journal of Statistical Software}, 51(11), 1-30.
#'
#' @seealso \code{\link{control_EM}}, \code{\link{control_parfm}},
#'   \code{\link[survival]{Surv}}
#'
#' @export
#'
#' @examples
#' set.seed(2026)
#' # Simulate clustered survival data with two well-separated components
#' dat_list <- list()
#' true_class <- integer()
#' for (g in 1:2) {
#'   for (j in 1:6) {
#'     # gamma frailty shared by the family, mean 1 and variance 0.5
#'     m_j <- stats::rgamma(4, shape = 2, scale = 0.5)
#'     t_ev <- ((-log(runif(4))) /
#'       (m_j * c(0.4, 0.12)[g] * exp(c(0.3, 0.5)[g] * rnorm(4))))^(1 / c(1.3, 1.7)[g])
#'     dat_list[[length(dat_list) + 1]] <- data.frame(
#'       time = pmin(t_ev, 8), status = as.integer(t_ev <= 8),
#'       x = rnorm(4),
#'       Z1 = rnorm(4, c(-1.5, 1.5)[g]),
#'       Z2 = rnorm(4, c(1, -1)[g]),
#'       family = paste0("f", g, j)
#'     )
#'     true_class <- c(true_class, rep(g, 4))
#'   }
#' }
#' dat <- do.call(rbind, dat_list)
#'
#' fit <- fit_mixparfm(
#'   formula = survival::Surv(time, status) ~ x,
#'   G = 2,
#'   class_init = true_class,
#'   grouping_variable = "family",
#'   X_gaussian_variables = c("Z1", "Z2"),
#'   data = dat,
#'   baseline = "weibull",
#'   frailty = "gamma",
#'   control_EM_algorithm = control_EM(itermax = 20, tol = 1e-3, n_start = 3),
#'   control_parfm_algorithm = control_parfm(maxit = 1000)
#' )
#'
#' # Highest classification log-likelihood across the three starts
#' fit$classification_loglik
#' fit$selected_start
#'
#' # Natural-scale parameters, one column per component
#' fit$parameters$AFT_parameters
#'
#' # Classification BIC
#' fit$bic
fit_mixparfm <-
  function(
    formula,
    G,
    class_init,
    grouping_variable = NULL,
    strata = NULL,
    X_gaussian_variables = NULL,
    X_multinomial_variables = NULL,
    data,
    baseline = c(
      "weibull",
      "inweibull",
      "frechet",
      "exponential",
      "gompertz",
      "loglogistic",
      "lognormal"
    ),
    frailty = c("none", "gamma", "ingau"),
    control_EM_algorithm = control_EM(),
    control_parfm_algorithm = control_parfm()
  ) {
    baseline <- match.arg(baseline)
    frailty <- match.arg(frailty)
    if (!is.null(strata)) {
      stop("Stratified models are currently not supported.", call. = FALSE)
    }
    if (is.null(grouping_variable) && !identical(frailty, "none")) {
      stop(
        "A grouping variable is required for a shared-frailty model.",
        call. = FALSE
      )
    }

    # CEM controls
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    tol_zero_var <- control_EM_algorithm$tol_zero_var
    objective_tolerance <- control_EM_algorithm$objective_tolerance
    move_tolerance <- control_EM_algorithm$move_tolerance
    max_c_sweeps <- control_EM_algorithm$max_c_sweeps
    minimum_component_size <- control_EM_algorithm$minimum_component_size
    minimum_component_events <- control_EM_algorithm$minimum_component_events
    runaway_parameter_limit <- control_EM_algorithm$runaway_parameter_limit
    n_start <- control_EM_algorithm$n_start

    # Component fit controls
    inip <- control_parfm_algorithm$inip
    iniFpar <- control_parfm_algorithm$iniFpar
    method <- control_parfm_algorithm$method
    maxit <- control_parfm_algorithm$maxit
    theta_floor <- control_parfm_algorithm$theta_floor

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

    N <- nrow(data)
    if (G < 1L || G != as.integer(G)) {
      stop("G must be a positive integer.", call. = FALSE)
    }
    if (N < G * minimum_component_size) {
      stop(
        "Not enough observations for ",
        G,
        " components of minimum size ",
        minimum_component_size,
        ".",
        call. = FALSE
      )
    }

    # Response: time and status ------------------------------------------
    lhs <- formula[[2]]
    response <- tryCatch(
      eval(lhs, envir = data),
      error = function(e) NULL
    )
    if (is.null(response)) {
      # Surv() not visible from the caller; evaluate with the survival
      # namespace as the enclosing environment
      response <- eval(lhs, envir = data, enclos = asNamespace("survival"))
    }
    if (!inherits(response, "Surv")) {
      stop(
        "The left-hand side of formula must be a Surv object.",
        call. = FALSE
      )
    }
    response_type <- attr(response, "type")
    if (!identical(response_type, "right")) {
      stop(
        "Only right-censored Surv responses are supported; left-truncated ",
        "(counting-process) responses are not.",
        call. = FALSE
      )
    }
    time <- as.numeric(response[, 1L])
    event <- as.numeric(response[, 2L])
    response_vars <- all.vars(lhs)
    if (any(vapply(data[response_vars], anyNA, TRUE))) {
      stop(
        "Missing values in the survival response are not supported.",
        call. = FALSE
      )
    }
    if (length(time) != N || length(event) != N) {
      stop("Survival response has a different length than data.", call. = FALSE)
    }
    if (any(!is.finite(time)) || any(time <= 0)) {
      stop(
        "Survival times must be finite and positive.",
        call. = FALSE
      )
    }
    if (!all(event %in% c(0, 1))) {
      stop("The event indicator must be 0 or 1.", call. = FALSE)
    }
    event <- as.integer(event)

    # Design matrix (built once from the full data so that full-data and
    # component-level evaluations can never disagree)
    rhs_vars <- all.vars(formula[[3]])
    if (any(vapply(data[rhs_vars], anyNA, TRUE))) {
      stop("Missing values in the covariates are not supported.", call. = FALSE)
    }
    full_design <- stats::model.matrix(formula, data)
    if (nrow(full_design) != N) {
      stop(
        "Failed to build the design matrix for all rows of data.",
        call. = FALSE
      )
    }
    design <- full_design[, -1L, drop = FALSE]
    beta_names <- colnames(design)

    # Minimum number of events per component (audit 07, finding N2):
    # degenerate zero-event or near-zero-event components have an
    # unbounded survival likelihood, and their parameter estimates run
    # away to the boundary. The default floor is the number of survival
    # parameters of one component.
    if (is.null(minimum_component_events)) {
      minimum_component_events <- .mixparfm_frailty_npar(frailty) +
        .mixparfm_baseline_npar(baseline) +
        ncol(design)
    }
    if (sum(event) < G * minimum_component_events) {
      stop(
        "Not enough events for ",
        G,
        " components of minimum ",
        minimum_component_events,
        " event(s) each.",
        call. = FALSE
      )
    }

    # Grouping structure
    if (is.null(grouping_variable)) {
      group_index <- rep(1L, N)
      group_levels <- "all"
    } else {
      if (!grouping_variable %in% names(data)) {
        stop("The grouping variable is not present in data.", call. = FALSE)
      }
      group_factor <- factor(data[[grouping_variable]])
      if (anyNA(group_factor)) {
        stop("The grouping variable contains missing values.", call. = FALSE)
      }
      group_index <- as.integer(group_factor)
      group_levels <- levels(group_factor)
    }

    # Cluster-weighted covariates --------------------------------------
    is_X_gaussian <- !is.null(X_gaussian_variables)
    is_X_multinomial <- !is.null(X_multinomial_variables)

    if (is_X_gaussian) {
      X_gaussian <- as.matrix(data[, X_gaussian_variables, drop = FALSE])
      p_gaussian <- ncol(X_gaussian)
      mclust_model_name <- if (p_gaussian == 1) "V" else "VVV"
      validate_gaussian_mstep <- function(gaussian_params, context) {
        sigma <- gaussian_params$parameters$variance$sigma
        sigmasq <- gaussian_params$parameters$variance$sigmasq
        min_eigenvalue <- if (!is.null(sigmasq)) {
          min(as.numeric(sigmasq))
        } else {
          min(vapply(
            seq_len(dim(sigma)[3L]),
            function(g) {
              min(
                eigen(sigma[,, g], symmetric = TRUE, only.values = TRUE)$values
              )
            },
            0
          ))
        }
        if (!is.finite(min_eigenvalue) || min_eigenvalue <= tol_zero_var) {
          stop(
            "Degenerate Gaussian covariance during ",
            context,
            " (minimum eigenvalue: ",
            if (is.finite(min_eigenvalue)) {
              format(min_eigenvalue)
            } else {
              "non-finite"
            },
            "; tolerance: ",
            tol_zero_var,
            ")",
            call. = FALSE
          )
        }
      }
    } else {
      X_gaussian <- NULL
      p_gaussian <- 0
      n_par_X_gaussian <- 0
    }

    multinomial_level_index <- list()
    multinomial_templates <- list()
    n_par_X_multinomial <- 0
    if (is_X_multinomial) {
      X_multinomial <- data[, X_multinomial_variables, drop = FALSE]
      p_multinomial <- ncol(X_multinomial)
      for (m in seq_len(p_multinomial)) {
        level_names <- levels(as.factor(X_multinomial[, m]))
        multinomial_templates[[m]] <- matrix(
          0,
          nrow = length(level_names),
          ncol = G,
          dimnames = list(level_names, NULL)
        )
        multinomial_level_index[[m]] <- match(
          as.character(X_multinomial[, m]),
          level_names
        )
        if (anyNA(multinomial_level_index[[m]])) {
          stop(
            "A categorical covariate contains values without a probability level.",
            call. = FALSE
          )
        }
      }
    } else {
      X_multinomial <- NULL
      p_multinomial <- 0
    }

    # Initial partition (validated before any fitting) --------------------
    class_init <- as.integer(class_init)
    if (
      length(class_init) != N ||
        anyNA(class_init) ||
        any(class_init < 1L | class_init > G)
    ) {
      stop(
        "class_init must be an integer vector of length nrow(data) with ",
        "values between 1 and G.",
        call. = FALSE
      )
    }
    init_counts <- tabulate(class_init, nbins = G)
    if (any(init_counts < minimum_component_size)) {
      stop(
        "class_init leaves component(s) ",
        paste(which(init_counts < minimum_component_size), collapse = ", "),
        " with fewer than ",
        minimum_component_size,
        " observation(s).",
        call. = FALSE
      )
    }
    init_event_counts <- tabulate(class_init[event == 1L], nbins = G)
    if (any(init_event_counts < minimum_component_events)) {
      stop(
        "class_init leaves component(s) ",
        paste(
          which(init_event_counts < minimum_component_events),
          collapse = ", "
        ),
        " with fewer than ",
        minimum_component_events,
        " event(s); components with too few events give degenerate ",
        "survival-parameter estimates. Change or disable the constraint ",
        "via control_EM(minimum_component_events).",
        call. = FALSE
      )
    }

    # Start partitions: the user-supplied one, plus random ones
    starts <- vector("list", n_start)
    starts[[1L]] <- class_init
    if (n_start > 1L) {
      for (s in seq_len(n_start - 1L)) {
        repeat {
          candidate <- sample.int(G, N, replace = TRUE)
          if (
            min(tabulate(candidate, nbins = G)) >= minimum_component_size &&
              min(tabulate(candidate[event == 1L], nbins = G)) >=
                minimum_component_events
          ) {
            break
          }
        }
        starts[[s + 1L]] <- candidate
      }
    }

    # Component survival log-likelihood at a subset of rows; the same
    # expression the optimizer inside .mixparfm_fit_component maximizes.
    survival_loglik_at <- function(p, index) {
      terms <- .mixparfm_component_terms(
        p,
        design[index, , drop = FALSE],
        time[index],
        event[index],
        baseline,
        frailty
      )
      .mixparfm_component_survival_loglik(
        terms,
        group_index[index],
        frailty
      )
    }

    # State: everything the ICM classification step and the
    # objective need at the current parameter values.
    build_state <- function(
      class,
      params,
      tau,
      gaussian_params,
      multinom_params
    ) {
      terms <- lapply(seq_len(G), function(g) {
        .mixparfm_component_terms(
          params[[g]],
          design,
          time,
          event,
          baseline,
          frailty
        )
      })
      event_log_hazard <- vapply(
        terms,
        function(component) component$event_log_hazard,
        numeric(N)
      )

      covariate_logdens <- NULL
      if (is_X_gaussian) {
        covariate_logdens <- as.matrix(do.call(
          mclust::cdens,
          c(list(data = X_gaussian, logarithm = TRUE), gaussian_params)
        ))
      }
      if (is_X_multinomial) {
        multinomial_logdens <- matrix(0, N, G)
        for (m in seq_len(p_multinomial)) {
          probs <- multinom_params[[m]]
          levels_m <- multinomial_level_index[[m]]
          for (g in seq_len(G)) {
            multinomial_logdens[, g] <- multinomial_logdens[, g] +
              log(probs[levels_m, g])
          }
        }
        covariate_logdens <- if (is.null(covariate_logdens)) {
          multinomial_logdens
        } else {
          covariate_logdens + multinomial_logdens
        }
      }

      list(
        G = G,
        class = as.integer(class),
        frailty = frailty,
        group_index = group_index,
        terms = terms,
        additive = .mixparfm_additive_log_matrix(
          tau,
          event_log_hazard,
          covariate_logdens
        )
      )
    }

    # One ICM-CEM run from a given initial partition.
    run_cem <- function(init_class) {
      class <- as.integer(init_class)

      # Initialization: component-wise fits on the initial partition,
      # with deterministic (or user-supplied) starting values.
      params <- vector("list", G)
      fit_list <- vector("list", G)
      for (g in seq_len(G)) {
        index <- class == g
        start_g <- .mixparfm_default_init(
          time[index],
          event[index],
          design[index, , drop = FALSE],
          baseline,
          frailty,
          inip,
          iniFpar
        )
        fit_list[[g]] <- .mixparfm_fit_component(
          design = design[index, , drop = FALSE],
          time = time[index],
          event = event[index],
          group_index = group_index[index],
          baseline = baseline,
          frailty = frailty,
          start = start_g,
          method = method,
          maxit = maxit,
          theta_floor = theta_floor
        )
        if (!fit_list[[g]]$ok) {
          stop(
            "Initial fit of component ",
            g,
            " failed: ",
            fit_list[[g]]$message,
            call. = FALSE
          )
        }
        params[[g]] <- fit_list[[g]]$estim_par
      }

      tau <- as.numeric(tabulate(class, nbins = G)) / N

      multinomial_params <- vector("list", p_multinomial)
      if (is_X_multinomial) {
        for (m in seq_len(p_multinomial)) {
          probs <- multinomial_templates[[m]]
          for (g in seq_len(G)) {
            tab <- table(factor(
              X_multinomial[class == g, m],
              levels = rownames(probs)
            ))
            probs[, g] <- as.numeric(tab) / sum(tab)
          }
          multinomial_params[[m]] <- probs
        }
      }

      gaussian_params <- NULL
      if (is_X_gaussian) {
        gaussian_params <- mclust::mstep(
          data = X_gaussian,
          z = mclust::unmap(class),
          modelName = mclust_model_name
        )
        validate_gaussian_mstep(gaussian_params, "initialization")
      }

      # Main loop -------------------------------------------------------
      iter <- 0L
      trace_records <- vector("list", itermax)
      n_m_step_rejections <- 0L
      err <- Inf

      repeat {
        state <- build_state(
          class,
          params,
          tau,
          gaussian_params,
          multinomial_params
        )
        objective_before <- .mixparfm_classification_loglik(state)

        # ICM classification step
        c_result <- .mixparfm_conditional_cstep(
          state,
          max_sweeps = max_c_sweeps,
          move_tolerance = move_tolerance,
          minimum_component_size = minimum_component_size,
          minimum_component_events = minimum_component_events,
          objective_tolerance = objective_tolerance
        )
        class <- c_result$class
        objective_after_c <- c_result$objective_after

        # M step -----------------------------------------------------
        tau <- as.numeric(tabulate(class, nbins = G)) / N
        if (is_X_multinomial) {
          for (m in seq_len(p_multinomial)) {
            for (g in seq_len(G)) {
              tab <- table(factor(
                X_multinomial[class == g, m],
                levels = rownames(multinomial_params[[m]])
              ))
              multinomial_params[[m]][, g] <- as.numeric(tab) / sum(tab)
            }
          }
        }
        if (is_X_gaussian) {
          gaussian_params <- mclust::mstep(
            data = X_gaussian,
            z = mclust::unmap(class),
            modelName = mclust_model_name
          )
          validate_gaussian_mstep(
            gaussian_params,
            paste("iteration", iter + 1L)
          )
        }

        for (g in seq_len(G)) {
          index <- class == g
          old_value <- survival_loglik_at(params[[g]], index)
          candidate <- .mixparfm_fit_component(
            design = design[index, , drop = FALSE],
            time = time[index],
            event = event[index],
            group_index = group_index[index],
            baseline = baseline,
            frailty = frailty,
            start = params[[g]], # warm start from the current estimates
            method = method,
            maxit = maxit,
            theta_floor = theta_floor
          )
          # strict non-decrease: a candidate that does not improve the
          # component's survival log-likelihood is rejected outright
          if (candidate$ok && candidate$loglik >= old_value) {
            params[[g]] <- candidate$estim_par
            fit_list[[g]] <- candidate
          } else {
            # keep the previous estimates: a partial M-step is still
            # non-decreasing for this component's contribution
            n_m_step_rejections <- n_m_step_rejections + 1L
          }
        }

        state <- build_state(
          class,
          params,
          tau,
          gaussian_params,
          multinomial_params
        )
        objective_after_m <- .mixparfm_classification_loglik(state)

        total_delta <- objective_after_m - objective_before
        err <- abs(total_delta) / (1 + abs(objective_before))

        iter <- iter + 1L
        trace_records[[iter]] <- data.frame(
          iteration = iter,
          objective_before = objective_before,
          objective_after_c = objective_after_c,
          objective_after_m = objective_after_m,
          c_step_delta = objective_after_c - objective_before,
          m_step_delta = objective_after_m - objective_after_c,
          total_delta = total_delta,
          c_moves = c_result$total_moves,
          c_sweeps = nrow(c_result$sweeps),
          c_step_converged = c_result$c_step_converged,
          m_step_rejections = n_m_step_rejections,
          relative_change = total_delta / (1 + abs(objective_before))
        )

        if (err <= tol || iter >= itermax) break
      }

      cem_trace <- do.call(
        rbind,
        trace_records[!vapply(trace_records, is.null, logical(1))]
      )
      converged <- is.finite(err) && err <= tol

      list(
        class = class,
        params = params,
        tau = tau,
        gaussian_params = gaussian_params,
        multinomial_params = multinomial_params,
        fit_list = fit_list,
        classification_loglik = tail(cem_trace$objective_after_m, 1L),
        cem_trace = cem_trace,
        n_iter = iter,
        converged = converged,
        final_relative_change = err,
        n_m_step_rejections = n_m_step_rejections
      )
    }

    # Degeneracy flags (audit 07, finding N2): a start whose final
    # partition leaves a component below the event floor, or whose
    # component survival parameters have run away (huge log frailty
    # variance, baseline or regression coefficients on the optimizer
    # scale), is flagged and excluded from the multi-start selection.
    flag_degenerate <- function(class, params) {
      reasons <- character(0)
      event_counts <- tabulate(class[event == 1L], nbins = G)
      short <- which(event_counts < minimum_component_events)
      if (length(short) > 0L) {
        reasons <- c(
          reasons,
          sprintf(
            "component(s) %s have fewer than %d event(s)",
            paste(short, collapse = ", "),
            minimum_component_events
          )
        )
      }
      if (is.finite(runaway_parameter_limit)) {
        for (g in seq_len(G)) {
          if (any(abs(params[[g]]) > runaway_parameter_limit)) {
            reasons <- c(
              reasons,
              sprintf(
                "component %d has runaway survival parameters",
                g
              )
            )
          }
        }
      }
      reasons
    }

    # Run the starts; failures are recorded, not thrown -------------------
    start_results <- lapply(starts, function(cl) {
      tryCatch(
        {
          r <- run_cem(cl)
          r$degeneracy_reasons <- flag_degenerate(r$class, r$params)
          r
        },
        error = function(e) e
      )
    })
    start_ok <- !vapply(start_results, inherits, TRUE, "error")
    if (!any(start_ok)) {
      stop(
        "All ",
        n_start,
        " start(s) failed. First error: ",
        conditionMessage(start_results[[1L]]),
        call. = FALSE
      )
    }
    # Degenerate starts (audit 07, finding N2) are excluded from the
    # selection: their objective values are spurious boundary solutions
    # of the classification log-likelihood.
    degenerate <- vapply(
      start_results,
      function(r) {
        !inherits(r, "error") && length(r$degeneracy_reasons) > 0L
      },
      TRUE
    )
    eligible <- which(start_ok & !degenerate)
    if (!length(eligible)) {
      stop(
        "All non-failed start(s) ended in degenerate solutions (runaway ",
        "survival parameters or components below the minimum event count). ",
        "Consider a larger n_start, a smaller G, or relax the checks via ",
        "control_EM(minimum_component_events, runaway_parameter_limit).",
        call. = FALSE
      )
    }
    logliks <- vapply(
      start_results[eligible],
      function(r) r$classification_loglik,
      0
    )
    selected <- eligible[which.max(logliks)]
    best <- start_results[[selected]]
    n_degenerate_starts <- sum(degenerate)
    if (n_degenerate_starts > 0L) {
      warning(
        n_degenerate_starts,
        " start(s) ended in degenerate solutions and were excluded from ",
        "the selection (see start_summary).",
        call. = FALSE
      )
    }

    class <- best$class
    params <- best$params
    tau <- best$tau
    loglik <- best$classification_loglik

    start_summary <- data.frame(
      start = seq_len(n_start),
      classification_loglik = vapply(
        start_results,
        function(r) {
          if (inherits(r, "error")) NA_real_ else r$classification_loglik
        },
        0
      ),
      n_iter = vapply(
        start_results,
        function(r) {
          if (inherits(r, "error")) NA_integer_ else r$n_iter
        },
        1L
      ),
      converged = vapply(
        start_results,
        function(r) {
          if (inherits(r, "error")) NA else r$converged
        },
        TRUE
      ),
      error = vapply(
        start_results,
        function(r) {
          if (inherits(r, "error")) conditionMessage(r) else ""
        },
        ""
      ),
      degenerate = vapply(
        start_results,
        function(r) {
          if (inherits(r, "error")) NA else length(r$degeneracy_reasons) > 0L
        },
        TRUE
      ),
      degeneracy_reason = vapply(
        start_results,
        function(r) {
          if (inherits(r, "error")) {
            ""
          } else {
            paste(r$degeneracy_reasons, collapse = "; ")
          }
        },
        ""
      )
    )

    # Natural-scale parameter table --------------------------------------
    transformed <- lapply(seq_len(G), function(g) {
      .mixparfm_transform_parameters(
        params[[g]],
        baseline,
        frailty,
        beta_names
      )
    })
    AFT_parameters <- t(do.call(rbind, transformed))
    colnames(AFT_parameters) <- seq_len(G)

    parameters <- list(
      tau = tau,
      AFT_parameters = AFT_parameters
    )
    if (is_X_gaussian) {
      parameters$X_gaussian_parameters <- list(
        mu = best$gaussian_params$parameters$mean,
        sigma = best$gaussian_params$parameters$variance$sigma
      )
    }
    if (is_X_multinomial) {
      parameters$X_multinomial_parameters <- best$multinomial_params
    }

    # Empirical-Bayes frailty effects -----------------------------------
    if (frailty == "none") {
      frailty_effect <- NULL
      frailty_var_effect <- NULL
    } else {
      terms_best <- lapply(seq_len(G), function(g) {
        .mixparfm_component_terms(
          params[[g]],
          design,
          time,
          event,
          baseline,
          frailty
        )
      })
      effect_frames <- vector("list", G)
      var_frames <- vector("list", G)
      for (g in seq_len(G)) {
        index <- class == g
        d_j <- tapply(event[index], group_index[index], sum)
        s_j <- tapply(terms_best[[g]]$risk[index], group_index[index], sum)
        j_index <- as.integer(names(d_j))
        posterior <- lapply(
          seq_along(d_j),
          function(k) {
            .mixparfm_frailty_posterior(
              frailty,
              d_j[k],
              s_j[k],
              terms_best[[g]]$theta
            )
          }
        )
        means <- vapply(posterior, function(p) p$mean, 0)
        variances <- vapply(posterior, function(p) p$var, 0)
        effect_frames[[g]] <- data.frame(
          group = group_levels[j_index],
          setNames(list(means), as.character(g)),
          check.names = FALSE
        )
        var_frames[[g]] <- data.frame(
          group = group_levels[j_index],
          setNames(list(variances), as.character(g)),
          check.names = FALSE
        )
      }
      join_frames <- function(frames) {
        Reduce(
          function(a, b) merge(a, b, by = "group", all = TRUE),
          frames
        )
      }
      frailty_effect <- join_frames(effect_frames)
      frailty_var_effect <- join_frames(var_frames)
      colnames(frailty_effect)[1] <- grouping_variable
      colnames(frailty_var_effect)[1] <- grouping_variable
    }

    # Frailty-variance floor warning -------------------------------------
    theta_floored <- vapply(
      best$fit_list,
      function(f) isTRUE(f$theta_floored),
      TRUE
    )
    if (frailty != "none" && any(theta_floored)) {
      warning(
        "The frailty variance estimate hit the floor theta_floor = ",
        theta_floor,
        " in component(s) ",
        paste(which(theta_floored), collapse = ", "),
        "; the corresponding frailty is effectively switched off. ",
        "Consider a larger floor or frailty = \"none\".",
        call. = FALSE
      )
    }

    # Classification BIC -------------------------------------------------
    n_par_tau <- G - 1
    n_par_parfm <- length(AFT_parameters)
    if (is_X_gaussian) {
      n_par_X_gaussian <- mclust::nMclustParams(
        modelName = mclust_model_name,
        d = p_gaussian,
        G = G
      ) -
        n_par_tau
    }
    if (is_X_multinomial) {
      n_par_X_multinomial <- sum(
        vapply(best$multinomial_params, length, 0L) - G
      )
    }
    bic_final <- 2 *
      loglik -
      (n_par_tau + n_par_parfm + n_par_X_gaussian + n_par_X_multinomial) *
        log(N)

    cem_trace <- best$cem_trace
    loglik_vec <- cem_trace$objective_after_m
    reached_itermax <- !best$converged && best$n_iter >= itermax
    stopping_reason <- if (best$converged) {
      "tolerance_reached"
    } else if (reached_itermax) {
      "iteration_limit"
    } else {
      "stopped_without_tolerance"
    }

    list(
      loglik = loglik,
      classification_loglik = loglik,
      parameters = parameters,
      frailty_effect = frailty_effect,
      frailty_var_effect = frailty_var_effect,
      z = mclust::unmap(class),
      class = as.integer(class),
      bic = bic_final,
      classification_bic = bic_final,
      baseline = baseline,
      frailty = frailty,
      fit_parfm = lapply(best$fit_list, function(f) {
        list(
          estim_par = f$estim_par,
          loglik = f$loglik,
          convergence = f$convergence,
          message = f$message,
          theta_floored = f$theta_floored
        )
      }),
      loglik_vec = loglik_vec,
      classification_loglik_vec = loglik_vec,
      cem_trace = cem_trace,
      objective_name = "classification_loglik",
      algorithm = "icm_cem",
      converged = best$converged,
      stopping_reason = stopping_reason,
      n_iter = best$n_iter,
      reached_itermax = reached_itermax,
      final_relative_change = best$final_relative_change,
      tol_zero_var = tol_zero_var,
      n_start = n_start,
      selected_start = selected,
      start_summary = start_summary,
      n_degenerate_starts = n_degenerate_starts,
      n_m_step_rejections = best$n_m_step_rejections,
      theta_floor = theta_floor
    )
  }
