#' Control Parameters for the ICM-CEM Algorithm
#'
#' @description
#' Auxiliary function for controlling the ICM-CEM algorithm (a
#' classification EM algorithm whose classification step is carried out by
#' iterated conditional modes) used in fitting finite mixtures of parametric
#' frailty models. This function specifies convergence criteria,
#' classification-step tolerances, and the number of random starts.
#'
#' @param itermax Integer. Maximum number of CEM iterations per start.
#'   Default is 1000.
#' @param tol Numeric. Convergence tolerance for the relative change in
#'   the classification log-likelihood. The algorithm stops when the
#'   relative change is smaller than \code{tol}. Default is 1e-08.
#' @param tol_zero_var Numeric. Tolerance for detecting degenerate
#'   covariance matrices of the Gaussian cluster-weighted covariates.
#'   Components whose covariance has an eigenvalue below this threshold
#'   produce an error for the corresponding start. Default is
#'   \code{sqrt(.Machine$double.eps)}.
#' @param objective_tolerance Numeric. Absolute tolerance used by the
#'   internal monotonicity checks of the classification log-likelihood.
#'   Default is 1e-7.
#' @param move_tolerance Numeric. Relative tolerance of the strict-
#'   improvement rule of the ICM classification step: a
#'   single-patient reassignment is accepted only if it improves the
#'   score by more than \code{move_tolerance * (1 + |score|)}. This
#'   keeps rounding-level "improvements" of the incrementally maintained
#'   cell sums from causing spurious moves or cycling. Default is 1e-9.
#' @param max_c_sweeps Integer. Maximum number of sequential
#'   classification sweeps per CEM iteration. Default is 100.
#' @param minimum_component_size Integer. Minimum number of patients per
#'   component enforced during the classification sweeps; the conditional
#'   step never moves a patient if that would leave the current component
#'   with fewer than this many patients. Default is 1, i.e. empty
#'   components are prevented.
#' @param minimum_component_events Integer or NULL. Minimum number of
#'   events per component, enforced like \code{minimum_component_size}:
#'   the conditional step never moves a patient if that would leave the
#'   current component with fewer than this many events. This prevents
#'   zero-event (or near-zero-event) components, whose survival-parameter
#'   estimates run away to degenerate boundary solutions. NULL (default)
#'   means the number of survival parameters of one component (frailty
#'   variance, baseline and regression coefficients); 0 disables the
#'   constraint. The initial partition (\code{class_init}) must satisfy
#'   the constraint.
#' @param runaway_parameter_limit Numeric. A component whose parameter
#'   vector on the optimizer scale (log frailty variance, transformed
#'   baseline parameters, regression coefficients) has any entry larger
#'   than this limit in absolute value is flagged as degenerate. Flagged
#'   starts are excluded from the multi-start selection and reported in
#'   \code{start_summary}. Default is 30; \code{Inf} disables the check.
#' @param n_start Integer. Number of starts. The first start uses the
#'   partition supplied via the \code{class_init} argument of
#'   \code{\link{fit_mixparfm}}; the remaining ones are random partitions
#'   respecting \code{minimum_component_size} and
#'   \code{minimum_component_events}. The start attaining the
#'   highest classification log-likelihood among the non-degenerate ones
#'   is returned. Set the seed of the random number generator before
#'   calling for reproducibility. Default is 5.
#' @param E_step_update Character. Only \code{"classification"} is
#'   accepted; the argument is retained for backward compatibility, and
#'   the stochastic (SEM) update of earlier versions is no longer
#'   available.
#'
#' @return A list containing the supplied algorithm controls.
#'
#' @seealso \code{\link{fit_mixparfm}}, \code{\link{control_parfm}}
#'
#' @export
#'
#' @examples
#' # Default control parameters
#' control_EM()
#'
#' # Faster convergence testing
#' control_EM(itermax = 100, tol = 1e-05)
#'
#' # Multiple random starts with best-objective selection
#' control_EM(itermax = 100, tol = 1e-05, n_start = 5)
control_EM <- function(
  itermax = 1000,
  tol = 1e-08,
  tol_zero_var = sqrt(.Machine$double.eps),
  objective_tolerance = 1e-7,
  move_tolerance = 1e-9,
  max_c_sweeps = 100L,
  minimum_component_size = 1L,
  minimum_component_events = NULL,
  runaway_parameter_limit = 30,
  n_start = 5L,
  E_step_update = "classification"
) {
  numeric_controls <- c(
    tol = tol,
    tol_zero_var = tol_zero_var,
    objective_tolerance = objective_tolerance,
    move_tolerance = move_tolerance
  )
  if (any(!is.finite(numeric_controls)) || any(numeric_controls < 0)) {
    stop("EM tolerances must be finite and non-negative.", call. = FALSE)
  }
  integer_controls <- c(
    itermax = itermax,
    max_c_sweeps = max_c_sweeps,
    minimum_component_size = minimum_component_size,
    n_start = n_start
  )
  if (
    any(!is.finite(integer_controls)) ||
      any(integer_controls < 1L) ||
      any(integer_controls != as.integer(integer_controls))
  ) {
    stop(
      "itermax, max_c_sweeps, minimum_component_size and n_start must be ",
      "positive integers.",
      call. = FALSE
    )
  }
  if (!is.null(minimum_component_events)) {
    if (
      length(minimum_component_events) != 1L ||
        !is.finite(minimum_component_events) ||
        minimum_component_events < 0 ||
        minimum_component_events != as.integer(minimum_component_events)
    ) {
      stop(
        "minimum_component_events must be a non-negative integer or NULL.",
        call. = FALSE
      )
    }
  }
  if (
    length(runaway_parameter_limit) != 1L ||
      is.na(runaway_parameter_limit) ||
      runaway_parameter_limit <= 0
  ) {
    stop(
      "runaway_parameter_limit must be a positive number (or Inf).",
      call. = FALSE
    )
  }
  if (
    length(E_step_update) != 1L ||
      !identical(E_step_update, "classification")
  ) {
    stop(
      "Only E_step_update = \"classification\" is available; the ",
      "stochastic (SEM) update is no longer implemented.",
      call. = FALSE
    )
  }

  list(
    itermax = as.integer(itermax),
    tol = tol,
    tol_zero_var = tol_zero_var,
    objective_tolerance = objective_tolerance,
    move_tolerance = move_tolerance,
    max_c_sweeps = as.integer(max_c_sweeps),
    minimum_component_size = as.integer(minimum_component_size),
    minimum_component_events = if (is.null(minimum_component_events)) {
      NULL
    } else {
      as.integer(minimum_component_events)
    },
    runaway_parameter_limit = runaway_parameter_limit,
    n_start = as.integer(n_start),
    E_step_update = E_step_update
  )
}

#' Control Parameters for the Component-Wise Frailty Model Fits
#'
#' @description
#' Auxiliary function for controlling the optimization of the
#' component-specific parametric frailty models fitted in the M-step of
#' the ICM-CEM algorithm.
#'
#' @param inip Numeric vector or NULL. Initial values on the optimizer
#'   scale for the baseline hazard parameters and the regression
#'   coefficients (the frailty parameter is supplied separately through
#'   \code{iniFpar}). If NULL (default), deterministic data-driven
#'   starting values are computed for each component.
#' @param iniFpar Numeric or NULL. Initial value for the frailty
#'   variance (natural scale). If NULL (default), 1 is used. Ignored when
#'   \code{frailty = "none"}.
#' @param method Character. Optimization method to use:
#'   \describe{
#'     \item{\code{"nlminb"}}{Port's \code{nlminb} optimizer with box
#'       constraints (default)}
#'     \item{\code{"Nelder-Mead"}}{Nelder-Mead simplex method}
#'     \item{\code{"BFGS"}}{Quasi-Newton BFGS algorithm}
#'   }
#' @param maxit Integer. Maximum number of iterations for the optimization
#'   algorithm. Default is 500.
#' @param theta_floor Numeric. Lower bound for the estimated frailty
#'   variance of each component. The bound protects the optimizer from
#'   the degenerate corner where the frailty variance collapses to zero;
#'   a warning is emitted when the bound is active in the final fit.
#'   Ignored when \code{frailty = "none"}. Default is 1e-4.
#' @param Fparscale,showtime,correct Legacy arguments of earlier versions
#'   of the package (which relied on \code{parfm}); they are accepted for
#'   backward compatibility and ignored.
#'
#' @return A list containing the supplied optimizer controls.
#'
#' @seealso \code{\link{fit_mixparfm}}, \code{\link{control_EM}}
#'
#' @export
#'
#' @examples
#' # Default control parameters
#' control_parfm()
#'
#' # Custom optimization settings with more iterations
#' control_parfm(maxit = 1000)
#'
#' # Tighter floor on the frailty variance
#' control_parfm(theta_floor = 1e-3)
control_parfm <- function(
  inip = NULL,
  iniFpar = NULL,
  method = c("nlminb", "Nelder-Mead", "BFGS"),
  maxit = 500,
  theta_floor = 1e-4,
  Fparscale = NULL,
  showtime = NULL,
  correct = NULL
) {
  method <- match.arg(method)
  if (
    length(maxit) != 1L ||
      !is.finite(maxit) ||
      maxit < 1L ||
      maxit != as.integer(maxit)
  ) {
    stop("maxit must be a positive integer.", call. = FALSE)
  }
  if (
    length(theta_floor) != 1L || !is.finite(theta_floor) || theta_floor <= 0
  ) {
    stop("theta_floor must be finite and positive.", call. = FALSE)
  }

  list(
    inip = inip,
    iniFpar = iniFpar,
    method = method,
    maxit = as.integer(maxit),
    theta_floor = theta_floor,
    Fparscale = Fparscale,
    showtime = showtime,
    correct = correct
  )
}
