#' Control Parameters for the EM Algorithm
#'
#' @description
#' Auxiliary function for controlling the EM algorithm used in fitting finite
#' mixtures of parametric frailty models. This function specifies convergence
#' criteria and the type of E-step update.
#'
#' @param itermax Integer. Maximum number of EM iterations. Default is 1000.
#' @param tol Numeric. Convergence tolerance for the relative change in
#'   log-likelihood. The algorithm stops when the relative change is smaller
#'   than \code{tol}. Default is 1e-08.
#' @param err Numeric. Initial error value for the EM algorithm. Default is
#'   \code{.Machine$double.xmax / 2}.
#' @param tol_zero_var Numeric. Tolerance for detecting zero variance components.
#'   Components with variance below this threshold may be considered degenerate.
#'   Default is \code{sqrt(.Machine$double.eps)}.
#' @param E_step_update Character. Type of E-step update to use:
#'   \describe{
#'     \item{\code{"classification"}}{Classification EM (CEM) - uses hard assignments
#'       based on maximum a posteriori probability (default)}
#'     \item{\code{"stochastic"}}{Stochastic EM (SEM) - uses stochastic assignments
#'       by sampling from the posterior probabilities}
#'   }
#'
#' @return A list with components corresponding to the control parameters.
#'
#' @seealso \code{\link{fit_mixparfm}}, \code{\link{control_parfm}}
#'
#' @export
#'
#' @examples
#' # Default control parameters
#' control_EM()
#'
#' # Custom parameters for faster convergence testing
#' control_EM(itermax = 100, tol = 1e-05)
#'
#' # Using stochastic EM
#' control_EM(E_step_update = "stochastic")
control_EM <- function(
  itermax = 1000,
  tol = 1e-08,
  err = .Machine$double.xmax / 2,
  tol_zero_var = sqrt(.Machine$double.eps),
  E_step_update = "classification"
) {
  # E_step_update can be either classification (CEM) or stochastic (SEM)
  list(
    itermax = itermax,
    tol = tol,
    err = err,
    E_step_update = E_step_update
  )
}

#' Control Parameters for Parametric Frailty Model Fitting
#'
#' @description
#' Auxiliary function for controlling the optimization of parametric frailty models
#' within each mixture component. These parameters are passed to the \code{parfm}
#' package functions for fitting individual cluster-specific models.
#'
#' @param inip Numeric vector or NULL. Initial values for the baseline hazard
#'   parameters. If NULL (default), initial values are computed automatically
#'   by the \code{parfm} package.
#' @param iniFpar Numeric or NULL. Initial value for the frailty parameter.
#'   If NULL (default), an initial value is computed automatically.
#' @param method Character. Optimization method to use:
#'   \describe{
#'     \item{\code{"nlminb"}}{Port's \code{nlminb} optimizer (default)}
#'     \item{\code{"Nelder-Mead"}}{Nelder-Mead simplex method}
#'     \item{\code{"BFGS"}}{Quasi-Newton BFGS algorithm}
#'   }
#' @param maxit Integer. Maximum number of iterations for the optimization
#'   algorithm. Default is 500.
#' @param Fparscale Numeric. Scaling parameter for the frailty parameter in
#'   the optimization routine. Default is 1.
#' @param showtime Logical. If TRUE, displays computation time for fitting
#'   each mixture component. Useful for diagnosing slow convergence. Default is FALSE.
#' @param correct Numeric. Correction factor for the optimization. Default is 0.
#'
#' @return A list with components corresponding to the control parameters.
#'
#' @seealso \code{\link{fit_mixparfm}}, \code{\link{control_EM}}, \code{\link[parfm]{parfm}}
#'
#' @export
#'
#' @examples
#' # Default control parameters
#' control_parfm()
#'
#' # Custom optimization settings with more iterations
#' control_parfm(method = "BFGS", maxit = 1000)
#'
#' # With custom initial values and timing
#' control_parfm(inip = c(0.5, 1.2), showtime = TRUE)
#'
#' # Using Nelder-Mead with increased iterations
#' control_parfm(method = "Nelder-Mead", maxit = 2000)
control_parfm <- function(
  inip = NULL,
  iniFpar = NULL,
  method = "nlminb",
  maxit = 500,
  Fparscale = 1,
  showtime = FALSE,
  correct = 0
) {
  list(
    inip = inip,
    iniFpar = iniFpar,
    method = method,
    maxit = maxit,
    Fparscale = Fparscale,
    showtime = showtime,
    correct = correct
  )
}
