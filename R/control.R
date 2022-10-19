#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2,
                       tol_zero_var = sqrt(.Machine$double.eps),
                       E_step_update = "classification") { # E_step_update can be either classification (CEM) or stochastic (SEM)
  list(
    itermax = itermax,
    tol = tol,
    err = err,
    E_step_update = E_step_update
  )
}

#' @export
control_parfm <- function(inip = NULL,
                          iniFpar = NULL,
                          method = "nlminb",
                          maxit = 500,
                          Fparscale = 1,
                          showtime = FALSE,
                          correct = 0){
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
