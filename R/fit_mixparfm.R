#' @export
fit_mixparfm <-
  function(formula,
           G,
           class_init,
           cluster = NULL,
           strata = NULL,
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
           frailty   = c("none", "gamma", "ingau", "possta",
                         "lognormal", "loglogistic"),
           control_EM_algorithm = control_EM(),
           control_parfm=control_parfm()) {

    if (missing(data)) {
      data <- eval(parse(text = paste(
        "data.frame(",
        paste(all.vars(formula), collapse = ", "),
        ")"
      )))
    }

    # EM controls
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err

    # parfm controls
    inip = control_parfm$inip
    iniFpar = control_parfm$iniFpar
    method = control_parfm$method
    maxit = control_parfm$maxit
    Fparscale = control_parfm$Fparscale
    showtime = control_parfm$showtime
    correct = control_parfm$correct

    # While loop controls
    iter <- 0
    loglik <- loglik_prev <- -.Machine$integer.max / 2
    loglik_vec <- NULL
    crit <- TRUE

    # Parameters initialization according to class_init

    tau <- colMeans(mclust::unmap(class_init))
    params <- vector("list", length = G)

    # obsdata object needed for log likelihood computation from the parfm package:
    # it will automatically account for different
    # dist and frailties
    obsdata <-
      parfm::obsdata_creator(
        formula = formula,
        data = data,
        cluster = cluster,
        strata = strata,
        frailty = frailty,
        dist = baseline
      )

    for (g in 1:G) {
      data_g <- data[class_init == g,]

      fit_parfm_g <-
        parfm::parfm(
          formula = formula,
          cluster = cluster,
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
      params[[g]] <- attributes(fit_parfm_g)$estim_par
    }

    while (crit) {

      # E step

      time_log_density <-
        sapply(1:G, function(g)
          parfm::loglikelihood_i(
            p = params[[g]],
            obs = obsdata,
            dist = baseline,
            frailty = frailty,
            correct = correct,
            transform = TRUE
          ))

      log_comp_mixture <-
        sweep(
          (time_log_density),
          MARGIN = 2,
          STATS = log(tau),
          FUN = "+"
        )

      z_max <- apply(log_comp_mixture, 1, max)
      log_density <-
        z_max + log(rowSums(exp(log_comp_mixture - z_max)))
      log_z <- log_comp_mixture - log_density
      z <- exp(log_z)

      # Classification step (alternatively, an S step can be performed here)
      map_z <- mclust::map(z)
      z <- mclust::unmap(map_z)

      # M step
      tau <- colMeans(z)

      for (g in 1:G) {
        data_g <- data[map_z == g,]

        fit_parfm_g <-
          parfm::parfm(
            formula = formula,
            cluster = cluster,
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
        params[[g]] <- attributes(fit_parfm_g)$estim_par
      }

      # Log likelihood

      for (g in 1:G) {
        log_density_g <-
          parfm::loglikelihood_i(
            p = params[[g]],
            obs = obsdata,
            dist = baseline,
            frailty = frailty,
            correct = correct
          )
        time_log_density[, g] <-  log_density_g
      }


      log_comp_mixture <-
        sweep(
          (time_log_density),
          MARGIN = 2,
          STATS = log(tau),
          FUN = "+"
        )

      z_max <- apply(log_comp_mixture, 1, max)
      loglik <- sum(z_max + log(rowSums(exp(
        log_comp_mixture - z_max
      ))))

      # Check convergence
      err <- abs(loglik - loglik_prev) / (1 + abs(loglik))
      loglik_prev <- loglik
      loglik_vec <- c(loglik_vec, loglik)
      iter <- iter + 1

      crit <- err > tol & iter < itermax
    }
    params_transformed <-
      lapply(1:G, function(g)
        parfm::paramaters_transformator(
          estim_par = params[[g]],
          frailty = frailty,
          dist = baseline,
          obsdata=obsdata
        ))
    if(G==1){
      AFT_parameters <- as.matrix(Reduce(x = params_transformed, "rbind"))
    } else{
      AFT_parameters <- t(Reduce(x = params_transformed, "rbind"))
    }

    colnames(AFT_parameters) <- 1:G
    parameters <- list(tau = tau, AFT_parameters = AFT_parameters)

    res <-
      list(
        loglik = loglik,
        parameters = parameters,
        z = z,
        class = mclust::map(z),
        baseline=baseline,
        frailty=frailty,
        loglik_vec = loglik_vec
      )

  }
