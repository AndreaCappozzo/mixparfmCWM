#' @export
fit_mixparfm <-
  function(formula,
           G,
           class_init,
           grouping_variable = NULL,
           strata = NULL,
           X_gaussian_variables = NULL,
           X_multinomial_variables = NULL, # FIXME not yet implemented
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

    # Retrieve marginal variables
    # Continous covariates
    if(!is.null(X_gaussian_variables)){
      X_gaussian <-   data[,X_gaussian_variables]
      mclust_model_name <- ifelse(length(X_gaussian_variables)==1,"V","VVV")
    }



    # EM controls
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err
    E_step_update <- control_EM_algorithm$E_step_update

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
      data_g <- data[class_init == g,]

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
    }

    # X_gaussian_params_flexCWM  <- flexCWM:::.PX_norm(colXn=ncol(X_gaussian),Xnorm=X_gaussian,modelXnorm="VVV",z=z_init,k=G,n=nrow(data),eps=10e-6)
    # mclust::covw(X=X_gaussian,Z=z_init, normalize = FALSE)
    X_gaussian_params  <- mclust::mstep(data = X_gaussian, z = z_init,modelName = mclust_model_name)

    while (crit) {

# E step ------------------------------------------------------------------


      time_log_density <-
        sapply(1:G, function(g)
          parfm::loglikelihood_i(
            p = parfm_params[[g]],
            obs = obsdata,
            dist = baseline,
            frailty = frailty,
            correct = correct,
            transform = TRUE
          ))

      X_gaussian_log_density <- tryCatch(
        do.call(mclust::cdens, c(
          list(data = X_gaussian,
               logarithm = TRUE), X_gaussian_params
        )),
        error = function(e)
          - Inf
      )

      log_comp_mixture <-
        sweep(
          (time_log_density+X_gaussian_log_density),
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
      if(E_step_update=="classification"){
        map_z <- mclust::map(z)
      } else if(E_step_update=="stochastic"){
        map_z <- mclust::map(t(apply(z, 1, function(prob)
          rmultinom(1, 1,
                    prob))))
        # FIXME check the "average" strategy of weibullRMM_SEM in mixtools \cite{Bordes2016} to see whether it may improve the estimates
      }

      z <- mclust::unmap(map_z)

      data_g_list <- split(data,map_z)


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

        # FIXME add MLE qualitative covariates
      }
      # MLE for continuous covariates
      X_gaussian_params  <-
        mclust::mstep(data = X_gaussian,
                      z = z,
                      modelName = mclust_model_name)


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
      log_density_X_gaussian_g <-
        colSums(do.call(mclust::cdens, c(
          list(data = X_gaussian,
               logarithm = TRUE), X_gaussian_params
        ))*z)

      loglik <- sum(log_density_parfm_g+log_density_X_gaussian_g+log(tau))

      # Check convergence
      err <- abs(loglik - loglik_prev) / (1 + abs(loglik))
      loglik_prev <- loglik
      loglik_vec <- c(loglik_vec, loglik)
      iter <- iter + 1

      crit <- err > tol & iter < itermax
    }

    parfm_params_transformed <-
      lapply(1:G, function(g)
        parfm::paramaters_transformator(
          estim_par = parfm_params[[g]],
          frailty = frailty,
          dist = baseline,
          obsdata=obsdata
        ))
    if(G==1){
      AFT_parameters <- as.matrix(Reduce(x = parfm_params_transformed, "rbind"))
    } else{
      AFT_parameters <- t(Reduce(x = parfm_params_transformed, "rbind"))
    }

    colnames(AFT_parameters) <- 1:G
    parameters <-
      list(
        tau = tau,
        AFT_parameters = AFT_parameters,
        X_gaussian_parameters = list(
          mu = X_gaussian_params$parameters$mean,
          sigma = X_gaussian_params$parameters$variance$sigma
        )
      )
    frailty_effect_df_list <- vector(mode = "list",length = G)

    for(g in 1:G){
      # all this mess is done as we may have some clusters for which no obs fall within a given group, and I want to keep them all
      frailty_g <- predict.parfm(fit_parfm_g_list[[g]])
      frailty_effect_df_list[[g]] <- data.frame(group = attributes(frailty_g)$names, frailty = c(frailty_g))
      colnames(frailty_effect_df_list[[g]])[2] <- g
    }

    frailty_effect=plyr::join_all(frailty_effect_df_list,by = "group",type = "full")
    # if no frailty_effect is present for a given group in a cluster (i.e., no obs from that group belong to the g-th cluster) it returns NA
    colnames(frailty_effect)[1] <- grouping_variable # I manually specify the grouping variable name


    res <-
      list(
        loglik = loglik,
        parameters = parameters,
        frailty_effect=frailty_effect,
        z = z,
        class = mclust::map(z),
        baseline=baseline,
        frailty=frailty,
        loglik_vec = loglik_vec
      )

  }
