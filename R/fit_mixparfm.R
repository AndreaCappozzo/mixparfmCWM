#' @export
fit_mixparfm <-
  function(formula,
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
    data <- as.data.frame(unclass(data),stringsAsFactors = TRUE) # converts all character variables to factors
    # (this automatically handles situations in which for some clusters not all the levels of the qualitative variable are observed)

    N <- nrow(data)
    # Retrieve marginal variables
    # Continous covariates
    is_X_gaussian <- !is.null(X_gaussian_variables)
    is_X_multinomial <- !is.null(X_multinomial_variables)
    if(is_X_gaussian){
      X_gaussian <-   data[,X_gaussian_variables]
      p_gaussian <- length(X_gaussian_variables)
      mclust_model_name <- ifelse(p_gaussian==1,"V","VVV")
    } else{
      # if no continuous random covariates are modeled, the associated loglik contribution is set to 0 and will stay like this
      # throughout the algorithm
      X_gaussian_log_density <- 0
      log_density_X_gaussian_g <- 0
      n_par_X_gaussian <- 0
    }

    if(is_X_multinomial){
      X_multinomial <-   data[,X_multinomial_variables,drop=FALSE]

      p_multinomial <- length(X_multinomial_variables) # the number of variables modeled independently via multinomial dist (also denoted with M hereafter)
      C_M_vec <- apply(X_multinomial, 2, function(X_mult_col) length(unique(X_mult_col)))
      X_multinomial_params <- vector(mode = "list",length = p_multinomial) # container for the C_m x G matrices of probabilities of variable m
      # C_m is the # of categories for variable m, m=1,...M aka p_multinomial (see notation in Section 3 of Variable selection for latent class analysis with application to low back pain diagnosis Fop et al)
      # I need to store them in a list since C_m (the total # of categories) may differ for the p_multinomial variables
      for (m in 1:p_multinomial) {
        X_multinomial_params[[m]] <- matrix(nrow = C_M_vec[m], ncol = G)
        rownames(X_multinomial_params[[m]]) <- names(table(X_multinomial[,m]))
      }
    } else{
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

      if (is_X_multinomial) {
        X_multinomial_g <- X_multinomial[class_init == g,,drop=FALSE]
        for (m in 1:p_multinomial) {
          tab_X_multinomial_gm <- table(factor(X_multinomial_g[,m],levels = rownames(X_multinomial_params[[m]])))
          # this handles the case in which for some levels
          # of the qualitative variable there are 0 entries in cluster g
          X_multinomial_params[[m]][,g] <- tab_X_multinomial_gm/sum(tab_X_multinomial_gm)
        }
      }
    }

    # 3 ways to compute the same quantities
    # X_gaussian_params  <- flexCWM:::.PX_norm(colXn=ncol(X_gaussian),Xnorm=X_gaussian,modelXnorm="VVV",z=z_init,k=G,n=nrow(data),eps=10e-6)
    # X_gaussian_params <- mclust::covw(X=X_gaussian,Z=z_init, normalize = FALSE)
    if (is_X_gaussian) {
      X_gaussian_params  <-
        mclust::mstep(data = X_gaussian,
                      z = z_init,
                      modelName = mclust_model_name)
    }

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
      if(is_X_gaussian) {
        X_gaussian_log_density <- tryCatch(
          do.call(mclust::cdens, c(
            list(data = X_gaussian,
                 logarithm = TRUE), X_gaussian_params
          )),
          error = function(e)
            - Inf
        )
      }

      if(is_X_multinomial) {
        X_multinomial_density_list <- lapply(1:p_multinomial, function(m)
          sapply(1:G, function(g)
            X_multinomial_params[[m]][, g][as.character(X_multinomial[, m])]))
        X_multinomial_log_density <-
          Reduce(f = "+",
                 x = lapply(X_multinomial_density_list, log))
      }

      log_comp_mixture <-
        sweep(
          (time_log_density+X_gaussian_log_density+X_multinomial_log_density),
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

        if (is_X_multinomial) {
          # MLE for categorical covariates
          X_multinomial_g_list <- split(X_multinomial,map_z)
          for (m in 1:p_multinomial) {
            tab_X_multinomial_gm <- table(factor(X_multinomial_g_list[[g]][,m],levels = rownames(X_multinomial_params[[m]])))
            X_multinomial_params[[m]][,g] <- tab_X_multinomial_gm/sum(tab_X_multinomial_gm)
          }
        }

      }

      if (is_X_gaussian) {
        # MLE for continuous covariates
        X_gaussian_params  <-
          mclust::mstep(data = X_gaussian,
                        z = z,
                        modelName = mclust_model_name)
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
          colSums(do.call(mclust::cdens, c(
            list(data = X_gaussian,
                 logarithm = TRUE), X_gaussian_params
          )) * z)
      }

      if (is_X_multinomial) {
        X_multinomial_density_list <- lapply(1:p_multinomial, function(m)
          sapply(1:G, function(g)
            X_multinomial_params[[m]][, g][as.character(X_multinomial[, m])]))
        X_multinomial_log_dens <-
          (Reduce(f = "+",
                 x = lapply(X_multinomial_density_list, log)))*z
        X_multinomial_log_dens_no_Nan <- ifelse(is.nan(X_multinomial_log_dens),no = X_multinomial_log_dens,yes=0)
        # this handles the case in which for some levels
        # of the qualitative variable there are 0 entries in cluster g
        log_density_X_multinomial_g <- colSums(X_multinomial_log_dens_no_Nan)
      }

      n_g <- colSums(z)

      loglik <-
        sum(
          log_density_parfm_g + log_density_X_gaussian_g + log_density_X_multinomial_g +
            log(tau)*n_g
        )

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
        AFT_parameters = AFT_parameters
      )

    if(is_X_gaussian) {
      parameters$X_gaussian_parameters = list(mu = X_gaussian_params$parameters$mean,
                                              sigma = X_gaussian_params$parameters$variance$sigma)
    }

    if(is_X_multinomial) {
      parameters$X_multinomial_parameters = X_multinomial_params
    }

    if(frailty=="none") {
      # in case no frailty is considered simply return NULL for frailty_effect (this is necessary to prevent errors in the call to predict.parfm)
      frailty_effect <- NULL
      frailty_var_effect <- NULL
    } else{
      frailty_effect_df_list <- vector(mode = "list", length = G)
      frailty_var_effect_df_list <- vector(mode = "list", length = G)

      for (g in 1:G) {
        # all this mess is done as we may have some clusters for which no obs fall within a given group, and I want to keep them all
        frailty_g <- predict.parfm(fit_parfm_g_list[[g]])[[1]]
        frailty_var <- predict.parfm(fit_parfm_g_list[[g]])[[2]]

        frailty_effect_df_list[[g]] <-
          data.frame(group = attributes(frailty_g)$names,
                     frailty = c(frailty_g))
        colnames(frailty_effect_df_list[[g]])[2] <- g

        frailty_var_effect_df_list[[g]] <-
          data.frame(group = attributes(frailty_var)$names,
                     frailty_var = c(frailty_var))
        colnames(frailty_var_effect_df_list[[g]])[2] <- g
      }

      frailty_effect = plyr::join_all(frailty_effect_df_list, by = "group", type = "full")
      frailty_var_effect = plyr::join_all(frailty_var_effect_df_list, by = "group", type = "full")

      # if no frailty_effect is present for a given group in a cluster (i.e., no obs from that group belong to the g-th cluster) it returns NA
      colnames(frailty_effect)[1] <-
        grouping_variable # I manually specify the grouping variable name

      colnames(frailty_var_effect)[1] <-
        grouping_variable
    }

    # Compute bic

    n_par_tau <- G - 1
    n_par_parfm <- length(parameters$AFT_parameters)
    if(is_X_gaussian) {
      n_par_X_gaussian <-
        mclust::nMclustParams(modelName = mclust_model_name, d = p_gaussian, G = G) -
        n_par_tau # nMclustParams counts also the G-1 mixture weights, I subtract them otherwise I would count them twice
    }
    if (is_X_multinomial) {
      n_par_X_multinomial <- sum(sapply(X_multinomial_params, length)-G)
    }
    bic_final <- 2*loglik-(n_par_tau+n_par_parfm+n_par_X_gaussian+n_par_X_multinomial)*log(N)

    res <-
      list(
        loglik = loglik,
        parameters = parameters,
        frailty_effect=frailty_effect,
        frailty_var_effect = frailty_var_effect,
        z = z,
        class = mclust::map(z),
        bic=bic_final,
        baseline=baseline,
        frailty=frailty,
        fit_parfm= fit_parfm_g_list, # #FIXME add droplevels + formula_g a list of dimension G containing the estimated parametric frailty models in the G clusters. Useful for assessing significance
        loglik_vec = loglik_vec
      )

  }
