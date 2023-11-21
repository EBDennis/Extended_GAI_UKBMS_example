# Functions to fit the extended GAI

library(data.table)
library(parallel)
library(optimParallel)

# Likelihood function
ll_func_eGAI <- function(parm,
                       count_df, # data
                       nT, # number of visits (weeks)
                       years, # vector of years
                       mu.type, # should mu be fixed, annual etc
                       sigma.type, # should sigma be fixed, annual etc
                       beta.type # should beta be linear or annual
){

  nyears <- length(years)

  avals <- data.table(YEAR = rep(years, each = nT),
                      WEEKNO = rep(1:nT, nyears),
                      NM = dnorm(rep(1:nT, nyears),
                                 switch(mu.type,
                                        "fixed"= {exp(parm["mu"])},
                                        "annual"={exp(rep(parm[which(names(parm)=="mu")], each = nT))}),
                                 switch(sigma.type,
                                        "fixed"= {exp(parm["sigma"])},
                                        "annual"={exp(rep(parm[which(names(parm)=="sigma")], each = nT))})))

  dfa <- count_df[avals, on = .(YEAR, WEEKNO)]

  betas <- data.table(YEAR = years, exp_beta =
                        switch(beta.type,
                               "annual" = {exp(c(0, parm[which(names(parm)=="beta")]))},
                               "linear" = {exp(parm["beta"]*(years - min(years) + 1))}))
  dfa <- dfa[betas, on = "YEAR"]

  # Concentrated likelihood formulation
  dfa[, exp_alpha := sum(COUNT)/sum(exp_beta*NM), by = SITENO]
  dfa[, N_est := exp_alpha*exp_beta]
  dfa[, lambda := NM*N_est]
  dfa[, loglik := dpois(COUNT, lambda = lambda, log = TRUE)]
  return(-sum(dfa$loglik, na.rm=TRUE))

}


# Starting values
start_val_func_eGAI <- function(mu.type,
                              mu.st = NULL, # optional starting point for mu
                              sigma.type,
                              sigma.st = NULL, # optional starting point for mu
                              beta.type,
                              beta.st = NULL,
                              years){

  nyears <- length(years)

  if(is.null(mu.st)){
    switch(mu.type,
           "fixed" = {mu_st <- sample(5:20, 1)},
           "annual" = {mu_st <- sample(5:20, nyears, replace = TRUE)})
  } else {
    mu_st <- mu.st
  }
  names(mu_st) <- rep("mu", length(mu_st))


  if(is.null(sigma.st)){
    switch(sigma.type,
           "fixed" = {sigma_st <- sample(1:5, 1)},
           "annual" = {sigma_st <- sample(1:5, nyears, replace = TRUE)})
  } else {
    sigma_st <- sigma.st
  }
  names(sigma_st) <- rep("sigma", length(sigma_st))


  if(is.null(beta.st)){
    switch(beta.type,
           "annual" = {beta_st <- log(rep(1, length(years) - 1))},
           "linear" = {beta_st <- 0})
  } else {
    beta_st <- log(beta.st)
  }
  names(beta_st)<- rep("beta", length(beta_st))

  parm <- c(log(mu_st), log(sigma_st), beta_st)

  return(parm)
}



# Wrapper for fitting extended GAI for multiple starts
fit_it_model_eGAI <- function(count_df, nT, years,
                            mu.type, mu.st = NULL,
                            sigma.type, sigma.st = NULL,
                            beta.type, beta.st = NULL,
                            nstart = 3,
                            parallel = FALSE,
                            ncpus = 2){
  if(nstart == 1){
    st <- proc.time()
    output <- try(fit_model_eGAI(count_df, nT, years,
                               mu.type, mu.st,
                               sigma.type, sigma.st,
                               beta.type, beta.st,
                               parallel, ncpus), silent=FALSE)
    et <- proc.time()
    output$time <- (et-st)[3]
  } else {
    fit_k <- list(); fit_k.ll <- rep(NA, nstart)
    for(k in 1:nstart){
      st <- proc.time()
      fit1 <- try(fit_model_eGAI(count_df, nT, years,
                               mu.type, mu.st,
                               sigma.type, sigma.st,
                               beta.type, beta.st,
                               parallel, ncpus), silent=FALSE)
      et <- proc.time()
      fit1$time <- (et-st)[3]
      fit_k[[k]] <- fit1
      if(!is.na(fit_k.ll[[k]][[1]]))fit_k.ll[k] <- fit_k[[k]]$ll.val
    }
    output <- list(fit_k[[min(c(1:nstart)[fit_k.ll==max(fit_k.ll,na.rm=T)], na.rm=TRUE)]],
                   fit_k,
                   fit_k.ll)
  }
  return(output)
}

# Fit the extended GAI model - optionally in parallel
fit_model_eGAI <- function(count_df, nT, years,
                         mu.type, mu.st,
                         sigma.type, sigma.st,
                         beta.type, beta.st,
                         parallel = FALSE,
                         ncpus = 2){

  parm <- start_val_func_eGAI(mu.type = mu.type,
                            mu.st = mu.st,
                            sigma.type = sigma.type,
                            sigma.st = sigma.st,
                            beta.type = beta.type,
                            beta.st = beta.st,
                            years = years
  )

  if(!parallel){
    this.fit <- try(optim(par = parm,
                          fn = ll_func_eGAI,
                          hessian = TRUE,
                          method = "L-BFGS-B",
                          control = list(maxit=100000),
                          count_df = count_df,
                          nT = nT,
                          years = years,
                          mu.type = mu.type,
                          sigma.type = sigma.type,
                          beta.type = beta.type),
                    silent=TRUE)
  } else {
    cl <- makeCluster(ncpus) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    clusterEvalQ(cl, library(data.table))

    this.fit <- try(optimParallel::optimParallel(par = parm,
                                                 fn = ll_func_eGAI,
                                                 hessian = TRUE,
                                                 control = list(maxit=100000),
                                                 count_df = count_df,
                                                 nT = nT,
                                                 years = years,
                                                 mu.type = mu.type,
                                                 sigma.type = sigma.type,
                                                 beta.type = beta.type),
                    silent=TRUE)
    setDefaultCluster(cl=NULL); stopCluster(cl)
  }

  if(is.list(this.fit) & class(try(solve(this.fit$hessian),silent=TRUE))[1] != "try-error"){
    # Model output
    nyears <- length(years)

    a.out <- data.table(YEAR = rep(years, each = nT),
                        WEEKNO = rep(1:nT, nyears),
                        NM = dnorm(rep(1:nT, nyears),
                                   switch(mu.type,
                                          "fixed"= {exp(this.fit$par["mu"])},
                                          "annual"={exp(rep(this.fit$par[which(names(this.fit$par)=="mu")], each = nT))}),
                                   switch(sigma.type,
                                          "fixed"= {exp(this.fit$par["sigma"])},
                                          "annual"={exp(rep(this.fit$par[which(names(this.fit$par)=="sigma")], each = nT))})))


    dfa.out <- count_df[a.out, on = .(WEEKNO, YEAR)]

    if(beta.type == "annual"){
      betas.out <- data.table(YEAR = years,
                              beta = c(0, this.fit$par[which(names(parm)=="beta")]),
                              beta_SE =  c(0, sqrt(diag(solve(this.fit$hessian))[which(names(parm)=="beta")])),
                              exp_beta = exp(c(0, this.fit$par[which(names(parm)=="beta")])),
                              exp_beta_SE = c(0, sapply(which(names(parm)=="beta"),
                                                        function(i)
                                                          msm::deltamethod(~exp(x1),
                                                                           this.fit$par[i],
                                                                           diag(solve(this.fit$hessian))[i]))))
    }

    if(beta.type == "linear"){
      syears <- years - min(years) + 1
      betas.out <- data.table(YEAR = years,
                              beta = this.fit$par[which(names(parm)=="beta")]*syears,
                              beta_SE =  sqrt(syears*syears*diag(solve(this.fit$hessian))[which(names(parm)=="beta")]),
                              exp_beta = exp(this.fit$par[which(names(parm)=="beta")]*syears))
    }


    dfa.out <- dfa.out[betas.out, on = "YEAR"]
    dfa.out[, exp_alpha := sum(COUNT)/sum(exp_beta*NM), by = SITENO]
    dfa.out[, N := exp_alpha*exp_beta]
    dfa.out[, Fitted := NM*N]

    output <- list(ll.val = -this.fit$value,
                   model_spec = list(mu.type = mu.type,
                                     sigma.type = sigma.type,
                                     beta.type = beta.type),
                   npar = length(this.fit$par),
                   pars = this.fit$par,
                   alpha.out = unique(dfa.out[, .(SITENO, exp_alpha)]),
                   N.out = unique(dfa.out[, .(SITENO, YEAR, N)]),
                   mu.out =  exp(this.fit$par[which(names(this.fit$par)=="mu")]),
                   sigma.out = exp(this.fit$par[which(names(this.fit$par)=="sigma")]),
                   betas.out = betas.out,
                   a.out = a.out,
                   modelfit = this.fit,
                   output_df = dfa.out,
                   starting_vals = parm)

    output$deviance <- 2*(dfa.out[COUNT != 0, sum(COUNT*log(COUNT/Fitted) - (COUNT-Fitted))] +
                            dfa.out[COUNT == 0, sum(Fitted)])

    output$D <- output$deviance/(sum(!is.na(dfa.out$COUNT)) - output$npar)

    return(output)
  } else {NA}
}

