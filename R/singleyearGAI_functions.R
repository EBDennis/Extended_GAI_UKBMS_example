# Functions to fit the univoltine phenomenological GAI (separately to each year)

library(data.table)

# Likelihood function
ll_func_df <- function(parm, df, nT){

  mu <- exp(parm[1])
  sigma = exp(parm[2])

  avals <- data.table(WEEKNO = 1:nT, NM = dnorm(1:nT, mu, sigma))

  dfa <- merge(df, avals, by = "WEEKNO")
  dfa[is.na(COUNT), NM := NA]

  # Concentrated likelihood formulation
  dfa[, N_est := sum(COUNT, na.rm = TRUE)/sum(NM, na.rm = TRUE), by = SITENO]
  dfa[, lambda := NM*N_est]
  dfa[, loglik := dpois(COUNT, lambda = lambda, log = TRUE)]
  return(-sum(dfa$loglik, na.rm=TRUE))
}


# Starting values
start_val_func_df <- function(){

  mu_st <- sample(5:15, 1)
  sigma_st <- sample(1:5, 1)
  parm <- c(log(mu_st), log(sigma_st))

  return(parm)
}


# Wrapper for fitting GAI the model for multiple starts
fit_it_model_df <- function(df, nT, nstart = 3){

  fit_k <- list(); fit_k.ll <- rep(NA, nstart)
  for(k in 1:nstart){
    st <- proc.time()
    fit1 <- try(fit_model_df(df, nT), silent=FALSE)
    et <- proc.time()
    fit1$time <- (et-st)[3]
    fit_k[[k]] <- fit1
    if(!is.na(fit_k[[k]][[1]]))fit_k.ll[k] <- fit_k[[k]]$ll.val
  }
  sapply(fit_k, function(x)x$ll)
  output <- list(fit_k[[min(c(1:nstart)[fit_k.ll==max(fit_k.ll,na.rm=T)], na.rm=TRUE)]],
                 fit_k,
                 fit_k.ll)
  return(output)
}

# Fit the GAI model
fit_model_df <- function(df, nT){

  parm <- start_val_func_df()

  this.fit <- try(optim(par=parm,
                        fn=ll_func_df,
                        hessian=TRUE,
                        method="Nelder-Mead",
                        control=list(maxit=100000),
                        df = df,
                        nT = nT),
                  silent=TRUE)

  if(is.list(this.fit) & class(try(solve(this.fit$hessian),silent=TRUE))[1] != "try-error"){
    # Model output
    N.out <- mu.out <- sigma.out <- NULL

    mu.out <- exp(this.fit$par[1])
    sigma.out <- exp(this.fit$par[2])

    a.out <- data.frame(WEEKNO = 1:nT, NM = dnorm(1:nT, mu.out, sigma.out))

    dfa.out <- merge(df, a.out)
    dfa.out[is.na(COUNT), NM := NA]

    # Concentrated likelihood formulation
    dfa.out[, N := sum(COUNT, na.rm = TRUE)/sum(NM, na.rm = TRUE), by = SITENO]
    dfa.out[, TOTALCOUNT := sum(COUNT, na.rm=TRUE), by = SITENO]
    dfa.out[, TOTALNM := sum(NM, na.rm=TRUE), by = SITENO]
    dfa.out[, Fitted := NM*N]


    output <- list(ll.val=-this.fit$value,
                   npar=length(this.fit$par),
                   pars = this.fit$par,
                   N.out = unique(dfa.out[, .(SITENO, TOTALCOUNT, TOTALNM, N)]),
                   mu.out = mu.out,
                   sigma.out = sigma.out,
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
