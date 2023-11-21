# Fit the extended GAI to UKBMS data

library(data.table)
library(ggplot2)
library(gridExtra)
library(speedglm)

# Functions for GAI fitting
source("R/extendedGAI_functions.R")
source("R/singleyearGAI_functions.R")

spp <- "Gatekeeper"
years <- 1976:2022

set.seed(123)

# Data and prep----
counts <- readRDS("Data/UKBMS_counts_2sp_1976-2022.rds")

# Data prep ----
counts_sp <- counts[SPECIES == spp]
# Number of visits per season (weeks)
nT <- uniqueN(counts_sp$WEEKNO)

# Fit separate (univoltine Normal) GAI to each year  ----
# Use these outputs as starting values to multiyear model fit
output_sep <- NULL
for(k in years){
  counts_sp_k <- counts_sp[YEAR == k]
  mdf_k <- fit_model_df(df = counts_sp_k, nT = nT)
  mdf_k$YEAR <- k
  output_sep[[paste0(k)]] <- mdf_k
}


# Fit extended GAI ----

# 1. Fit model with fixed mu and sigma to get beta values for use as starting point to fully annual model
mdffa <- fit_it_model_eGAI(count_df = counts_sp, nT = nT, years = years,
                         mu.type = "fixed",
                         sigma.type = "fixed",
                         beta.type = "annual",
                         mu.st = mean(sapply(output_sep, function(x)x$mu.out)),
                         sigma.st = mean(sapply(output_sep, function(x)x$sigma.out)),
                         nstart = 1,
                         parallel = TRUE,
                         ncpus = 4)


# 2. Fit full extended GAI with annual parameters
# Use estimates from separate GAI as starting values for mu and sigma
# Use estimates from simpler extended GAI (above) for starting values for beta
mda <- fit_it_model_eGAI(count_df = counts_sp, nT = nT, years = years,
                       mu.type = "annual",
                       sigma.type = "annual",
                       beta.type = "annual",
                       mu.st = sapply(output_sep, function(x)x$mu.out),
                       sigma.st = sapply(output_sep, function(x)x$sigma.out),
                       beta.st = tail(mdffa$betas.out$exp_beta, -1),
                       nstart = 1,
                       parallel = TRUE,
                       ncpus = 4)

saveRDS(mda, paste0("Output/", spp, "_myGAI_annual",
                    "_", min(years),"to", max(years),".rds"))


# Fit extended model with linear beta
mdaal <- fit_it_model_eGAI(count_df = counts_sp, nT = nT, years = years,
                         mu.type = "annual",
                         sigma.type = "annual",
                         beta.type = "linear",
                         mu.st = sapply(output_sep, function(x)x$mu.out),
                         sigma.st = sapply(output_sep, function(x)x$sigma.out),
                         nstart = 1,
                         parallel = TRUE,
                         ncpus = 4)

saveRDS(mdaal, paste0("Output/", spp, "_myGAI_aal",
                      "_", min(years),"to", max(years),".rds"))







