# This script is mildly computationally intensive
# It was run on a supercomputer with 3 cores and 10GB of memory per core, and took ~1 hour
# It would probably be able to run on a desktop in a pinch

library(tidyverse)
library(parallel)
library(nimble)

setwd(here::here("data"))
final <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv") %>% 
  dplyr::group_by(scientificName) %>% 
  dplyr::mutate(sp_disease = cur_group_id())

data <- list(
  y = final$positive, 
  sr_site_mean = final$sr_site_mean, 
  sr_site_sd = final$sr_site_sd)

constants <- list(
  nsp = length(unique(final$sp_disease)),
  nsite = length(unique(final$site)), 
  nind = nrow(final), 
  site = final$site,
  sp = final$sp_disease)

code <- nimbleCode({
  mu_gamma0 ~ dnorm(0, sd = 1)
  sd_gamma0 ~ dexp(1)
  sd_epsilon ~ dexp(1)
  gamma1 ~ dnorm(0, sd = 1)
  
  for(i in 1:nsp){
    gamma0[i] ~ dnorm( mu_gamma0, sd = sd_gamma0 )
  }
  for( i in 1:nsite){
    epsilon[i] ~ dnorm( 0, sd = sd_epsilon )
  }
  mean_sr <- mean( sr_site[1:nind] )
  sd_sr <- sd( sr_site[1:nind] )
  for( i in 1:nind ) {
    sr_site[i] ~ T( dnorm( sr_site_mean[i], sd = sr_site_sd[i] ), 0, 36 )
    sr_site_scaled[i] <- ( sr_site[i] - mean_sr ) / sd_sr
    y[i] ~ dbern( kappa[i] )
    logit(kappa[i]) <- gamma0[sp[i]] + gamma1*sr_site_scaled[i] + epsilon[site[i]]
  }
})

inits <- function(){
  list(
    mu_gamma0 = rnorm(1, 0, 0.1),
    sr_site = data$sr_site_mean,
    mean_sr = mean( data$sr_site_mean), 
    sd_sr = mean( data$sr_site_sd),
    sd_gamma0 = runif(1, 0, 1), 
    gamma1 = rnorm(1, 0, 0.25), 
    gamma0 = rnorm(constants$nsp, 0, 1),
    sd_epsilon = rexp(1),
    epsilon = rnorm(constants$nsite, 0, 1)
  )
}

params <- c("mu_gamma0", "sd_gamma0", "gamma1", "gamma0", "sd_epsilon", "epsilon", "mean_sr", "sd_sr")

# model <- nimbleModel(code = code,
#                      # name = "code",
#                      constants = constants,
#                      data = data,
#                      inits = inits())

nc <- 3
nb <- 10000
ni <- nb + 5000
nt <- 5

library(parallel)
start <- Sys.time()
cl <- makeCluster(nc)

parallel::clusterExport(cl, c("code",
                              "inits", 
                              "data", 
                              "constants", 
                              "params", 
                              "nb", 
                              "ni", 
                              "nt"))

for(j in seq_along(cl)) {
  set.seed(j)
  init <- inits()
  clusterExport(cl[j], "init")
}

out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  model <- nimbleModel(code = code, 
                       name = "code", 
                       constants = constants, 
                       data = data, 
                       inits = init)
  
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  out1 <- runMCMC(CmodelMCMC, 
                  nburnin = nb, 
                  niter = ni, 
                  thin = nt)
  
  return(as.mcmc(out1))
})
end <- Sys.time()
print(end - start)
stopCluster(cl)

save( code, data, constants, out, 
      file = paste0("rodent_pathogen_sr_site_", Sys.Date(), ".RData"))