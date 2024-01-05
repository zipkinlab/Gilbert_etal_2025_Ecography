library(tidyverse)
library(parallel)
library(nimble)

# setwd(here::here("data"))
final <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv")

data <- list(
  y = final$positive, 
  shan_site_mean = final$shan_site_mean, 
  shan_site_sd = final$shan_site_sd)

constants <- list(
  nsite = length(unique(final$site)), 
  nind = nrow(final), 
  site = final$site)

code <- nimbleCode({
  mu_gamma0 ~ dnorm(0, sd = 1)
  sd_gamma0 ~ dexp(1)
  gamma1 ~ dnorm(0, sd = 0.5)
  for( i in 1:nsite){
    gamma0[i] ~ dnorm( mu_gamma0, sd = sd_gamma0 )
  }
  mean_shan <- mean( shan_site[1:nind] )
  sd_shan <- sd( shan_site[1:nind] )
  for( i in 1:nind ) {
    shan_site[i] ~ T( dnorm( shan_site_mean[i], sd = shan_site_sd[i] ), 0, 36 )
    shan_site_scaled[i] <- ( shan_site[i] - mean_shan ) / sd_shan
    y[i] ~ dbern( kappa[i] )
    logit( kappa[i] ) <- gamma0[ site[i] ] + gamma1 * shan_site_scaled[i]
  }
})

inits <- function(){
  list(
    mu_gamma0 = rnorm(1, 0, 0.1),
    shan_site = data$shan_site_mean,
    mean_shan = mean( data$shan_site_mean), 
    sd_shan = mean( data$shan_site_sd),
    sd_gamma0 = runif(1, 0, 1), 
    gamma1 = rnorm(1, 0, 0.25), 
    gamma0 = rnorm(constants$nsite, 0, 1)
  )
}

params <- c("mu_gamma0", "sd_gamma0", "gamma1", "gamma0", "mean_shan", "sd_shan")

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
      file = paste0("rodent_pathogen_shan_site_", Sys.Date(), ".RData"))