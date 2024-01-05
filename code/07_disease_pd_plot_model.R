library(tidyverse)
library(parallel)
library(nimble)

# setwd(here::here("data"))
final <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv")

data <- list(
  y = final$positive, 
  pd_plot_mean = final$pd_plot_mean, 
  pd_plot_sd = final$pd_plot_sd)

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
  mean_pd <- mean( pd_plot[1:nind] )
  sd_pd <- sd( pd_plot[1:nind] )
  for( i in 1:nind ) {
    pd_plot[i] ~ T( dnorm( pd_plot_mean[i], sd = pd_plot_sd[i] ), 0, 36 )
    pd_plot_scaled[i] <- ( pd_plot[i] - mean_pd ) / sd_pd
    y[i] ~ dbern( kappa[i] )
    logit( kappa[i] ) <- gamma0[ site[i] ] + gamma1 * pd_plot_scaled[i]
  }
})

inits <- function(){
  list(
    mu_gamma0 = rnorm(1, 0, 0.1),
    pd_plot = data$pd_plot_mean,
    mean_pd = mean( data$pd_plot_mean), 
    sd_pd = mean( data$pd_plot_sd),
    sd_gamma0 = runif(1, 0, 1), 
    gamma1 = rnorm(1, 0, 0.25), 
    gamma0 = rnorm(constants$nsite, 0, 1)
  )
}

params <- c("mu_gamma0", "sd_gamma0", "gamma1", "gamma0", "mean_pd", "sd_pd")

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
      file = paste0("rodent_pathogen_pd_plot_", Sys.Date(), ".RData"))