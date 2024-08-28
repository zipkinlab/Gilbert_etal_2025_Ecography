# Run multispecies, spatially-stratified capture-recapture model
# NOTE: This is a computationally intensive script
# It was run on a supercomputer using 3 cores and 20GB of memory per core, and took ~1.5 days to run
# This script will not be feasible to run on a laptop or typical desktop

library(tidyverse)
library(nimble)
library(parallel)
library(here)

setwd(here::here("data"))
load("neon_cr_data_2024-08-27.RData")

zst <- final %>%
  dplyr::group_by(siteID, site, plotID, plot, period, scientificName, sp, ind) %>% 
  dplyr::summarise( zst = max(y)) %>% 
  dplyr::arrange(siteID, period, sp, ind) %>% 
  pull(zst)

plotst <- get_z_index %>% 
  pull(plotst)

mu_beta1_st <- rnorm(1, 0, 0.05)
sd_beta1_st <- rexp(1, 10)
mu_beta2_st <- rnorm(1, 0, 0.1)
sd_beta2_st <- rexp(1, 10)
mu_beta3_st <- rnorm(1, 0, 0.1)
sd_beta3_st <- rexp(1, 10)
mu_beta1_sp_st <- rnorm(constants$nsp, mu_beta1_st, sd_beta1_st )
sd_beta1_sp_st <- rexp(constants$nsp, 5)
mu_gamma1_st <- rnorm(1, -5, 0.05)
sd_gamma1_st <- rexp(1, 20)
gamma2_st <- rnorm(1, 10, 0.1)
beta1_st <- array(NA,
                  dim = c(constants$nsp, constants$nsite, max(constants$max_period), max(constants$max_plot)))
for(h in 1:constants$nsp){
  for(s in 1:constants$nsite){
    for(t in 1:max(constants$max_period)){
      for(j in 1:max(constants$max_plot)){
        beta1_st[h,s,t,j] <- rnorm(1, mu_beta1_sp_st[h], sd_beta1_sp_st[h])
      }
    }
  }
}
gamma1_st <- rnorm(constants$nsp, mu_gamma1_st, sd_gamma1_st)
beta2_st <- rnorm(constants$nsp, mu_beta2_st, sd_beta2_st)
beta3_st <- rnorm(constants$nsp, mu_beta3_st, sd_beta3_st)
mu_alpha1_st <- rnorm(1, 0, 0.5)
sd_alpha1_st <- rexp(1, 10)
mu_alpha2_st <- rnorm(1, 0, 0.1)
sd_alpha2_st <- rexp(1, 10)
mu_alpha3_st <- rnorm(1, 0, 0.1)
sd_alpha3_st <- rexp(1, 10)
mu_alpha4_st <- rnorm(1, 0, 0.1)
sd_alpha4_st <- rexp(1, 10)
alpha1_st <- rnorm(constants$nsp, mu_alpha1_st, sd_alpha1_st)
alpha2_st <- rnorm(constants$nsp, mu_alpha2_st, sd_alpha2_st)  
alpha3_st <- rnorm(constants$nsp, mu_alpha3_st, sd_alpha3_st)
alpha4_st <- rnorm(constants$nsp, mu_alpha4_st, sd_alpha4_st)
lambda_st <- N_st <- probs_st <- array(NA,
                                       dim = c(constants$nsp,
                                               constants$nsite,
                                               max(constants$max_period, na.rm = TRUE),
                                               max(constants$max_plot, na.rm = TRUE)))
phi_st <- theta_st <- array(NA, dim = c(constants$nsp, constants$nsite))
psi_st <- array(NA, dim = c(constants$nsp, constants$nsite, max(constants$max_period)))
for(h in 1:constants$nsp){
  for(s in 1:constants$nsite){ 
    phi_st[h,s] <- ilogit( gamma1_st[h] + gamma2_st*data$in_range[h,s])
    theta_st[h,s] <- rbinom(1, 1, phi_st[h,s])
    for( t in 1:max(constants$max_period, na.rm = TRUE)){ # t = trapping period (3-day burst)
      for( j in 1:max(constants$max_plot, na.rm = TRUE)) { # j = NEON plotID within site
        lambda_st[h, s, t, j] <- exp( beta1_st[h,s,t,j] + beta2_st[h] * data$CANOPY[s,j] + beta3_st[h] * data$CANOPY[s,j]* data$CANOPY[s,j]) * data$in_range[h,s] # CONSTRAIN BY WHETHER OR NOT NEON SITE IS IN SPECIES' RANGE 
        N_st[h, s, t, j] <- rpois(1, lambda_st[h, s, t, j] )
      }
      psi_st[h,s,t] <- sum( lambda_st[h,s,t,1:constants$max_plot[s,t]]) / ( constants$M[h,s])
    }
  }
}

inits <- function(){
  list(
    mu_beta1 = mu_beta1_st, 
    sd_beta1 = sd_beta1_st, 
    mu_beta2 = mu_beta2_st, 
    sd_beta2 = sd_beta2_st,
    mu_beta3 = mu_beta3_st, 
    sd_beta3 = sd_beta3_st, 
    mu_alpha1 = mu_alpha1_st, 
    sd_alpha1 = sd_alpha1_st, 
    mu_alpha2 = mu_alpha2_st, 
    sd_alpha2 = sd_alpha2_st, 
    mu_alpha3 = mu_alpha3_st, 
    sd_alpha3 = sd_alpha3_st, 
    mu_alpha4 = mu_alpha4_st, 
    sd_alpha4 = sd_alpha4_st, 
    mu_beta1_sp = mu_beta1_sp_st,
    sd_beta1_sp = sd_beta1_sp_st,
    beta1 = beta1_st, 
    beta2 = beta2_st, 
    beta3 = beta3_st,
    alpha1 = alpha1_st, 
    alpha2 = alpha2_st, 
    alpha3 = alpha3_st, 
    alpha4 = alpha4_st,
    lambda = lambda_st,
    N = N_st, 
    z = zst,
    plot = plotst
  )
}

model.code <- nimbleCode({
  mu_beta1 ~ dnorm(0, sd = 2)
  sd_beta1 ~ dexp(1)
  mu_beta2 ~ dnorm(0, sd = 0.5)
  sd_beta2 ~ dexp(1)
  mu_beta3 ~ dnorm(0, sd = 0.5)
  sd_beta3 ~ dexp(1)
  mu_alpha1 ~ dnorm(0, sd = 2)
  sd_alpha1 ~ dexp(1)
  mu_alpha2 ~ dnorm(0, sd = 0.5)
  sd_alpha2 ~ dexp(1)
  mu_alpha3 ~ dnorm(0, sd = 0.5)
  sd_alpha3 ~ dexp(1)
  mu_alpha4 ~ dnorm(0, sd = 0.5)
  sd_alpha4 ~ dexp(1)
  for(h in 1:nsp){
    mu_beta1_sp[h] ~ dnorm( mu_beta1, sd = sd_beta1 )
    sd_beta1_sp[h] ~ dexp(1)
    beta2[h]  ~ dnorm( mu_beta2,  sd = sd_beta2  )
    beta3[h]  ~ dnorm( mu_beta3,  sd = sd_beta3  )
    alpha1[h] ~ dnorm( mu_alpha1, sd = sd_alpha1 )
    alpha2[h] ~ dnorm( mu_alpha2, sd = sd_alpha2 )
    alpha3[h] ~ dnorm( mu_alpha3, sd = sd_alpha3 )
    alpha4[h] ~ dnorm( mu_alpha4, sd = sd_alpha4 )
    for(s in 1:nsite){ 
      for(t in 1:max_period[s]){
        psi[h,s,t] <- sum( lambda[h, s, t, 1:max_plot[s,t]]) / M[h,s]
        for( j in 1:max_plot[s, t] ) { # j = NEON plotID within site
          beta1[h,s,t,j] ~ dnorm( mu_beta1_sp[h], sd = sd_beta1_sp[h] )
          lambda[h, s, t, j] <- exp( beta1[h,s,t,j] + beta2[h] * CANOPY[s,j] + beta3[h]*CANOPY[s,j]*CANOPY[s,j] ) * in_range[h,s] # CONSTRAIN BY WHETHER OR NOT NEON SITE IS IN SPECIES' RANGE
          N[h, s, t, j] ~ dpois( lambda[h, s, t, j] )
          probs[h, s, t, j] <- lambda[h, s, t, j] / sum( lambda[h, s, t, 1:max_plot[s, t]] )
        }
      }
    }
  }
  for( i in 1:nz ) {
    z[i] ~ dbern( psi[ species_z[i], site_z[i], period_z[i]] )
    plot[i] ~ dcat ( probs[ species_z[i], site_z[i], period_z[i], 1:max_plot_z[i]])
  }
  for(i in 1:ny) {
    logit(p[i]) <- alpha1[species_y[i]] + alpha2[species_y[i]]*TMIN[i] + alpha3[species_y[i]]*TMIN[i]*TMIN[i] + alpha4[species_y[i]]*PRCP[i]
    y[i] ~ dbern( p[i] * z[z_index[i]])
  }
})

params <- c(
  "mu_beta1", "sd_beta1", 
  "mu_beta2", "sd_beta2", 
  "mu_beta3", "sd_beta3",
  "mu_beta1_sp", "sd_beta1_sp",
  "beta1", "beta2", "beta3",
  "mu_alpha1", "sd_alpha1", 
  "mu_alpha2", "sd_alpha2", 
  "mu_alpha3", "sd_alpha3", 
  "mu_alpha4", "sd_alpha4",
  "alpha1", "alpha2", 
  "alpha3", "alpha4",
  "N", "psi")

nc <- 3
nburn <- 50000
ni <- nburn + 100000
nt <- 100

# model <- nimbleModel(code = model.code, 
#                      name = "model.code", 
#                      constants = constants, 
#                      data = data, 
#                      inits = inits())

start <- Sys.time()
cl <- makeCluster(nc)

parallel::clusterExport(cl, c("model.code",
                              "inits", 
                              "data", 
                              "constants", 
                              "params", 
                              "nburn", 
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
  
  model <- nimbleModel(code = model.code, 
                       name = "model.code", 
                       constants = constants, 
                       data = data, 
                       inits = init)
  
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  out1 <- runMCMC(CmodelMCMC, 
                  nburnin = nburn, 
                  niter = ni, 
                  thin = nt)
  
  return(as.mcmc(out1))
})
end <- Sys.time()
print(end - start)
stopCluster(cl)

settings <- list(
  n.chains = nc,
  n.iterations = ni,
  n.burnin = nburn,
  n.thin = nt
)

save( model.code,
      data,
      constants,
      final, 
      get_z_index,
      out, 
      settings,
      file = paste0("neon_capture_recapture_results_", Sys.Date(), ".RData"))