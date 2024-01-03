library(tidyverse)
library(nimble)
library(parallel)

load("neon_cr_data_2024-01-03.RData")

zst <- final %>%
  dplyr::group_by(siteID, site, plotID, plot, period, scientificName, sp, ind) %>% 
  dplyr::summarise( zst = max(y)) %>% 
  dplyr::arrange(siteID, period, sp, ind) %>% 
  pull(zst)

plotst <- get_z_index %>% 
  pull(plotst)

mu_beta1_st <- rnorm(1, 0, 1)
sd_beta1_st <- rexp(1, 10)
mu_beta2_st <- rnorm(1, 0, 0.2)
sd_beta2_st <- rexp(1, 10)
beta1_st <- array( rnorm(constants$nsp * constants$nsite, mu_beta1_st, sd_beta1_st), 
                   dim = c(constants$nsp, constants$nsite))
beta2_st <- rnorm(constants$nsp, mu_beta2_st, sd_beta2_st)
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
sd_theta_st <- rexp(1, 25)
theta_st <- array( rnorm(constants$nsp * constants$nsite * (max(constants$max_period, na.rm = TRUE) - 1), 1, sd_theta_st), 
                   dim = c(constants$nsp, constants$nsite, (max(constants$max_period, na.rm = TRUE) - 1)))

lambda_st <- N_st <- probs_st <- array(NA,
                                       dim = c(constants$nsp,
                                               constants$nsite,
                                               max(constants$max_period, na.rm = TRUE),
                                               max(constants$max_plot, na.rm = TRUE)))
for(h in 1:constants$nsp){
  for(s in 1:constants$nsite){ 
    for( j in 1:max(constants$max_plot, na.rm = TRUE)) { # j = NEON plotID within site
      lambda_st[h, s, 1, j] <- exp( beta1_st[h,s] + beta2_st[h] * data$CANOPY[s,j]) * data$in_range[h,s] # CONSTRAIN BY WHETHER OR NOT NEON SITE IS IN SPECIES' RANGE 
      N_st[h, s, 1, j] <- rpois(1, lambda_st[h, s, 1, j] )
    }
  }
}
for(h in 1:constants$nsp){
  for(s in 1:constants$nsite){ 
    for( t in 2:max(constants$max_period, na.rm = TRUE)){ # t = trapping period (3-day burst)
      for( j in 1:max(constants$max_plot, na.rm = TRUE)) { # j = NEON plotID within site
        lambda_st[h, s, t, j] <- lambda_st[h,s,t-1,j] * theta_st[h,s,t-1]
        N_st[h, s, t, j] <- rpois(1, lambda_st[h, s, t, j] )
      }
    }
  }
}

inits <- function(){
  list(
    mu_beta1 = mu_beta1_st, 
    sd_beta1 = sd_beta1_st, 
    mu_beta2 = mu_beta2_st, 
    sd_beta2 = sd_beta2_st, 
    mu_alpha1 = mu_alpha1_st, 
    sd_alpha1 = sd_alpha1_st, 
    mu_alpha2 = mu_alpha2_st, 
    sd_alpha2 = sd_alpha2_st, 
    mu_alpha3 = mu_alpha3_st, 
    sd_alpha3 = sd_alpha3_st, 
    mu_alpha4 = mu_alpha4_st, 
    sd_alpha4 = sd_alpha4_st, 
    beta2 = beta2_st, 
    alpha1 = alpha1_st, 
    alpha2 = alpha2_st, 
    alpha3 = alpha3_st, 
    alpha4 = alpha4_st,
    sd_theta = sd_theta_st, 
    theta = theta_st, 
    beta1 = beta1_st, 
    lambda = lambda_st, 
    N = N_st, 
    z = zst,
    plot = plotst
  )
}

model.code <- nimbleCode({
  sd_theta ~ dexp(25)
  mu_beta1 ~ dnorm(0, sd = 2)
  sd_beta1 ~ dexp(1)
  mu_beta2 ~ dnorm(0, sd = 0.5)
  sd_beta2 ~ dexp(1)
  mu_alpha1 ~ dnorm(0, sd = 2)
  sd_alpha1 ~ dexp(1)
  mu_alpha2 ~ dnorm(0, sd = 0.5)
  sd_alpha2 ~ dexp(1)
  mu_alpha3 ~ dnorm(0, sd = 0.5)
  sd_alpha3 ~ dexp(1)
  mu_alpha4 ~ dnorm(0, sd = 0.5)
  sd_alpha4 ~ dexp(1)
  for(h in 1:nsp){
    beta2[h]  ~ dnorm( mu_beta2,  sd = sd_beta2  )
    alpha1[h] ~ dnorm( mu_alpha1, sd = sd_alpha1 )
    alpha2[h] ~ dnorm( mu_alpha2, sd = sd_alpha2 )
    alpha3[h] ~ dnorm( mu_alpha3, sd = sd_alpha3 )
    alpha4[h] ~ dnorm( mu_alpha4, sd = sd_alpha4 )
    for(s in 1:nsite){ 
      beta1[h,s] ~ dnorm( mu_beta1, sd = sd_beta1 )
      psi[h, s, 1] <- sum(lambda[h, s, 1, 1:max_plot[s, 1]]) / M[h, s]
      for( j in 1:max_plot[s, 1] ) { # j = NEON plotID within site
        lambda[h, s, 1, j] <- exp( beta1[h,s] + beta2[h] * CANOPY[s,j] ) * in_range[h,s] # CONSTRAIN BY WHETHER OR NOT NEON SITE IS IN SPECIES' RANGE 
        N[h, s, 1, j] ~ dpois( lambda[h, s, 1, j] )
        probs[h, s, 1, j] <- lambda[h, s, 1, j] / sum( lambda[h, s, 1, 1:max_plot[s, 1]] )
      }
    }
  }
  for(h in 1:nsp){
    for(s in 1:nsite){ 
      for( t in 2:max_period[s] ){ # t = trapping period (3-day burst)
        psi[h, s, t] <- sum(lambda[h, s, t, 1:max_plot[s, t]]) / M[h, s]
        theta[h,s,t-1] ~ dnorm( 1, sd = sd_theta ) 
        for( j in 1:max_plot[s, t] ) { # j = NEON plotID within site
          lambda[h, s, t, j] <- lambda[h,s,t-1,j] * theta[h,s,t-1]
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
  "beta1", "beta2", 
  "mu_alpha1", "sd_alpha1", 
  "mu_alpha2", "sd_alpha2", 
  "mu_alpha3", "sd_alpha3", 
  "mu_alpha4", "sd_alpha4",
  "alpha1", "alpha2", 
  "alpha3", "alpha4", 
  "sd_theta", 
  "theta",
  "N")

nc <- 3
nburn <- 75000
ni <- nburn + 100000
nt <- 175

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
