# Quick script - ran in ~2 minutes on basic desktop
# Check influence of sample type (blood only, ear tissue only, or both)
library(tidyverse)
library(parallel)
library(nimble)
library(here)

setwd(here::here("data"))
final <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv") %>% 
  dplyr::group_by(scientificName) %>% 
  dplyr::mutate(sp_disease = cur_group_id()) |> 
  dplyr::group_by(type) |> 
  dplyr::mutate(method = cur_group_id())

data <- list(
  y = final$positive)

constants <- list(
  nsp = length(unique(final$sp_disease)),
  nsite = length(unique(final$site)), 
  nind = nrow(final), 
  nmethod = length(unique(final$method)),
  site = final$site,
  sp = final$sp_disease,
  method = final$method)

code <- nimbleCode({
  mu_gamma0 ~ dnorm(0, sd = 1)
  sd_gamma0 ~ dexp(1)
  sd_epsilon ~ dexp(1)
  for(i in 1:nmethod){
    zeta[i] ~ dnorm(0, sd = 0.5)
  }
  for(i in 1:nsp){
    gamma0[i] ~ dnorm( mu_gamma0, sd = sd_gamma0 )
  }
  for( i in 1:nsite){
    epsilon[i] ~ dnorm( 0, sd = sd_epsilon )
  }
  for( i in 1:nind ) {
    y[i] ~ dbern( kappa[i] )
    logit( kappa[i] ) <- gamma0[sp[i]] + epsilon[site[i]] + zeta[method[i]]
  }
})

inits <- function(){
  list(
    mu_gamma0 = rnorm(1, 0, 0.1),
    sd_gamma0 = runif(1, 0, 1),
    gamma0 = rnorm(constants$nsp, 0, 1),
    sd_epsilon = rexp(1),
    epsilon = rnorm(constants$nsite, 0, 1),
    zeta = rnorm(constants$nmethod, 0, 0.5)
  )
}

params <- c("mu_gamma0", "sd_gamma0", "gamma0",
            "sd_epsilon", "epsilon",  "zeta")

nc <- 3
nb <- 20000
ni <- nb + 10000
nt <- 10

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

MCMCvis::MCMCsummary(out, params = "zeta") |>
  tibble::as_tibble(rownames = "params") |> 
  dplyr::mutate(method = parse_number(params)) |> 
  dplyr::select(method, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
  dplyr::full_join( final |>
                      dplyr::select(method, type) |>
                      dplyr::distinct()
  ) |>
  dplyr::mutate(type = factor(type, levels = c("blood", "ear", "both"))) |> 
  ggplot2::ggplot(aes(x = mean, y = type)) + 
  ggplot2::geom_vline(xintercept = 0,
                      color = MetBrewer::MetPalettes$Hiroshige[[1]][2],
                      linetype = "dashed") +
  ggplot2::geom_errorbar(aes(xmin = l95, xmax = u95), width = 0,
                         color = MetBrewer::MetPalettes$Hiroshige[[1]][9],
                         linewidth = 1) +
  ggplot2::geom_point(size = 3, color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Estimated effect", 
                y = "Sample type")

setwd(here::here("figures"))
ggsave(
  filename = "figure_s03.png", 
  width = 3.5, 
  height = 3, 
  units = "in", 
  dpi = 600
)

setwd(here::here("results"))
save( code, data, constants, out, 
      file = paste0("rodent_pathogen_method_", Sys.Date(), ".RData"))