# Quick script - ran in ~2 minutes on basic desktop
# Check influence of sample type (blood only, ear tissue only, or both)
library(tidyverse)
library(parallel)
library(nimble)
library(here)
library(MetBrewer)
library(gtable)
library(lemon)
library(MCMCvis)

setwd(here::here("data"))

#dist <- readr::read_csv("site_distance_area.csv")
dist <- readr::read_csv("nplots.csv")

final <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv") %>% 
  dplyr::group_by(scientificName) %>% 
  dplyr::mutate(sp_disease = cur_group_id()) |>
  dplyr::ungroup() |> 
  dplyr::left_join(dist) |> 
  dplyr::filter(nplots > 1) |> 
  dplyr::group_by(siteID) |> 
  dplyr::mutate(site = cur_group_id()) |> 
  dplyr::ungroup()

run_model <- function( df, metric){
  
  positive <- df$positive
  metric <- as.numeric(scale(df[[metric]]))
  dist <- as.numeric(scale(df$nplots))
  
  data <- list(
    y = positive, 
    metric = metric, 
    dist = dist
  )
  
  constants <- list(
    nsp = length(unique(df$sp_disease)), 
    nsite = length(unique(df$site)), 
    nind = nrow(df),
    site = df$site,
    sp = df$sp_disease)
  
  code <- nimbleCode({
    mu_gamma0 ~ dnorm(0, sd = 1)
    sd_gamma0 ~ dexp(1)
    sd_epsilon ~ dexp(1)
    gamma1 ~ dnorm(0, sd = 1)
    gamma2 ~ dnorm(0, sd = 1)
    gamma3 ~ dnorm(0, sd = 1)
    for(i in 1:nsp){
      gamma0[i] ~ dnorm( mu_gamma0, sd = sd_gamma0 )
    }
    for( i in 1:nsite){
      epsilon[i] ~ dnorm( 0, sd = sd_epsilon )
    }
    for( i in 1:nind ) {
      y[i] ~ dbern( kappa[i] )
      logit( kappa[i] ) <- gamma0[sp[i]] + gamma1*metric[i] + gamma2*dist[i] + gamma3*metric[i]*dist[i] + epsilon[site[i]]
    }
  })
  
  inits <- function(){
    list(
      mu_gamma0 = rnorm(1, 0, 0.1),
      sd_gamma0 = runif(1, 0, 1),
      gamma0 = rnorm(constants$nsp, 0, 1),
      sd_epsilon = rexp(1),
      epsilon = rnorm(constants$nsite, 0, 1),
      gamma1 = rnorm(1, 0, 0.1), 
      gamma2 = rnorm(1, 0, 0.1), 
      gamma3 = rnorm(1, 0, 0.1)
    )
  }
  
  params <- c("mu_gamma0", "sd_gamma0", "gamma0", "gamma1", "gamma2", "gamma3",
              "sd_epsilon", "epsilon")
  
  nc <- 3
  nb <- 20000
  ni <- nb + 10000
  nt <- 10
  
  out <- nimble::nimbleMCMC(
    code = code, 
    constants = constants, 
    data = data, 
    inits = inits(),
    monitors = params, 
    thin = nt, 
    niter = ni, 
    nburnin = nb, 
    nchains = nc)
  
  settings <- list(ni = ni, nb = nb, nc = nc, nt = nt)
  
  return(
    list(
      code, 
      data, 
      constants, 
      settings, 
      out
    )
  )
  
}

sr <- run_model( final, "sr_site_mean" )
pd <- run_model( final, "pd_site_mean" )
totN <- run_model( final, "n_site_mean" )
shan <- run_model( final, "shan_site_mean" )
faith <- run_model( final, "faith_site_mean" )

comp_marginal <- function( df, output, metric){
  
  pdat <- tidyr::expand_grid(
    metric = seq(from = min(output[[2]]$metric),
                 to = max(output[[2]]$metric), 
                 by = 0.25),
    dist = c(min(output[[2]]$dist), 
             0, 
             max(output[[2]]$dist)))
  
  metric_sc <- scale(df[[metric]])
  dist_sc <- scale(df$nplots)
  
  pred <- MCMCvis::MCMCpstr( output[[5]], params = c("mu_gamma0"), type = "chains")[[1]] |> 
    tibble::as_tibble() |>
    tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma0") |>
    dplyr::full_join(
      MCMCvis::MCMCpstr( output[[5]], params = c("gamma1"), type = "chains")[[1]] |>
        tibble::as_tibble() |>
        tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma1")) |>
    dplyr::full_join(
      MCMCvis::MCMCpstr( output[[5]], params = c("gamma2"), type = "chains")[[1]] |>
        tibble::as_tibble() |>
        tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma2")) |>
    dplyr::full_join(
      MCMCvis::MCMCpstr( output[[5]], params = c("gamma3"), type = "chains")[[1]] |>
        tibble::as_tibble() |>
        tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma3")) |>
    dplyr::cross_join( pdat ) |>
    dplyr::mutate( p = ilogit( gamma0 + gamma1*metric + gamma2*dist + gamma3*metric*dist)) |>
    dplyr::group_by(metric, dist) |>
    dplyr::summarise( mean = mean(p),
                      l95 = quantile(p, 0.025),
                      u95 = quantile(p, 0.975)) |>
    dplyr::mutate(metric_raw = metric*attr(metric_sc, "scaled:scale") + attr(metric_sc, "scaled:center"),
                  dist_raw = round(dist*attr(dist_sc, "scaled:scale") +
                                     attr(dist_sc, "scaled:center"), 2)) |> 
    tibble::add_column(metric_name = metric)
  
  return(pred)
}

sr_marginal <- comp_marginal(df = final, output = sr, metric = "sr_site_mean")
pd_marginal <- comp_marginal(df = final, output = pd, metric = "pd_site_mean")
n_marginal <- comp_marginal(df = final, output = totN, metric = "n_site_mean")
shan_marginal <- comp_marginal(df = final, output = shan, metric = "shan_site_mean")
faith_marginal <- comp_marginal(df = final, output = faith, metric = "faith_site_mean")

key <- tibble::tribble(
  ~metric_name, ~label, 
  "sr_site_mean", "Species richness", 
  "pd_site_mean", "Peromyscus dominance", 
  "n_site_mean", "Total abundance", 
  "shan_site_mean", "Shannon diversity",
  "faith_site_mean", "Phylogenetic diversity")

plot_df <- dplyr::full_join(sr_marginal, pd_marginal) |> 
  dplyr::full_join(n_marginal) |> 
  dplyr::full_join(shan_marginal) |> 
  dplyr::full_join(faith_marginal) |> 
  dplyr::full_join(key) |> 
  dplyr::mutate(label  = factor(label, levels = c(
    "Species richness", 
    "Peromyscus dominance", 
    "Total abundance", 
    "Shannon diversity", 
    "Phylogenetic diversity")))

p <- ggplot2::ggplot( plot_df, aes(x = metric_raw, 
                                   y = mean, 
                                   color = factor(dist_raw), 
                                   fill = factor(dist_raw))) +
  ggplot2::facet_wrap(~label, scales = "free_x") +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) +
  ggplot2::geom_line(size = 1.5) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Biodiversity metric value",
                y = "Probability of infection") +
  ggplot2::scale_color_manual(
    "Number of plots",
    values = MetBrewer::MetPalettes$Archambault[[1]][c(1, 3, 6)]) +
  ggplot2::scale_fill_manual(
    "Number of plots",
    values = MetBrewer::MetPalettes$Archambault[[1]][c(1, 3, 6)]) +
  ggplot2::theme( plot.background = element_rect(color = NA, fill = "white"), 
                  panel.background = element_rect(color = NA, fill = "white"), 
                  axis.title = element_text(size = 10, color = "black"), 
                  axis.text = element_text(size = 9, color = "black"), 
                  strip.text = element_text(size = 10, color = "black"), 
                  legend.title = element_text(size = 10, color = "black"), 
                  legend.text = element_text(size = 9, color = "black"))

shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

p_shift <- shift_legend2(p)

get_summary <- function( output, metric){
  sum <- MCMCvis::MCMCsummary( output[[5]],
                               params = c("gamma1", "gamma2", "gamma3"),
                               probs = c(0.025, 0.160, 0.840, 0.975)) |> 
    tibble::as_tibble(rownames = "param") |> 
    tibble::add_column(metric = metric) |> 
    dplyr::select(metric, param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`)
  
  return(sum)
}

sr_sum <- get_summary(sr, "sr_site_mean")
pd_sum <- get_summary(pd, "pd_site_mean")
n_sum <- get_summary(totN, "n_site_mean")
shan_sum <- get_summary(shan, "shan_site_mean")
faith_sum <- get_summary(faith, "faith_site_mean")

allsum <- dplyr::full_join(sr_sum, pd_sum) |> 
  dplyr::full_join(n_sum) |> 
  dplyr::full_join(shan_sum) |> 
  dplyr::full_join(faith_sum) |> 
  dplyr::full_join(key |>
                     dplyr::rename(metric = metric_name)) 

setwd(here::here("results"))
save(
  sr, 
  pd, 
  totN, 
  shan, 
  faith, 
  plot_df, 
  allsum,
  file = "nplots_metacommunity_results.RData")

setwd(here::here("figures"))
ggsave("figure_s07.png",
       plot = p_shift, 
       width = 6.5, height = 4.5, units = "in", dpi = 600)

allsum |> 
  filter(param == "gamma3") |> 
  dplyr::mutate(label  = factor(label, levels = c(
    "Species richness", 
    "Peromyscus dominance", 
    "Total abundance", 
    "Shannon diversity", 
    "Phylogenetic diversity"))) |> 
  
  ggplot(aes(x = mean, y = label)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1.5) + 
  geom_point(size = 3) +
  labs(x = "N_plots x Biodiversity interaction") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = 10), 
        axis.title = element_text(color = "black", size = 11),
        plot.background = element_rect(color = NA, fill = "white"),
        panel.background = element_rect(color = NA, fill = "white"))

ggsave("figure_s06.png",
       width = 4, height = 3, units = "in", dpi = 600)
