library(ape)
library(here)
library(brms)
library(tidyverse)

setwd(here::here("data"))

d <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv") |> 
  dplyr::mutate(tipLabel = stringr::str_replace_all(scientificName, 
                                                    " ", "_")) |> 
  dplyr::mutate( faith_plot_mean = as.numeric(scale(faith_plot_mean)), 
                 faith_site_mean = as.numeric(scale(faith_site_mean)))

dtips <- d |>
  dplyr::select(tipLabel) |>
  dplyr::distinct() |>
  dplyr::arrange() |> 
  dplyr::pull()

tree <- ape::read.nexus( "output.nex" )[[1]]

dropsy <- tree$tip[!tree$tip.label %in% dtips]

pruned_tree <- ape::drop.tip(tree, dropsy)

A <- ape::vcv.phylo( pruned_tree )

options(mc.cores = 4)

setwd(here::here("results"))
m1 <- brms::brm(positive ~ 1 + sr_plot_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "sr_plot_phylo")

m2 <- brms::brm(positive ~ 1 + sr_site_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "sr_site_phylo")

m3 <- brms::brm(positive ~ 1 + pd_plot_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "pd_plot_phylo")

m4 <- brms::brm(positive ~ 1 + pd_site_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "pd_site_phylo")

m5 <- brms::brm(positive ~ 1 + n_plot_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "n_plot_phylo")

m6 <- brms::brm(positive ~ 1 + n_site_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "n_site_phylo")

m7 <- brms::brm(positive ~ 1 + shan_plot_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "shan_plot_phylo")

m8 <- brms::brm(positive ~ 1 + shan_site_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "shan_site_phylo")

m9 <- brms::brm(positive ~ 1 + faith_plot_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                data = d, 
                family = bernoulli(link = "logit"),
                data2 = list(A = A), 
                prior = c(
                  prior(normal(0, 1), "b"), 
                  prior(normal(0, 2), "Intercept"),
                  prior(exponential(1), "sd")),
                chains = 4,
                iter = 4000,
                file = "faith_plot_phylo")

m10 <- brms::brm(positive ~ 1 + faith_site_mean + (1|gr(tipLabel, cov = A)) + (1|site),
                 data = d, 
                 family = bernoulli(link = "logit"),
                 data2 = list(A = A), 
                 prior = c(
                   prior(normal(0, 1), "b"), 
                   prior(normal(0, 2), "Intercept"),
                   prior(exponential(1), "sd")),
                 chains = 4,
                 iter = 4000,
                 file = "faith_site_phylo")