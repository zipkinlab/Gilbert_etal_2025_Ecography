library(brms)
library(tidyverse)
library(here)
library(here)

setwd(here::here("results"))

load("rodent_pathogen_sr_plot_2024-04-18.RData")
sr_plot <- out
load("rodent_pathogen_sr_site_2024-04-18.RData")
sr_site <- out
load("rodent_pathogen_pd_plot_2024-04-19.RData")
pd_plot <- out
load("rodent_pathogen_pd_site_2024-04-19.RData")
pd_site <- out
load("rodent_pathogen_shan_plot_2024-04-18.RData")
shan_plot <- out
load("rodent_pathogen_shan_site_2024-04-18.RData")
shan_site <- out
load("rodent_pathogen_n_plot_2024-04-18.RData")
n_plot <- out
load("rodent_pathogen_n_site_2024-04-18.RData")
n_site <- out
load("rodent_pathogen_faith_plot_2024-08-21.RData")
faith_plot <- out
load("rodent_pathogen_faith_site_2024-08-21.RData")
faith_site <- out

sr_site_phylo <- readRDS("sr_site_phylo.rds")
sr_plot_phylo <- readRDS("sr_plot_phylo.rds")
pd_site_phylo <- readRDS("pd_site_phylo.rds")
pd_plot_phylo <- readRDS("pd_plot_phylo.rds")
n_site_phylo <- readRDS("n_site_phylo.rds")
n_plot_phylo <- readRDS("n_plot_phylo.rds")
shan_site_phylo <- readRDS("shan_site_phylo.rds")
shan_plot_phylo <- readRDS("shan_plot_phylo.rds")
faith_site_phylo <- readRDS("faith_site_phylo.rds")
faith_plot_phylo <- readRDS("faith_plot_phylo.rds")

phylo_ests <- brms::fixef( sr_site_phylo, 
       probs = c(0.025, 0.160, 0.840, 0.975)) |> 
  tibble::as_tibble(rownames = "param") |> 
  dplyr::filter(param == "sr_site_mean") |> 
  tibble::add_column(metric = "Species richness", 
                     scale = "Metacommunity") |> 
  dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5) |> 
  dplyr::full_join(
    brms::fixef( sr_plot_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "sr_plot_mean") |> 
      tibble::add_column(metric = "Species richness", 
                         scale = "Local community") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( pd_site_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "pd_site_mean") |> 
      tibble::add_column(metric = "Peromyscus dominance", 
                         scale = "Metacommunity") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( pd_plot_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "pd_plot_mean") |> 
      tibble::add_column(metric = "Peromyscus dominance", 
                         scale = "Local community") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( n_plot_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "n_plot_mean") |> 
      tibble::add_column(metric = "Total abundance", 
                         scale = "Local community") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( n_site_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "n_site_mean") |> 
      tibble::add_column(metric = "Total abundance", 
                         scale = "Metacommunity") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( shan_plot_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "shan_plot_mean") |> 
      tibble::add_column(metric = "Shannon diversity", 
                         scale = "Local community") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( shan_site_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "shan_site_mean") |> 
      tibble::add_column(metric = "Shannon diversity", 
                         scale = "Metacommunity") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( faith_plot_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "faith_plot_mean") |> 
      tibble::add_column(metric = "Phylogenetic diversity", 
                         scale = "Local community") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::full_join(
    brms::fixef( faith_site_phylo, 
                 probs = c(0.025, 0.160, 0.840, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::filter(param == "faith_site_mean") |> 
      tibble::add_column(metric = "Phylogenetic diversity", 
                         scale = "Metacommunity") |> 
      dplyr::select(metric, scale, mean = Estimate, Q2.5:Q97.5)
  ) |> 
  dplyr::mutate(metric_styled = ifelse(grepl("Pero", metric), "<i>Peromyscus</i> dominance", metric)) |> 
  dplyr::mutate(metric_styled = factor(metric_styled, levels = c(
    "Phylogenetic diversity",
    "Shannon diversity",
    "Total abundance", 
    "<i>Peromyscus</i> dominance", 
    "Species richness"))) |> 
  dplyr::rename(l95 = Q2.5, 
                l68 = Q16, 
                u68 = Q84, 
                u95 = Q97.5) |> 
  tibble::add_column(method = "Phylogenetic")

non_phylo_ests <- MCMCvis::MCMCsummary( pd_site, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
  tibble::as_tibble(rownames = "param") |> 
  dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
  tibble::add_column(metric = "PD", scale = "site") |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( pd_plot, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "PD", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( sr_site, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "SR", scale = "site")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( sr_plot, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "SR", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( n_plot, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Total abundance", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( n_site, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Total abundance", scale = "site")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( shan_plot, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Shannon diversity", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( shan_site, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Shannon diversity", scale = "site"))  |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( faith_site, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Phylogenetic diversity", scale = "site")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( faith_plot, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, l68 = `16%`, u68 = `84%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Phylogenetic diversity", scale = "plot")) |> 
  dplyr::mutate( metric = ifelse(metric == "SR", "Species richness",
                                 ifelse(metric == "PD", "Peromyscus dominance", metric))) |> 
  dplyr::mutate(metric_styled = ifelse(grepl("Pero", metric), "<i>Peromyscus</i> dominance", metric)) |> 
  dplyr::mutate(metric_styled = factor(metric_styled, levels = c(
    "Phylogenetic diversity",
    "Shannon diversity",
    "Total abundance", 
    "<i>Peromyscus</i> dominance", 
    "Species richness"))) |> 
  dplyr::select(-param) |> 
  dplyr::mutate(scale = ifelse(scale == "site", "Metacommunity", 
                               "Local community")) |> 
  tibble::add_column(method = "No phylogeny")

plot.d <- dplyr::full_join( phylo_ests, 
                  non_phylo_ests)
# ggplot2::geom_errorbar( aes(xmin = l68, xmax = u68, color = scale),
#                         width = 0, size = 1.75, position = position_dodge(width = 0.4)) + 
#   ggplot2::geom_errorbar( aes(xmin = l95, xmax = u95, color = scale),
#                           width = 0, size = 0.75, position = position_dodge(width = 0.4)) + 
#   ggplot2::geom_point(aes(color = scale), size = 3.5, position = position_dodge(width = 0.4)) + 


ggplot( plot.d, aes(x = mean, y = metric_styled, 
                    color = scale, shape = method, linetype = method)) +
  ggplot2::geom_vline(xintercept = 0, color = "gray30", linetype = "dashed") +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0,
                position = position_dodge(width = 0.4),
                size = 1) +
  # geom_errorbar(aes(xmin = l68, xmax = u68), width = 0,
  #               position = position_dodge(width = 0.4),
  #               size = 1.25) +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  ggplot2::labs( x = "Biodiversity effect (logit scale)", 
                 color = "Scale") +
  ggplot2::theme_classic() + 
  ggplot2::scale_color_manual(values = MetBrewer::MetPalettes$Hiroshige[[1]][c(1, 8)]) +
  ggplot2::guides( color = guide_legend(reverse = TRUE)) +
  ggplot2::theme(axis.title.y = element_blank(),
                 axis.title = element_text(color = "black", size = 11), 
                 axis.text.y = ggtext::element_markdown(color = "black", size = 10),
                 axis.text.x = element_text(color = "black", size = 10), 
                 legend.text = element_text(color = "black", size = 8), 
                 legend.title = element_blank(),
                 axis.line = element_line(color = "black", size = 0.1), 
                 axis.ticks.x = element_line(color = "black", size = 0.1),
                 axis.ticks.y = element_blank(),
                 legend.position = "bottom",
                 legend.box = "vertical",
                 legend.margin = margin(0, 0, 0, 0),
                 legend.box.margin = margin(0, 0, 0, 0)) +
  ggplot2::guides(color = guide_legend(ncol = 2, reverse = T))

setwd(here::here("figures"))
ggsave(
  filename = "figure_s02.png", 
  width = 5, 
  height = 6, 
  units = "in", 
  dpi = 600
)