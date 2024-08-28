library(here)
library(tidyverse)
library(MCMCvis)
library(nimble)

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

setwd(here::here("data"))
d <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv")
load("neon_cr_data_2024-03-29.RData")
load("neon_mammal_box_trapping_v01.RData")

MCMCvis::MCMCsummary( pd_site, params = c("gamma1"), probs = c(0.025, 0.16, 0.84, 0.975)) |> 
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
  dplyr::mutate( scale = ifelse(scale == "plot", "Local community", "Metacommunity")) |> 
  dplyr::mutate(scale = factor(scale, levels = c("Local community", "Metacommunity"))) |>
  ggplot2::ggplot( aes( x = mean, y = metric_styled)) + 
  ggplot2::geom_vline(xintercept = 0, color = "gray30", linetype = "dashed") +
  ggplot2::geom_errorbar( aes(xmin = l68, xmax = u68, color = scale),
                          width = 0, size = 1.75, position = position_dodge(width = 0.4)) + 
  ggplot2::geom_errorbar( aes(xmin = l95, xmax = u95, color = scale),
                          width = 0, size = 0.75, position = position_dodge(width = 0.4)) + 
  ggplot2::geom_point(aes(color = scale), size = 3.5, position = position_dodge(width = 0.4)) + 
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
                 legend.position = "bottom") +
  ggplot2::guides(color = guide_legend(ncol = 2, reverse = T)) 

setwd(here::here("figures"))
ggsave(
  "figure_03.png", 
  width = 4.85, 
  height = 4.25, 
  units = "in", 
  dpi = 600
)  

MCMCpstr(sr_plot, params = c("gamma1"), type = "chains")[[1]] |>
  as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
  summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
  add_column( metric = "SR", 
              scale = "plot") |> 
  dplyr::select(metric, scale, mean, l95, u95, p_positive) |> 
  
  full_join(
    MCMCpstr(sr_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "SR", 
                  scale = "site") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  
  full_join(
    MCMCpstr(pd_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "PD", 
                  scale = "site") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(pd_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "PD", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(n_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "N", 
                  scale = "site") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(n_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "N", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(shan_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "Shannon", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(shan_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "Shannon", 
                  scale = "site") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(faith_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "Faith", 
                  scale = "site") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  full_join(
    MCMCpstr(faith_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value)), mean = mean(value), l95 = quantile(value, c(0.025)), u95 = quantile(value, c(0.975))) |> 
      add_column( metric = "Faith", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, mean, l95, u95, p_positive)
  ) |> 
  mutate(p_positive = round(p_positive, 2),
         p_negative = round(1 - p_positive, 2))

site_key <-
  final |>
  dplyr::select(site, siteID) |>
  dplyr::distinct() |> 
  dplyr::arrange(siteID) |>
  tibble::add_column(
    full_name = c("Bartlett Experimental Forest",
                  "Blandy Experimental Farm",
                  "Dead Lake",
                  "Disney Wilderness Preserve",
                  "Great Smoky Mountains National Park",
                  "Harvard Forest",
                  "The Jones Center",
                  "Konza Prairie Agroecosystem",
                  "Konza Prairie Biological Station",
                  "Mountain Lake Biological Station",
                  "Northern Great Plains Research Laboratory",
                  "Marvin Klemme Range Research Station",
                  "Oak Ridge National Laboratory",
                  "Ordway-Swisher Biological Station",
                  "Smithsonian Conservation Biology Institute",
                  "Smithsonian Environmental Research Center",
                  "Steigerwald-Chequamegon",
                  "Talladega National Forest",
                  "Treehaven",
                  "Unversity of Kansas Field Station",
                  "University of Notre Dame Environmental Research Center",
                  "Chase Lake National Wildlife Refuge")) |> 
  dplyr::left_join(
    boxtrap$mam_pertrapnight |> 
      dplyr::select(siteID, plotID, decimalLatitude, decimalLongitude) |> 
      dplyr::distinct() |> 
      dplyr::group_by(siteID) |> 
      dplyr::summarise(x = mean(decimalLongitude), 
                       y = mean(decimalLatitude)))

sr_plot_d <- d |>
  dplyr::select(site, sr_plot_mean, sr_plot_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(sr_plot_mean == min(sr_plot_mean) | sr_plot_mean == max(sr_plot_mean)) |> 
  dplyr::mutate( sr_l95 = sr_plot_mean - 2*sr_plot_sd,
                 sr_u95 = sr_plot_mean + 2*sr_plot_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(sr_plot_mean == min(sr_plot_mean), "min", "max")) |>
  dplyr::select(-sr_plot_mean, -sr_plot_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( sr_plot = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.25))) |>
  tidyr::unnest(cols = sr_plot) |>
  dplyr::ungroup() |>
  dplyr::mutate(sr_plot_sc = as.numeric(scale(sr_plot))) |>
  dplyr::select(site, sr_plot = sr_plot_sc, sr_plot_unscaled = sr_plot)

sr_plot_com <-
  sr_plot_d |>
  dplyr::select(sr = sr_plot_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Species richness", 
                     scale = "Local community")

sr_site_d <- d |>
  dplyr::select(site, sr_site_mean, sr_site_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(sr_site_mean == min(sr_site_mean) | sr_site_mean == max(sr_site_mean)) |> 
  dplyr::mutate( sr_l95 = sr_site_mean - 2*sr_site_sd,
                 sr_u95 = sr_site_mean + 2*sr_site_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(sr_site_mean == min(sr_site_mean), "min", "max")) |>
  dplyr::select(-sr_site_mean, -sr_site_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( sr_site = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.25))) |>
  tidyr::unnest(cols = sr_site) |>
  dplyr::ungroup() |>
  dplyr::mutate(sr_site_sc = as.numeric(scale(sr_site))) |>
  dplyr::select(site, sr_site = sr_site_sc, sr_site_unscaled = sr_site)

sr_site_com <- sr_site_d |>
  dplyr::select(sr = sr_site_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_site = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_site, x_unscaled = sr) |> 
  tibble::add_column(metric = "Species richness", 
                     scale = "Metacommunity")

pd_plot_d <- d |>
  dplyr::select(site, pd_plot_mean, pd_plot_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(pd_plot_mean == min(pd_plot_mean) | pd_plot_mean == max(pd_plot_mean)) |> 
  dplyr::mutate( sr_l95 = pd_plot_mean - 2*pd_plot_sd,
                 sr_u95 = pd_plot_mean + 2*pd_plot_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(pd_plot_mean == min(pd_plot_mean), "min", "max")) |>
  dplyr::select(-pd_plot_mean, -pd_plot_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( pd_plot = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = pd_plot) |>
  dplyr::ungroup() |>
  dplyr::mutate(pd_plot_sc = as.numeric(scale(pd_plot))) |>
  dplyr::select(site, pd_plot = pd_plot_sc, pd_plot_unscaled = pd_plot)

pd_plot_com <- tibble(pd_plot = seq(from = 0, to = 1, by = 0.1)) |>
  dplyr::mutate(x = as.numeric(scale(pd_plot))) |> 
  dplyr::rename(x_unscaled = pd_plot) |> 
  tibble::add_column(metric = "Peromyscus dominance") |> 
  tibble::add_column(scale = "Local community") 

pd_site_d <- d |>
  dplyr::select(site, pd_site_mean, pd_site_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(pd_site_mean == min(pd_site_mean) | pd_site_mean == max(pd_site_mean)) |> 
  dplyr::mutate( sr_l95 = pd_site_mean - 2*pd_site_sd,
                 sr_u95 = pd_site_mean + 2*pd_site_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(pd_site_mean == min(pd_site_mean), "min", "max")) |>
  dplyr::select(-pd_site_mean, -pd_site_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( pd_site = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = pd_site) |>
  dplyr::ungroup() |>
  dplyr::mutate(pd_site_sc = as.numeric(scale(pd_site))) |>
  dplyr::select(site, pd_site = pd_site_sc, pd_site_unscaled = pd_site)

pd_site_com <- tibble(pd_site = seq(from = 0, to = 1, by = 0.1)) |>
  dplyr::mutate(x = as.numeric(scale(pd_site))) |> 
  dplyr::rename(x_unscaled = pd_site) |> 
  tibble::add_column(metric = "Peromyscus dominance") |> 
  tibble::add_column(scale = "Metacommunity") 

n_plot_d <- d |>
  dplyr::select(site, n_plot_mean, n_plot_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(n_plot_mean == min(n_plot_mean) | n_plot_mean == max(n_plot_mean)) |> 
  dplyr::mutate( sr_l95 = n_plot_mean - 2*n_plot_sd,
                 sr_u95 = n_plot_mean + 2*n_plot_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(n_plot_mean == min(n_plot_mean), "min", "max")) |>
  dplyr::select(-n_plot_mean, -n_plot_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( n_plot = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = n_plot) |>
  dplyr::ungroup() |>
  dplyr::mutate(n_plot_sc = as.numeric(scale(n_plot)),
                nraw = expm1(n_plot)) |>
  dplyr::select(site, n_plot = n_plot_sc, n_plot_unscaled = n_plot)

n_plot_com <- n_plot_d |>
  dplyr::select(sr = n_plot_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Log( Total abundance )", 
                     scale = "Local community")

n_site_d <- d |>
  dplyr::select(site, n_site_mean, n_site_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(n_site_mean == min(n_site_mean) | n_site_mean == max(n_site_mean)) |> 
  dplyr::mutate( sr_l95 = n_site_mean - 2*n_site_sd,
                 sr_u95 = n_site_mean + 2*n_site_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(n_site_mean == min(n_site_mean), "min", "max")) |>
  dplyr::select(-n_site_mean, -n_site_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( n_site = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = n_site) |>
  dplyr::ungroup() |>
  dplyr::mutate(n_site_sc = as.numeric(scale(n_site)),
                nraw = expm1(n_site)) |>
  dplyr::select(site, n_site = n_site_sc, n_site_unscaled = n_site)

n_site_com <- n_site_d |>
  dplyr::select(sr = n_site_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Log( Total abundance )", 
                     scale = "Metacommunity")

shan_plot_d <- d |>
  dplyr::select(site, shan_plot_mean, shan_plot_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(shan_plot_mean == min(shan_plot_mean) | shan_plot_mean == max(shan_plot_mean)) |> 
  dplyr::mutate( sr_l95 = shan_plot_mean - 2*shan_plot_sd,
                 sr_u95 = shan_plot_mean + 2*shan_plot_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(shan_plot_mean == min(shan_plot_mean), "min", "max")) |>
  dplyr::select(-shan_plot_mean, -shan_plot_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( shan_plot = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = shan_plot) |>
  dplyr::ungroup() |>
  dplyr::mutate(shan_plot_sc = as.numeric(scale(shan_plot))) |>
  dplyr::select(site, shan_plot = shan_plot_sc, shan_plot_unscaled = shan_plot)

shan_plot_com <- shan_plot_d |>
  dplyr::select(sr = shan_plot_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Shannon diversity", 
                     scale = "Local community")

shan_site_d <- d |>
  dplyr::select(site, shan_site_mean, shan_site_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(shan_site_mean == min(shan_site_mean) | shan_site_mean == max(shan_site_mean)) |> 
  dplyr::mutate( sr_l95 = shan_site_mean - 2*shan_site_sd,
                 sr_u95 = shan_site_mean + 2*shan_site_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(shan_site_mean == min(shan_site_mean), "min", "max")) |>
  dplyr::select(-shan_site_mean, -shan_site_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( shan_site = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = shan_site) |>
  dplyr::ungroup() |>
  dplyr::mutate(shan_site_sc = as.numeric(scale(shan_site))) |>
  dplyr::select(site, shan_site = shan_site_sc, shan_site_unscaled = shan_site)

shan_site_com <- shan_site_d |>
  dplyr::select(sr = shan_site_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Shannon diversity", 
                     scale = "Metacommunity")

faith_plot_d <- d |>
  dplyr::select(site, faith_plot_mean, faith_plot_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(faith_plot_mean == min(faith_plot_mean) | faith_plot_mean == max(faith_plot_mean)) |> 
  dplyr::mutate( sr_l95 = faith_plot_mean - 2*faith_plot_sd,
                 sr_u95 = faith_plot_mean + 2*faith_plot_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(faith_plot_mean == min(faith_plot_mean), "min", "max")) |>
  dplyr::select(-faith_plot_mean, -faith_plot_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( faith_plot = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = faith_plot) |>
  dplyr::ungroup() |>
  dplyr::mutate(faith_plot_sc = as.numeric(scale(faith_plot))) |>
  dplyr::select(site, faith_plot = faith_plot_sc, faith_plot_unscaled = faith_plot)

faith_plot_com <- faith_plot_d |>
  dplyr::select(sr = faith_plot_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Phylogenetic diversity", 
                     scale = "Local community")

faith_site_d <- d |>
  dplyr::select(site, faith_site_mean, faith_site_sd) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(faith_site_mean == min(faith_site_mean) | faith_site_mean == max(faith_site_mean)) |> 
  dplyr::mutate( sr_l95 = faith_site_mean - 2*faith_site_sd,
                 sr_u95 = faith_site_mean + 2*faith_site_sd) |>
  dplyr::mutate(sr_l95 = ifelse(sr_l95 < 0, 0, sr_l95)) |>
  dplyr::mutate( type = ifelse(faith_site_mean == min(faith_site_mean), "min", "max")) |>
  dplyr::select(-faith_site_mean, -faith_site_sd) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_l95:sr_u95) |>
  dplyr::mutate(sr_u95_max = ifelse(is.na(sr_u95_max), sr_u95_min, sr_u95_max)) |> 
  dplyr::select(-sr_l95_max, -sr_u95_min) |>  
  dplyr::mutate( faith_site = list(seq(from = sr_l95_min, to = sr_u95_max, by = 0.1))) |>
  tidyr::unnest(cols = faith_site) |>
  dplyr::ungroup() |>
  dplyr::mutate(faith_site_sc = as.numeric(scale(faith_site))) |>
  dplyr::select(site, faith_site = faith_site_sc, faith_site_unscaled = faith_site)

faith_site_com <- faith_site_d |>
  dplyr::select(sr = faith_site_unscaled) |> 
  dplyr::filter( sr == min(sr) | sr == max(sr)) |> 
  dplyr::distinct() |> 
  tibble::add_column(type = c("min", "max")) |> 
  tidyr::pivot_wider(names_from = type, values_from = sr) |> 
  dplyr::mutate( sr = list(seq(from = min, to = max, by = 0.1))) |> 
  tidyr::unnest(sr) |> 
  dplyr::mutate(sr_plot = as.numeric(scale(sr))) |> 
  dplyr::select(x = sr_plot, x_unscaled = sr) |> 
  tibble::add_column(metric = "Phylogenetic diversity", 
                     scale = "Metacommunity")

sr_com_marginal <- MCMCvis::MCMCpstr( sr_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::select(iter, mu_gamma0, gamma1) |>
  dplyr::cross_join( sr_plot_com) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
  dplyr::group_by(x, x_unscaled, metric, scale) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(
    
    MCMCvis::MCMCpstr( sr_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( sr_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::select(iter, mu_gamma0, gamma1) |>
      dplyr::cross_join( sr_site_com) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
      dplyr::group_by(x, x_unscaled, metric, scale) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))))

sr_marginal <-
MCMCvis::MCMCpstr( sr_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_plot, params = c("epsilon"), type = "chains")[[1]] |> 
      tibble::as_tibble(rownames = "site") |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
      dplyr::mutate(site = parse_number(site),
                    iter = parse_number(iter))) |> 
  dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
  dplyr::full_join( sr_plot_d) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * sr_plot + epsilon)) |> 
  dplyr::group_by(site, sr_plot, sr_plot_unscaled) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(site_key) |> 
  tibble::add_column(scale = "Local community") |> 
  dplyr::rename(xcoord = x, 
                ycoord = y,
                x = sr_plot,
                x_unscaled = sr_plot_unscaled) |> 
  
  dplyr::full_join(
    
    MCMCvis::MCMCpstr( sr_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( sr_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( sr_site, params = c("epsilon"), type = "chains")[[1]] |> 
          tibble::as_tibble(rownames = "site") |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
          dplyr::mutate(site = parse_number(site),
                        iter = parse_number(iter))) |> 
      dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
      dplyr::full_join( sr_site_d) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * sr_site + epsilon)) |> 
      dplyr::group_by(site, sr_site, sr_site_unscaled) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))) |> 
      dplyr::full_join(site_key) |> 
      
      tibble::add_column(scale = "Metacommunity") |> 
      dplyr::rename(xcoord = x, 
                    ycoord = y,
                    x = sr_site,
                    x_unscaled = sr_site_unscaled)) |> 
  tibble::add_column(metric = "Species richness")

pd_com_marginal <- MCMCvis::MCMCpstr( pd_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( pd_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::select(iter, mu_gamma0, gamma1) |>
  dplyr::cross_join( pd_plot_com) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
  dplyr::group_by(x, x_unscaled, metric, scale) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(
    
    MCMCvis::MCMCpstr( pd_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( pd_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::select(iter, mu_gamma0, gamma1) |>
      dplyr::cross_join( pd_site_com) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
      dplyr::group_by(x, x_unscaled, metric, scale) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))))

pd_marginal <-
  MCMCvis::MCMCpstr( pd_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( pd_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( pd_plot, params = c("epsilon"), type = "chains")[[1]] |> 
      tibble::as_tibble(rownames = "site") |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
      dplyr::mutate(site = parse_number(site),
                    iter = parse_number(iter))) |> 
  dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
  dplyr::full_join( pd_plot_d) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * pd_plot + epsilon)) |> 
  dplyr::group_by(site, pd_plot, pd_plot_unscaled) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(site_key) |> 
  tibble::add_column(scale = "Local community") |> 
  dplyr::rename(
    xcoord = x, 
    ycoord = y,
    x= pd_plot, x_unscaled = pd_plot_unscaled) |> 
  
  dplyr::full_join(
    MCMCvis::MCMCpstr( pd_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( pd_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( pd_site, params = c("epsilon"), type = "chains")[[1]] |> 
          tibble::as_tibble(rownames = "site") |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
          dplyr::mutate(site = parse_number(site),
                        iter = parse_number(iter))) |> 
      dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
      dplyr::full_join( pd_site_d) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * pd_site + epsilon)) |> 
      dplyr::group_by(site, pd_site, pd_site_unscaled) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))) |> 
      dplyr::full_join(site_key) |> 
      tibble::add_column(scale = "Metacommunity") |> 
      dplyr::rename(xcoord = x,
                    ycoord = y,
                    x = pd_site, x_unscaled = pd_site_unscaled)
  ) |> 
  tibble::add_column(metric = "Peromyscus dominance") |> 
  dplyr::filter(x_unscaled >=0 & x_unscaled <= 1)

n_com_marginal <- MCMCvis::MCMCpstr( n_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( n_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::select(iter, mu_gamma0, gamma1) |>
  dplyr::cross_join( n_plot_com) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
  dplyr::group_by(x, x_unscaled, metric, scale) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(
    
    MCMCvis::MCMCpstr( n_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( n_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::select(iter, mu_gamma0, gamma1) |>
      dplyr::cross_join( n_site_com) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
      dplyr::group_by(x, x_unscaled, metric, scale) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))))

n_marginal <-
  MCMCvis::MCMCpstr( n_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( n_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( n_plot, params = c("epsilon"), type = "chains")[[1]] |> 
      tibble::as_tibble(rownames = "site") |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
      dplyr::mutate(site = parse_number(site),
                    iter = parse_number(iter))) |> 
  dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
  dplyr::full_join( n_plot_d) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * n_plot + epsilon)) |> 
  dplyr::group_by(site, n_plot, n_plot_unscaled) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(site_key) |> 
  tibble::add_column(scale = "Local community") |> 
  dplyr::rename(xcoord = x, 
                ycoord = y,
                x= n_plot, x_unscaled = n_plot_unscaled) |> 
  
  dplyr::full_join(
    MCMCvis::MCMCpstr( n_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( n_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( n_site, params = c("epsilon"), type = "chains")[[1]] |> 
          tibble::as_tibble(rownames = "site") |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
          dplyr::mutate(site = parse_number(site),
                        iter = parse_number(iter))) |> 
      dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
      dplyr::full_join( n_site_d) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * n_site + epsilon)) |> 
      dplyr::group_by(site, n_site, n_site_unscaled) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))) |> 
      dplyr::full_join(site_key) |> 
      tibble::add_column(scale = "Metacommunity") |> 
      dplyr::rename(xcoord = x, 
                    ycoord = y,
                    x = n_site, x_unscaled = n_site_unscaled)
  ) |> 
  tibble::add_column(metric = "Log( Total abundance )")

shan_com_marginal <- MCMCvis::MCMCpstr( shan_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( shan_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::select(iter, mu_gamma0, gamma1) |>
  dplyr::cross_join( shan_plot_com) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
  dplyr::group_by(x, x_unscaled, metric, scale) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(
    
    MCMCvis::MCMCpstr( shan_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( shan_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::select(iter, mu_gamma0, gamma1) |>
      dplyr::cross_join( shan_site_com) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
      dplyr::group_by(x, x_unscaled, metric, scale) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))))

shan_marginal <-
  MCMCvis::MCMCpstr( shan_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( shan_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( shan_plot, params = c("epsilon"), type = "chains")[[1]] |> 
      tibble::as_tibble(rownames = "site") |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
      dplyr::mutate(site = parse_number(site),
                    iter = parse_number(iter))) |> 
  dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
  dplyr::full_join( shan_plot_d) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * shan_plot + epsilon)) |> 
  dplyr::group_by(site, shan_plot, shan_plot_unscaled) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(site_key) |> 
  tibble::add_column(scale = "Local community") |> 
  dplyr::rename(xcoord = x, 
                ycoord = y,
                x= shan_plot, x_unscaled = shan_plot_unscaled) |> 
  
  dplyr::full_join(
    MCMCvis::MCMCpstr( shan_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( shan_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( shan_site, params = c("epsilon"), type = "chains")[[1]] |> 
          tibble::as_tibble(rownames = "site") |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
          dplyr::mutate(site = parse_number(site),
                        iter = parse_number(iter))) |> 
      dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
      dplyr::full_join( shan_site_d) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * shan_site + epsilon)) |> 
      dplyr::group_by(site, shan_site, shan_site_unscaled) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))) |> 
      dplyr::full_join(site_key) |> 
      tibble::add_column(scale = "Metacommunity") |> 
      dplyr::rename(xcoord = x, 
                    ycoord = y,
                    x = shan_site, x_unscaled = shan_site_unscaled)
  ) |> 
  tibble::add_column(metric = "Shannon diversity")

faith_com_marginal <- MCMCvis::MCMCpstr( faith_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( faith_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::select(iter, mu_gamma0, gamma1) |>
  dplyr::cross_join( faith_plot_com) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
  dplyr::group_by(x, x_unscaled, metric, scale) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(
    
    MCMCvis::MCMCpstr( faith_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( faith_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::select(iter, mu_gamma0, gamma1) |>
      dplyr::cross_join( faith_site_com) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * x)) |> 
      dplyr::group_by(x, x_unscaled, metric, scale) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))))

faith_marginal <-
  MCMCvis::MCMCpstr( faith_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
  dplyr::mutate(iter = parse_number(iter)) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( faith_plot, params = c("gamma1"), type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
      dplyr::mutate(iter = parse_number(iter))) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( faith_plot, params = c("epsilon"), type = "chains")[[1]] |> 
      tibble::as_tibble(rownames = "site") |> 
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
      dplyr::mutate(site = parse_number(site),
                    iter = parse_number(iter))) |> 
  dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
  dplyr::full_join( faith_plot_d) |> 
  dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * faith_plot + epsilon)) |> 
  dplyr::group_by(site, faith_plot, faith_plot_unscaled) |> 
  dplyr::summarise(mean = mean(p), 
                   l95 = quantile(p, c(0.025)), 
                   u95 = quantile(p, c(0.975))) |> 
  dplyr::full_join(site_key) |> 
  tibble::add_column(scale = "Local community") |> 
  dplyr::rename(xcoord = x, 
                ycoord = y,
                x= faith_plot, x_unscaled = faith_plot_unscaled) |> 
  
  dplyr::full_join(
    MCMCvis::MCMCpstr( faith_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "mu_gamma0") |> 
      dplyr::mutate(iter = parse_number(iter)) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( faith_site, params = c("gamma1"), type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "gamma1") |> 
          dplyr::mutate(iter = parse_number(iter))) |> 
      dplyr::full_join(
        MCMCvis::MCMCpstr( faith_site, params = c("epsilon"), type = "chains")[[1]] |> 
          tibble::as_tibble(rownames = "site") |> 
          tidyr::pivot_longer(starts_with("V"), names_to = "iter", values_to = "epsilon") |> 
          dplyr::mutate(site = parse_number(site),
                        iter = parse_number(iter))) |> 
      dplyr::select(site, iter, mu_gamma0, gamma1, epsilon) |>
      dplyr::full_join( faith_site_d) |> 
      dplyr::mutate( p = plogis( mu_gamma0 + gamma1 * faith_site + epsilon)) |> 
      dplyr::group_by(site, faith_site, faith_site_unscaled) |> 
      dplyr::summarise(mean = mean(p), 
                       l95 = quantile(p, c(0.025)), 
                       u95 = quantile(p, c(0.975))) |> 
      dplyr::full_join(site_key) |> 
      tibble::add_column(scale = "Metacommunity") |> 
      dplyr::rename(xcoord = x, 
                    ycoord = y,
                    x = faith_site, x_unscaled = faith_site_unscaled)
  ) |> 
  tibble::add_column(metric = "Phylogenetic diversity")

site_marginal <- full_join(sr_marginal, pd_marginal) |> 
  full_join(n_marginal) |> 
  full_join(shan_marginal) |>
  full_join(faith_marginal) |> 
  dplyr::mutate(metric = ifelse(metric == "Log( Total abundance )", 
                                "Total abundance", metric)) |> 
  dplyr::mutate(metric = factor(metric, levels = c("Species richness", 
                                                   "Peromyscus dominance", 
                                                   "Total abundance",
                                                   "Shannon diversity",
                                                   "Phylogenetic diversity"))) |> 
  dplyr::mutate(siteID = factor(siteID, levels = c("DSNY", "OSBS", "JERC", "DELA", "TALL", "OAES", "GRSM", "ORNL", 
                                                   "MLBS", "SERC", "SCBI", "UKFS", "BLAN", "KONZ", "KONA", "HARV", 
                                                   "BART", "TREE", "STEI", "UNDE", "NOGP", "WOOD"))) 
  
com_marginal <- dplyr::full_join(sr_com_marginal, pd_com_marginal) |>
  dplyr::full_join(n_com_marginal) |> 
  dplyr::full_join(shan_com_marginal) |> 
  dplyr::full_join(faith_com_marginal) |> 
  dplyr::mutate(metric = ifelse(metric == "Log( Total abundance )", 
                                "Total abundance", metric)) |> 
  dplyr::mutate(metric = factor(metric, levels = c("Species richness", 
                                                   "Peromyscus dominance", 
                                                   "Total abundance",
                                                   "Shannon diversity",
                                                   "Phylogenetic diversity")))
ggplot() +
  geom_ribbon( data = com_marginal, aes(x = x_unscaled, ymin = l95, ymax = u95), 
               alpha = 0.3) +
  geom_line( data = com_marginal,
             aes(x = x_unscaled, y = mean),
             size = 1.5) +
  geom_line(data = site_marginal, 
            aes(x = x_unscaled, y = mean, color = siteID), 
            size = 0.5,
            alpha = 0.75) +
  scale_color_viridis_d(begin = 0.5, end = 1) +
  facet_grid(scale~metric, scales = "free_x") +
  labs( x = "Biodiversity metric value",
        y = expression('P( ' ~italic(Borrelia)~ 'infection )'),
        color = "NEON site") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(color = "black",
                                  size = 9.4),
        axis.title = element_text(color = "black", 
                                  size = 10),
        axis.text = element_text(color = "black", size = 9), 
        legend.title = element_text(color = "black", size = 10), 
        legend.text = element_text(color = "black", size = 9),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray80"),
        panel.grid.minor = element_line(linewidth = 0.05, color = "gray80"),
        plot.background = element_rect(color = NA, fill = "white"), 
        panel.background = element_rect(color = NA, fill = "white")) +
  guides(color = guide_legend(nrow = 3))

setwd(here::here("figures"))
ggsave(
  "figure_04.png", 
  width = 8.5, 
  height = 5.5, 
  units = "in", 
  dpi = 600
)
