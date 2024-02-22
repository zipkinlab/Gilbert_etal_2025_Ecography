library(here)
library(tidyverse)
library(MCMCvis)
library(nimble)

setwd(here::here("results"))

load("rodent_pathogen_sr_plot_2024-01-21.RData")
sr_plot <- out
load("rodent_pathogen_sr_site_2024-01-21.RData")
sr_site <- out
load("rodent_pathogen_pd_plot_2024-01-21.RData")
pd_plot <- out
load("rodent_pathogen_pd_site_2024-01-21.RData")
pd_site <- out
load("rodent_pathogen_shan_plot_2024-01-21.RData")
shannon_plot <- out
load("rodent_pathogen_shan_site_2024-01-21.RData")
shannon_site <- out
load("rodent_pathogen_n_plot_2024-01-21.RData")
n_plot <- out
load("rodent_pathogen_n_site_2024-01-21.RData")
n_site <- out

setwd(here::here("data"))

d <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv")
load("neon_cr_data_2024-01-17.RData")

site_key <- final |>
  dplyr::select(site, siteID) |>
  dplyr::distinct()

sr_plot_d <- d |>
  dplyr::select(site, sr_plot_mean) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(sr_plot_mean == min(sr_plot_mean) | sr_plot_mean == max(sr_plot_mean)) |>
  dplyr::mutate( type = ifelse(sr_plot_mean == min(sr_plot_mean), "min", "max")) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_plot_mean) |>
  dplyr::mutate( max = ifelse(is.na(max), min, max)) |>
  dplyr::mutate( sr_plot = list(seq(from = min, to = max, by = 0.05))) |>
  tidyr::unnest() |>
  dplyr::ungroup() |>
  dplyr::mutate(sr_plot = as.numeric(scale(sr_plot))) |>
  dplyr::select(-min, -max)

sr_site_d <- d |>
  dplyr::select(site, sr_site_mean) |>
  dplyr::distinct() |>
  dplyr::group_by(site) |>
  dplyr::filter(sr_site_mean == min(sr_site_mean) | sr_site_mean == max(sr_site_mean)) |>
  dplyr::mutate( type = ifelse(sr_site_mean == min(sr_site_mean), "min", "max")) |>
  tidyr::pivot_wider(names_from = type, values_from = sr_site_mean) |>
  dplyr::mutate( max = ifelse(is.na(max), min, max)) |>
  dplyr::mutate( sr_site = list(seq(from = min, to = max, by = 0.05))) |>
  tidyr::unnest() |>
  dplyr::ungroup() |>
  dplyr::mutate(sr_site = as.numeric(scale(sr_site)))

sr_plot_sc <- MCMCvis::MCMCsummary(sr_plot, params = c("sd_sr"))$mean
sr_plot_me <- MCMCvis::MCMCsummary(sr_plot, params = c("mean_sr"))$mean

sr_plot_p <- MCMCvis::MCMCpstr( sr_plot, params = c("gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble(rownames = "site") |>
  tidyr::pivot_longer(2:3001, names_to = "iter", values_to = "gamma0") |>
  dplyr::mutate(site = stringr::str_remove(site, "gamma0")) |>
  dplyr::mutate(site = readr::parse_number(site)) |>
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_plot, params = c("gamma1"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma1")) |>
  dplyr::full_join(sr_plot_d) |>
  dplyr::mutate( p = nimble::ilogit(gamma0 + gamma1 * sr_plot )) |>
  dplyr::group_by( site, sr_plot ) |>
  dplyr::summarise( mean = mean(p),
                    l95 = quantile(p, c(0.025)),
                    u95 = quantile(p, c(0.975))) |>
  dplyr::ungroup() |>
  dplyr::mutate( sr_plot_unscaled = sr_plot * sr_plot_sc  + sr_plot_me )

sr_site_sc <- MCMCvis::MCMCsummary(sr_site, params = c("sd_sr"))$mean
sr_site_me <- MCMCvis::MCMCsummary(sr_site, params = c("mean_sr"))$mean

sr_site_p <- MCMCvis::MCMCpstr( sr_site, params = c("gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble(rownames = "site") |>
  tidyr::pivot_longer(2:3001, names_to = "iter", values_to = "gamma0") |>
  dplyr::mutate(site = stringr::str_remove(site, "gamma0")) |>
  dplyr::mutate(site = readr::parse_number(site)) |>
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_site, params = c("gamma1"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma1")) |>
  dplyr::full_join(sr_site_d) |>
  dplyr::mutate( p = nimble::ilogit(gamma0 + gamma1 * sr_site )) |>
  dplyr::group_by( site, sr_site ) |>
  dplyr::summarise( mean = mean(p),
                    l95 = quantile(p, c(0.025)),
                    u95 = quantile(p, c(0.975))) |>
  dplyr::ungroup() |>
  dplyr::mutate( sr_site_unscaled = sr_site * sr_site_sc  + sr_site_me )

sr_site_preds <- sr_site_p |>
  tibble::add_column(scale = "Metacommunity") |>
  dplyr::rename(sr = sr_site_unscaled) |>
  dplyr::full_join(
    sr_plot_p |>
      tibble::add_column(scale = "Local community") |>
      dplyr::rename(sr = sr_plot_unscaled)) |>
  dplyr::full_join(site_key) |>
  dplyr::mutate(siteID = factor(siteID,
                                levels = c("TREE", "SCBI", "SERC", "ORNL", "BLAN", "JERC",
                                           "DELA", "DSNY", "BART", "MLBS", "STEI", "OAES",
                                           "UKFS", "HARV", "TALL", "NOGP", "UNDE", "WOOD",
                                           "GRSM", "OSBS", "KONA", "KONZ")))

sr_site_overall <- tibble::tibble(
  sr_site = seq(from = min(d$sr_site_mean), to = max(sr_site_p$sr_site_unscaled), by = 0.05)) |> 
  dplyr::mutate(sr_site_sc = as.numeric(scale(sr_site)))

sr_plot_overall <- tibble::tibble(
  sr_plot = seq(from = min(d$sr_plot_mean), to = max(sr_plot_p$sr_plot_unscaled), by = 0.05)) |> 
  dplyr::mutate(sr_plot_sc = as.numeric(scale(sr_plot)))

overall <- MCMCvis::MCMCpstr( sr_site, params = c("mu_gamma0"), type = "chains")[[1]] |>
  tibble::as_tibble(rownames = "site") |>
  tidyr::pivot_longer(2:3001, names_to = "iter", values_to = "gamma0") |>
  dplyr::select(-site) |> 
  # dplyr::mutate(site = stringr::str_remove(site, "gamma0")) |>
  # dplyr::mutate(site = readr::parse_number(site)) |>
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_site, params = c("gamma1"), type = "chains")[[1]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma1")) |>
  dplyr::cross_join(sr_site_overall) |>
  dplyr::mutate( p = nimble::ilogit(gamma0 + gamma1 * sr_site_sc )) |>
  dplyr::group_by( sr_site, sr_site_sc) |>
  dplyr::summarise( mean = mean(p),
                    l95 = quantile(p, c(0.025)),
                    u95 = quantile(p, c(0.975))) |>
  dplyr::ungroup() |> 
  dplyr::rename(sr = sr_site, 
                sr_sc = sr_site_sc) |> 
  tibble::add_column(scale = "Metacommunity") |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr( sr_plot, params = c("mu_gamma0"), type = "chains")[[1]] |>
      tibble::as_tibble(rownames = "site") |>
      tidyr::pivot_longer(2:3001, names_to = "iter", values_to = "gamma0") |>
      dplyr::select(-site) |> 
      # dplyr::mutate(site = stringr::str_remove(site, "gamma0")) |>
      # dplyr::mutate(site = readr::parse_number(site)) |>
      dplyr::full_join(
        MCMCvis::MCMCpstr( sr_plot, params = c("gamma1"), type = "chains")[[1]] |>
          tibble::as_tibble() |>
          tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "gamma1")) |>
      dplyr::cross_join(sr_plot_overall) |>
      dplyr::mutate( p = nimble::ilogit(gamma0 + gamma1 * sr_plot_sc )) |>
      dplyr::group_by( sr_plot, sr_plot_sc) |>
      dplyr::summarise( mean = mean(p),
                        l95 = quantile(p, c(0.025)),
                        u95 = quantile(p, c(0.975))) |>
      dplyr::ungroup() |> 
      dplyr::rename(sr = sr_plot, 
                    sr_sc = sr_plot_sc) |> 
      tibble::add_column(scale = "Local community"))

ggplot() + 
  facet_wrap(~scale) +
  geom_ribbon(data = overall, aes(x = sr, ymin = l95, ymax = u95), color = NA, alpha = 0.4) + 
  geom_line(data = overall, aes(x = sr, y = mean), size = 2) +
  geom_ribbon( data = sr_site_preds, aes(x = sr, ymin = l95, ymax = u95, fill = siteID),
               color = NA, alpha = 0.1) +
  geom_line( data = sr_site_preds, aes(x = sr, y = mean, color = siteID), size = 1) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ggplot2::theme_classic() +
  ggplot2::labs( x = "Species richness",
                 y = expression('P( ' ~italic(Borrelia)~ 'infection )'),
                 color = "NEON site",
                 fill = "NEON site") +
  ggplot2::theme( strip.background = element_blank(),
                  legend.position = "bottom",
                  strip.text = element_text(color = "black", size = 10),
                  axis.title = element_text(color = "black", size = 11),
                  axis.text = element_text(color = "black", size = 10),
                  legend.text = element_text(color = "black", size = 9),
                  legend.title = element_text(color = "black", size = 8),
                  axis.line = element_line(color = "black", size = 0.1),
                  axis.ticks = element_line(color = "black", size = 0.1))  +
  ggplot2::guides(color = guide_legend(nrow = 4, title.position = "top", title.hjust = 0.5),
                  fill = guide_legend(nrow = 4, title.position = "top", title.hjust = 0.5))

setwd(here::here("figures"))
ggsave(
  filename = "species_richness_effect_v01.png",
  width = 5,
  height = 4.5,
  units = "in",
  dpi = 300
)


MCMCvis::MCMCsummary( pd_site, params = c("gamma1")) |> 
  tibble::as_tibble(rownames = "param") |> 
  dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
  tibble::add_column(metric = "PD", scale = "site") |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( pd_plot, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "PD", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( sr_site, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "SR", scale = "site")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( sr_plot, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "SR", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( n_plot, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Total abundance", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( n_site, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Total abundance", scale = "site")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( shannon_plot, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Shannon diversity", scale = "plot")) |> 
  dplyr::full_join(
    MCMCvis::MCMCsummary( shannon_site, params = c("gamma1")) |> 
      tibble::as_tibble(rownames = "param") |> 
      dplyr::select(param, mean, l95 = `2.5%`, u95 = `97.5%`) |> 
      tibble::add_column(metric = "Shannon diversity", scale = "site"))  |> 
  dplyr::mutate( metric = ifelse(metric == "SR", "Species richness",
                                 ifelse(metric == "PD", "Peromyscus dominance", metric))) |> 
  dplyr::mutate(metric_styled = ifelse(grepl("Pero", metric), "<i>Peromyscus</i> dominance", metric)) |> 
  dplyr::mutate(metric_styled = factor(metric_styled, levels = c(
    "Shannon diversity",
    "Total abundance", 
    "<i>Peromyscus</i> dominance", 
    "Species richness"))) |>
  dplyr::mutate( scale = ifelse(scale == "plot", "Local community", "Metacommunity")) |> 
  dplyr::mutate(scale = factor(scale, levels = c("Local community", "Metacommunity"))) |>
  ggplot2::ggplot( aes( x = mean, y = metric_styled)) + 
  ggplot2::geom_vline(xintercept = 0, color = "gray30", linetype = "dashed") +
  ggplot2::geom_errorbar( aes(xmin = l95, xmax = u95, color = scale),
                          width = 0, size = 1.5, position = position_dodge(width = 0.25)) + 
  ggplot2::geom_point(aes(color = scale), size = 3, position = position_dodge(width = 0.25)) + 
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
  "coefficient_plot_v01.png", 
  width = 4.85, 
  height = 4, 
  units = "in", 
  dpi = 300
)  

MCMCpstr(sr_plot, params = c("gamma1"), type = "chains")[[1]] |>
  as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
  summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
  add_column( metric = "SR", 
              scale = "plot") |> 
  dplyr::select(metric, scale, p_positive) |> 
  
  full_join(
    MCMCpstr(sr_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "SR", 
                  scale = "site") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  
  full_join(
    MCMCpstr(pd_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "PD", 
                  scale = "site") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  full_join(
    MCMCpstr(pd_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "PD", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  full_join(
    MCMCpstr(n_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "N", 
                  scale = "site") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  full_join(
    MCMCpstr(n_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "N", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  full_join(
    MCMCpstr(shannon_plot, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "Shannon", 
                  scale = "plot") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  full_join(
    MCMCpstr(shannon_site, params = c("gamma1"), type = "chains")[[1]] |>
      as_tibble() |> pivot_longer(1:3000, names_to = "iter", values_to = "value") |>
      summarise( p_positive = sum(value >= 0) / sum(!is.na(value))) |> 
      add_column( metric = "Shannon", 
                  scale = "site") |> 
      dplyr::select(metric, scale, p_positive)
  ) |> 
  mutate(p_positive = round(p_positive, 2),
         p_negative = round(1 - p_positive, 2))

