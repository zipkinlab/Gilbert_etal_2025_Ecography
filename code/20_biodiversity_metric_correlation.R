# 19 August 2024
# Script to explore collinearity of biodiversity metrics
# Objective 1: Evaluate correlation between plot- and site-level scales for each metric
# Objective 2: Evaluate how spacing of plots within a site affects plot-site metric correlations
# Objective 3: Within each scale, evluate correlation between difference metrics
library(here)
library(tidyverse)
library(corrr)
library(MetBrewer)
library(corrplot)

setwd(here::here("data"))

d <- readr::read_csv( "disease_with_biodiversity_metrics_v01.csv" )
site_size <- readr::read_csv("site_distance_area.csv")

metrics <- d |> 
  dplyr::select(siteID, plotID, period, sr_site_mean:faith_plot_sd ) |> 
  dplyr::distinct() |> 
  dplyr::left_join(site_size) |> 
  dplyr::select(
    siteID:period, 
    sr_site_mean, 
    sr_plot_mean, 
    pd_site_mean, 
    pd_plot_mean, 
    n_site_mean, 
    n_plot_mean, 
    shan_site_mean, 
    shan_plot_mean,
    faith_site_mean, 
    faith_plot_mean, 
    mean, med, area_km2)

# Correlation between plot-level and site-level biodiversity metrics
key <- tibble::tribble(
  ~metric1, ~name, 
  "sr", "Species richness", 
  "pd", "Peromyscus dominance", 
  "n", "Total abundance", 
  "shan", "Shannon diversity", 
  "faith", "Phylogenetic diversity")

scale_cor <- metrics |> 
  dplyr::select(sr_site_mean:faith_plot_mean) |> 
  cor() |> 
  tibble::as_tibble(rownames = "metric1") |> 
  tidyr::pivot_longer(sr_site_mean:faith_plot_mean,
                      names_to = "metric2") |> 
  dplyr::filter(!metric1 == metric2) |> 
  tidyr::separate(metric1, into = c("metric1", "scale1", "junk1")) |> 
  tidyr::separate(metric2, into = c("metric2", "scale2", "junk2")) |> 
  dplyr::filter(metric1 == metric2) |> 
  dplyr::select(-junk1, -junk2) |> 
  dplyr::select(metric1, value) |> 
  dplyr::distinct() |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(value = round(value, 2)) |> 
  dplyr::select(`Biodiversity metric` = name, 
                `Pearson correlation between scales` = value)

flextable::set_flextable_defaults(font.size = 10)
ft <- flextable::flextable( data = scale_cor, cwidth = 2)  

tmp <- tempfile(fileext = ".docx")

officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)

# Looking at how area/size of NEON sites relates to 
# how similar plot-level and site-level biodiversity metrics are
metrics |> 
  dplyr::select(siteID, plotID, period, 
                sr_site_mean:faith_plot_mean) |> 
  dplyr::mutate(sr = sr_site_mean - sr_plot_mean, 
                pd = pd_site_mean - pd_plot_mean, 
                n = n_site_mean - n_plot_mean,
                shan = shan_site_mean - shan_plot_mean, 
                faith = faith_site_mean - faith_plot_mean) |> 
  dplyr::select(siteID:period, 
                sr:faith) |> 
  dplyr::left_join(site_size) |> 
  tidyr::pivot_longer(sr:faith, 
                      names_to = "metric", 
                      values_to = "diff") |> 
  dplyr::left_join(key |> 
                     dplyr::rename(metric = metric1)) |> 
  
  ggplot2::ggplot(aes(x = log(med), y = diff)) + 
  ggplot2::facet_wrap(~name, scales = "free") + 
  ggplot2::geom_point(alpha = 0.2, color = MetBrewer::MetPalettes$Cassatt2[[1]][2]) +
  ggplot2::geom_smooth(method = "lm", se = FALSE,
                       color = MetBrewer::MetPalettes$Cassatt2[[1]][8]) +
  ggplot2::labs( x = "Log ( median distance between plots )",
                 y = "Difference in biodiversity metric value (site - plot)") +
  ggplot2::theme_classic() +
  ggplot2::theme(strip.background = element_rect(color = NA, 
                                                 fill = "white"))

setwd(here::here("figures"))
ggsave(
  filename = "figure_s08.png", 
  width = 5.5, 
  height = 4, 
  units = "in", 
  dpi = 600)

# calculate correlation between biodiversity metric DIFFERENCE (between scales)
# and spacing of plots within a site
space_cor <- metrics |> 
  dplyr::select(siteID, plotID, period, 
                sr_site_mean:faith_plot_mean) |> 
  dplyr::mutate(sr = sr_site_mean - sr_plot_mean, 
                pd = pd_site_mean - pd_plot_mean, 
                n = n_site_mean - n_plot_mean,
                shan = shan_site_mean - shan_plot_mean, 
                faith = faith_site_mean - faith_plot_mean) |> 
  dplyr::select(siteID:period, 
                sr:faith) |> 
  dplyr::left_join(site_size) |> 
  tidyr::pivot_longer(sr:faith, 
                      names_to = "metric", 
                      values_to = "diff") |> 
  dplyr::left_join(key |> 
                     dplyr::rename(metric = metric1)) |> 
  dplyr::group_by(metric) |> 
  dplyr::summarise( cor = round(cor(diff, med), 2)) |> 
  dplyr::left_join(key |> 
                     dplyr::rename(metric = metric1)) |> 
  dplyr::select(`Biodiversity metric` = name, 
                `Correlation between median distance between plots and scale difference in metric` = cor)


flextable::set_flextable_defaults(font.size = 10)
ft <- flextable::flextable( data = space_cor, cwidth = 2)  

tmp <- tempfile(fileext = ".docx")

officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)

# correlation between site-level metrics
metrics |> 
  dplyr::select(matches("site_")) |>
  dplyr::rename_with(~str_replace(., "_site_mean", "")) |> 
  corrr::correlate() |> 
  corrr::shave() |> 
  tidyr::pivot_longer(sr:faith, names_to = "term2") |>
  dplyr::left_join(
    key |> 
      dplyr::rename(term = metric1, 
             metric1 = name)
  ) |> 
  dplyr::left_join(
    key |> 
      dplyr::rename(
        term2 = metric1, 
        metric2 = name)
      ) |> 
  dplyr::select(metric1, metric2, value) |> 
  dplyr::filter(!is.na(value)) # |> 
# dplyr::filter(abs(value) > 0.6) # 6/10 have correlation > 0.6

# correlation between plot-level metrics
metrics |> 
  dplyr::select(matches("plot_")) |>
  dplyr::rename_with(~str_replace(., "_plot_mean", "")) |> 
  # dplyr::mutate_all(function(x) stringr::str_remove_all(x, "_site_mean"))
  corrr::correlate() |> 
  corrr::shave() |> 
  tidyr::pivot_longer(sr:faith, names_to = "term2") |> 
  dplyr::left_join(
    key |> 
      dplyr::rename(term = metric1, 
                    metric1 = name)
  ) |> 
  dplyr::left_join(
    key |> 
      dplyr::rename(
        term2 = metric1, 
        metric2 = name)
  ) |> 
  dplyr::select(metric1, metric2, value) |> 
  dplyr::filter(!is.na(value))  #|> 
# dplyr::filter(abs(value) > 0.6) # 4/10 have correlation > 0.6

p <- metrics |> 
  dplyr::select(matches("plot_")) |>
  dplyr::rename_with(~str_replace(., "_plot_mean", "")) |>
  setNames(c("SR", "PD",
             "N", "Shan", 
             "Phylo")) |>
  cor()

corrplot(p, method = "number", type = "upper", diag = FALSE)

s <- metrics |> 
  dplyr::select(matches("site_")) |>
  dplyr::rename_with(~str_replace(., "_site_mean", "")) |>
  setNames(c("SR", "PD",
             "N", "Shan", 
             "Phylo")) |>
  cor()

corrplot(s, method = "number", type = "upper", diag = FALSE)