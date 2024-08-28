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

all_diffs <- MCMCvis::MCMCpstr(sr_site, params = "gamma1", type = "chains")[[1]] |> 
  tibble::as_tibble() |> 
  tidyr::pivot_longer(1:3000, names_to = "iter") |> 
  tibble::add_column( scale = "site") |> 
  
  dplyr::full_join(
    
    MCMCvis::MCMCpstr(sr_plot, params = "gamma1", type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(1:3000, names_to = "iter") |> 
      tibble::add_column( scale = "plot")) |> 
  tidyr::pivot_wider(names_from = scale, values_from = value) |> 
  dplyr::mutate(diff = site - plot) |> 
  tibble::add_column( metric = "Species richness") |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr(pd_site, params = "gamma1", type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(1:3000, names_to = "iter") |> 
      tibble::add_column( scale = "site") |> 
      
      dplyr::full_join(
        
        MCMCvis::MCMCpstr(pd_plot, params = "gamma1", type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(1:3000, names_to = "iter") |> 
          tibble::add_column( scale = "plot")
        
      ) |> 
      
      tidyr::pivot_wider(names_from = scale, values_from = value) |> 
      dplyr::mutate(diff = site - plot) |> 
      tibble::add_column( metric = "Peromyscus dominance")
  ) |>  dplyr::full_join(
    
    
    MCMCvis::MCMCpstr(n_site, params = "gamma1", type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(1:3000, names_to = "iter") |> 
      tibble::add_column( scale = "site") |> 
      
      dplyr::full_join(
        
        MCMCvis::MCMCpstr(n_plot, params = "gamma1", type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(1:3000, names_to = "iter") |> 
          tibble::add_column( scale = "plot")) |> 
      tidyr::pivot_wider(names_from = scale, values_from = value) |> 
      dplyr::mutate(diff = site - plot) |> 
      tibble::add_column( metric = "Total abundance") ) |> 
  dplyr::full_join(
    MCMCvis::MCMCpstr(shan_site, params = "gamma1", type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(1:3000, names_to = "iter") |> 
      tibble::add_column( scale = "site") |> 
      
      dplyr::full_join(
        
        MCMCvis::MCMCpstr(shan_plot, params = "gamma1", type = "chains")[[1]] |> 
          tibble::as_tibble() |> 
          tidyr::pivot_longer(1:3000, names_to = "iter") |> 
          tibble::add_column( scale = "plot")) |> 
      tidyr::pivot_wider(names_from = scale, values_from = value) |> 
      dplyr::mutate(diff = site - plot) |> 
      tibble::add_column( metric = "Shannon diversity")) |> 
  dplyr::full_join(

MCMCvis::MCMCpstr(faith_site, params = "gamma1", type = "chains")[[1]] |> 
  tibble::as_tibble() |> 
  tidyr::pivot_longer(1:3000, names_to = "iter") |> 
  tibble::add_column( scale = "site") |> 
  
  dplyr::full_join(
    
    MCMCvis::MCMCpstr(faith_plot, params = "gamma1", type = "chains")[[1]] |> 
      tibble::as_tibble() |> 
      tidyr::pivot_longer(1:3000, names_to = "iter") |> 
      tibble::add_column( scale = "plot")) |> 
  tidyr::pivot_wider(names_from = scale, values_from = value) |> 
  dplyr::mutate(diff = site - plot) |> 
  tibble::add_column( metric = "Phylogenetic diversity")
  ) |> 
  dplyr::mutate( metric = factor(metric, levels = c(
    "Species richness", 
    "Peromyscus dominance", 
    "Total abundance", 
    "Shannon diversity", 
    "Phylogenetic diversity"
  )))

props <- all_diffs |> 
  group_by(metric) |> 
  summarise( prop = round( sum(diff > 0) / sum(!is.na(diff)), 3)) |> 
  tibble::add_column(diffs = 2.5,
                     y = 600)

setwd(here::here("figures"))

ggplot() +
  facet_wrap(~metric) +
  geom_histogram(data = all_diffs, aes(x = diff),
                 fill = MetBrewer::MetPalettes$Hiroshige[[1]][8],
                 color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  geom_label(data = props, 
             aes(x = diffs, y = y, label = prop),
             size = 3) +
  geom_vline(xintercept = 0, color = MetBrewer::MetPalettes$Hiroshige[[1]][1], 
             linewidth = 1) +
  labs(x = "Difference in biodiversity effect (site - plot)",
       y = "Count") +
  theme_classic() +
  theme(strip.background = element_rect(color = NA),
        strip.text = element_text(size = 10, color = "black"), 
        axis.text = element_text(size = 9, color = "black"), 
        axis.title = element_text(size = 10, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.ticks = element_line(color = "black", linewidth = 0.2))

ggsave(
  filename = "figure_s01.png", 
  width = 5.2, 
  height = 3, 
  units = "in", 
  dpi = 600
) 