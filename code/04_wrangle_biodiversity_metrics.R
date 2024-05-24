library(tidyverse)
library(here)
library(MCMCvis)
library(reshape2)
library(vegan)

# NOTE: this results file is too large to push to GitHub
# User can either run script 03_run_capture_recapture_model.R to produce a similar result file
# OR the file can be downloaded from this OneDrive link: https://1drv.ms/u/s!AtvYBfNq7AMkhIIjVUn6z4ooDwt1Kw?e=iUlju0
setwd(here::here("results"))
load("neon_capture_recapture_results_2024-04-11.RData")

N <- MCMCvis::MCMCpstr( out, params = "N", type = "chains" )

# select 500 random samples from the posterior to summarize biodiversity metrics
# ( to save computation time )
# setting a seed for reproducibility
# if you select different samples from the posterior, some of the reported values may vary a little
set.seed(42)
iters_to_select <- sort(sample(1:3000, 500, replace = FALSE))

N_df <- reshape2::melt( N[[1]], c("sp", "site", "period", "plot", "iter") ) |> 
  dplyr::filter(!is.na(value)) |> 
  dplyr::filter(iter %in% iters_to_select)

sr_site <- N_df |> 
  dplyr::group_by(sp, site, period, iter) |> 
  dplyr::summarise( value = sum(value, na.rm = TRUE)) |> 
  dplyr::group_by(site, period, iter) |> 
  dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |> 
  dplyr::group_by(site, period) |> 
  dplyr::summarise(sr_site_mean = mean(sr), 
                   sr_site_sd = sd(sr))

sr_plot <- N_df |> 
  dplyr::group_by(site, period, plot, iter) |> 
  dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |> 
  dplyr::group_by(site, period, plot) |> 
  dplyr::summarise(sr_plot_mean = mean(sr), 
                   sr_plot_sd = sd(sr))

pd_site <- N_df |> 
  dplyr::group_by(site, period, iter) |> 
  dplyr::mutate(totN = sum(value, na.rm = TRUE), 
                sp_prop = value / totN) |> 
  dplyr::filter(sp %in% c(16, 17, 18, 19)) |> 
  dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |> 
  dplyr::group_by(site, period) |> 
  dplyr::summarise(pd_site_mean = mean(pd), 
                   pd_site_sd = sd(pd))

pd_plot <- N_df |> 
  dplyr::group_by(site, period, plot, iter) |> 
  dplyr::mutate(totN = sum(value, na.rm = TRUE), 
                sp_prop = value / totN) |> 
  dplyr::filter(sp %in% c(16, 17, 18, 19)) |> 
  dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |> 
  dplyr::group_by(site, period, plot) |> 
  dplyr::summarise(pd_plot_mean = mean(pd), 
                   pd_plot_sd = sd(pd))

totN_site <- N_df |> 
  dplyr::group_by( site, period, iter) |> 
  dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
  dplyr::group_by(site, period) |> 
  dplyr::summarise(n_site_mean = mean(log1p(totN)), 
                   n_site_sd = sd(log1p(totN)))

totN_plot <- N_df |> 
  dplyr::group_by(site, period, plot, iter) |> 
  dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
  dplyr::group_by(site, period, plot) |> 
  dplyr::summarise( n_plot_mean = mean(log1p(totN)), 
                    n_plot_sd = sd(log1p(totN)))

shannon_site <- N_df |> 
  dplyr::group_by(sp, site, period, iter) |> 
  dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
  dplyr::group_by(site, period, iter) |> 
  dplyr::summarise(shannon = vegan::diversity(totN)) |> 
  dplyr::group_by(site, period) |> 
  dplyr::summarise(shan_site_mean = mean(shannon), 
                   shan_site_sd = sd(shannon))

shannon_plot <- N_df |> 
  dplyr::group_by(site, period, plot, iter) |> 
  dplyr::summarise( shannon = vegan::diversity(value)) |> 
  dplyr::group_by(site, period, plot) |> 
  dplyr::summarise( shan_plot_mean = mean(shannon), 
                    shan_plot_sd = sd(shannon))

setwd(here::here("data"))
final <- readr::read_csv("neon_model_data_update.csv")

disease_biodiversity <- final |> 
  dplyr::ungroup() |> 
  dplyr::select(siteID, site, plotID, plot, period, scientificName, sp, tagID, positive) |> 
  dplyr::filter(!is.na(positive)) |> 
  dplyr::distinct() |> 
  dplyr::left_join(sr_site) |> 
  dplyr::left_join(sr_plot) |> 
  dplyr::left_join(pd_site) |> 
  dplyr::left_join(pd_plot) |> 
  dplyr::left_join(totN_site) |> 
  dplyr::left_join(totN_plot) |> 
  dplyr::left_join(shannon_site) |> 
  dplyr::left_join(shannon_plot) 

setwd(here::here("data"))
write_csv(disease_biodiversity, "disease_with_biodiversity_metrics_v01.csv")
