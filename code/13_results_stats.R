library(tidyverse)
library(here)
library(MCMCvis)
library(reshape2)
library(vegan)

setwd(here::here("results"))
load("neon_capture_recapture_results_2024-01-19.RData")

N <- MCMCvis::MCMCpstr( out, params = "N", type = "chains" )

# select 500 random samples from the posterior to summarize biodiversity metrics
# ( to save computation time )
# setting a seed for reproducibility
# if you select different samples from the posterior, some of the reported values may vary a little
set.seed(42)
iters_to_select <- sort(sample(1:4500, 500, replace = FALSE))

N_df <- reshape2::melt( N[[1]], c("sp", "site", "period", "plot", "iter") ) |> 
  dplyr::filter(!is.na(value)) |> 
  dplyr::filter(iter %in% iters_to_select)

# report the sites with lowest and highest metacommunity species richness
( sr_site <- N_df |> 
    dplyr::group_by(sp, site, period, iter) |> 
    dplyr::summarise( value = sum(value, na.rm = TRUE)) |> 
    dplyr::group_by(site, period, iter) |> 
    dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |> 
    dplyr::group_by(site, period) |> 
    dplyr::summarise(sr_site_median = median(sr),
                     sr_site_l95 = quantile(sr, c(0.025)), 
                     sr_site_u95 = quantile(sr, c(0.975))) |> 
    dplyr::group_by(site) |> 
    dplyr::filter(sr_site_median == max(sr_site_median)) |> 
    dplyr::slice(1) |> 
    dplyr::ungroup() |> 
    dplyr::filter(sr_site_median == min(sr_site_median) | sr_site_median == max(sr_site_median)) |> 
    dplyr::left_join(
      final |>
        dplyr::select(siteID, site) |>
        dplyr::distinct()) )  

( pd_site <- N_df |> 
    dplyr::group_by(site, period, iter) |> 
    dplyr::mutate(totN = sum(value, na.rm = TRUE), 
                  sp_prop = value / totN) |> 
    dplyr::filter(sp %in% c(16, 17, 18, 19)) |> 
    dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |> 
    dplyr::group_by(site, period) |> 
    dplyr::summarise(pd_site_median = median(pd), 
                     pd_site_l95 = quantile(pd, c(0.025)), 
                     pd_site_u95 = quantile(pd, c(0.975))) |> 
    dplyr::group_by(site) |> 
    dplyr::filter(pd_site_median == max(pd_site_median)) |>
    dplyr::slice(1) |> 
    dplyr::ungroup() |> 
    dplyr::filter(pd_site_median == min(pd_site_median) | pd_site_median == max(pd_site_median) ) |> 
    dplyr::left_join(
      final |>
        dplyr::select(siteID, site) |>
        dplyr::distinct()) |> 
    dplyr::mutate(across(pd_site_median:pd_site_u95, function(x) round(x, 2))) )  

( totN_site <- N_df |> 
    dplyr::group_by( site, period, iter) |> 
    dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
    dplyr::group_by(site, period) |> 
    dplyr::summarise(n_site_median = median(totN), 
                     n_site_l95 = quantile(totN, c(0.025)),
                     n_site_u95 = quantile(totN, c(0.975))) |> 
    dplyr::group_by(site) |> 
    dplyr::filter(n_site_median == max(n_site_median)) |> 
    dplyr::slice(1) |> 
    dplyr::ungroup() |> 
    dplyr::filter(n_site_median == min(n_site_median) | n_site_median == max(n_site_median)) |> 
    dplyr::left_join(
      final |>
        dplyr::select(siteID, site) |>
        dplyr::distinct()) )

( shannon_site <- N_df |> 
    dplyr::group_by(sp, site, period, iter) |> 
    dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
    dplyr::group_by(site, period, iter) |> 
    dplyr::summarise(shannon = vegan::diversity(totN)) |> 
    dplyr::group_by(site, period) |> 
    dplyr::summarise(shan_site_median = median(shannon), 
                     shan_site_l95 = quantile(shannon, c(0.025)),
                     shan_site_u95 = quantile(shannon, c(0.975))) |> 
    dplyr::group_by(site) |> 
    dplyr::filter(shan_site_median == max(shan_site_median)) |> 
    dplyr::slice(1) |> 
    dplyr::ungroup() |> 
    dplyr::filter(shan_site_median == min(shan_site_median) | shan_site_median == max(shan_site_median)) |> 
    dplyr::left_join(
      final |>
        dplyr::select(siteID, site) |>
        dplyr::distinct()) |> 
    dplyr::mutate(across(shan_site_median:shan_site_u95, function(x) round(x, 2))) )

( sr_sd <- N_df |> 
    dplyr::group_by(site, period, plot, iter) |> 
    dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |> 
    dplyr::group_by(site, period) |> 
    dplyr::summarise(sr_sd = sd(sr)) |> 
    dplyr::ungroup() |> 
    dplyr::filter(sr_sd > 0) |> 
    dplyr::filter(sr_sd == min(sr_sd) | sr_sd == max(sr_sd)) |> 
    dplyr::left_join(
      final |>
        dplyr::select(siteID, site) |>
        dplyr::distinct()) |> 
    dplyr::mutate(sr_sd = round(sr_sd, 1)) )