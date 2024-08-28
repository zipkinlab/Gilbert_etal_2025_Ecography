library(tidyverse)
library(here)
library(MCMCvis)
library(reshape2)
library(vegan)
library(furrr)
library(ape)
library(caper)

# # NOTE: this results file is too large to push to GitHub
# # User can either run script 03_run_capture_recapture_model.R to produce a similar result file
# # OR the file can be downloaded from this OneDrive link: https://1drv.ms/u/s!AtvYBfNq7AMkhIIjVUn6z4ooDwt1Kw?e=iUlju0
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

setwd(here::here("results"))

# option to save posterior distributions of biodiversity metrics
# N_df |>
#   dplyr::group_by(sp, site, period, iter) |>
#   dplyr::summarise( value = sum(value, na.rm = TRUE)) |>
#   dplyr::group_by(site, period, iter) |>
#   dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |>
#   readr::write_csv("sr_site_posterior.csv")

sr_site <- N_df |>
  dplyr::group_by(sp, site, period, iter) |>
  dplyr::summarise( value = sum(value, na.rm = TRUE)) |>
  dplyr::group_by(site, period, iter) |>
  dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |>
  dplyr::group_by(site, period) |>
  dplyr::summarise(sr_site_mean = mean(sr),
                   sr_site_sd = sd(sr))

# N_df |> 
#   dplyr::group_by(site, period, plot, iter) |> 
#   dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |> 
#   readr::write_csv("sr_plot_posterior.csv")

sr_plot <- N_df |>
  dplyr::group_by(site, period, plot, iter) |>
  dplyr::summarise(sr = sum(value > 0, na.rm = TRUE)) |>
  dplyr::group_by(site, period, plot) |>
  dplyr::summarise(sr_plot_mean = mean(sr),
                   sr_plot_sd = sd(sr))

# N_df |> 
#   dplyr::group_by(site, period, iter) |> 
#   dplyr::mutate(totN = sum(value, na.rm = TRUE), 
#                 sp_prop = value / totN) |> 
#   dplyr::filter(sp %in% c(16, 17, 18, 19)) |> 
#   dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |> 
#   readr::write_csv("pd_site_posterior.csv")

pd_site <- N_df |>
  dplyr::group_by(site, period, iter) |>
  dplyr::mutate(totN = sum(value, na.rm = TRUE),
                sp_prop = value / totN) |>
  dplyr::filter(sp %in% c(16, 17, 18, 19)) |>
  dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |>
  dplyr::group_by(site, period) |>
  dplyr::summarise(pd_site_mean = mean(pd),
                   pd_site_sd = sd(pd))

# N_df |> 
#   dplyr::group_by(site, period, plot, iter) |> 
#   dplyr::mutate(totN = sum(value, na.rm = TRUE), 
#                 sp_prop = value / totN) |> 
#   dplyr::filter(sp %in% c(16, 17, 18, 19)) |> 
#   dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |> 
#   readr::write_csv("pd_plot_posterior.csv")

pd_plot <- N_df |>
  dplyr::group_by(site, period, plot, iter) |>
  dplyr::mutate(totN = sum(value, na.rm = TRUE),
                sp_prop = value / totN) |>
  dplyr::filter(sp %in% c(16, 17, 18, 19)) |>
  dplyr::summarise(pd = sum(sp_prop, na.rm = TRUE)) |>
  dplyr::group_by(site, period, plot) |>
  dplyr::summarise(pd_plot_mean = mean(pd),
                   pd_plot_sd = sd(pd))

# N_df |> 
#   dplyr::group_by( site, period, iter) |> 
#   dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
#   readr::write_csv("n_site_posterior.csv")

totN_site <- N_df |>
  dplyr::group_by( site, period, iter) |>
  dplyr::summarise(totN = sum(value, na.rm = TRUE)) |>
  dplyr::group_by(site, period) |>
  dplyr::summarise(n_site_mean = mean(log1p(totN)),
                   n_site_sd = sd(log1p(totN)))
# N_df |> 
#   dplyr::group_by(site, period, plot, iter) |> 
#   dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
#   readr::write_csv("n_plot_posterior.csv")

totN_plot <- N_df |>
  dplyr::group_by(site, period, plot, iter) |>
  dplyr::summarise(totN = sum(value, na.rm = TRUE)) |>
  dplyr::group_by(site, period, plot) |>
  dplyr::summarise( n_plot_mean = mean(log1p(totN)),
                    n_plot_sd = sd(log1p(totN)))

# N_df |> 
#   dplyr::group_by(sp, site, period, iter) |> 
#   dplyr::summarise(totN = sum(value, na.rm = TRUE)) |> 
#   dplyr::group_by(site, period, iter) |> 
#   dplyr::summarise(shannon = vegan::diversity(totN)) |> 
#   readr::write_csv("shan_site_posterior.csv")

shannon_site <- N_df |>
  dplyr::group_by(sp, site, period, iter) |>
  dplyr::summarise(totN = sum(value, na.rm = TRUE)) |>
  dplyr::group_by(site, period, iter) |>
  dplyr::summarise(shannon = vegan::diversity(totN)) |>
  dplyr::group_by(site, period) |>
  dplyr::summarise(shan_site_mean = mean(shannon),
                   shan_site_sd = sd(shannon))

# N_df |> 
#   dplyr::group_by(site, period, plot, iter) |> 
#   dplyr::summarise( shannon = vegan::diversity(value)) |> 
#   readr::write_csv("shan_plot_posterior.csv")

shannon_plot <- N_df |>
  dplyr::group_by(site, period, plot, iter) |>
  dplyr::summarise( shannon = vegan::diversity(value)) |>
  dplyr::group_by(site, period, plot) |>
  dplyr::summarise( shan_plot_mean = mean(shannon),
                    shan_plot_sd = sd(shannon))

# calculate abundance-weighted Faith index (phylogenetic diversity)
sp_key <- final |>
  dplyr::select(sp, scientificName) |>
  dplyr::distinct() |>
  dplyr::mutate(tipLabel = str_replace(scientificName, " ", "_")) |>
  dplyr::select(sp, tipLabel)

site_n <- N_df |>
  dplyr::group_by(sp, site, period, iter) |>
  dplyr::summarise(n = sum(value)) |>
  dplyr::left_join(sp_key) |>
  dplyr::group_by(site, period, iter) |>
  dplyr::mutate(index = dplyr::cur_group_id())

setwd(here::here("data"))
trees <- ape::read.nexus("output.nex")
tree <- trees[[1]]

plot_n <- N_df |>
  dplyr::left_join(sp_key) |>
  dplyr::group_by(site, period, plot, iter) |>
  dplyr::mutate(index = dplyr::cur_group_id()) |>
  dplyr::select(index, tipLabel, sp, site, period, plot, iter, n = value)

# Function to compute abundance-weighted Faith index for a single site
compute_phylo_diversity <- function(site_data, tree) {

# Species with zero-counts
  zeros <- site_data |>
    dplyr::filter(n == 0) |>
    dplyr::pull(tipLabel)

  # Drop species that don't occur at a site
  pruned_tree <- ape::drop.tip(tree, zeros)

  # Check if the pruned tree is empty (i.e., no tips left)
  if (length(pruned_tree$tip.label) == 0) {
    return(tibble::tibble(index = site_data$index[1], n_pd = NA_real_))
  }

  branches <- pruned_tree$edge
  branch_lengths <- pruned_tree$edge.length

  # Calculate mean abundance for subtending tips for a given node
  branch_abundance <- purrr::map_dbl(1:nrow(branches), function(j) {
    leaves.node <- caper::clade.members(x = branches[j, 2], phy = pruned_tree, tip.labels = TRUE)
    mean_abundance <- site_data |>
      dplyr::filter(n > 0) |>
      dplyr::filter(tipLabel %in% leaves.node) |>
      dplyr::pull(n) |>
      mean(na.rm = TRUE)  # Ensure mean handles NA values
    return(mean_abundance)
  })

  # Calculate abundance-weighted Faith index
  denom <- sum(branch_abundance, na.rm = TRUE)
  numerator <- sum(branch_lengths * branch_abundance, na.rm = TRUE)
  n_branches <- nrow(branches)

  # Handle cases where denom might be NA or zero to avoid division by zero
  if (is.na(denom) || denom == 0) {
    n_pd <- NA_real_
  } else {
    n_pd <- n_branches * (numerator / denom)
  }

  tibble::tibble(index = site_data$index[1], n_pd = n_pd)
}

plan(multisession, workers = 8)

phylo_plot <- plot_n |>
  dplyr::ungroup() |>
  dplyr::group_split(index) |>
  furrr::future_map_dfr(~compute_phylo_diversity(.x, tree),
                        .options = furrr::furrr_options(seed = TRUE))

# phylo_plot |>
#   dplyr::left_join(
#     plot_n |>
#       dplyr::select(index, site, period, plot, iter) |>
#       dplyr::distinct()) |>
#   readr::write_csv("faith_plot_posterior.csv")

phylo_plot_summary <- phylo_plot |>
  dplyr::left_join(
    plot_n |>
      dplyr::select(index, site, period, plot, iter) |>
      dplyr::distinct()) |>
  dplyr::group_by(site, period, plot) |>
  dplyr::summarise( faith_plot_mean = mean(n_pd, na.rm = TRUE),
                    faith_plot_sd = sd(n_pd, na.rm = TRUE) )

phylo_site <- site_n |>
  dplyr::ungroup() |>
  dplyr::group_split(index) |>
  furrr::future_map_dfr(~compute_phylo_diversity(.x, tree),
                        .options = furrr::furrr_options(seed = TRUE))

# phylo_site |> 
#   dplyr::left_join(
#     site_n |> 
#       dplyr::select(index, site, period, iter) |> 
#       dplyr::distinct()) |> 
#   readr::write_csv("faith_site_posterior.csv")

phylo_site_summary <- phylo_site |>
  dplyr::left_join(
    site_n |>
      dplyr::select(index, site, period, iter) |>
      dplyr::distinct()) |>
  dplyr::group_by(site, period ) |>
  dplyr::summarise( faith_site_mean = mean(n_pd, na.rm = TRUE),
                    faith_site_sd = sd(n_pd, na.rm = TRUE) )

# option to read in the files with the posteriors
setwd(here::here("results"))

sr_plot <- read_csv("sr_plot_posterior.csv") |> 
  dplyr::group_by(site, plot, period) |> 
  summarise(sr_plot_mean = mean(sr, na.rm = TRUE), 
            sr_plot_sd = sd(sr, na.rm = TRUE))

sr_site <- read_csv("sr_site_posterior.csv") |> 
  dplyr::group_by(site, period) |> 
  summarise(sr_site_mean = mean(sr, na.rm = TRUE), 
            sr_site_sd = sd(sr, na.rm = TRUE))

pd_plot <- read_csv("pd_plot_posterior.csv") |> 
  dplyr::group_by(site, plot, period) |> 
  summarise(pd_plot_mean = mean(pd, na.rm = TRUE), 
            pd_plot_sd = sd(pd, na.rm = TRUE))

pd_site <- read_csv("pd_site_posterior.csv") |> 
  dplyr::group_by(site, period) |> 
  summarise(pd_site_mean = mean(pd, na.rm = TRUE), 
            pd_site_sd = sd(pd, na.rm = TRUE))

n_plot <- read_csv("n_plot_posterior.csv") |> 
  dplyr::group_by(site, plot, period) |> 
  summarise(n_plot_mean = mean(totN, na.rm = TRUE), 
            n_plot_sd = sd(totN, na.rm = TRUE))

n_site <- read_csv("n_site_posterior.csv") |> 
  dplyr::group_by(site, period) |> 
  summarise(n_site_mean = mean(totN, na.rm = TRUE), 
            n_site_sd = sd(totN, na.rm = TRUE))

shan_plot <- read_csv("shan_plot_posterior.csv") |> 
  dplyr::group_by(site, plot, period) |> 
  summarise(shan_plot_mean = mean(shannon, na.rm = TRUE), 
            shan_plot_sd = sd(shannon, na.rm = TRUE))

shan_site <- read_csv("shan_site_posterior.csv") |> 
  dplyr::group_by(site, period) |> 
  summarise(shan_site_mean = mean(shannon, na.rm = TRUE), 
            shan_site_sd = sd(shannon, na.rm = TRUE))

faith_plot <- read_csv("faith_plot_posterior.csv") |> 
  dplyr::group_by(site, plot, period) |> 
  summarise(faith_plot_mean = mean(n_pd, na.rm = TRUE), 
            faith_plot_sd = sd(n_pd, na.rm = TRUE))

faith_site <- read_csv("faith_site_posterior.csv") |> 
  dplyr::group_by(site, period) |> 
  summarise(faith_site_mean = mean(n_pd, na.rm = TRUE), 
            faith_site_sd = sd(n_pd, na.rm = TRUE))

setwd(here::here("data"))

load("neon_cr_data_2024-08-27.RData")
load("neon_mammal_box_trapping_v01.RData")

sample_key <- boxtrap[[5]] |>
  tibble::as_tibble() |> 
  dplyr::select(tagID, bloodSampleID, earSampleID) |> 
  dplyr::distinct() |> 
  dplyr::filter(!is.na(tagID)) |> 
  dplyr::filter(!(is.na(earSampleID) & is.na(bloodSampleID))) |> 
  dplyr::group_by(tagID) |> 
  dplyr::summarise( n_blood = sum(!is.na(bloodSampleID)), 
                    n_ear = sum(!is.na(earSampleID))) |> 
  dplyr::mutate(type = ifelse(n_blood > 0 & n_ear == 0, "blood",
                              ifelse(n_ear > 0 & n_blood == 0, "ear", 
                                     ifelse(n_blood > 0 & n_ear > 0, "both", "whoops")))) |> 
  dplyr::select(tagID, type)

disease_biodiversity <-
  final |>
  dplyr::ungroup() |> 
  dplyr::select(siteID, site, plotID, plot, period, scientificName, sp, tagID, positive, type, path) |> 
  dplyr::filter(!is.na(positive)) |> 
  dplyr::distinct() |> 
  dplyr::left_join(sample_key) |> 
  dplyr::left_join(sr_site) |> 
  dplyr::left_join(sr_plot) |> 
  dplyr::left_join(pd_site) |> 
  dplyr::left_join(pd_plot) |> 
  dplyr::left_join(n_site) |> 
  dplyr::left_join(n_plot) |> 
  dplyr::left_join(shan_site) |> 
  dplyr::left_join(shan_plot) |> 
  dplyr::left_join(faith_site) |> 
  dplyr::left_join(faith_plot) |> 
  dplyr::filter(across(sr_site_mean:faith_plot_sd, function(x)!is.na(x)))
    
setwd(here::here("data"))
write_csv(disease_biodiversity, "disease_with_biodiversity_metrics_v01.csv")