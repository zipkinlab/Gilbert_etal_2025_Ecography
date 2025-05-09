## Repository Directory

### [code](./code): Contains code for downloading, formatting, analyzing, and visualizing data
 * [01_download_neon_data.R](./code/01_download_neon_data.R) Script to download relevant NEON data products
 * [02_format_capture_recapture_data.R](./code/02_format_capture_recapture_data.R) Prepare small mammal trapping data for model
 * [03_run_capture_recapture_model.R](./code/03_run_capture_recapture_model.R) Run the multispecies, spatially-stratified capture-recapture (MSSCR) model
 * [04_wrangle_biodiversity_metrics.R](./code/04_wrangle_biodiversity_metrics.R) Script to derive biodiversity metrics from MSSCR model posterior
 * [05_disease_sr_site_model.R](./code/05_disease_sr_site_model.R) Infection prevalence predicted by species richness at metacommunity level
 * [06_disease_sr_plot_model.R](./code/06_disease_sr_plot_model.R) Infection prevalence predicted by species richness at local community level
 * [07_disease_pd_plot_model.R](./code/07_disease_pd_plot_model.R) Infection prevalence predicted by Peromyscus dominance at local community level
 * [08_disease_pd_site_model.R](./code/08_disease_pd_site_model.R) Infection prevalence predicted by Peromyscus dominance at metacommunity level
 * [09_disease_n_site_model.R](./code/09_disease_n_site_model.R) Infection prevalence predicted by total abundance at the metacommunity level
 * [10_disease_n_plot_model.R](./code/10_disease_n_plot_model.R) Infection prevalence predicted by total abundance at the local community level
 * [11_disease_shan_plot_model.R](./code/11_disease_shan_plot_model.R) Infection prevalence predicted by Shannon diversity at local community level
 * [12_disease_shan_site_model.R](./code/12_disease_shan_site_model.R) Infection prevalence predicted by Shannon diversity at metacommunity level
 * [13_disease_faith_site_model.R](./code/13_disease_faith_site_model.R) Infection prevalence predicted by phylo diversity at metacommunity level
 * [14_disease_faith_plot_model.R](./code/14_disease_faith_plot_model.R) Infection prevalence predicted by phylo diversity at local community level
 * [15_phylogenetic_regressions.R](./code/15_phylogenetic_regressions.R) Fit versions of disease models that account for phylogenetic dependence
 * [16_check_influence_of_method.R](./code/16_check_influence_of_method.R) Check if sample type (blood, ear, or both) influences results
 * [17_calculate_distances_areas_sites.R](./code/17_calculate_distances_areas_sites.R) Calculate plot spacing within sites
 * [18_check_influence_of_plot_spacing.R](./code/18_check_influence_of_plot_spacing.R) Run model that includes variable for plot spacing
 * [19_check_influence_nplot.R](./code/19_check_influence_nplot.R) Run model that includes variable for number of plots sampled
 * [20_biodiversity_metric_correlation.R](./code/20_biodiversity_metric_correlation.R) Investigate collinearity of biodiversity metrics
