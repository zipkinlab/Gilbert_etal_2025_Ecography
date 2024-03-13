# Spatial scaling of biodiversity–disease relationships 

### [Neil A. Gilbert](https://gilbertecology.com), [Graziella V. DiRenzo](https://grazielladirenzo.weebly.com/), [Elise F. Zipkin](https://zipkinlab.org/)

### Data/code DOI: To be added upon acceptance

#### Please contact the first author for questions about the code or data: Neil A. Gilbert (neil.allen.gilbert@gmail.com)

__________________________________________________________________________________________________________________________________________

## Abstract

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
 * [13_results_stats.R](./code/13_results_stats.R) Script to grab statistics that are reported in the results

### [data](./data): Contains data for analyses 
 * [range_maps](./data/range_maps) Folder with IUCN range maps for the species included in the analysis
 * [NEON_canopy_all_sites.csv](./data/NEON_canopy_all_sites.csv) Mean percent canopy cover (NLCD) within a 500 m buffer of NEON plot (calculated from Google Earth Engine). Has columns for NEON siteID, plotID, and "mean", which is the mean percent canopy cover
 * [PanTHERIA_1-0_WR05_Aug2008.txt](./data/PanTHERIA_1-0_WR05_Aug2008.txt) PanTHERIA database
 * [disease_with_biodiversity_metrics_v01.csv](./data/disease_with_biodiversity_metrics_v01.csv) Formatted infection presence / host biodiversity data. Metadata information below
   
   | Variable name | Meaning |
    |---------------|---------|
    | siteID | NEON site (text) |
    | site | NEON site (numeric) |
    | plotID | NEON plot (text) |
    | plot | NEON plot (numeric) |
    | period | Trapping bout |
    | scientificName | Scientific name of individual screened for Lyme infection |
    | sp | Species numeric code for individual screened for Lyme infection |
    | tagID | Unique identifier for individual small mammal |
    | positive | Indicates whether or not individual tested positive for Lyme (1) or not (0) |
    | sr_site_mean | Posterior mean of species richness for the entire NEON site |
    | sr_site_sd | Posterior standard deviation of species richness for the entire NEON site |
    | sr_plot_mean | Posterior mean of species richness for NEON plot |
    | sr_plot_sd | Posterior standard deviation of species richness for NEON plot |
    | pd_site_mean | Posterior mean of <i>Peromyscus</i> dominance for the entire NEON site |
    | pd_site_sd | Posterior standard deviation of <i>Peromyscus</i> dominance for the entire NEON site |
    | pd_plot_mean | Posterior mean of <i>Peromyscus</i> dominance for NEON plot |
    | pd_plot_sd | Posterior standard deviation of <i>Peromyscus</i> dominance for NEON plot |
    | n_site_mean | Posterior mean of total abundance for the entire NEON site |
    | n_site_sd | Posterior standard deviation of total abundance for the entire NEON site |
    | n_plot_mean | Posterior mean of total abundance for NEON plot |
    | n_plot_sd | Posterior standard deviation of total abundance for NEON plot |
    | shan_site_mean | Posterior mean of Shannon diversity for the entire NEON site |
    | shan_site_sd | Posterior standard deviation of Shannon diversity for the entire NEON site |
    | shan_plot_mean | Posterior mean of Shannon diversity for NEON plot |
    | shan_plot_sd | Posterior standard deviation of Shannon diversity for NEON plot |
   
