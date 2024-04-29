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

### [data](./data): Contains data for analyses 
 * [range_maps](./data/range_maps) Folder with IUCN range maps for the species included in the analysis
 * [NEON_canopy_all_sites.csv](./data/NEON_canopy_all_sites.csv) Mean % canopy cover (NLCD) within a 500 m buffer of NEON plot (calculated from Google Earth Engine). Has columns for NEON siteID, plotID, and "mean", which is the mean % canopy cover
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

* [neon_cr_data_2024-03-29.RData](./data/neon_cr_data_2024-03-29.RData) Formatted data ready to go into model; .RData object with four items
   * **constants** Constants for model
     
    | Variable name | Meaning |
    |---------------|---------|
    | nsp | Number of species |
    | nsite | Number of NEON sites |
    | max_plot | Matrix reporting the maximum number of plots trapped for a site-bout combination |
    | M | Total number of individuals (real and augmented) for each site-species combination |
    | max_period | Vector reporting the maximum trapping bout trapped for each site |
    | nz | Number of "latent existence states" |
    | species_z | Species identifier used for indexing in z loop |
    | site_z | Site identifier used for indexing in z loop |
    | period_z | Period identifier used for indexing in z loop |
    | max_plot_z | Identifier for maximum plot for indexing in z loop |
    | ny | Length of detection/nondetection data loop |
    | species_y | Species identifier used for indexing in y loop |
    | z_index | Indexes latent Z state within the y loop |
  
   * **data** Data for model

    | Variable name | Meaning |
    |---------------|---------|
    | y | Vector of detection-nondetection data (capture history) |
    | TMIN | Minimum temperature on day of trapping (z-standardized) |
    | PRCP | Precipitation on the day of trapping (z-standardized) |
    | in_range | Matrix reporting whether (1) or not (0) a NEON site falls within the range of a given species |
    | plot | Plot identifier for where an individual was trapped (NA for augmented individuals) |
    | CANOPY | Matrix with mean percent canopy cover within a 500 m radius of the center of the NEON plot (z-standardized) |

   * **final** Dataframe from which some of the data and constants are pulled
 
    | Variable name | Meaning |
    |---------------|---------|
    | siteID |  NEON site (character) |
    | period | Trapping bout identifier (numeric) |
    | rep | Replicate trapping night within bout (numeric |
    | scientificName | Scientific name of individual |
    | ind | Identifier for individual (numeric) |
    | plotID | NEON plot (character) |
    | y | Vector of detection-nondetection data (capture history) |
    | year | Year of trapping |
    | month | Month of trapping |
    | day | Day of trapping |
    | prcp | Precipitation (mm) on the day of trapping |
    | tmin | Minimum temperature (C) on day of trapping |
    | plots_sampled | list-column reporting which plots were sampled during the given period |
    | plot | Plot identifier for where individual was trapped (NA for augmented individuals |
    | site | Site identifier (numeric) |
    | sp | Species index (numeric) |
    | plotst | Starting value for plot (known for trapped individuals; randomly assigned for augmented |
    | z_index | Index for the latent-z state for an individual-bout combination |
    | max_plot | Maximum plot trapped during a given bout |
    | in_range | Binary variable reporting whether (1) or not (0) a NEON site falls within the range of a given species |

   * **get_z_index** Dataframe from which ssome of the data and constants are pulled. Simplified version of **final** without replicate information

    | Variable name | Meaning |
    |---------------|---------|
    | siteID |  NEON site (character) |
    | site | NEON site index (numeric) |
    | period | Trapping bout identifier (numeric) |
    | sp | Species index (numeric) |
    | scientificName | Scientific name of individual |
    | ind | Identifier for individual (numeric) |
    | plotID | NEON plot (character) |
    | plot | Plot identifier for where individual was trapped (NA for augmented individuals |
    | plots_sampled | list-column reporting which plots were sampled during the given period |
    | plotst | Starting value for plot (known for trapped individuals; randomly assigned for augmented |
    | z_index | Index for the latent-z state for an individual-bout combination |
    | max_plot | Maximum plot trapped during a given bout |

  * [neon_mammal_box_trapping_v01.RData](./data/neon_mammal_box_trapping_v01.RData) Small mammal box trapping data; see [NEON documentation](https://data.neonscience.org/data-products/DP1.10072.001) for detail
  * [neon_mammal_sequences_v01.RData](./data/neon_mammal_sequences_v01.RData) Small mammal DNA sequence data; see [NEON documention](https://data.neonscience.org/data-products/DP1.10076.001) for detail
  * [neon_mammal_tick_pathogen_v01.RData](./data/neon_mammal_tick_pathogen_v01.RData) Small mammal tick pathogen screening data; see [NEON documentation](https://data.neonscience.org/data-products/DP1.10064.002) for detail

### [figures](./figures): Contains figures and code for creating figures

* [code_for_figures](./figures/code_for_figures) Folder with scripts to create figures
   * [01_figure_01.R](./figures/code_for_figures/01_figure_01.R) Create figure 1 (Map of NEON site/plot sampling structure)
   * [02_figure_02.R](./figures/code_for_figures/02_figure_02.R) Create Figure 2 (Map of biodiversity metrics) 
   * [03_figure_03_04.R](./figures/code_for_figures/03_figure_03_04.R) Create Figures 3 (coefficient plot) and 4 (marginal effects)
   * [04_figure_05.R](./figures/code_for_figures/04_figure_05.R) Create Figure 5 (infection % map)
   * [05_table_s1.R](./figures/code_for_figures/05_table_s1.R) Create Table S1      
* [figure_01.png](./figures/figure_01.png) Figure 1 (PNG)
* [figure_01.pptx](./figures/figure_01.pptx) Figure 1 (PPT file for some post hoc editing)
* [figure_02.png](./figures/figure_02.png) Figure 2
* [figure_03.png](./figures/figure_03.png) Figure 3
* [figure_04.png](./figures/figure_04.png) Figure 4
* [figure_05.png](./figures/figure_05.png) Figure 5

### [results](./results): Contains results files
* neon_capture_recapture_results_2024-04-11.RData. Output from MSSCR model (Step 1). File too large for GitHub: [download link](https://1drv.ms/u/s!AtvYBfNq7AMkhIIjVUn6z4ooDwt1Kw?e=iUlju0)
* [rodent_pathogen_n_plot_2024-04-18.RData](./results/rodent_pathogen_n_plot_2024-04-18.RData) Results from Step 2 model (total abundance, local community)
* [rodent_pathogen_n_site_2024-04-18.RData](./results/rodent_pathogen_n_site_2024-04-18.RData) Results from Step 2 model (total abundance, metacommunity)
* [rodent_pathogen_pd_plot_2024-04-19.RData](./results/rodent_pathogen_pd_plot_2024-04-19.RData) Results from Step 2 model (<i>Peromyscus</i> dominance, local community)
* [rodent_pathogen_pd_site_2024-04-19.RData](./results/rodent_pathogen_pd_site_2024-04-19.RData) Results from Step 2 model (<i>Peromyscus</i> dominance, metacommunity)
* [rodent_pathogen_shan_plot_2024-04-18.RData](./results/rodent_pathogen_shan_plot_2024-04-18.RData) Results from Step 2 model (Shannon diversity, local community)
* [rodent_pathogen_shan_site_2024-04-18.RData](./results/rodent_pathogen_shan_site_2024-04-18.RData) Results from Step 2 model (Shannon diversity, metacommunity)
* [rodent_pathogen_sr_plot_2024-04-18.RData](./results/rodent_pathogen_sr_plot_2024-04-18.RData) Results from Step 2 model (species richness, local community)
* [rodent_pathogen_sr_site_2024-04-18.RData](./results/rodent_pathogen_sr_site_2024-04-18.RData) Resutls from Step 2 model (species richness, metacommunity)
