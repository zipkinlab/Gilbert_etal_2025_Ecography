# Code for: Idiosyncratic spatial scaling of biodiversity–disease relationships 
__________________________________________________________________________________________________________________________________________

### Authors: Neil A. Gilbert, Graziella V. DiRenzo, Elise F. Zipkin
#### Please contact the first author for questions about the code or data: Neil A. Gilbert (neil.allen.gilbert@gmail.com)
__________________________________________________________________________________________________________________________________________

### Information

Repository Type: Program R scripts

Year of Origin: 2024

Year of Version: 2024

Version: 1.0.0

Digital Object Identifier (DOI): https://doi.org/10.5066/P13P2VWG

USGS Information Product Data System (IPDS) no.: IP-165614

__________________________________________________________________________________________________________________________________________
### Suggested Citation for Software

Gilbert, N.A., G. V. DiRenzo, & E. F. Zipkin. Code for idiosyncratic spatial scaling of biodiversity–disease relationships. Version 1.0.0; U.S. Geological Survey software release. Reston, VA. https://doi.org/10.5066/P13P2VWG
__________________________________________________________________________________________________________________________________________

### Repository description

This repository contains data and code supporting the Gilbert et. al (in press). The objective of the paper is to evaluate how estimated relationships between Lyme disease risk and host biodiversity change when host biodiversity is quantified at different spatial scales. Data come from the National Ecological Observatory Network's small mammal monitoring program. Data are freely available at the [NEON Data Portal](https://data.neonscience.org/home) and are provided in this repository for ease of reproducing analyses. The code (see Repository Directory below) includes scripts to format data, fit models, interpret results, and create figures. 

**Manuscript citation**.
Gilbert, N.A., G.V. DiRenzo, and E.F. Zipkin. In press. Idiosyncratic spatial scaling of biodiversity–disease relationships. _Ecography_
__________________________________________________________________________________________________________________________________________

## Manuscript Abstract
Host biodiversity is hypothesized to dilute the risk of vector-borne diseases because many host species are “dead ends” which cannot effectively transmit the disease to other individuals. However, many studies on biodiversity–disease relationships characterize host biodiversity at single, local spatial scales, which complicates efforts to forecast disease risk if associations between host biodiversity and disease change with spatial scale. Here, our objective was to evaluate the spatial scaling of relationships between host biodiversity and Borrelia infection prevalence in small mammals; Borrelia is a bacterial taxon which causes Lyme disease, a zoonotic disease affecting >400,000 people in the United States annually. We compared the associations between infection prevalence and small mammal host diversity for local communities (individual plots) and metacommunities (multiple plots aggregated within a landscape) sampled by the National Ecological Observatory Network, an emerging continental-scale environmental monitoring program with a hierarchical sampling design. We applied a multispecies, spatially-stratified capture–recapture model to a trapping dataset to estimate five small mammal biodiversity metrics, which we used to predict infection status for a subset of trapped individuals. We found that relationships between Borrelia infection prevalence and biodiversity did indeed vary when biodiversity was quantified at the different spatial scales, and that these scaling behaviors were idiosyncratic among the five biodiversity metrics. For example, species richness of local communities showed a negative (dilution) effect on infection prevalence, while species richness of the small mammal metacommunity showed a positive (amplification) effect on infection prevalence. In contrast, total abundance of the local community showed a negative (dilution) effect on infection prevalence, while total abundance of the metacommunity showed no association with infection prevalence. Our modeling approach can inform future analyses as data from similar monitoring programs accumulate and become increasingly available through time. Our results indicate that a focus on single spatial scales when assessing the influence of biodiversity on disease risk provides an incomplete picture of the complexity of disease dynamics in ecosystems. 

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
    | scientificName | Scientific name of individual screened for <i>Borrelia</i> infection |
    | sp | Species numeric code for individual screened for <i>Borrelia</i> infection |
    | tagID | Unique identifier for individual small mammal |
    | positive | Indicates whether or not individual tested positive for <i>Borrelia</i> (1) or not (0) |
    | type | Indicates sample type (ear tissue, blood, or both) |
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
    | faith_site_mean | Posterior mean of phylogenetic diversity (abundance-weighted Faith index) for the entire NEON site |
    | faith_site_sd | Posterior standard deviation of phylogenetic diversity (abundance-weighted Faith index) for the entire NEON site |
    | faith_plot_mean | Posterior mean of phylogenetic diversity (abundance-weighted Faith index)  for NEON plot |
    | faith_plot_sd | Posterior standard deviation of phylogenetic diversity (abundance-weighted Faith index) for NEON plot |

* [neon_cr_data_2024-08-27.RData](./data/neon_cr_data_2024-08-27.RData) Formatted data ready to go into model; .RData object with four items
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
    | tagID | Identifier for individual (NA for augmented individuals) |
    | y | Vector of detection-nondetection data (capture history) |
    | positive | Indicates whether <i>Borrelia</i> infection was detected |
    | type | Indicates sample type (ear, blood, or both) |
    | path | Indicates which <i>Borrelia</i> taxon was detected (<i>Borrelia</i> spp. only or <i>Borrelia</i> spp. + <i>Borrelia burgdorferi</i>) |
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
  * [nplots.csv](./data/nplots.csv) CSV file with number of plots sampled at each site for each period
  * [output.nex](./data/output.nex) Phylogeny for focal species; downloaded from [VertLife](https://vertlife.org/data/mammals/)
  * [site_distance_area.csv](./data/site_distance_area.csv) Information on distances between plots within sites

### [figures](./figures): Contains figures and code for creating figures

* [code_for_figures](./figures/code_for_figures) Folder with scripts to create figures
   * [01_figure_01.R](./figures/code_for_figures/01_figure_01.R) Create figure 1 (Map of NEON site/plot sampling structure)
   * [02_figure_02.R](./figures/code_for_figures/02_figure_02.R) Create Figure 2 (Map of biodiversity metrics) 
   * [03_figure_03_04.R](./figures/code_for_figures/03_figure_03_04.R) Create Figures 3 (coefficient plot) and 4 (marginal effects)
   * [04_figure_05.R](./figures/code_for_figures/04_figure_05.R) Create Figure 5 (infection % map)
   * [05_table_s1.R](./figures/code_for_figures/05_table_s1.R) Create Table S1
   * [06_video_01.R](./figures/code_for_figures/06_video_01.R) Create Supplemental Video
   * [07_figure_s01.R](./figures/code_for_figures/07_figure_s01.R) Create Fig. S1
   * [08_figure_s02.R](./figures/code_for_figures/08_figure_s02.R) Creat Fig. S2   
* [figure_01.png](./figures/figure_01.png) Figure 1 (PNG)
* [figure_01.pptx](./figures/figure_01.pptx) Figure 1 (PPT file for some post hoc editing)
* [figure_02.png](./figures/figure_02.png) Figure 2
* [figure_03.png](./figures/figure_03.png) Figure 3
* [figure_04.png](./figures/figure_04.png) Figure 4
* [figure_05.png](./figures/figure_05.png) Figure 5
* [figure_s01.png](./figures/figure_s01.png) Figure S1
* [figure_s01.pptx](./figures/figure_s01.pptx) Figure S1 (PPT file for post hoc editing)
* [figure_s02.png](./figures/figure_s02.png) Figure S2
* [figure_s03.png](./figures/figure_s03.png) Figure S3
* [figure_s04.png](./figures/figure_s04.png) Figure S4
* [figure_s05.png](./figures/figure_s05.png) Figure S5
* [figure_s06.png](./figures/figure_s06.png) Figure S6
* [figure_s07.png](./figures/figure_s07.png) Figure S7
* [figure_s08.png](./figures/figure_s08.png) Figure S8
* [figure_s09.png](./figures/figure_s09.png) Figure S9
* [figure_s10.png](./figures/figure_s10.png) Figure S10
* [video_s01.mp4](./figures/video_s01.mp4) Supplemental video of biodiversity metrics

### [results](./results): Contains results files
* neon_capture_recapture_results_2024-04-11.RData. Output from MSSCR model (Step 1). File too large for GitHub: [download link](https://1drv.ms/u/s!AtvYBfNq7AMkhIIjVUn6z4ooDwt1Kw?e=iUlju0)
* [faith_plot_phylo.rds](./results/faith_plot_phylo.rds) Results from alternative Step 2 model (phylogenetic diversity, local community) that includes phylogenetic dependence
* [faith_plot_posterior.csv](./results/faith_plot_posterior.csv) CSV with posterior distribution of plot-level phylogenetic diversity
* [faith_site_phylo.rds](./results/faith_site_phylo.rds) Results from alternative Step 2 model (phylogenetic diversity, metacommunity) that includes phylogenetic dependence
* [faith_site_posterior.csv](./results/faith_site_posterior.csv) CSV with posterior distribution of site-level phylogenetic diversity
* [n_plot_phylo.rds](./results/n_plot_phylo.rds) Results from alternative Step 2 model (total abundance, local community) that includes phylogenetic dependence
* [n_plot_posterior.csv](./results/n_plot_posterior.csv) CSV with posterior distribution of plot-level total abundance
* [n_site_phylo.rds](./results/n_site_phylo.rds) Results from alternative Step 2 model (total abundance, metacommunity) that includes phylogenetic dependence
* [n_site_posterior.csv](./results/n_site_posterior.csv) CSV with posterior distribution of site-level total abundance
* [pd_plot_phylo.rds](./results/pd_plot_phylo.rds)  Results from alternative Step 2 model (<i>Peromyscus</i> dominance, local community) that includes phylogenetic dependence
* [pd_plot_posterior.csv](./results/pd_plot_posterior.csv) CSV with posterior distribution of plot-level <i>Peromyscus</i> dominance
* [pd_site_phylo.rds](./results/pd_site_phylo.rds)  Results from alternative Step 2 model (<i>Peromyscus</i> dominance, metacommunity) that includes phylogenetic dependence
* [pd_site_posterior.csv](./results/pd_site_posterior.csv) CSV with posterior distribution of site-level <i>Peromyscus</i> dominance
* [rodent_pathogen_faith_plot_2024-08-21.RData](./results/rodent_pathogen_faith_plot_2024-08-21.RData) Results from Step 2 model (phylogenetic diversity, local community)
* [rodent_pathogen_faith_site_2024-09-16.RData](./results/rodent_pathogen_faith_site_2024-09-16.RData) Results from Step 2 model (phylogenetic diversity, metacommunity)
* [rodent_pathogen_n_plot_2024-04-18.RData](./results/rodent_pathogen_n_plot_2024-04-18.RData) Results from Step 2 model (total abundance, local community)
* [rodent_pathogen_n_site_2024-09-16.RData](./results/rodent_pathogen_n_site_2024-09-16.RData) Results from Step 2 model (total abundance, metacommunity)
* [rodent_pathogen_pd_plot_2024-04-19.RData](./results/rodent_pathogen_pd_plot_2024-04-19.RData) Results from Step 2 model (<i>Peromyscus</i> dominance, local community)
* [rodent_pathogen_pd_site_2024-09-16.RData](./results/rodent_pathogen_pd_site_2024-09-16.RData) Results from Step 2 model (<i>Peromyscus</i> dominance, metacommunity)
* [rodent_pathogen_shan_plot_2024-04-18.RData](./results/rodent_pathogen_shan_plot_2024-04-18.RData) Results from Step 2 model (Shannon diversity, local community)
* [rodent_pathogen_shan_site_2024-09-16.RData](./results/rodent_pathogen_shan_site_2024-09-16.RData) Results from Step 2 model (Shannon diversity, metacommunity)
* [rodent_pathogen_sr_plot_2024-04-18.RData](./results/rodent_pathogen_sr_plot_2024-04-18.RData) Results from Step 2 model (species richness, local community)
* [rodent_pathogen_sr_site_2024-09-16.RData](./results/rodent_pathogen_sr_site_2024-09-16.RData) Resutls from Step 2 model (species richness, metacommunity)
* [shan_plot_phylo.rds](./results/shan_plot_phylo.rds)  Results from alternative Step 2 model (Shannon diversity, local community) that includes phylogenetic dependence
* [shan_plot_posterior.csv](./results/shan_plot_posterior.csv) CSV with posterior distribution of plot-level Shannon diversity
* [shan_site_phylo.rds](./results/shan_site_phylo.rds)  Results from alternative Step 2 model (Shannon diversity, metacommunity) that includes phylogenetic dependence
* [shan_site_posterior.csv](./results/shan_site_posterior.csv) CSV with posterior distribution of site-level Shannon diversity
* [sr_plot_phylo.rds](./results/sr_plot_phylo.rds)  Results from alternative Step 2 model (species richness, local community) that includes phylogenetic dependence
* [sr_plot_posterior.csv](./results/sr_plot_posterior.csv) CSV with posterior distribution of plot-level species richness
* [sr_site_phylo.rds](./results/sr_site_phylo.rds)  Results from alternative Step 2 model (species richness, metacommunity) that includes phylogenetic dependence
* [sr_site_posterior.csv](./results/sr_site_posterior.csv) CSV with posterior distribution of site-level species richness
