## Repository Directory

### data: Contains data for analyses 
 * [range_maps](./range_maps) Folder with IUCN range maps for the species included in the analysis
 * [NEON_canopy_all_sites.csv](./NEON_canopy_all_sites.csv) Mean % canopy cover (NLCD) within a 500 m buffer of NEON plot (calculated from Google Earth Engine). Has columns for NEON siteID, plotID, and "mean", which is the mean % canopy cover
 * [PanTHERIA_1-0_WR05_Aug2008.txt](./PanTHERIA_1-0_WR05_Aug2008.txt) PanTHERIA database
 * [disease_with_biodiversity_metrics_v01.csv](./disease_with_biodiversity_metrics_v01.csv) Formatted infection presence / host biodiversity data. Metadata information below
   
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

* [neon_cr_data_2024-08-27.RData](./neon_cr_data_2024-08-27.RData) Formatted data ready to go into model; .RData object with four items
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

  * [neon_mammal_box_trapping_v01.RData](./neon_mammal_box_trapping_v01.RData) Small mammal box trapping data; see [NEON documentation](https://data.neonscience.org/data-products/DP1.10072.001) for detail
  * [neon_mammal_sequences_v01.RData](./neon_mammal_sequences_v01.RData) Small mammal DNA sequence data; see [NEON documention](https://data.neonscience.org/data-products/DP1.10076.001) for detail
  * [neon_mammal_tick_pathogen_v01.RData](./neon_mammal_tick_pathogen_v01.RData) Small mammal tick pathogen screening data; see [NEON documentation](https://data.neonscience.org/data-products/DP1.10064.002) for detail
  * [nplots.csv](./nplots.csv) CSV file with number of plots sampled at each site for each period
  * [output.nex](./output.nex) Phylogeny for focal species; downloaded from [VertLife](https://vertlife.org/mammals/)
  * [site_distance_area.csv](./site_distance_area.csv) Information on distances between plots within sites