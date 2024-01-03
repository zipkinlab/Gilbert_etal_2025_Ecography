library(here)
library(tidyverse)
library(lubridate)
library(janitor)
library(daymetr)

setwd(here::here("data"))

canopy <- readr::read_csv("NEON_canopy_all_sites.csv") |> 
  dplyr::select(siteID, plotID, pcanopy_500m = mean)

load("neon_mammal_tick_pathogen_v01.RData")
load("neon_mammal_box_trapping_v01.RData")
load("neon_mammal_sequences_v01.RData")

pantheria <- utils::read.delim("PanTHERIA_1-0_WR05_Aug2008.txt")  |>  
  dplyr::rename( genus = MSW05_Genus,
                 order = MSW05_Order,
                 family = MSW05_Family,
                 sp = MSW05_Species,
                 mass = X5.1_AdultBodyMass_g ) |> 
  tibble::as_tibble() |> 
  dplyr::select( order, family, genus, sp, mass ) |> 
  dplyr::mutate( scientificName = base::paste( genus, sp ) )

# which month-year combinations for sites have associated pathogen screening?
path_surveys <- rodent_pathogen$rpt2_pathogentesting |> 
  tibble::as_tibble() |> 
  dplyr::filter(siteID %in% c("BART", "BLAN", "CLBJ", "DCFS", "DELA", "DSNY", "GRSM", "HARV", 
                              "JERC", "KONA", "KONZ", "LENO", "MLBS", "NOGP", "OAES", "ORNL", 
                              "OSBS", "SCBI", "SERC", "STEI", "TALL", "TREE", "UKFS", "UNDE", 
                              "WOOD")) |>
  dplyr::select(siteID, plotID, collectDate) |> 
  dplyr::mutate(collectDate = as.Date(collectDate)) |> 
  dplyr::distinct() |> 
  dplyr::rename( path_collect_date = collectDate )

# which month-year combinations for sites have box trapping data
mam_surveys <- boxtrap$mam_pertrapnight |> 
  tibble::as_tibble() |> 
  dplyr::filter(siteID %in% c("BART", "BLAN", "CLBJ", "DCFS", "DELA", "DSNY", "GRSM", "HARV", 
                              "JERC", "KONA", "KONZ", "LENO", "MLBS", "NOGP", "OAES", "ORNL", 
                              "OSBS", "SCBI", "SERC", "STEI", "TALL", "TREE", "UKFS", "UNDE", 
                              "WOOD")) |> 
  dplyr::select(siteID, plotID, collectDate) |> 
  dplyr::mutate(collectDate = as.Date(collectDate)) |> 
  dplyr::distinct() |> 
  dplyr::rename(box_collect_date = collectDate)

# master list of site x plot x sampling dates assigned to site-level periods
survey_key <- full_join( path_surveys, mam_surveys ) |> 
  dplyr::mutate( diff = abs(as.numeric(box_collect_date - path_collect_date))) |> 
  dplyr::filter( diff <= 10) |>
  dplyr::select(siteID, box_collect_date ) |> 
  dplyr::distinct() |> 
  dplyr::arrange(siteID, box_collect_date ) |> 
  dplyr::group_by( siteID ) |> 
  dplyr::mutate( dlag = base::as.numeric( box_collect_date - dplyr::lag(box_collect_date)),
                 dlag = tidyr::replace_na(dlag, 0), 
                 period = base::cumsum(dlag > 7) + 1) |> 
  dplyr::select(-dlag) |> 
  dplyr::full_join(mam_surveys) |> 
  dplyr::full_join( path_surveys ) |> 
  dplyr::mutate( diff = abs(as.numeric(box_collect_date - path_collect_date))) |> 
  dplyr::filter( diff <= 10 ) |> 
  dplyr::select( siteID, plotID, period, box_collect_date, path_collect_date )

site_dates <- survey_key |> 
  tidyr::pivot_longer( box_collect_date:path_collect_date, names_to = "name", values_to = "collectDate") |> 
  dplyr::select(-name) |> 
  dplyr::distinct() |> 
  dplyr::ungroup() |> 
  dplyr::summarise(min_date = lubridate::year(min(collectDate)),
                   max_date = lubridate::year(max(collectDate)))

site_coords <- survey_key |> 
  dplyr::select(siteID, plotID, collectDate = box_collect_date) |> 
  dplyr::distinct() |> 
  dplyr::left_join(boxtrap[[5]]) |> 
  dplyr::select(siteID, plotID, xcoord = decimalLongitude, ycoord = decimalLatitude) |> 
  dplyr::distinct() |> 
  dplyr::summarise( xcoord = mean(xcoord), 
                    ycoord = mean(ycoord))

df <- list(list())
for(i in 1:nrow(site_coords)){
  df[[i]] <- daymetr::download_daymet(
    site = site_coords[[i, "siteID"]],
    lat = site_coords[[i, "ycoord"]],
    lon = site_coords[[i, "xcoord"]],
    start = site_dates[[1, "min_date"]],
    end = site_dates[[1, "max_date"]],
    path = tempdir(),
    internal = TRUE, 
    silent = FALSE, 
    force = FALSE, 
    simplify = FALSE)$data |> 
    tibble::as_tibble() |>
    janitor::clean_names() |> 
    tibble::add_column( siteID = site_coords[[i, "siteID"]] ) |> 
    dplyr::mutate( date = as.Date( paste( yday, year ), format = "%j %Y")) |> 
    dplyr::select( siteID, year, yday, date, dayl_s:vp_pa )
}

detection_variables <- dplyr::bind_rows(df) |> 
  dplyr::right_join(
    survey_key |> 
      dplyr::select(siteID, date = box_collect_date) |> 
      dplyr::distinct()) |> 
  dplyr::select(siteID, collectDate = date, prcp = prcp_mm_day, tmin = tmin_deg_c)

setwd(here::here("data/range_maps"))

range_filenames <- list.files(pattern = ".shp")
res <- list(list())

site_shp <- site_coords |> 
  sf::st_as_sf(coords = c("xcoord", "ycoord"), 
               crs = 4326)

for( i in 1:length(range_filenames)){
  
  if(range_filenames[i] == "mus_musculus.shp") next() # issue with house mouse shapefile...skip
  
  shp <- sf::st_read(range_filenames[i])  
  
  res[[i]] <- sf::st_filter(site_shp, shp) |> 
    sf::st_drop_geometry() |> 
    tibble::add_column(in_range = 1) |> 
    dplyr::full_join(
      sf::st_drop_geometry(site_shp)) |>
    dplyr::mutate(in_range = ifelse(is.na(in_range), 0, in_range)) |> 
    tibble::add_column(scientific_name = range_filenames[i]) |> 
    dplyr::mutate(scientific_name = str_remove_all(scientific_name, ".shp")) |> 
    dplyr::mutate(scientific_name = str_replace_all(scientific_name, "_", " ")) |> 
    dplyr::mutate(scientific_name = str_to_sentence(scientific_name)) |> 
    dplyr::select(scientific_name, siteID, in_range)
  
  print(paste("Finished", i, "of", length(range_filenames)))
}

# Assume that Mus musculus and Rattus rattus ranges overlap with all NEON sites
site_in_range <- dplyr::bind_rows(res) |> 
  dplyr::full_join(
    tibble::tibble(
      scientific_name = "Mus musculus",
      siteID = site_shp$siteID, 
      in_range = 1)) |> 
  dplyr::full_join(
    tibble::tibble(
      scientific_name = "Rattus rattus",
      siteID = site_shp$siteID, 
      in_range = 1)) |> 
  dplyr::arrange(scientific_name, siteID) |> 
  # update a few cases where the species was detected at sites just outside IUCN range polygon
  # also correct a couple types / taxonomic mismatches
  dplyr::mutate( in_range = ifelse( scientific_name == "Zapus princeps" & siteID == "NOGP", 1, in_range)) |> 
  dplyr::mutate( in_range = ifelse( scientific_name == "Glaucomys volans" & siteID %in% c("STEI", "TREE", "UNDE"), 1, in_range)) |> 
  dplyr::mutate( in_range = ifelse( scientific_name == "Peromyscus leucopus" & siteID == "UNDE", 1, in_range),
                 in_range = ifelse(scientific_name == "Peromyscus maniculatus" & siteID == "SERC", 1, in_range),
                 scientific_name = ifelse(scientific_name == "Microtus pennslyvanicus", "Microtus pennsylvanicus", scientific_name),
                 scientific_name = ifelse(scientific_name == "Podonmys floridanus", "Podomys floridanus", scientific_name),
                 scientific_name = ifelse(scientific_name == "Neotamias minimus", "Tamias minimus", scientific_name)) |> 
  dplyr::rename(scientificName = scientific_name)

lyme_status <- survey_key |> 
  dplyr::select(siteID, plotID, period, collectDate = path_collect_date ) |> 
  dplyr::distinct() |> 
  dplyr::left_join( rodent_pathogen[[5]] ) |> 
  filter( grepl("Borrelia", testPathogenName)) |> 
  dplyr::select(siteID, plotID, collectDate, sampleID, testPathogenName, testResult) |> 
  dplyr::mutate( sample_type = ifelse( grepl(".E$", sampleID), "earSampleID", 
                                       ifelse(grepl(".B$", sampleID), "bloodSampleID", "whoops"))) |>
  dplyr::group_by(siteID, plotID, collectDate, sampleID, sample_type) |>
  dplyr::summarise( positive = sum(testResult == "Positive")) |> 
  dplyr::mutate( positive = ifelse(positive > 0, 1, 0)) |> 
  dplyr::mutate( individual = stringr::str_remove(sampleID, ".B$")) |> 
  dplyr::mutate( individual = stringr::str_remove(individual, ".E$")) |> 
  tidyr::pivot_wider(names_from = sample_type, 
                     values_from = sampleID)

captures <- survey_key |> 
  dplyr::select(siteID, plotID, period, collectDate = box_collect_date) |> 
  dplyr::distinct() |>
  dplyr::group_by(siteID, plotID, period) |> 
  dplyr::mutate(rep = dplyr::row_number(),
                nreps = base::max(rep)) |> 
  dplyr::left_join(boxtrap$mam_pertrapnight) |> 
  dplyr::filter( trapStatus == "5 - capture") |>
  dplyr::full_join(lyme_status) |> 
  dplyr::left_join(
    dplyr::select(mamseq[[6]], earSampleID = sampleID, sp_seq = species)) |> 
  dplyr::mutate( scientificName = ifelse( !is.na(sp_seq), sp_seq, scientificName ) ) |> 
  dplyr::left_join( pantheria )  |>  
  dplyr::filter( ! base::grepl("sp.", scientificName ) ) |>
  dplyr::filter( ! base::grepl("/", scientificName ) ) |>
  dplyr::filter( order == "Rodentia" | order == "Soricomorpha" ) |>
  dplyr::select(siteID, plotID, period, rep, nreps, scientificName, tagID, positive ) |> 
  tibble::add_column( y = 1 ) |> 
  tidyr::pivot_wider( names_from = rep, values_from = c(y, positive) ) |> 
  dplyr::mutate( y_1 = base::ifelse( is.na( y_1 ) & ( ! is.na( y_2 ) | ! is.na( y_3 ) ), 0, y_1 ),
                 y_2 = base::ifelse( is.na( y_2 ) & nreps > 1, 0, y_2 ),
                 y_3 = base::ifelse( is.na( y_3 ) & nreps > 2, 0, y_3 ),
                 positive = ifelse( !is.na(positive_1), positive_1, 
                                    ifelse( !is.na(positive_2), positive_2, 
                                            ifelse(!is.na(positive_3), positive_3, NA)))) |> 
  dplyr::select(-positive_1, -positive_2, -positive_3) |> 
  dplyr::arrange( siteID, plotID, period, scientificName ) |> 
  dplyr::group_by( siteID, period, scientificName ) |> 
  dplyr::mutate( ind = dplyr::row_number()) |> 
  tidyr::pivot_longer(c(y_1, y_2, y_3), 
                      names_to = "rep", 
                      values_to = "y") |> 
  dplyr::mutate(rep = parse_number(rep)) |> 
  dplyr::filter(rep <= nreps ) |> 
  dplyr::left_join(
    survey_key |>
      dplyr::select(siteID, plotID, period, collectDate = box_collect_date) |>
      dplyr::distinct() |>
      dplyr::group_by(siteID, plotID, period) |>
      dplyr::mutate(rep = dplyr::row_number()))

# figure out "super population size" for data augmentation
figure_m <- captures |>
  dplyr::filter(!is.na(y)) |> 
  dplyr::group_by(siteID, scientificName) |> 
  dplyr::summarise( max_ind = max(ind))   |> 
  dplyr::mutate( M = max_ind + 100 )

data_expand <- list( list() )
for( i in 1:base::nrow(figure_m)){
  data_expand[[i]] <- tibble::tibble(
    scientificName = base::rep( figure_m[[i, "scientificName"]], figure_m[[i, "M"]]),
    siteID = base::rep( figure_m[[i, "siteID"]], figure_m[[i, "M"]]),
    ind = 1:figure_m[[i, "M"]] )
}  

da_captures <- dplyr::bind_rows( data_expand ) |> 
  dplyr::full_join(
    captures |> 
      dplyr::filter(!is.na(y)) |> 
      dplyr::select(siteID, period, rep, collectDate, scientificName) |> 
      dplyr::distinct() |> 
      dplyr::group_by(siteID, period, rep, scientificName) |> 
      dplyr::summarise(collectDate = median(collectDate))) |>
  dplyr::arrange( siteID, period, rep, scientificName, ind ) |> 
  dplyr::full_join( captures |> 
                      dplyr::select(-collectDate)) |> 
  dplyr::mutate( y = tidyr::replace_na( y, 0 ) ) |> 
  dplyr::left_join(detection_variables)

get_z_index <-
  da_captures |> 
  dplyr::select(siteID, plotID, period, scientificName, ind, positive ) |> 
  dplyr::distinct() |> 
  ( function(x) tibble::add_column(x, z_index = 1:base::nrow(x)))() |> 
  dplyr::full_join(
    captures |> 
      dplyr::ungroup() |> 
      dplyr::select(siteID, plotID) |> 
      dplyr::distinct() |> 
      dplyr::arrange(plotID) |> 
      dplyr::group_by(siteID) |> 
      dplyr::mutate(plot = dplyr::row_number()) |> 
      dplyr::full_join(
        captures |> 
          dplyr::ungroup() |> 
          dplyr::select(siteID, plotID, period) |> 
          dplyr::distinct()
      ) |> 
      dplyr::group_by(siteID, period) |> 
      dplyr::mutate(plots_sampled = list(unique(plot))) |> 
      dplyr::select(-plotID, -plot) |> 
      dplyr::distinct()) |> 
  dplyr::ungroup() |> 
  dplyr::full_join(
    captures |> 
      dplyr::ungroup() |> 
      dplyr::select(siteID, plotID) |> 
      dplyr::distinct() |> 
      dplyr::arrange(plotID) |> 
      dplyr::group_by(siteID) |> 
      dplyr::mutate(plot = dplyr::row_number())) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(plotst = ifelse(!is.na(plot),
                                plot, 
                                sample(plots_sampled, 1, replace = TRUE))) |> 
  dplyr::select(siteID, period, plotID, plot, plotst, scientificName, ind, z_index) |> 
  dplyr::group_by(scientificName) |> 
  dplyr::mutate(sp = dplyr::cur_group_id()) |> 
  dplyr::group_by(siteID) |> 
  dplyr::mutate(site = dplyr::cur_group_id()) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(siteID, period) |> 
  dplyr::mutate(max_plot = max(plot, na.rm = TRUE)) |> 
  dplyr::ungroup() |> 
  dplyr::arrange( siteID, period, sp, ind )

plot_canopy <- da_captures |>
  dplyr::full_join(get_z_index) |> 
  dplyr::left_join(site_in_range) |>
  dplyr::select(siteID, plotID, plot) |> 
  dplyr::distinct() |> 
  dplyr::filter(!is.na(plotID)) |> 
  dplyr::arrange(siteID, plot) |> 
  dplyr::left_join(canopy) |>
  dplyr::group_by(siteID, plotID, plot) |> 
  dplyr::summarise(pcanopy_500m = mean(pcanopy_500m)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(canopy = base::as.numeric( base::scale(pcanopy_500m))) |>
  dplyr::select(-plotID, -pcanopy_500m) |> 
  tidyr::pivot_wider(names_from = plot, values_from = canopy) |> 
  dplyr::mutate( across(`1`:`8`, function(x) tidyr::replace_na( x, 0 ) ) ) |>
  dplyr::select(-siteID) |> 
  base::as.matrix() |>
  base::unname()

range_mat <- site_in_range |> 
  dplyr::filter( siteID %in% unique(da_captures$siteID)) |> 
  dplyr::filter( scientificName %in% unique(da_captures$scientificName)) |> 
  dplyr::arrange(scientificName, siteID) |> 
  dplyr::group_by(scientificName) |> 
  dplyr::mutate(sp = cur_group_id()) |> 
  dplyr::group_by(siteID) |> 
  dplyr::mutate(site = cur_group_id()) |> 
  dplyr::ungroup() |> 
  dplyr::select(sp, site, in_range) |> 
  tidyr::pivot_wider(names_from = site, values_from = in_range) |> 
  dplyr::select(-sp) |> 
  as.matrix() |> 
  unname()

final <- da_captures |>
  dplyr::full_join(get_z_index) |> 
  dplyr::left_join(site_in_range) |> 
  dplyr::group_by(scientificName) |> 
  dplyr::mutate(sp = cur_group_id()) |> 
  dplyr::group_by(siteID) |> 
  dplyr::mutate(site = cur_group_id()) |> 
  dplyr::arrange(siteID, period, sp, ind)

max_plot <- get_z_index |> 
  dplyr::select(siteID, period, plot) |>
  dplyr::distinct() |> 
  dplyr::group_by(siteID, period) |> 
  dplyr::summarise(max_plot = max(plot, na.rm = TRUE)) |> 
  tidyr::pivot_wider(names_from = period, values_from = max_plot) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(siteID) |>
  dplyr::select(-siteID) |> 
  dplyr::mutate(across(1:max(get_z_index$period, na.rm = TRUE), function(x) ifelse(is.na(x), 1, x))) |> 
  base::as.matrix() |> 
  base::unname()

max_period <- get_z_index |> 
  dplyr::select(siteID, period) |> 
  dplyr::distinct() |> 
  dplyr::group_by(siteID) |> 
  dplyr::summarise(max_period = max(period, na.rm = TRUE)) |> 
  dplyr::arrange(siteID) |>
  # two sites that have only 1 period. Change it to 2 to protect looping; shouldn't affect things?
  dplyr::mutate(max_period = ifelse(max_period == 1, 2, max_period)) |> 
  dplyr::pull(max_period)

sp_site_M <- figure_m |> 
  dplyr::arrange(siteID, scientificName) |>
  dplyr::group_by(siteID) |> 
  dplyr::mutate(site = cur_group_id()) |>
  dplyr::ungroup() |> 
  dplyr::select(-siteID, -max_ind) |> 
  tidyr::pivot_wider(names_from = site, values_from = M) |> 
  dplyr::arrange(scientificName) |> 
  dplyr::select(-scientificName) |> 
  dplyr::mutate(across(1:length(unique(figure_m$siteID)), function(x) replace_na(x, 100))) |> 
  as.matrix() |> 
  unname()

data <- list(
  y = final$y,
  TMIN = as.numeric(scale(final$tmin)), 
  PRCP = as.numeric(scale(final$prcp)),
  in_range = range_mat, 
  plot = get_z_index$plot,
  CANOPY = plot_canopy)

constants <- list(
  nsp = length(unique(final$sp)),
  nsite = length(unique(final$site)),
  max_plot = max_plot,
  M = sp_site_M,
  max_period = max_period,
  nz = nrow(get_z_index),
  species_z = get_z_index$sp, 
  site_z = get_z_index$site, 
  period_z = get_z_index$period,
  max_plot_z = get_z_index$max_plot,
  ny = nrow(final), 
  species_y = final$sp, 
  z_index = final$z_index)

setwd(here::here("data"))

save( data, 
      constants, 
      final,
      get_z_index,
      file = paste0("neon_cr_data_", Sys.Date(), ".RData"))