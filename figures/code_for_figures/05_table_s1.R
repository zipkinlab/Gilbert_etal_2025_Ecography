library(here)
library(tidyverse)
library(lubridate)
library(janitor)
library(sf)
library(officer)
library(flextable)
library(magrittr)

setwd(here::here("data"))

d <- read_csv("disease_with_biodiversity_metrics_v01.csv")

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
  dplyr::mutate(collectDate = lubridate::as_date(collectDate)) |> 
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
  dplyr::filter(!siteID == "LENO") |> 
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
  dplyr::left_join(mam_surveys) |> 
  dplyr::left_join( path_surveys ) |> 
  dplyr::mutate( diff = abs(as.numeric(box_collect_date - path_collect_date))) |> 
  dplyr::filter( diff <= 10 ) |> 
  dplyr::select( siteID, plotID, period, box_collect_date, path_collect_date ) |> 
  dplyr::select(siteID, plotID) |>
  dplyr::distinct() |>
  dplyr::arrange(siteID, plotID) |>
  dplyr::group_by(siteID) |>
  dplyr::mutate(plot = row_number()) |> 
  dplyr::full_join(
    full_join( path_surveys, mam_surveys ) |>
      dplyr::filter(!siteID == "LENO") |> 
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
      dplyr::left_join(mam_surveys) |> 
      dplyr::left_join( path_surveys ) |> 
      dplyr::mutate( diff = abs(as.numeric(box_collect_date - path_collect_date))) |> 
      dplyr::filter( diff <= 10 ) |> 
      dplyr::select( siteID, plotID, period, box_collect_date, path_collect_date ))

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

plot_coords <- survey_key |> 
  dplyr::select(siteID, plotID, collectDate = box_collect_date) |> 
  dplyr::distinct() |> 
  dplyr::left_join(boxtrap[[5]]) |> 
  dplyr::select(siteID, plotID, xcoord = decimalLongitude, ycoord = decimalLatitude) |> 
  dplyr::distinct()

# calculate average area of NEON site (polygone around plot centroids)
plot_coords |> 
  sf::st_as_sf(coords = c("xcoord", "ycoord"), 
           crs = 4326) |> 
  dplyr::group_by(siteID) |> 
  dplyr::summarise(geometry = st_union(geometry)) |> 
  sf::st_convex_hull() |> 
  sf::st_area() |> 
  mean() / 1e6

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

glimpse(captures)

n_table <- captures |> 
  dplyr::ungroup() |> 
  dplyr::select(tagID, scientificName) |> 
  dplyr::distinct() |> 
  count(scientificName) |> 
  dplyr::full_join(d |>
                     count(scientificName) |>
                     rename(ntest = n)) |> 
  dplyr::arrange(scientificName) |> 
  dplyr::rename( Species = scientificName, 
                 `Number captured` = n, 
                 `Number screened` = ntest)

flextable::set_flextable_defaults(font.size = 10)
ft <- flextable::flextable( data = n_table, cwidth = 2)  

tmp <- tempfile(fileext = ".docx")

officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)