# create animations of biodiversity metrics through time. this script 
# generates individual gifs, and then I patched them together in powerpoint post hoc

library(tidyverse)
library(here)
library(patchwork)
library(sf)

# NOTE: this results file is too large to push to GitHub
# User can either run script 03_run_capture_recapture_model.R to produce a similar result file
# OR the file can be downloaded from this OneDrive link: https://1drv.ms/u/s!AtvYBfNq7AMkhIIjVUn6z4ooDwt1Kw?e=iUlju0
setwd(here::here("results"))
load("neon_capture_recapture_results_2024-04-11.RData")

read_file_with_name <- function(file_name) {
  # Read the CSV file
  data <- read_csv(file_name)
  
  # Add a new column with the filename
  data <- data %>% mutate(filename = file_name)
  
  return(data)
}

file_names <- list.files(pattern = "*posterior.csv")

combined_data <- file_names |> 
  map_dfr(read_file_with_name) |> 
  dplyr::mutate(filename = stringr::str_remove_all(filename, "_posterior.csv")) |> 
  dplyr::mutate( value = ifelse(!is.na(n_pd), n_pd,
                                 ifelse(!is.na(totN), totN, 
                                        ifelse(!is.na(pd), pd, 
                                               ifelse(!is.na(shannon), shannon, sr))))) |> 
  tidyr::separate( filename, into = c("metric", "scale"), sep = "_") |> 
  dplyr::select(metric, scale, index, site, period, plot, iter, value)


summary <- combined_data |> 
  dplyr::filter(scale == "site") |> 
  dplyr::group_by(metric, scale, site, period) |> 
  dplyr::summarise( mean = mean(value, na.rm = TRUE),
                    sd = sd(value, na.rm = TRUE))


usa <- sf::st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) |> 
  dplyr::filter(! ID %in% c("washington", "oregon", "california", "nevada", "idaho", "utah", "arizona", "new mexico", 
                            "colorado", "wyoming", "montana") )

sites <- tibble::tibble(
  siteID = c("BART", "BLAN", "CLBJ", "DCFS", "DELA", "DSNY", "GRSM", "HARV", 
             "JERC", "KONA", "KONZ", "LENO", "MLBS", "NOGP", "OAES", "ORNL", 
             "OSBS", "SCBI", "SERC", "STEI", "TALL", "TREE", "UKFS", "UNDE", 
             "WOOD"), 
  xcoord = c(-71.3001021666667, -77.9643916, -97.61411925, -99.096803, -87.8097465, 
             -81.4102095555556, -83.4818927777778, -72.2440151111111, -84.457995, 
             -96.6270666666667, -96.570225125, -88.1918962857143, -80.5240751666667, 
             -100.911671285714, -99.0628615, -84.269079, -81.980435375, -78.1503481428571, 
             -76.542701125, -90.0730375, -87.4199015, -89.5628562857143, -95.196369, 
             -89.520750625, -99.2389639090909),
  ycoord = c(44.049036, 39.0864602, 33.3610965, 47.1750165, 32.53394875, 
             28.082952, 35.6918571111111, 42.4386822222222, 31.193419, 39.1550028333333, 
             39.093005125, 31.8143305714286, 37.400306, 46.7946964285714, 
             35.412609875, 35.9345959, 29.69836225, 38.89244, 38.8898185, 
             45.8091743333333, 32.90788275, 45.4896738571429, 39.0467195, 
             46.236238875, 47.1379026363636)) |> 
  dplyr::left_join(
    final |> 
      dplyr::select(siteID, site) |> 
      dplyr::distinct()
  ) |> 
  dplyr::filter(!is.na(site))

for_plot <- final |> 
  dplyr::select(siteID, site) |> 
  dplyr::distinct() |> 
  dplyr::full_join(
    sites
  ) |> 
  dplyr::full_join(
    summary
  ) |> 
  sf::st_as_sf(coords = c("xcoord", "ycoord"), 
               crs = 4326) |> 
  dplyr::mutate(across(mean:sd, function(x) ifelse(x == 0, NA, x)))

sr_p <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
  ggplot2::geom_sf(data = filter(for_plot, metric == "sr"),
                   aes(geometry = geometry, fill = mean),
                   pch = 21, 
                   color = "black", 
                   size =5) +  
  ggplot2::scale_fill_viridis_c("Mean") +
  transition_states(period) +
  ggplot2::theme_void() +
  ggplot2::ggtitle("Species richness, period: {closest_state}") +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                           size = 10,
                                           color = "black"),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 9, color = "black"),
                 legend.box.margin = margin(0, 5, 0, 0),
                 legend.frame = element_rect(color = "black"),
                 legend.ticks = element_blank())  + 
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

sr_anim <- animate(sr_p, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "sr_mean_animation.gif",
          animation = sr_anim)

sr_sd_p <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
  ggplot2::geom_sf(data = filter(for_plot, metric == "sr"),
                   aes(geometry = geometry, fill = sd),
                   pch = 21, 
                   color = "black", 
                   size =5) +  
  ggplot2::scale_fill_viridis_c("SD", option = "A") + 
  transition_states(period) +
  ggplot2::theme_void() +
  ggplot2::ggtitle("Species richness, period: {closest_state}") +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                           size = 10,
                                           color = "black"),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 9, color = "black"),
                 legend.box.margin = margin(0, 5, 0, 0),
                 legend.frame = element_rect(color = "black"),
                 legend.ticks = element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

sr_sd_anim <- animate(sr_sd_p, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "sr_sd_animation.gif",
          animation = sr_sd_anim)

pd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "pd"),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black", 
                     size = 5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    transition_states(period) +
    ggplot2::ggtitle( "Peromyscus dominance, period: {closest_state}") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

pd_mean_anim <- animate(pd, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "pd_mean_animation.gif",
          animation = pd_mean_anim)

pd_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "pd"),
                     aes(geometry = geometry, fill = sd),
                     pch = 21, 
                     color = "black", 
                     size = 5) +  
    ggplot2::scale_fill_viridis_c("SD", 
                                  option = "A") + 
    ggplot2::theme_void() +
    transition_states(period) +
    ggplot2::ggtitle( "Peromyscus dominance, period: {closest_state}") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

pd_sd_anim <- animate(pd_sd, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "pd_sd_animation.gif",
          animation = pd_sd_anim)

n <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "n"),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black",
                     size = 5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    transition_states(period) +
    ggplot2::ggtitle("Log( Total abundance ), period: {closest_state}") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

n_mean_anim <- animate(n, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "n_mean_animation.gif",
          animation = n_mean_anim)

n_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "n"),
                     aes(geometry = geometry, fill = sd),
                     pch = 21, 
                     color = "black",
                     size = 5) +  
    ggplot2::scale_fill_viridis_c("SD",
                                  option = "A") + 
    ggplot2::theme_void() +
    transition_states(period) +
    ggplot2::ggtitle("Log( Total abundance ), period: {closest_state}") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

n_sd_anim <- animate(n_sd, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "n_sd_animation.gif",
          animation = n_sd_anim)

shan <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "shan"),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black",
                     size = 5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Shannon diversity, period: {closest_state}") +
    transition_states(period) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

shan_mean_anim <- animate(shan, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "shan_mean_animation.gif",
          animation = shan_mean_anim)

shan_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "shan"),
                     aes(geometry = geometry, fill = sd),
                     size = 5,
                     pch = 21, 
                     color = "black") +  
    ggplot2::scale_fill_viridis_c("SD", option = "A") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Shannon diversity, period: {closest_state}") +
    transition_states(period) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

shan_sd_anim <- animate(shan_sd, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "shan_sd_animation.gif",
          animation = shan_sd_anim)

faith <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
  ggplot2::geom_sf(data = filter(for_plot, metric == "faith"),
                   aes(geometry = geometry, fill = mean),
                   pch = 21, 
                   color = "black",
                   size = 5) +  
  ggplot2::scale_fill_viridis_c("Mean") + 
  ggplot2::theme_void() +
  ggplot2::ggtitle("Phylogenetic diversity, period: {closest_state}") +
  transition_states(period) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                           size = 10,
                                           color = "black"),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 9, color = "black"),
                 legend.box.margin = margin(0, 5, 0, 0),
                 legend.frame = element_rect(color = "black"),
                 legend.ticks = element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

faith_mean_anim <- animate(faith, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "faith_mean_animation.gif",
          animation = faith_mean_anim)

faith_sd <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
  ggplot2::geom_sf(data = filter(for_plot, metric == "faith"),
                   aes(geometry = geometry, fill = sd),
                   size = 5,
                   pch = 21, 
                   color = "black") +  
  ggplot2::scale_fill_viridis_c("SD", option = "A") + 
  ggplot2::theme_void() +
  ggplot2::ggtitle("Phylogenetic diversity, period: {closest_state}") +
  transition_states(period) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                           size = 10,
                                           color = "black"),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 9, color = "black"),
                 legend.box.margin = margin(0, 5, 0, 0),
                 legend.frame = element_rect(color = "black"),
                 legend.ticks = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3))

faith_sd_anim <- animate(faith_sd, fps = 10, height = 250, width = 300)
setwd(here::here("figures"))
anim_save(filename = "faith_sd_animation.gif",
          animation = faith_sd_anim)
