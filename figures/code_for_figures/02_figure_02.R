library(tidyverse)
library(here)
library(patchwork)
library(sf)

setwd(here::here("data"))
d <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv") |> 
  dplyr::select(siteID, site) |> 
  dplyr::distinct()

setwd(here::here("results"))

df <- read_csv("sr_site_posterior.csv") |> 
  dplyr::group_by(site, period) |> 
  dplyr::summarise(mean = mean(sr), 
                   sd = sd(sr)) |> 
  dplyr::filter(mean > 0) |> 
  tibble::add_column(metric = "Species richness") |> 
  dplyr::full_join(
    read_csv("pd_site_posterior.csv") |> 
      dplyr::group_by(site, period) |> 
      dplyr::summarise(mean = mean(pd), 
                       sd = sd(pd)) |> 
      dplyr::filter(mean > 0) |> 
      tibble::add_column(metric = "Peromyscus dominance")
    ) |> 
  dplyr::full_join(
    read_csv("n_site_posterior.csv") |> 
      dplyr::group_by(site, period) |> 
      dplyr::summarise(mean = mean(log1p(totN)), 
                       sd = sd(log1p(totN))) |> 
      dplyr::filter(mean > 0) |> 
      tibble::add_column(metric = "Total abundance")
  ) |> 
  dplyr::full_join(
    read_csv("shan_site_posterior.csv") |> 
      dplyr::group_by(site, period) |> 
      dplyr::summarise(mean = mean(shannon), 
                       sd = sd(shannon)) |> 
      dplyr::filter(mean > 0) |> 
      tibble::add_column(metric = "Shannon diversity")
  ) |> 
  dplyr::full_join(
    read_csv("faith_site_posterior.csv") |> 
      dplyr::group_by(site, period) |> 
      dplyr::summarise(mean = mean(n_pd), 
                       sd = sd(n_pd)) |> 
      dplyr::filter(mean > 0) |> 
      tibble::add_column(metric = "Phylogenetic diversity")
  )

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
             46.236238875, 47.1379026363636)) 

for_plot <- df |> 
  dplyr::left_join(d) |> 
  dplyr::left_join(sites)|>
  sf::st_as_sf(coords = c("xcoord", "ycoord"), 
               crs = 4326) 

( sr <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Species richness" & period == 1),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black", 
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Species richness") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( sr_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Species richness" & period == 1),
                     aes(geometry = geometry, fill = sd),
                     pch = 21, 
                     color = "black", 
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("SD", option = "A") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Species richness") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))
  

( pd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Peromyscus dominance" & period == 1),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black", 
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle( expression( paste( italic("Peromyscus"), " dominance"))) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( pd_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Peromyscus dominance" & period == 1),
                     aes(geometry = geometry, fill = sd),
                     pch = 21, 
                     color = "black", 
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("SD", 
                                   option = "A") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle( expression( paste( italic("Peromyscus"), " dominance"))) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( n <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Total abundance" & period == 1),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black",
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Log( Total abundance )") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( n_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Total abundance" & period == 1),
                     aes(geometry = geometry, fill = sd),
                     pch = 21, 
                     color = "black",
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("SD",
                                   option = "A") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Log( Total abundance )") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank())  + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( shan <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Shannon diversity" & period == 1),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black",
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Shannon diversity") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))


( shan_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Shannon diversity" & period == 1),
                     aes(geometry = geometry, fill = sd),
                     size = 2.5,
                     pch = 21, 
                     color = "black") +  
    ggplot2::scale_fill_viridis_c("SD", option = "A") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Shannon diversity") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( faith <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Phylogenetic diversity" & period == 1),
                     aes(geometry = geometry, fill = mean),
                     pch = 21, 
                     color = "black",
                     size = 2.5) +  
    ggplot2::scale_fill_viridis_c("Mean") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Phylogenetic diversity") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) + 
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( faith_sd <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = usa, aes(geometry = geom)) +
    ggplot2::geom_sf(data = filter(for_plot, metric == "Phylogenetic diversity" & period == 1),
                     aes(geometry = geometry, fill = sd),
                     size = 2.5,
                     pch = 21, 
                     color = "black") +  
    ggplot2::scale_fill_viridis_c("SD", option = "A") + 
    ggplot2::theme_void() +
    ggplot2::ggtitle("Phylogenetic diversity") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             size = 10,
                                             color = "black"),
                   legend.title = element_text(size = 9, color = "black"),
                   legend.text = element_text(size = 9, color = "black"),
                   legend.box.margin = margin(0, 5, 0, 0),
                   legend.frame = element_rect(color = "black"),
                   legend.ticks = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 3)))

( sr / sr_sd) | ( pd / pd_sd ) | ( n / n_sd ) | (shan / shan_sd) | ( faith / faith_sd)

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_02.png", 
  width = 9.5, 
  height = 3, 
  units = "in", 
  dpi = 600)
