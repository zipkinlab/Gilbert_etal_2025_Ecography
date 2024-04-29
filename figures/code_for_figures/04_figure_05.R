library(here)
library(tidyverse)
library(MCMCvis)
library(nimble)

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

setwd(here::here("data"))
d <- read_csv("disease_with_biodiversity_metrics_v01.csv")

dis <- d |> 
  group_by(siteID) |> 
  summarise(npos = sum(positive), 
            n = sum(!is.na(positive))) |> 
  dplyr::mutate(prop = npos / n) |> 
  dplyr::left_join(sites) |> 
  sf::st_as_sf(coords = c("xcoord", "ycoord"), 
               crs = 4326) |> 
  dplyr::rename(`Number tested` = n)

ggplot() + 
  geom_sf( data = usa, aes(geometry = geom)) +
  geom_sf( data = dis,
           aes(geometry = geometry, fill = prop, size = `Number tested`),
           pch = 21, 
           color = "black") +
  scale_size(breaks = c(1, 30, 70, 100, 150, 180),
             range = c(2, 8)) +
  scale_fill_viridis_c("Proportion infected ", 
                         limits = c(0, 0.12),
                       breaks = c(0, 0.05, 0.10)) +
# breaks = c(0, 0.05, 0.10)) +
theme_void() +
  theme(legend.title = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 10, color = "black"), 
        legend.frame = element_rect(color = "black"),
        legend.ticks = element_blank(),
        plot.background = element_rect(color = NA, fill = "white"), 
        panel.background = element_rect(color = NA, fill = "white")) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4))

setwd(here::here("figures"))  
ggsave(filename = "figure_05.png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 600)