library(tidyverse)
library(sf)
library(basemaps)
library(ggsn)
library(cowplot)
library(MetBrewer)

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
  sf::st_as_sf(coords = c("xcoord", "ycoord"), 
               crs = 4326, 
               agr = "constant")

usa <- sf::st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) |> 
  dplyr::filter(! ID %in% c("washington", "oregon", "california", "nevada", "idaho", "utah", "arizona", "new mexico", 
                            "colorado", "wyoming", "montana") )

( neon_site_map <- ggplot2::ggplot() +
    ggplot2::geom_sf( data = usa, aes( geometry = geom),
                      fill = "gray95", 
                      color = "gray60",
                      size = 0.01) +
    ggplot2::geom_sf( data = sites, aes( geometry = geometry ),
                      color =  MetBrewer::MetPalettes$Hiroshige[[1]][10],
                      size = 1.5) +
    ggplot2::geom_sf( data = filter(sites, siteID == "BART"), aes( geometry = geometry ),
                      color =  MetBrewer::MetPalettes$Hiroshige[[1]][1],
                      size = 3) +
    ggplot2::theme_void() +
    ggplot2::theme( plot.background = element_rect(fill = "white", color = "black")))

bart <- tibble::tibble( plotID = c("BART_012", "BART_001", "BART_015", "BART_007", "BART_062", 
                                   "BART_084"),
                        xcoord = c(-71.315451, -71.305115, -71.275059, -71.298377, -71.308646, 
                                   -71.297965),
                        ycoord = c(44.047478, 44.046877, 44.042898, 44.051112, 44.05735, 44.048501)) |> 
  sf::st_as_sf( coords = c("xcoord", "ycoord"), crs = 4326, agr = "constant") |> 
  sf::st_transform(crs = 3857)

bart_hull <- bart |> 
  dplyr::summarise(geometry = sf::st_union(geometry)) |> 
  sf::st_convex_hull()

bart_hull <- sf::st_simplify( sf::st_buffer( 
  sf::st_convex_hull( sf::st_union( sf::st_geometry(bart))), dist = 300))

bart_ext <- sf::st_bbox(bart) |>  
  sf::st_as_sfc() |>  
  sf::st_buffer(500)

basemaps::set_defaults(map_service = "mapbox",
                       map_type = "satellite",
                       map_token = "pk.eyJ1IjoibmFnaWxiZXJ0IiwiYSI6ImNrd3ppdGF2bjB2aHQyem5zcnVyc214djEifQ.UMRkcBYCkKnn4nAsnbq_Dw")

bart_gg <- basemaps::basemap_ggplot(bart_ext)

( bart_plot_map <- bart_gg +
    ggplot2::geom_sf( data = bart_hull, aes(geometry = geometry), 
                      fill = MetBrewer::MetPalettes$Hiroshige[[1]][7],
                      alpha = 0.6, color = NA) +
    ggplot2::geom_sf( data = bart, aes(geometry = geometry),
                      size = 9,
                      color = MetBrewer::MetPalettes$Hiroshige[[1]][5],
                      alpha = 0.8) +
    ggplot2::theme_void() +
    ggsn::scalebar(data = bart,
                   location = "bottomleft",
                   dist = 0.5,
                   dist_unit = "km",
                   transform = FALSE,
                   st.color = "white",
                   st.dist = 0.07,
                   border.size = 0.2,
                   height = 0.02,
                   st.bottom = TRUE) )

cowplot::ggdraw() +
  cowplot::draw_plot( bart_plot_map ) +
  cowplot::draw_plot( neon_site_map, 
                      x = 0.605, 
                      y = 0.56, 
                      width = 0.4, 
                      height = 0.4)

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_01.png", 
  width = 5.95, 
  height = 3.5, 
  units = "in", 
  dpi = 600)
