# 18 December 2023
# Download the 2 relevant data prducts from NEON
# Tick pathogen status in rodent tissue and small mammal box trapping data

library(neonUtilities)
library(here)

setwd(here::here("data"))

# DP1.10064.002 - rodent pathogen status, tick-borne
# Only available from 2020 onward
rodent_pathogen <- neonUtilities::loadByProduct(dpID = "DP1.10064.002", check.size = FALSE)
save(rodent_pathogen, file = "neon_mammal_tick_pathogen_v01.RData")

# DP1.10072.001 - small mammal box trapping
boxtrap <- neonUtilities::loadByProduct(dpID = "DP1.10072.001",
                                        check.size = FALSE,
                                        startdate = "2020-05", # tick pathogen screening begins in June 2020
                                        # this omits island (HI, PR) and western sites (xcoord > -103)
                                        site = c("BART", "BLAN", "CLBJ", "DCFS", "DELA", "DSNY", "GRSM", "HARV", 
                                                 "JERC", "KONA", "KONZ", "LENO", "MLBS", "NOGP", "OAES", "ORNL", 
                                                 "OSBS", "SCBI", "SERC", "STEI", "TALL", "TREE", "UKFS", "UNDE", 
                                                 "WOOD"))

save(boxtrap, file = "neon_mammal_box_trapping_v01.RData")


mamseq <- neonUtilities::loadByProduct(dpID = "DP1.10076.001",
                                       site = c("BART", "BLAN", "CLBJ", "DCFS", "DELA", "DSNY", "GRSM", "HARV", 
                                                "JERC", "KONA", "KONZ", "LENO", "MLBS", "NOGP", "OAES", "ORNL", 
                                                "OSBS", "SCBI", "SERC", "STEI", "TALL", "TREE", "UKFS", "UNDE", 
                                                "WOOD"),
                                     startdate = "2020-05",
                                     check.size = FALSE)
save(mamseq, file = "neon_mammal_sequences_v01.RData")
