##%######################################################%##
#                                                          #
####             Pesticide data processing              ####
#                                                          #
##%######################################################%##

# This script pulls together the PEST CHEM-GRID data into one raster of
# total pesticide application per grid cell.


# The data:
# 1200 maps for top 20 active ingredients in pesticides
# estimated 2015
# projected for 2020 and 2025

# file names:
# APR = application rate
# QI = qality index
# HarvestedArea = slightly different to Monfreda versions. 

# Resolution = 5 arc-minute (10km at the equator)


# files are per crop, per active ingredient, per year, for a high or low estimate


rm(list = ls())

# load libraries
library(raster)

# directories - on RDS
datadir <- "Z:/Datasets/PEST CHEMGRIDS/ApplicationRates/PEST-CHEMGRIDS_v1_01_APR/PEST-CHEMGRIDS_v1_01_APR/GEOTIFF"
outdir <- "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/8. Forest cover - insects, pollinators, pest controllers/PREDICTS_InsectsLandscapes/Data"

# list the files
files <- list.files(datadir) # 1200 files

# just select the 2015 maps
files <- list.files(datadir, pattern = "2015") # 400 files, high and low estimates

# take a look at one of the maps
test_map <- raster(paste0(datadir, "/", files[1]))

plot(test_map)
summary(test_map)

test_map2 <- round(test_map, digits = 3)
plot(test_map2)
summary(test_map2)
unique(test_map2)


test_map3 <- test_map
test_map3[test_map3 == -1] <- "B_NA"
test_map3[test_map3 == -1.5] <- "No_Data"
test_map3[test_map3 == -2] <- NA
plot(test_map3)





# from looking at the pngs and data ranges, I think all values 0 + are the rates
# convert all that is below 0 to NA
test_map[test_map < 0 ] <- NA
plot(test_map)

### Summarise the data across crops and active ingredients ###


# separate out the low and high estimates
files_H <- files[grep(files, pattern = "5_H")] # 200 files
files_L <- files[grep(files, pattern = "5_L")] # 200 files

# stack the raster files and summarise the information
pest_H <- stack(paste0(datadir, "/", files_H))

# first convert all negative values to NA
pest_H[pest_H < 0] <- NA

# get the sum of app rates across the crops/pesticides
pest_H_total <- calc(x = pest_H, fun = sum)

# remove large raster stack
rm(pest_H)

# same for the low estimates
pest_L <- stack(paste0(datadir, "/", files_L))

#  convert all negative values to NA
pest_H[pest_H < 0] <- NA

# get the sum of app rates across the crops/pesticides
pest_L_total <- calc(x = pest_L, fun = sum)

# remove large raster stack
rm(pest_L)

plot(pest_H_total)
plot(pest_L_total)

# save these maps

