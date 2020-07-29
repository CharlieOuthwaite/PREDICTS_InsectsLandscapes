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

# list the files
files <- list.files(datadir) # 1200 files

# just select the 2015 maps
files <- list.files(datadir, pattern = "2015") # 400 files, high and low estimates

# separate out the low and high estimates
files_H <- files[grep(files, pattern = "5_H")]
files_L <- files[grep(files, pattern = "5_L")]

# stack the raster files
pest_H <- stack(paste0(datadir, "/", files_H))
pest_L <- stack(paste0(datadir, "/", files_L))


pest_H_total <- calc(x = pest_H, fun = sum)
pest_L_total <- calc(x = pest_L, fun = sum)




