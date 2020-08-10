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
library(ggplot2)
library(maps)
library(viridis)

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
#pest_H[pest_H < 0] <- NA


# try the reclassify function
#test_map <- reclassify(test_map, matrix(data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)
pest_H <- reclassify(pest_H, matrix(data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# get the sum of app rates across the crops/pesticides
pest_H_total <- calc(x = pest_H, fun = sum, na.rm = TRUE)

# take a look at the total application rates
plot(pest_H_total$layer)
pest_H_total_data <- as.data.frame(pest_H_total, xy = TRUE)

# save the raster of high estimate totals
writeRaster(pest_H_total, filename = paste0(outdir, "/Pesticide_totalAPR_High.tif"))

# remove large raster stack
#rm(pest_H)



# same for the low estimates
pest_L <- stack(paste0(datadir, "/", files_L))

#  convert all negative values to NA
pest_L <- reclassify(pest_L, matrix(data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)
# started 17:16 07/08/2020

# get the sum of app rates across the crops/pesticides
pest_L_total <- calc(x = pest_L, fun = sum, na.rm = TRUE)

# take a look at the total application rates
plot(pest_L_total)
range(pest_L_total$layer)

# remove large raster stack
rm(pest_L)

# save the raster of low application rates
writeRaster(pest_L_total, filename = paste0(outdir, "/Pesticide_totalAPR_Low.tif"))



##### plot of the datasets for supplementary info  #####

pest_L_total <- raster(paste0(outdir, "/Pesticide_totalAPR_Low.tif"))
pest_H_total <- raster(paste0(outdir, "/Pesticide_totalAPR_High.tif"))  


# mask out the sea area that is set to 0

data("wrld_simpl", package = 'maptools')
pest_L_total_crop <- mask(pest_L_total$Pesticide_totalAPR_Low, wrld_simpl)
pest_L_total_crop <- crop(pest_L_total_crop, wrld_simpl)
pest_H_total_crop <- mask(pest_H_total$Pesticide_totalAPR_High, wrld_simpl)
pest_H_total_crop <- crop(pest_H_total_crop, wrld_simpl)


# convert for use in ggplot
pest_H_data <- as.data.frame(pest_H_total_crop, xy = TRUE)
pest_L_data <- as.data.frame(pest_L_total_crop, xy = TRUE)

pest_H_data <- pest_H_data[!is.na(pest_H_data$Pesticide_totalAPR_High), ]
pest_L_data <- pest_L_data[!is.na(pest_L_data$Pesticide_totalAPR_Low), ]


pest_H_data$est <- "High"
pest_L_data$est <- "Low"

names(pest_H_data)[3] <- "APR"
names(pest_L_data)[3] <- "APR"


pest_data <- rbind(pest_H_data, pest_L_data)
pest_data$est <- factor(pest_data$est, levels = c("Low", "High"))

# set breaks
brk <- c(0, 1, 3, 10, 30, 100, 300)

# plot
ggplot(data = pest_data) +
  geom_raster(aes(x = x, y = y, fill = APR)) +
  facet_grid(~ est) +
  scale_fill_gradientn(name = "Total APR\nkg/ha", breaks = brk, trans = "log1p", colours = viridis(10), labels = brk) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = "bottom", 
        legend.key.width = unit(2,"cm")) 

# save
ggsave(filename = paste0(outdir, "/SuppFigs/Pesticide_data_maps.pdf"), width = 8, height = 4)

