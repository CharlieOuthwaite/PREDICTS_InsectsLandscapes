##%######################################################%##
#                                                          #
####             1. Organise PREDICTS data              ####
#                                                          #
##%######################################################%##

# This script organises the PREDICTS data for all insects, but subsets to just 
# those can be considered pollinators or pest controllers. 


# clear environment
rm(list = ls())

# load libraries
library(predictsFunctions)
library(StatisticalModels)
library(sf)
library(ggplot2)
library(raster)
library(BioStatR)
library(ggthemes)


# set the data folder
datadir <- "Data/"

# where to save the final dataset
outdir <- "1_PREDICTS_PLUS_VARIABLES"
if(!dir.exists(outdir)) dir.create(outdir)


# read in the complete PREDICTS dataset
pred.data <- readRDS(paste0(datadir, "/database.rds")) # 3250404 rows

# subset to the insects only
pred.data <- pred.data[(pred.data$Class=="Insecta"),] # 935078

# read in the predicts subset for pollinators
pred.data.pols <- readRDS(paste0(datadir, "/PREDICTS_pollinators_8_exp.rds")) # 357930 rows

# subset to insect pollinators only
pred.data.pols <- pred.data.pols[pred.data.pols$Class == "Insecta", ] # 279647 rows

# read in the species lists  pest controllers
pcs <- read.csv(paste0(datadir, "/pc_database.csv"))  # 4301 rows

# select pest controllers that have a classification of 2c and above.
pcs <- pcs[grep("2", pcs$Pest_control), ] # 1163 rows

# create subset pest controllers (selecting from insect only subset)
pred.data.pc <- pred.data[pred.data$Genus %in% pcs$Genus, ]  # 133988 rows

# correct sampling effort 
pred.data.pols <- CorrectSamplingEffort(pred.data.pols)
pred.data.pc <- CorrectSamplingEffort(pred.data.pc)

# merge sites: this combines potential subsamples within one site
pred.data.pols <- MergeSites(pred.data.pols[, 1:67]) # 267112 rows, doesn't work with additional columns
pred.data.pc <- MergeSites(pred.data.pc) # 90177 rows

# Calculate site level metrics
pred.sites.metrics.pols <- SiteMetrics(pred.data.pols, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 5559 rows
pred.sites.metrics.pc <- SiteMetrics(pred.data.pc, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 3497 rows

# which land uses are we interested in, I have kept plantation in here too
sites.sub.pols <- pred.sites.metrics.pols[!pred.sites.metrics.pols$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide"), ]
sites.sub.pc <- pred.sites.metrics.pc[!pred.sites.metrics.pc$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide"), ]

# 4049 rows
# 2721 rows

# remove sites with NA in lat/long columns
sites.sub.pols <- sites.sub.pols[!is.na(sites.sub.pols$Longitude),  ] # 4049 rows
sites.sub.pc <- sites.sub.pc[!is.na(sites.sub.pc$Longitude),  ] # 2715 rows

####  Now the landscape data ####

# Changes made here, just interested in ncrop, field size, pesticide app and
# then percentage of NH, either as a variable to compare, or to look at 
# buffering effect.


############################################################
#                                                          #
#       Step 2: Determine the number of crops              #
#                in each grid square                       #  
#                                                          #
############################################################



## pollinators dataset ##

sites.sub_xy_pols <- sites.sub.pols[, c("Longitude", "Latitude")]
sites.sub_xy_pols <- SpatialPoints(sites.sub_xy_pols, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))

# get the intersection values, what are the production values for each site?
sites.sub.pols$ncrop <- extract(ncrops, sites.sub_xy_pols, na.rm = FALSE)

# replace NAs with 0
sites.sub.pols[is.na(sites.sub.pols$ncrop), "ncrop"] <- 0 # 4049 rows



## pest controllers dataset ##

sites.sub_xy_pc <- sites.sub.pc[, c("Longitude", "Latitude")]
sites.sub_xy_pc <- SpatialPoints(sites.sub_xy_pc, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))

# get the intersection values, what are the production values for each site?
sites.sub.pc$ncrop <- extract(ncrops, sites.sub_xy_pc, na.rm = FALSE)

# replace NAs with 0
sites.sub.pc[is.na(sites.sub.pc$ncrop), "ncrop"] <- 0 # 2715 rows


##%######################################################%##
#                                                          #
####             Percentage natural habitat             ####
#                                                          #
##%######################################################%##

# using reprojected percentage natural habitat raster based on Hoskins et al 2016:

# Hoskins, A.J., Bush, A., Gilmore, J., Harwood, T., Hudson, L.N., Ware, C., Williams, 
# K.J., and Ferrier, S. (2016). Downscaling land-use data to provide global estimates 
# of five land-use classes. Ecol. Evol.

# read in the raster
percNH <- raster(paste0(datadir,"/PercentNatural.tif"))

# 
percNH <- percNH/10

# extract the dataset info for the PREDICTS sites
sites.sub.pols$percNH <- extract(percNH, sites.sub_xy_pols, na.rm = FALSE)
sites.sub.pc$percNH <- extract(percNH, sites.sub_xy_pc, na.rm = FALSE)

# how many NAs
nrow(sites.sub.pols[is.na(sites.sub.pols$percNH),]) #3
nrow(sites.sub.pc[is.na(sites.sub.pc$percNH),]) #3




##%######################################################%##
#                                                          #
####             Pesticide application rate             ####
#                                                          #
##%######################################################%##


# read in the raster
pest_H <- raster(paste0(datadir,"/Pesticide_totalAPR_High.tif"))
pest_L <- raster(paste0(datadir,"/Pesticide_totalAPR_Low.tif"))

# extract the dataset info for the PREDICTS sites
sites.sub$pest_H <- extract(pest_H, sites.sub_xy, na.rm = FALSE)
sites.sub$pest_L <- extract(pest_L, sites.sub_xy, na.rm = FALSE)

# how many NAs
nrow(sites.sub[is.na(sites.sub$pest_H),]) #0
nrow(sites.sub[is.na(sites.sub$pest_L),]) #0



##%######################################################%##
#                                                          #
####                     Field Size                     ####
#                                                          #
##%######################################################%##
# see file "legend items" in data folder for info on coding of data in the raster.

# read in the raster
fields <- raster(paste0(datadir,"/Global Field Sizes/dominant_field_size_categories.tif"))

# extract the info for the predicts sites
sites.sub$fields <- extract(fields, sites.sub_xy, na.rm = FALSE)

# how many NAs
nrow(sites.sub[is.na(sites.sub$fields),]) #0

sites.sub$fields <- as.factor(sites.sub$fields)
plot(sites.sub$fields)



##%######################################################%##
#                                                          #
####                   Tropical sites                   ####
#                                                          #
##%######################################################%##


# get the tropical values
sites.sub$Tropical <- NA

sites.sub[sites.sub$Latitude > -23.44 & sites.sub$Latitude < 23.44, 'Tropical'] <- "Tropical"

# label the remaining as temperate
sites.sub[is.na(sites.sub$Tropical), 'Tropical'] <- "Temperate"

# set as a factor
sites.sub$Tropical <- as.factor(sites.sub$Tropical)
# levels: Temperate Tropical

table(sites.sub$Tropical)

# Temperate  Tropical 
# 3481           2125  





##%######################################################%##
#                                                          #
####                 Assess the dataset                 ####
#                                                          #
##%######################################################%##


# remove any rows that have NA in the variable columns
summary(is.na(sites.sub))


# remove rows with NAs for any variable of interest
sites.sub <- sites.sub[!is.na(sites.sub$Hansen_mindist), ] # 133 NAs 
sites.sub <- sites.sub[!is.na(sites.sub$homogen), ] # 45 NAs 
sites.sub <- sites.sub[!is.na(sites.sub$percNH), ] # 3 NAs 4170 rows

sites.sub.pols <- sites.sub.pols[!is.na(sites.sub.pols$Hansen_mindist), ] # 133 NAs 
sites.sub.pols <- sites.sub.pols[!is.na(sites.sub.pols$homogen), ] # 44 NAs 
sites.sub.pols <- sites.sub.pols[!is.na(sites.sub.pols$percNH), ] # 3 NAs 3337 rows

sites.sub.pc <- sites.sub.pc[!is.na(sites.sub.pc$Hansen_mindist), ] # 133 NAs 
sites.sub.pc <- sites.sub.pc[!is.na(sites.sub.pc$homogen), ] # 18 NAs 



# remove those sites that have "Cannot decide" as a use intensity
nrow(sites.sub[sites.sub$Use_intensity == "Cannot decide", ]) # 549
sites.sub <- sites.sub[!sites.sub$Use_intensity == "Cannot decide", ] 
nrow(sites.sub.pols[sites.sub.pols$Use_intensity == "Cannot decide", ]) # 435
sites.sub.pols <- sites.sub.pols[!sites.sub.pols$Use_intensity == "Cannot decide", ] 
nrow(sites.sub.pc[sites.sub.pc$Use_intensity == "Cannot decide", ]) # 402
sites.sub.pc <- sites.sub.pc[!sites.sub.pc$Use_intensity == "Cannot decide", ] 

# nrows of dataset
nrow(sites.sub) # 3621
nrow(sites.sub.pols) # 2902
nrow(sites.sub.pc) # 1454

