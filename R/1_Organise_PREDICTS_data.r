##%######################################################%##
#                                                          #
####             1. Organise PREDICTS data              ####
#                                                          #
##%######################################################%##


# This script takes the PREDICTS dataset and subsets it into insects only.
# Then add detail of those species that are pollinators and/or pest controllers.


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
datadir <- "DATA"

# where to save the final dataset
outdir <- "1_PREDICTS_PLUS_VARIABLES"
dir.create(outdir)

# read in the complete PREDICTS dataset
pred.data <- readRDS(paste0(datadir, "/database.rds")) # 3250404 rows

# subset to the insects only
pred.data <- pred.data[(pred.data$Class=="Insecta"),] # 935078

# read in the predicts subset for pollinators
pred.data.pols <- readRDS(paste0(datadir, "/PREDICTS_pollinators_5.rds")) # 372646 rows

# subset to insect pollinators only
pred.data.pols <- pred.data.pols[pred.data.pols$Class == "Insecta", ] # 292916 rows

# read in the species lists  pest controllers
pcs <- read.csv(paste0(datadir, "/pc_database.csv"))  # 4301 rows

# select pest controllers that have a classification of 2c and above.
pcs <- pcs[grep("2", pcs$Pest_control), ] # 1163 rows

# create subset pest controllers (selecting from insect only subset)
pred.data.pc <- pred.data[pred.data$Genus %in% pcs$Genus, ]  # 133988 rows


### organise datasets using functions from predictsFunctions package ###

# correct sampling effort 
pred.data <- CorrectSamplingEffort(pred.data)
pred.data.pols <- CorrectSamplingEffort(pred.data.pols)
pred.data.pc <- CorrectSamplingEffort(pred.data.pc)

# merge sites: this combines potential subsamples within one site
pred.data <- MergeSites(pred.data) # 826525 rows
pred.data.pols <- MergeSites(pred.data.pols[, 1:67]) # 278351 rows, doesn't work with additional columns
pred.data.pc <- MergeSites(pred.data.pc) # 90177 rows

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(pred.data, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 7800 rows
pred.sites.metrics.pols <- SiteMetrics(pred.data.pols, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 6087 rows
pred.sites.metrics.pc <- SiteMetrics(pred.data.pc, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 3497 rows



# only interested in natural habitats plus cropland, drop other land uses #
# site level data cropland only
sites.sub <- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]
sites.sub.pols <- pred.sites.metrics.pols[!pred.sites.metrics.pols$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]
sites.sub.pc <- pred.sites.metrics.pc[!pred.sites.metrics.pc$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]

# 5120 rows
# 4176 rows
# 2437 rows

# remove sites with NA in lat/long columns
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ] # 5115 rows
sites.sub.pols <- sites.sub.pols[!is.na(sites.sub.pols$Longitude),  ] # 4176 rows
sites.sub.pc <- sites.sub.pc[!is.na(sites.sub.pc$Longitude),  ] # 2432 rows



############################################################
#                                                          #
#      Step 2: Subset data to those in forest biomes       #
#                                                          #
############################################################

## determine the forest biome for the all insecte sites, then match with the subsets

# load the WWF ecoregions shapefile
ecobio <- st_read(dsn = paste0(datadir, "/", "wwf_terr_ecos.shp"))

# convert the biome info to a factor so its not plotted as a continuous variable
ecobio$BIOME <- as.factor(ecobio$BIOME)

# get the lat longs of the PREDICTS sites
sites.xy <- sites.sub[, c("Longitude", "Latitude")] # 15612

# take only the biome attribute
wwfbiomes <- ecobio['BIOME']


# transform the lat/longs into spatial points
# set the CRS to the same as that of the forest biome polygons
pointsT <- st_as_sf(sites.xy, coords = c("Longitude", "Latitude"))
st_crs(pointsT) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84"
pointsT <- st_transform(pointsT, crs = st_crs(wwfbiomes))


# The polygons seem to be slightly mismatched with the land area, for the NAs try the within distance function to fill in the gaps.
biomes <- sf::st_intersects(pointsT, wwfbiomes)

# replace any that are NULL with NA
biomes[sapply(biomes, function(x) length(x)==0L)] <- NA

# Get the values for the biomes
biomes <- as.character(wwfbiomes$BIOME[unlist(biomes)])

# add info to the data table
sites.sub$Forest_biome <- biomes

# check for the NAs
nrow(sites.sub[is.na(sites.sub$Forest_biome), ]) # 29 sites


#### now try and fill in the sites that are NAs ####


# subset the points to just those with NA
pointsST_NA <- st_as_sf(sites.sub[is.na(sites.sub$Forest_biome), ], coords = c("Longitude", "Latitude"))
st_crs(pointsST_NA) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84"
pointsST_NA <- st_transform(pointsST_NA, crs = st_crs(wwfbiomes))

# take a look at why they are NAs
mapview::mapview(pointsST_NA) + mapview::mapview(wwfbiomes)
# just outside of polygons

# try and get the biomes for the NAs, using the within distance - these were checked by eye
# in some cases, more than one biome will be selected
biomes_NA <- sf::st_join(pointsST_NA, wwfbiomes, st_is_within_distance, dist = 10000)

pointsST_NA$Forest_biome <- NA

# list the sites with the NAs
vals <- unique(biomes_NA$SSBS)

# go through each site, get the biome if only one found, if not add NA
for(i in 1:length(vals)){
  
  x <- biomes_NA[biomes_NA$SSBS == vals[i], ]
  nrow(x)
  
  # if the multiple values are the same, just keep one
  if(length(unique(x$BIOME)) == 1){
    
    pointsST_NA[pointsST_NA$SSBS == vals[i], 'Forest_biome'] <- unique(x$BIOME)
  }
  else{pointsST_NA[pointsST_NA$SSBS == vals[i], 'Forest_biome'] <- NA}
}


# check the remaining NAs on the map
pointsST_NA[is.na(pointsST_NA$Forest_biome), ]
mapview::mapview(pointsST_NA) + mapview::mapview(wwfbiomes)

# just outside the boundary of biome 14
pointsST_NA[rownames(pointsST_NA) == 3485, "Forest_biome"] <- 14
pointsST_NA[rownames(pointsST_NA) == 3486, "Forest_biome"] <- 14
pointsST_NA[rownames(pointsST_NA) == 3493, "Forest_biome"] <- 14
pointsST_NA[rownames(pointsST_NA) == 3494, "Forest_biome"] <- 14


# add these in to the sites.sub table, theyre in the same order
sites.sub[is.na(sites.sub$Forest_biome), "Forest_biome"] <- pointsST_NA$Forest_biome


### add in the forest biome to the data subsets 



# subset to forest biomes only
# forest biomes include 1:6 and 12
sites.sub <- sites.sub[sites.sub$Forest_biome %in% c(1:6, 12), ] # 4345

# site per forest biome
table(sites.sub$Forest_biome)

# rename the forest biomes, taken from original dataset info
sites.sub$Forest_biome <- sub(12, "Mediterranean Forests, Woodlands & Scrub", sites.sub$Forest_biome)
sites.sub$Forest_biome <- sub(1, "Tropical & Subtropical Moist Broadleaf Forests", sites.sub$Forest_biome)
sites.sub$Forest_biome <- sub(2, "Tropical & Subtropical Dry Broadleaf Forests", sites.sub$Forest_biome)
sites.sub$Forest_biome <- sub(3, "Tropical & Subtropical Coniferous Forests", sites.sub$Forest_biome)
sites.sub$Forest_biome <- sub(4, "Temperate Broadleaf & Mixed Forests", sites.sub$Forest_biome)
sites.sub$Forest_biome <- sub(5, "Temperate Conifer Forests", sites.sub$Forest_biome)
sites.sub$Forest_biome <- sub(6, "Boreal Forests/Taiga", sites.sub$Forest_biome)

sites.sub$Forest_biome <- factor(sites.sub$Forest_biome,
                                 levels=c("Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests",
                                          "Boreal Forests/Taiga", "Mediterranean Forests, Woodlands & Scrub",  
                                          "Tropical & Subtropical Coniferous Forests", "Tropical & Subtropical Dry Broadleaf Forests",
                                          "Tropical & Subtropical Moist Broadleaf Forests"))

sites.sub$SSBS <- as.factor(sites.sub$SSBS)


############################################################
#                                                          #
#     Step 3: Input the production and fertiliser          #
#                 intensity information                    #
#                                                          #
############################################################


# total production per grid cell and production per hectare per grid cell (across all crops) 
# have been estimated globally from the EarthStat data.


# read in the relevent rasters and extract values to the predicts data tables
total.prod <- raster(paste0(datadir, "/Earthstat_prod_total_newmethod.tif"))
prodperhec <- raster(paste0(datadir, "/Earthstat_prod_perhec_newmethod.tif"))
total.fert <- raster(paste0(datadir, "/Earthstat_fert_total_newmethod.tif"))
fertperhec <- raster(paste0(datadir, "/Earthstat_fert_perhec_newmethod_17crops.tif"))


# convert the PREDICTS lat/longs into spatial points
sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]
sites.sub_xy <- SpatialPoints(sites.sub_xy, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))


# get the intersection values, what are the production values for each site?
sites.sub$prod.total <- extract(total.prod, sites.sub_xy, na.rm = FALSE)
sites.sub$prod.per.hec <- extract(prodperhec, sites.sub_xy, na.rm = FALSE)
sites.sub$fert.total <- extract(total.fert, sites.sub_xy, na.rm = FALSE)
sites.sub$fert.per.hec <- extract(fertperhec, sites.sub_xy, na.rm = FALSE)

# where there are NAs, change to 0
sites.sub[is.na(sites.sub$prod.total), "prod.total"] <- 0
sites.sub[is.na(sites.sub$prod.per.hec), "prod.per.hec"] <- 0
sites.sub[is.na(sites.sub$fert.total), "fert.total"] <- 0
sites.sub[is.na(sites.sub$fert.per.hec), "fert.per.hec"] <- 0



############################################################
#                                                          #
#       Step 4: Determine the number of crops              #
#                in each grid square                       #  
#                                                          #
############################################################

# add in the number of crops for each grid square using the raster generated from the earthstat data

# read in the ncrop raster that was made using the earthstat data
ncrops <- raster(paste0(datadir, "/Earthstat_NCrops_newmethod.tif"))

# get the intersection values, what are the production values for each site?
sites.sub$ncrop <- extract(ncrops, sites.sub_xy, na.rm = FALSE)

# replace NAs with 0
sites.sub[is.na(sites.sub$ncrop), "ncrop"] <- 0


############################################################
#                                                          #
#       Step 5: Organise EarthStat data on fraction        #
#                   of area harvested                      #
#                                                          #
############################################################


# use the raster layer creasted in the 1. Earthstat data organisation script
fract.totals <- raster(paste0(datadir, "/Earthstat_areaharv_total_frac_newmethod.tif"))

# get the intersection values, what are the production values for each site?
sites.sub$frac.harvested <- extract(fract.totals, sites.sub_xy, na.rm = FALSE)

# replace NAs with 0
sites.sub[is.na(sites.sub$frac.harvested), "frac.harvested"] <- 0


############################################################
#                                                          #
####   Step 6: Calculate distance to forest habitat     ####
#                                                          #
############################################################


# Distance of site to areas of 80% dense forest or greater has been determined elsewhere.
# Due to the size of the dataset and processing time, this could not be done here.

# read in the hansen distance
hans <- read.csv(paste0(datadir,"/hans_min_dist_80.csv"))

# change hansen column name
colnames(hans)[27] <- "Hansen_mindist"

# check for missing info
summary(sites.sub$SSBS %in% hans$SSBS) # 133 rows with no hansen info 

#cut down the hansen dataset
hans <- hans[, c("SSBS", "Hansen_mindist")]

# merge the datasets by site ID
sites.sub <- merge(sites.sub, hans, by = "SSBS", all.x = T)




############################################################
#                                                          #
####    Step 7: Calculate number of land covers in      ####
#             in buffer around predicts site               #
#                                                          #
############################################################

# this currently uses the Kehoe land use dataset, land cover classifications have been 
# merged to ignore the livestock and suitability elements.  

# read in the kehoe raster
landuse <- raster(paste0(datadir, "/", "GLS_v02.bj_RECLASSIFIED.tif"))

# convert the 0 to NA - this is water. But it shouldnt be counted as a landuse
landuse[landuse$GLS_v02.bj_RECLASSIFIED == 0] <- NA


# create the buffer options around the data points
sites.buffer.100 <- buffer(sites.sub_xy, width = 100, dissolve = FALSE)
sites.buffer.500 <- buffer(sites.sub_xy, width = 500, dissolve = FALSE)
sites.buffer.1k <- buffer(sites.sub_xy, width = 1000, dissolve = FALSE)
sites.buffer.3k <- buffer(sites.sub_xy, width = 3000, dissolve = FALSE)
sites.buffer.5k <- buffer(sites.sub_xy, width = 5000, dissolve = FALSE)


#mapview::mapview(sites.sub_xy)
#mapview::mapview(sites.buffer.5k)


# extract the raster values for the buffer sites
landcovers_extract.100 <- extract(landuse,  sites.buffer.100, df = TRUE)
landcovers_extract.500 <- extract(landuse, sites.buffer.500, df = TRUE)
landcovers_extract.1k <- extract(landuse, sites.buffer.1k, df = TRUE)
landcovers_extract.3k <- extract(landuse, sites.buffer.3k, df = TRUE)
landcovers_extract.5k <- extract(landuse, sites.buffer.5k, df = TRUE)




# Function to look at the number of unique landcovers for each ID. 

calc_landcovers <- function(x){
  
  results <- NULL
  
  for(i in unique(x$ID)){
    
    subtab <- x[x$ID == i, ]
    
    nlandcovers <- length(unique(subtab$layer))
    
    result <- c(i, nlandcovers)
    
    results <- rbind(results, result)
    
  }
  
  colnames(results) <- c("ID", "Nlandcovers")
  
  return(results)
  
}

# create the table for each set with n landcovers
res.100 <- calc_landcovers(landcovers_extract.100)
res.500 <- calc_landcovers(landcovers_extract.500)
res.1k <- calc_landcovers(landcovers_extract.1k)
res.3k <- calc_landcovers(landcovers_extract.3k)
res.5k <- calc_landcovers(landcovers_extract.5k)

# add these values to the sites dataset
sites.sub$landcovers.100 <- res.100[,2]
sites.sub$landcovers.500 <- res.500[,2]
sites.sub$landcovers.1k <- res.1k[,2]
sites.sub$landcovers.3k <- res.3k[,2]
sites.sub$landcovers.5k <- res.5k[,2]




##%######################################################%##
#                                                          #
####       Step 9: determine the homogeneity (and       ####
#        potentially other metrics) for each site          #
#                                                          #
##%######################################################%##

# This code organises the homogeneity metric from the Tuanmu 2015 paper:

# Tuanmu, M. N. & Jetz, W. A global, remote sensing-based characterization of terrestrial
# habitat heterogeneity for biodiversity and ecosystem modelling. 
# Glob. Ecol. Biogeogr. 24, 1329-1339 (2015).

# read in the required datasets
homogen <- raster(paste0(datadir, "/", "Homogeneity_01_05_5km_uint16.tif"))

# info from website hosting dataset:
# values of the data layers should be mulitplied by 0.0001 to obtain the actual values of the metrics.
homogen <- homogen*0.0001

# extract the dataset info for the PREDICTS sites
sites.sub$homogen <- extract(homogen, sites.sub_xy, na.rm = FALSE)

# how many NAs
nrow(sites.sub[is.na(sites.sub$homogen), ]) #45



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
sites.sub$percNH <- extract(percNH, sites.sub_xy, na.rm = FALSE)

# how many NAs
nrow(sites.sub[is.na(sites.sub$percNH),]) #3



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
# 2992      1353 



##%######################################################%##
#                                                          #
####      Merge data subsets with site information      ####
#                                                          #
##%######################################################%##

# all the data has now been added to the sites.sub dataset, now to merge
# this with the pols and pc subsets

# merge datasets to add in forest biomes for each site
sites.sub.pols <- merge(sites.sub.pols, sites.sub, by = "SSBS", all.x = TRUE)
sites.sub.pc <- merge(sites.sub.pc, sites.sub, by = "SSBS", all.x = TRUE)

# remove duplicate columns
sites.sub.pols <- sites.sub.pols[ , c(1:19, 38:55)]
colnames(sites.sub.pols) <- sub(".x", "", colnames(sites.sub.pols) )
sites.sub.pc <- sites.sub.pc[ , c(1:19, 38:55)]
colnames(sites.sub.pc) <- sub(".x", "", colnames(sites.sub.pc) )






##%######################################################%##
#                                                          #
####                 Assess the dataset                 ####
#                                                          #
##%######################################################%##


# remove any rows that have NA in the variable columns
summary(is.na(sites.sub))
summary(is.na(sites.sub.pols))
summary(is.na(sites.sub.pc))


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






##%######################################################%##
#                                                          #
####      checking correlations between variables       ####
#                                                          #
##%######################################################%##



# function to create plot
panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 2)
}

# save figure
pdf(file = paste0(outdir, "/Correlations_all_variables.pdf"), width =14, height = 9)

# correlations, including all nlandcovers buffers
pairs(sites.sub[ , c(21:36)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

pairs(sites.sub.pols[ , c(21:36)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

pairs(sites.sub.pc[ , c(21:36)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

dev.off()



# figure with reduced number of variables
pdf(file = paste0(outdir, "/Correlations_reduced_variables.pdf"), width =11, height = 9)

# correlations subset of variables, 5km landcovers buffer only, high pesticide estimates

pairs(sites.sub[ , c(21:27, 32:35)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

pairs(sites.sub.pols[ , c(21:27, 32:35)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

pairs(sites.sub.pc[ , c(21:27, 32:35)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

dev.off()


# save the untransformed datasets including all variables
save(sites.sub, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_INSECTS.rdata"))
save(sites.sub.pols, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_POLS.rdata"))
save(sites.sub.pc, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_PCS.rdata"))





##%######################################################%##
#                                                          #
####                Data transformations                ####
#                                                          #
##%######################################################%##

# subset columns
final.data <- sites.sub[, c(1,8:10, 14:37)]
final.data.pols <- sites.sub.pols[, c(1,8:10, 14:37)]
final.data.pc <- sites.sub.pc[, c(1,8:10, 14:37)]

final.data.trans <- final.data
final.data.trans.pols <- final.data.pols
final.data.trans.pc <- final.data.pc


# log transform some of the continuous variables where they are skewed
final.data.trans$fert.total_log <-log(final.data.trans$fert.total+1) # there are 0s so +1
final.data.trans$Hansen_mindist_log <-log(final.data.trans$Hansen_mindist+1) # there are 0s so +1
final.data.trans$pest_H_log <-log(final.data.trans$pest_H+1) # there are 0s so +1

final.data.trans.pols$fert.total_log <-log(final.data.trans.pols$fert.total+1) # there are 0s so +1
final.data.trans.pols$Hansen_mindist_log <-log(final.data.trans.pols$Hansen_mindist+1) # there are 0s so +1
final.data.trans.pols$pest_H_log <-log(final.data.trans.pols$pest_H+1) # there are 0s so +1

final.data.trans.pc$fert.total_log <-log(final.data.trans.pc$fert.total+1) # there are 0s so +1
final.data.trans.pc$Hansen_mindist_log <-log(final.data.trans.pc$Hansen_mindist+1) # there are 0s so +1
final.data.trans.pc$pest_H_log <-log(final.data.trans.pc$pest_H+1) # there are 0s so +1


# standardise all continuous variables
final.data.trans$landcovers.5k <- scale(final.data.trans$landcovers.5k)
final.data.trans$homogen <- scale(final.data.trans$homogen)
final.data.trans$fert.total_log <- scale(final.data.trans$fert.total_log)
final.data.trans$Hansen_mindist_log <-scale(final.data.trans$Hansen_mindist_log)
final.data.trans$percNH <-scale(final.data.trans$percNH)
final.data.trans$ncrop <-scale(final.data.trans$ncrop)
final.data.trans$pest_H_log <-scale(final.data.trans$pest_H_log)

final.data.trans.pols$landcovers.5k <- scale(final.data.trans.pols$landcovers.5k)
final.data.trans.pols$homogen <- scale(final.data.trans.pols$homogen)
final.data.trans.pols$fert.total_log <- scale(final.data.trans.pols$fert.total_log)
final.data.trans.pols$Hansen_mindist_log <-scale(final.data.trans.pols$Hansen_mindist_log)
final.data.trans.pols$percNH <-scale(final.data.trans.pols$percNH)
final.data.trans.pols$ncrop <-scale(final.data.trans.pols$ncrop)
final.data.trans.pols$pest_H_log <-scale(final.data.trans.pols$pest_H_log)

final.data.trans.pc$landcovers.5k <- scale(final.data.trans.pc$landcovers.5k)
final.data.trans.pc$homogen <- scale(final.data.trans.pc$homogen)
final.data.trans.pc$fert.total_log <- scale(final.data.trans.pc$fert.total_log)
final.data.trans.pc$Hansen_mindist_log <-scale(final.data.trans.pc$Hansen_mindist_log)
final.data.trans.pc$percNH <-scale(final.data.trans.pc$percNH)
final.data.trans.pc$ncrop <-scale(final.data.trans.pc$ncrop)
final.data.trans.pc$pest_H_log <-scale(final.data.trans.pc$pest_H_log)



# get data sections for the scaling info for plotting later
Hansen_mindist_log <-final.data.trans$Hansen_mindist_log
landcovers.5k <- final.data.trans$landcovers.5k
fert.total_log <- final.data.trans$fert.total_log
homogen <- final.data.trans$homogen
percNH <- final.data.trans$percNH
ncrop <- final.data.trans$ncrop
pest_H_log <- final.data.trans$pest_H_log

# save the scaling values for projections
h <- c("homogen", attr(homogen, "scaled:scale"), attr(homogen, "scaled:center"))
l <- c("landcovers.5k", attr(landcovers.5k, "scaled:scale"), attr(landcovers.5k, "scaled:center"))
f <- c("fert.total_log", attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center"))
d <- c("Hansen_mindist_log", attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center"))
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))
n <- c("ncrop", attr(ncrop, "scaled:scale"), attr(ncrop, "scaled:center"))
c <- c("pest_H", attr(pest_H_log, "scaled:scale"), attr(pest_H_log, "scaled:center"))

values <- rbind(d, f, l, h, p, n, c)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values_INSECTS.csv"), row.names = FALSE)



# get data sections for the scaling info for plotting later
Hansen_mindist_log <-final.data.trans.pols$Hansen_mindist_log
landcovers.5k <- final.data.trans.pols$landcovers.5k
fert.total_log <- final.data.trans.pols$fert.total_log
homogen <- final.data.trans.pols$homogen
percNH <- final.data.trans.pols$percNH
ncrop <- final.data.trans.pols$ncrop
pest_H_log <- final.data.trans.pols$pest_H_log

# save the scaling values for projections
h <- c("homogen", attr(homogen, "scaled:scale"), attr(homogen, "scaled:center"))
l <- c("landcovers.5k", attr(landcovers.5k, "scaled:scale"), attr(landcovers.5k, "scaled:center"))
f <- c("fert.total_log", attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center"))
d <- c("Hansen_mindist_log", attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center"))
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))
n <- c("ncrop", attr(ncrop, "scaled:scale"), attr(ncrop, "scaled:center"))
c <- c("pest_H", attr(pest_H_log, "scaled:scale"), attr(pest_H_log, "scaled:center"))

values <- rbind(d, f, l, h, p, n, c)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values_POLS.csv"), row.names = FALSE)



# get data sections for the scaling info for plotting later
Hansen_mindist_log <-final.data.trans.pc$Hansen_mindist_log
landcovers.5k <- final.data.trans.pc$landcovers.5k
fert.total_log <- final.data.trans.pc$fert.total_log
homogen <- final.data.trans.pc$homogen
percNH <- final.data.trans.pc$percNH
ncrop <- final.data.trans.pc$ncrop
pest_H_log <- final.data.trans.pc$pest_H_log

# save the scaling values for projections
h <- c("homogen", attr(homogen, "scaled:scale"), attr(homogen, "scaled:center"))
l <- c("landcovers.5k", attr(landcovers.5k, "scaled:scale"), attr(landcovers.5k, "scaled:center"))
f <- c("fert.total_log", attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center"))
d <- c("Hansen_mindist_log", attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center"))
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))
n <- c("ncrop", attr(ncrop, "scaled:scale"), attr(ncrop, "scaled:center"))
c <- c("pest_H", attr(pest_H_log, "scaled:scale"), attr(pest_H_log, "scaled:center"))

values <- rbind(d, f, l, h, p, n, c)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values_PC.csv"), row.names = FALSE)





### organise factor levels ###

# drop unused levels of factors
final.data.trans <- droplevels(final.data.trans)
final.data.trans.pols <- droplevels(final.data.trans.pols)
final.data.trans.pc <- droplevels(final.data.trans.pc)

# set factor levels in the best order
final.data.trans$Use_intensity <- relevel(final.data.trans$Use_intensity, ref = "Minimal use")
final.data.trans.pols$Use_intensity <- relevel(final.data.trans.pols$Use_intensity, ref = "Minimal use")
final.data.trans.pc$Use_intensity <- relevel(final.data.trans.pc$Use_intensity, ref = "Minimal use")

# nsites per use intensity
table(final.data.trans$Use_intensity)
table(final.data.trans.pols$Use_intensity)
table(final.data.trans.pc$Use_intensity)


# nsites per biome
table(final.data.trans$Forest_biome)
table(final.data.trans.pols$Forest_biome)
table(final.data.trans.pc$Forest_biome)

# set land use as character variable
final.data.trans$Predominant_land_use <- as.character(final.data.trans$Predominant_land_use)
final.data.trans.pols$Predominant_land_use <- as.character(final.data.trans.pols$Predominant_land_use)
final.data.trans.pc$Predominant_land_use <- as.character(final.data.trans.pc$Predominant_land_use)


# combine secondary land uses
final.data.trans$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans[final.data.trans$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"

final.data.trans.pols$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans.pols$Predominant_land_use)
final.data.trans.pols$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans.pols$Predominant_land_use)
final.data.trans.pols$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans.pols$Predominant_land_use)
final.data.trans.pols[final.data.trans.pols$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"

final.data.trans.pc$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans.pc$Predominant_land_use)
final.data.trans.pc$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans.pc$Predominant_land_use)
final.data.trans.pc$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans.pc$Predominant_land_use)
final.data.trans.pc[final.data.trans.pc$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"


table(final.data.trans$Predominant_land_use)
table(final.data.trans.pols$Predominant_land_use)
table(final.data.trans.pc$Predominant_land_use)


# set factor levels of predominant land use
final.data.trans$Predominant_land_use <- factor(final.data.trans$Predominant_land_use,
                                                levels=c("Primary vegetation","Secondary vegetation", "Cropland"))
final.data.trans.pols$Predominant_land_use <- factor(final.data.trans.pols$Predominant_land_use,
                                                levels=c("Primary vegetation","Secondary vegetation", "Cropland"))
final.data.trans.pc$Predominant_land_use <- factor(final.data.trans.pc$Predominant_land_use,
                                                levels=c("Primary vegetation","Secondary vegetation", "Cropland"))


# nsites per land use
table(final.data.trans$Predominant_land_use)
table(final.data.trans.pols$Predominant_land_use)
table(final.data.trans.pc$Predominant_land_use)


pdf(file = paste0(outdir, "/Correlations_final_variables.pdf"), width =9, height = 9)
# correlations of final set of variables - transformed
pairs(final.data.trans[ , c(16, 23:25, 29:31)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

pairs(final.data.trans.pols[ , c(16, 23:25, 29:31)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

pairs(final.data.trans.pc[ , c(16, 23:25, 29:31)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

dev.off()



# save transformed dataset
save(final.data.trans, file = paste0(outdir, "/PREDICTS_dataset_TRANS_INSECTS.rdata"))
save(final.data.trans.pols, file = paste0(outdir, "/PREDICTS_dataset_TRANS_POLS.rdata"))
save(final.data.trans.pc, file = paste0(outdir, "/PREDICTS_dataset_TRANS_INSECTS_PC.rdata"))





##%######################################################%##
#                                                          #
####                 Get data summaries                 ####
#                                                          #
##%######################################################%##

# get the number of studies, sites, etc for the paper

length(unique(final.data.trans$SS)) # 209 studies
nrow(final.data.trans) # 3621 sites

length(unique(final.data.trans.pols$SS)) # 142 studies
nrow(final.data.trans.pols) # 2902 sites

length(unique(final.data.pc$SS)) # 99 studies
nrow(final.data.pc) # 1454 sites



##%######################################################%##
#                                                          #
####   Plot of site distribution across forest biomes   ####
#                                                          #
##%######################################################%##



# extract the PREDICTS points
plot_data <- final.data.trans[, c("SS", "SSBS", "Longitude", "Latitude")]

# look at nsites per study
nsites <- as.matrix(table(plot_data$SS))
nsites <- cbind(rownames(nsites), nsites[,1])

# Get the lat/long per study
lon <- unlist(apply(nsites, MARGIN = 1, FUN = function(x){plot_data[plot_data$SS == x[1], 3][1]}))
lat <- unlist(apply(nsites, MARGIN = 1, FUN = function(x){plot_data[plot_data$SS == x[1], 4][1]}))

nsites <- cbind(nsites, lon, lat)

nsites <- as.data.frame(nsites)

colnames(nsites) <- c("SS", "nsites", "lon", "lat")

nsites$nsites <- as.numeric(as.character(nsites$nsites))
nsites$lon <- as.numeric(as.character(nsites$lon))
nsites$lat <- as.numeric(as.character(nsites$lat))

# select the forest biomes only
ecobio <- ecobio[ecobio$BIOME %in% c(1:6, 12), 'BIOME']

ecobio$BIOME <- sub(12, "Mediterranean Forests, Woodlands & Scrub", ecobio$BIOME)
ecobio$BIOME <- sub(1, "Tropical & Subtropical Moist Broadleaf Forests", ecobio$BIOME)
ecobio$BIOME <- sub(2, "Tropical & Subtropical Dry Broadleaf Forests", ecobio$BIOME)
ecobio$BIOME <- sub(3, "Tropical & Subtropical Coniferous Forests", ecobio$BIOME)
ecobio$BIOME <- sub(4, "Temperate Broadleaf & Mixed Forests", ecobio$BIOME)
ecobio$BIOME<- sub(5, "Temperate Conifer Forests", ecobio$BIOME)
ecobio$BIOME <- sub(6, "Boreal Forests/Taiga", ecobio$BIOME)


# plot the raster in ggplot
map.world <- map_data('world')

# for colourblind pallette work around
n = length(unique(ecobio$BIOME))

# plot of predicts sites across biomes including size per n sites
ggplot()+
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "lightgrey", colour="lightgrey", size=0.2) +
  geom_sf(data = ecobio, aes(fill = BIOME), alpha = 0.7, col = NA) +
  geom_point(data = nsites, aes(x = lon, y = lat, size = nsites), col = c("#0F0F0F"), alpha = 0.5) +
  theme(panel.grid.major = element_line(colour = "transparent"), 
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size = 16)) +
  xlab("") +
  ylab("") +
  scale_size_continuous(range = c(0.2, 5), breaks = c(1, 10, 50, 100, 200)) +
  guides(fill=guide_legend(nrow=7,byrow=TRUE), size = guide_legend(nrow = 7, byrow = T)) +
  scale_fill_manual(breaks = ecobio$BIOME, values = colorblind_pal()(n + 1)[-1])


ggsave(filename = paste0(outdir, "/MAP_Predicts_points_biome_INSECTS.pdf"),
       plot = last_plot(),
       width = 8,
       height = 6)






