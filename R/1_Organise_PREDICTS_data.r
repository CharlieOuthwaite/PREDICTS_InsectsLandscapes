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



# read in the species lists for pollinators and pest controllers

polls <- 
  
pcs <- read.csv(paste0(datadir, "/pc_database.csv"))  


# create subsets for pollinators and pestcontrollers
pred.data.pol <- 
  
pred.data.pc <-   


### organise using functions from predictsFunctions package ###

#### 1. Insect subset  ####
  
# correct sampling effort 
pred.data <- CorrectSamplingEffort(pred.data)

# merge sites: this combines potential subsamples within one site
pred.data <- MergeSites(pred.data) # 826525 rows

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(pred.data, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 7800 rows

# only interested in natural habitats plus cropland, drop other land uses #

# site level data cropland only
sites.sub <- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]

# remove sites with NA in lat/long columns
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ] # 5115 rows


#### 2. Pollinators subset ####






#### 3. Pest controllers subset ####









############################################################
#                                                          #
#      Step 2: Subset data to those in forest biomes       #
#                                                          #
############################################################


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

pointsST_NA

# check the NA ones manually on the map
pointsST_NA[rownames(pointsST_NA) == 3485, "Forest_biome"] <- 14
pointsST_NA[rownames(pointsST_NA) == 3486, "Forest_biome"] <- 14
pointsST_NA[rownames(pointsST_NA) == 3493, "Forest_biome"] <- 14
pointsST_NA[rownames(pointsST_NA) == 3494, "Forest_biome"] <- 14


# add these in to the sites.sub table, theyre in the same order
sites.sub[is.na(sites.sub$Forest_biome), "Forest_biome"] <- pointsST_NA$Forest_biome


# sites per biome
table(sites.sub$Forest_biome)

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


