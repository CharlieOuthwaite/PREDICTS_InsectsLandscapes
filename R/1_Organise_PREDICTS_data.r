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
####                 Assess the dataset                 ####
#                                                          #
##%######################################################%##


# remove any rows that have NA in the variable columns
summary(is.na(sites.sub))

# remove rows with NAs for any variable of interest
sites.sub <- sites.sub[!is.na(sites.sub$Hansen_mindist), ] # 133 NAs 4212 rows
sites.sub <- sites.sub[!is.na(sites.sub$homogen), ] # 45 NAs 4170 rows
sites.sub <- sites.sub[!is.na(sites.sub$percNH), ] # 3 NAs 4170 rows

# nrows of dataset
nrow(sites.sub) # 4170

# remove those sites that have "Cannot decide" as a use intensity
nrow(sites.sub[sites.sub$Use_intensity == "Cannot decide", ]) # 550
sites.sub <- sites.sub[!sites.sub$Use_intensity == "Cannot decide", ] 

nrow(sites.sub) # 3620


### checking correlations between variables ###


# function to create plot
panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 2)
}


pdf(file = paste0(outdir, "/Correlations_all_insects.pdf"), width =14, height = 9)
# correlations, including all nlandcovers buffers
pairs(sites.sub[ , c(21:34)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

dev.off()

pdf(file = paste0(outdir, "/Correlations_subset_insects.pdf"), width =11, height = 9)
# correlations subset of variables, 5km landcpvers buffer only
pairs(sites.sub[ , c(21:27, 32:34)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)
dev.off()


# save the untransformed dataset including all variables
save(sites.sub, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_INSECTS.rdata"))





##%######################################################%##
#                                                          #
####                Data transformations                ####
#                                                          #
##%######################################################%##

# subset columns
final.data <- sites.sub[, c(1,8:10, 14:35)]

final.data.trans <- final.data


# log transform some of the continuous variables where they are skewed
final.data.trans$fert.total_log <-log(final.data.trans$fert.total+1) # there are 0s so +1
final.data.trans$Hansen_mindist_log <-log(final.data.trans$Hansen_mindist+1) # there are 0s so +1

# standardise all continuous variables
final.data.trans$landcovers.5k <- scale(final.data.trans$landcovers.5k)
final.data.trans$homogen <- scale(final.data.trans$homogen)
final.data.trans$fert.total_log <- scale(final.data.trans$fert.total_log)
final.data.trans$Hansen_mindist_log <-scale(final.data.trans$Hansen_mindist_log)
final.data.trans$percNH <-scale(final.data.trans$percNH)


# get data sections for the scaling info for plotting later
Hansen_mindist_log <-final.data.trans$Hansen_mindist_log
landcovers.5k <- final.data.trans$landcovers.5k
fert.total_log <- final.data.trans$fert.total_log
homogen <- final.data.trans$homogen
percNH <- final.data.trans$percNH


# save the scaling values for projections
h <- c("homogen", attr(homogen, "scaled:scale"), attr(homogen, "scaled:center"))
l <- c("landcovers.5k", attr(landcovers.5k, "scaled:scale"), attr(landcovers.5k, "scaled:center"))
f <- c("fert.total_log", attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center"))
d <- c("Hansen_mindist_log", attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center"))
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))

values <- rbind(d, f, l, h, p)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values.csv"), row.names = FALSE)


### organise factor levels ###

# drop unused levels of factors
final.data.trans <- droplevels(final.data.trans)

# set factor levels in the best order
final.data.trans$Use_intensity <- relevel(final.data.trans$Use_intensity, ref = "Minimal use")

# nsites per use intensity
table(final.data.trans$Use_intensity)


# set the forest biome to names rather than numbers
final.data.trans$Forest_biome <- sub(12, "Mediterranean Forests, Woodlands & Scrub", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(1, "Tropical & Subtropical Moist Broadleaf Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(2, "Tropical & Subtropical Dry Broadleaf Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(3, "Tropical & Subtropical Coniferous Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(4, "Temperate Broadleaf & Mixed Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(5, "Temperate Conifer Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(6, "Boreal Forests/Taiga", final.data.trans$Forest_biome)



final.data.trans$Forest_biome <- factor(final.data.trans$Forest_biome,
                                        levels=c("Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests",
                                                 "Boreal Forests/Taiga", "Mediterranean Forests, Woodlands & Scrub",
                                                 "Tropical & Subtropical Coniferous Forests", 
                                                 "Tropical & Subtropical Dry Broadleaf Forests",
                                                 "Tropical & Subtropical Moist Broadleaf Forests"))
# nsites per biome
table(final.data.trans$Forest_biome)

# set land use as character variable
final.data.trans$Predominant_land_use <- as.character(final.data.trans$Predominant_land_use)


# combine secondary land uses
final.data.trans$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans[final.data.trans$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"

table(final.data.trans$Predominant_land_use)


# set factor levels of predominant land use
final.data.trans$Predominant_land_use <- factor(final.data.trans$Predominant_land_use,
                                                levels=c("Primary vegetation","Secondary vegetation", "Cropland"))


# nsites per land use
table(final.data.trans$Predominant_land_use)


pdf(file = paste0(outdir, "/Correlations_final_variables_INSECTS.pdf"), width =9, height = 9)
# correlations of final set of variables - transformed
pairs(final.data.trans[ , c(23:25, 27:28)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)
dev.off()



# save transformed dataset
save(final.data.trans, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_TRANS_INSECTS.rdata"))




##%######################################################%##
#                                                          #
####                  Dataset subsets                   ####
#                                                          #
##%######################################################%##

# Subsets of the data required for the abundance models and the range size models

## abundance subset
final.data.abun <- final.data.trans[!is.na(final.data.trans$Total_abundance), ] # 3385 rows

# log the abundance values
final.data.abun$logAbun <- log(final.data.abun$Total_abundance+1)


# Incorporate the RCAR data

# Set the path to your local copy of the database
predicts.path <- paste0(datadir, "/PREDICTS_Sites_Mean_RangeSizes.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)

# Use unique SSBS values to match up the two datasets.
final.data.rcar<- merge(final.data.trans, predicts, by = "SSBS", all.x= T)

# how many of these have values
sum(!is.na(final.data.rcar$RCAR_110km)) # 2404

# remove duplicated columns
final.data.rcar <- final.data.rcar[ , c(1:28, 37)]

# remove NAs
final.data.rcar <- final.data.rcar[!is.na(final.data.rcar$RCAR_110km), ] # 2404 rows

# remove the .x from column names
colnames(final.data.rcar) <- sub(".x", "", colnames(final.data.rcar))

final.data.abun <- droplevels(final.data.abun)
final.data.rcar <- droplevels(final.data.rcar)

# save
save(final.data.abun, file = paste0(outdir, "/PREDICTS_abun_subset_INSECTS.rdata"))
save(final.data.rcar, file = paste0(outdir, "/PREDICTS_rcar_subset_INSECTS.rdata"))



##%######################################################%##
#                                                          #
####                 Get data summaries                 ####
#                                                          #
##%######################################################%##

# get the number of studies, sites, etc for the paper

length(unique(final.data.trans$SS)) # 209 studies
nrow(final.data.trans) # 3620 sites

length(unique(final.data.abun$SS)) # 192 studies
nrow(final.data.abun) # 3385 sites

length(unique(final.data.rcar$SS)) # 133 studies
nrow(final.data.rcar) # 2404 sites



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






