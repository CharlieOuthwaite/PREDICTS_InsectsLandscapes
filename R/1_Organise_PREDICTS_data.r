##%######################################################%##
#                                                          #
####             1. Organise PREDICTS data              ####
#                                                          #
##%######################################################%##


# This script takes the PREDICTS dataset and subsets it into insects only.



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
outdir <- "/1_PREDICTS_PLUS_VARIABLES/"
if(!dir.exists(outdir)) dir.create(outdir)

# read in the complete PREDICTS dataset
pred.data <- readRDS(paste0(datadir, "/database.rds")) # 3250404 rows

# subset to the insects only
pred.data <- pred.data[(pred.data$Class=="Insecta"),] # 935078

### organise datasets using functions from predictsFunctions package ###

# correct sampling effort 
pred.data <- CorrectSamplingEffort(pred.data)

# merge sites: this combines potential subsamples within one site
pred.data <- MergeSites(pred.data) # 826525 rows

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(pred.data, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 7800 rows


# only interested in natural habitats plus cropland, drop other land uses #
# site level data cropland only
# sites.sub <- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]
# sites.sub.pols <- pred.sites.metrics.pols[!pred.sites.metrics.pols$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]
# sites.sub.pc <- pred.sites.metrics.pc[!pred.sites.metrics.pc$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]

# which land uses are we interested in, I have kept plantation in here too
sites.sub <- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide"), ]

# 5612 rows


# remove sites with NA in lat/long columns
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ] # 5606 rows


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

# add in the number of crops for each grid square using the raster generated from the earthstat data
# convert the PREDICTS lat/longs into spatial points

## all insects dataset ##

# get coordinates for sites
sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]
sites.sub_xy <- SpatialPoints(sites.sub_xy, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))
 
# read in the ncrop raster that was made using the earthstat data
ncrops <- raster(paste0(datadir, "/Earthstat_NCrops_newmethod.tif"))

# get the intersection values, what are the production values for each site?
sites.sub$ncrop <- extract(ncrops, sites.sub_xy, na.rm = FALSE)

# replace NAs with 0
sites.sub[is.na(sites.sub$ncrop), "ncrop"] <- 0 # 5606 rows



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
sites.sub <- sites.sub[!is.na(sites.sub$percNH), ] # 3 NAs 5603 rows

# remove those sites that have "Cannot decide" as a use intensity
nrow(sites.sub[sites.sub$Use_intensity == "Cannot decide", ]) # 616
sites.sub <- sites.sub[!sites.sub$Use_intensity == "Cannot decide", ] 

# nrows of dataset
nrow(sites.sub) # 4987



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
pairs(sites.sub[ , c(22:26)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

dev.off()

# only high ones are the high and low pesticide (as expected) and NH with field size (-0.42).



# save the untransformed datasets including all variables
save(sites.sub, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_INSECTS.rdata"))



##%######################################################%##
#                                                          #
####                Data transformations                ####
#                                                          #
##%######################################################%##

# subset columns
final.data.trans <- sites.sub

hist(final.data.trans$ncrop)
hist(final.data.trans$pest_H)
hist(final.data.trans$pest_L)
hist(final.data.trans$percNH)

hist(log(final.data.trans$pest_H + 1))
hist(log(final.data.trans$pest_L+ 1))

# log transform some of the continuous variables where they are skewed
# final.data.trans$pest_H_log <-log(final.data.trans$pest_H+1) # there are 0s so +1
# final.data.trans$pest_L_log <-log(final.data.trans$pest_L+1) # there are 0s so +1



# standardise all continuous variables

final.data.trans$percNH_RS <-scale(final.data.trans$percNH)
final.data.trans$ncrop_RS <-scale(final.data.trans$ncrop)
final.data.trans$pest_L_RS <-scale(final.data.trans$pest_L)
final.data.trans$pest_H_RS <-scale(final.data.trans$pest_H)



# get data sections for the scaling info for plotting later
percNH <- final.data.trans$percNH_RS
ncrop <- final.data.trans$ncrop_RS
pest_L<- final.data.trans$pest_L_RS
pest_H <- final.data.trans$pest_H_RS

# save the scaling values for projections
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))
n <- c("ncrop", attr(ncrop, "scaled:scale"), attr(ncrop, "scaled:center"))
cl <- c("pest_L", attr(pest_L, "scaled:scale"), attr(pest_L, "scaled:center"))
ch <- c("pest_H", attr(pest_H, "scaled:scale"), attr(pest_H, "scaled:center"))

values <- rbind(p, n, cl, ch)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values_INSECTS.csv"), row.names = FALSE)


### organise factor levels ###

# drop unused levels of factors
final.data.trans <- droplevels(final.data.trans)

# set factor levels in the best order
final.data.trans$Use_intensity <- relevel(final.data.trans$Use_intensity, ref = "Minimal use")

# nsites per use intensity
table(final.data.trans$Use_intensity)

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
                                                levels=c("Primary vegetation","Secondary vegetation", "Plantation forest", "Cropland"))

# nsites per land use
table(final.data.trans$Predominant_land_use)


pdf(file = paste0(outdir, "/Correlations_final_variables.pdf"), width =9, height = 9)
# correlations of final set of variables - transformed
pairs(final.data.trans[ , c(26:31)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

dev.off()



# save transformed dataset
save(final.data.trans, file = paste0(outdir, "/PREDICTS_dataset_TRANS_INSECTS.rdata"))




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
ecobio$BIOME <- sub(5, "Temperate Conifer Forests", ecobio$BIOME)
ecobio$BIOME <- sub(6, "Boreal Forests/Taiga", ecobio$BIOME)


# plot the raster in ggplot
map.world <- map_data('world')

# for colourblind pallette work around
n = length(unique(ecobio$BIOME))

# plot of predicts sites across biomes including size per n sites
ggplot()+
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= c("#7CCD7C"), colour= "transparent", size=0.2, alpha = 0.5) +
  geom_sf(data = ecobio,  aes(fill = BIOME), alpha = 1, col = NA) +
  geom_point(data = nsites, aes(x = lon, y = lat, size = nsites), col = c("#0F0F0F"), alpha = 0.5) +
  #geom_polygon(data = wwfbiomes, x = )
  theme(panel.grid.major = element_line(colour = "transparent"), 
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
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






