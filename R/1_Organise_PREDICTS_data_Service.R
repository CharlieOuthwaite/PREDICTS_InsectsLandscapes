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

# read in the ncrop raster that was made using the earthstat data
ncrops <- raster(paste0(datadir, "/Earthstat_NCrops_newmethod.tif"))


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
nrow(sites.sub.pc[is.na(sites.sub.pc$percNH),]) #0




##%######################################################%##
#                                                          #
####             Pesticide application rate             ####
#                                                          #
##%######################################################%##


# read in the raster
pest_H <- raster(paste0(datadir,"/Pesticide_totalAPR_High.tif"))
pest_L <- raster(paste0(datadir,"/Pesticide_totalAPR_Low.tif"))

# extract the dataset info for the PREDICTS sites
sites.sub.pols$pest_H <- extract(pest_H, sites.sub_xy_pols, na.rm = FALSE)
sites.sub.pc$pest_H <- extract(pest_H, sites.sub_xy_pc, na.rm = FALSE)

sites.sub.pols$pest_L <- extract(pest_L, sites.sub_xy_pols, na.rm = FALSE)
sites.sub.pc$pest_L <- extract(pest_L, sites.sub_xy_pc, na.rm = FALSE)

# how many NAs
nrow(sites.sub.pols[is.na(sites.sub.pols$pest_H),]) #0
nrow(sites.sub.pc[is.na(sites.sub.pc$pest_H),]) #0

nrow(sites.sub.pols[is.na(sites.sub.pols$pest_L),]) #0
nrow(sites.sub.pc[is.na(sites.sub.pc$pest_L),]) #0



##%######################################################%##
#                                                          #
####                     Field Size                     ####
#                                                          #
##%######################################################%##
# see file "legend items" in data folder for info on coding of data in the raster.

# read in the raster
fields <- raster(paste0(datadir,"/Global Field Sizes/dominant_field_size_categories.tif"))

# extract the info for the predicts sites
sites.sub.pols$fields <- extract(fields, sites.sub_xy_pols, na.rm = FALSE)
sites.sub.pc$fields <- extract(fields, sites.sub_xy_pc, na.rm = FALSE)

# how many NAs
nrow(sites.sub.pols[is.na(sites.sub.pols$fields),]) #0
nrow(sites.sub.pc[is.na(sites.sub.pc$fields),]) #0

# sites.sub$fields <- as.factor(sites.sub$fields)
# plot(sites.sub$fields)



##%######################################################%##
#                                                          #
####                   Tropical sites                   ####
#                                                          #
##%######################################################%##


# get the tropical values
sites.sub.pols$Tropical <- NA
sites.sub.pc$Tropical <- NA

sites.sub.pols[sites.sub.pols$Latitude > -23.44 & sites.sub.pols$Latitude < 23.44, 'Tropical'] <- "Tropical"
sites.sub.pc[sites.sub.pc$Latitude > -23.44 & sites.sub.pc$Latitude < 23.44, 'Tropical'] <- "Tropical"

# label the remaining as temperate
sites.sub.pols[is.na(sites.sub.pols$Tropical), 'Tropical'] <- "Temperate"
sites.sub.pc[is.na(sites.sub.pc$Tropical), 'Tropical'] <- "Temperate"

# set as a factor
sites.sub.pols$Tropical <- as.factor(sites.sub.pols$Tropical)
sites.sub.pc$Tropical <- as.factor(sites.sub.pc$Tropical)
# levels: Temperate Tropical

table(sites.sub.pols$Tropical)
table(sites.sub.pc$Tropical)

# Temperate  Tropical 
# 2986      1063 
# 
# Temperate  Tropical 
# 1299      1416 





##%######################################################%##
#                                                          #
####                 Assess the dataset                 ####
#                                                          #
##%######################################################%##


# remove any rows that have NA in the variable columns
summary(is.na(sites.sub.pols))
summary(is.na(sites.sub.pc))


# remove rows with NAs for any variable of interest
sites.sub.pols <- sites.sub.pols[!is.na(sites.sub.pols$percNH), ] 

# remove those sites that have "Cannot decide" as a use intensity
sites.sub.pols <- sites.sub.pols[!sites.sub.pols$Use_intensity == "Cannot decide", ] 
sites.sub.pc <- sites.sub.pc[!sites.sub.pc$Use_intensity == "Cannot decide", ] 

# nrows of dataset
nrow(sites.sub.pols) # 3577
nrow(sites.sub.pc) # 2312



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
pdf(file = paste0(outdir, "/Correlations_all_variables_POLS.pdf"), width =14, height = 9)

# correlations, including all nlandcovers buffers
pairs(sites.sub.pols[ , c(22:26)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

dev.off()

# save figure
pdf(file = paste0(outdir, "/Correlations_all_variables_PCS.pdf"), width =14, height = 9)

# correlations, including all nlandcovers buffers
pairs(sites.sub.pc[ , c(22:26)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex = 2)

dev.off()



# save the untransformed datasets including all variables
save(sites.sub.pols, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_POLS.rdata"))
save(sites.sub.pc, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_PCS.rdata"))

# load(file = paste0(outdir, "/PREDICTS_dataset_inc_variables_POLS.rdata"))
# load(file = paste0(outdir, "/PREDICTS_dataset_inc_variables_PCS.rdata"))


##%######################################################%##
#                                                          #
####                Data transformations                ####
#                                                          #
##%######################################################%##

# subset columns
final.data.trans.pols <- sites.sub.pols
final.data.trans.pc <- sites.sub.pc

hist(final.data.trans.pols$ncrop)
hist(final.data.trans.pols$pest_H)
hist(final.data.trans.pols$pest_L)
hist(final.data.trans.pols$percNH)

hist(log(final.data.trans.pols$pest_H + 1))
hist(log(final.data.trans.pols$pest_L+ 1))

# log transform some of the continuous variables where they are skewed
# final.data.trans$pest_H_log <-log(final.data.trans$pest_H+1) # there are 0s so +1
# final.data.trans$pest_L_log <-log(final.data.trans$pest_L+1) # there are 0s so +1


hist(final.data.trans.pc$ncrop)
hist(final.data.trans.pc$pest_H)
hist(final.data.trans.pc$pest_L)
hist(final.data.trans.pc$percNH)

hist(log(final.data.trans.pc$pest_H + 1))
hist(log(final.data.trans.pc$pest_L+ 1))

# log transform some of the continuous variables where they are skewed
# final.data.trans$pest_H_log <-log(final.data.trans$pest_H+1) # there are 0s so +1
# final.data.trans$pest_L_log <-log(final.data.trans$pest_L+1) # there are 0s so +1

# standardise all continuous variables

final.data.trans.pols$percNH_RS <-scale(final.data.trans.pols$percNH)
final.data.trans.pols$ncrop_RS <-scale(final.data.trans.pols$ncrop)
final.data.trans.pols$pest_L_RS <-scale(final.data.trans.pols$pest_L)
final.data.trans.pols$pest_H_RS <-scale(final.data.trans.pols$pest_H)

hist(final.data.trans.pc$ncrop_RS)
hist(final.data.trans.pc$pest_H_RS)
hist(final.data.trans.pc$pest_L_RS)
hist(final.data.trans.pc$percNH_RS)

final.data.trans.pc$percNH_RS <-scale(final.data.trans.pc$percNH)
final.data.trans.pc$ncrop_RS <-scale(final.data.trans.pc$ncrop)
final.data.trans.pc$pest_L_RS <-scale(final.data.trans.pc$pest_L)
final.data.trans.pc$pest_H_RS <-scale(final.data.trans.pc$pest_H)



# get data sections for the scaling info for plotting later
percNH <- final.data.trans.pols$percNH_RS
ncrop <- final.data.trans.pols$ncrop_RS
pest_L<- final.data.trans.pols$pest_L_RS
pest_H <- final.data.trans.pols$pest_H_RS

# save the scaling values for projections
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))
n <- c("ncrop", attr(ncrop, "scaled:scale"), attr(ncrop, "scaled:center"))
cl <- c("pest_L", attr(pest_L, "scaled:scale"), attr(pest_L, "scaled:center"))
ch <- c("pest_H", attr(pest_H, "scaled:scale"), attr(pest_H, "scaled:center"))

values <- rbind(p, n, cl, ch)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values_POLS.csv"), row.names = FALSE)


# get data sections for the scaling info for plotting later
percNH <- final.data.trans.pc$percNH_RS
ncrop <- final.data.trans.pc$ncrop_RS
pest_L<- final.data.trans.pc$pest_L_RS
pest_H <- final.data.trans.pc$pest_H_RS

# save the scaling values for projections
p <- c("percNH", attr(percNH, "scaled:scale"), attr(percNH, "scaled:center"))
n <- c("ncrop", attr(ncrop, "scaled:scale"), attr(ncrop, "scaled:center"))
cl <- c("pest_L", attr(pest_L, "scaled:scale"), attr(pest_L, "scaled:center"))
ch <- c("pest_H", attr(pest_H, "scaled:scale"), attr(pest_H, "scaled:center"))

values <- rbind(p, n, cl, ch)
colnames(values) <- c("variable", "scale", "centre")

#save
write.csv(values, paste0(outdir, "/Scaling_values_PCS.csv"), row.names = FALSE)



### organise factor levels ###

# drop unused levels of factors
final.data.trans.pols <- droplevels(final.data.trans.pols)
final.data.trans.pc <- droplevels(final.data.trans.pc)

# set factor levels in the best order
final.data.trans.pols$Use_intensity <- relevel(final.data.trans.pols$Use_intensity, ref = "Minimal use")
final.data.trans.pc$Use_intensity <- relevel(final.data.trans.pc$Use_intensity, ref = "Minimal use")

# nsites per use intensity
table(final.data.trans.pols$Use_intensity)
table(final.data.trans.pc$Use_intensity)

# set land use as character variable
final.data.trans.pols$Predominant_land_use <- as.character(final.data.trans.pols$Predominant_land_use)
final.data.trans.pc$Predominant_land_use <- as.character(final.data.trans.pc$Predominant_land_use)


# combine secondary land uses
final.data.trans.pols$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans.pols$Predominant_land_use)
final.data.trans.pols$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans.pols$Predominant_land_use)
final.data.trans.pols$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans.pols$Predominant_land_use)
final.data.trans.pols[final.data.trans.pols$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"

final.data.trans.pc$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans.pc$Predominant_land_use)
final.data.trans.pc$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans.pc$Predominant_land_use)
final.data.trans.pc$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans.pc$Predominant_land_use)
final.data.trans.pc[final.data.trans.pc$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"


table(final.data.trans.pols$Predominant_land_use)
table(final.data.trans.pc$Predominant_land_use)



# combine plantation and cropland sites
final.data.trans.pols$LU <- final.data.trans.pols$Predominant_land_use
final.data.trans.pc$LU <- final.data.trans.pc$Predominant_land_use

# set cropland and plantation to "Agriculture"
final.data.trans.pols$LU[final.data.trans.pols$LU %in% c("Cropland", "Plantation")] <- "Agriculture"
final.data.trans.pc$LU[final.data.trans.pc$LU %in% c("Cropland", "Plantation")] <- "Agriculture"

# set factor levels of predominant land use
final.data.trans.pols$LU <- factor(final.data.trans.pols$LU,
                              levels=c("Primary vegetation","Secondary vegetation", "Agriculture"))
final.data.trans.pc$LU <- factor(final.data.trans.pc$LU,
                                   levels=c("Primary vegetation","Secondary vegetation", "Agriculture"))

# nsites per land use
table(final.data.trans.pols$LU)
table(final.data.trans.pc$LU)

# Primary vegetation Secondary vegetation          Agriculture 
#                959                 1055                 1395
# 
# Primary vegetation Secondary vegetation          Agriculture 
#                978                  758                  294 

pdf(file = paste0(outdir, "/Correlations_final_variables_POLS.pdf"), width =9, height = 9)
# correlations of final set of variables - transformed
pairs(final.data.trans.pols[ , c(26:31)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

dev.off()
pdf(file = paste0(outdir, "/Correlations_final_variables_PCS.pdf"), width =9, height = 9)
# correlations of final set of variables - transformed
pairs(final.data.trans.pols[ , c(26:31)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "",
      cex.labels = 1.3)

dev.off()

# edit field size names
final.data.trans.pols$fields <- as.character(final.data.trans.pols$fields)
final.data.trans.pols$fields[final.data.trans.pols$fields == "0"] <- "No fields"
final.data.trans.pols$fields[final.data.trans.pols$fields == "3502"] <- "Very large"
final.data.trans.pols$fields[final.data.trans.pols$fields == "3503"] <- "Large"
final.data.trans.pols$fields[final.data.trans.pols$fields == "3504"] <- "Medium"
final.data.trans.pols$fields[final.data.trans.pols$fields == "3505"] <- "Small"
final.data.trans.pols$fields[final.data.trans.pols$fields == "3506"] <- "Very small"

table(final.data.trans.pols$fields)

# Large     Medium  No fields      Small Very small 
# 187        955       1653        572        210

final.data.trans.pols$fields <- factor(final.data.trans.pols$fields, levels = c("No fields", "Very small", "Small", "Medium", "Large"))

# edit field size names
final.data.trans.pc$fields <- as.character(final.data.trans.pc$fields)
final.data.trans.pc$fields[final.data.trans.pc$fields == "0"] <- "No fields"
final.data.trans.pc$fields[final.data.trans.pc$fields == "3502"] <- "Very large"
final.data.trans.pc$fields[final.data.trans.pc$fields == "3503"] <- "Large"
final.data.trans.pc$fields[final.data.trans.pc$fields == "3504"] <- "Medium"
final.data.trans.pc$fields[final.data.trans.pc$fields == "3505"] <- "Small"
final.data.trans.pc$fields[final.data.trans.pc$fields == "3506"] <- "Very small"

table(final.data.trans.pc$fields)

# Large     Medium  No fields      Small Very small 
# 101        339       1547        213        112

final.data.trans.pc$fields <- factor(final.data.trans.pc$fields, levels = c("No fields", "Very small", "Small", "Medium", "Large"))




# transform abundance
final.data.trans.pols$logAbun <- log(final.data.trans.pols$Total_abundance + 1)
final.data.trans.pc$logAbun <- log(final.data.trans.pc$Total_abundance + 1)


# save transformed dataset
save(final.data.trans.pols, file = paste0(outdir, "/PREDICTS_dataset_TRANS_POLS.rdata"))
save(final.data.trans.pc, file = paste0(outdir, "/PREDICTS_dataset_TRANS_PCS.rdata"))




##%######################################################%##
#                                                          #
####                 Get data summaries                 ####
#                                                          #
##%######################################################%##

# get the number of studies, sites, etc for the paper

length(unique(final.data.trans.pols$SS)) # 149 studies
nrow(final.data.trans.pols) # 3577 sites

length(unique(final.data.trans.pc$SS)) # 125 studies
nrow(final.data.trans.pc) # 2312 sites




##%######################################################%##
#                                                          #
####            Plot points on a global map             ####
#                                                          #
##%######################################################%##


# extract the PREDICTS points
plot_data <- final.data.trans.pols[, c("SS", "SSBS", "Longitude", "Latitude")]

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


# plot the raster in ggplot
map.world <- map_data('world')


# plot of predicts sites across biomes including size per n sites
ggplot()+
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= c("#7CCD7C"), colour= "transparent", size=0.2, alpha = 0.5) +
  #geom_sf(data = ecobio,  aes(fill = BIOME), alpha = 1, col = NA) +
  geom_point(data = nsites, aes(x = lon, y = lat, size = nsites), col = c("#0F0F0F"), alpha = 0.5) +
  #geom_polygon(data = wwfbiomes, x = )
  theme(panel.grid.major = element_line(colour = "transparent"), 
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size = 16), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("") +
  scale_size_continuous(range = c(0.2, 5), breaks = c(1, 10, 50, 100, 200)) +
  guides(fill=guide_legend(nrow=7,byrow=TRUE), size = guide_legend(nrow = 7, byrow = T)) 


ggsave(filename = paste0(outdir, "/MAP_Predicts_points_biome_POLS.pdf"),
       plot = last_plot(),
       width = 8,
       height = 6)





# extract the PREDICTS points
plot_data <- final.data.trans.pc[, c("SS", "SSBS", "Longitude", "Latitude")]

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



# plot of predicts sites across biomes including size per n sites
ggplot()+
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= c("#7CCD7C"), colour= "transparent", size=0.2, alpha = 0.5) +
  #geom_sf(data = ecobio,  aes(fill = BIOME), alpha = 1, col = NA) +
  geom_point(data = nsites, aes(x = lon, y = lat, size = nsites), col = c("#0F0F0F"), alpha = 0.5) +
  #geom_polygon(data = wwfbiomes, x = )
  theme(panel.grid.major = element_line(colour = "transparent"), 
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size = 16), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("") +
  scale_size_continuous(range = c(0.2, 5), breaks = c(1, 10, 50, 100, 200)) +
  guides(fill=guide_legend(nrow=7,byrow=TRUE), size = guide_legend(nrow = 7, byrow = T)) 


ggsave(filename = paste0(outdir, "/MAP_Predicts_points_biome_PCS.pdf"),
       plot = last_plot(),
       width = 8,
       height = 6)

