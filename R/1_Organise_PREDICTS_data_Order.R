##%######################################################%##
#                                                          #
####             1. Organise PREDICTS data              ####
#                                                          #
##%######################################################%##

# This script organises the PREDICTS data for all insects, but including 
# insect order so that the models can be run with Order as a factor level. 


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

pred.data <- droplevels(pred.data)
## split the data by order

length(unique(pred.data$Order)) # 25 orders

# remove rows where order is not known
pred.data <- pred.data[!pred.data$Order == "", ] # 934446 rows

# table(pred.data$Order)
# some orders will need to be combined where few records. 
# look at numbers of sites after data aggregated then come back to combine some groups. 


pred.data$Order2 <-as.character(pred.data$Order)

#Dictyoptera is a superorder including Blaattodea and Mantodea
pred.data$Order2[pred.data$Order2 == "Blattodea" | pred.data$Order2 == "Mantodea"] <- "Dictyoptera"

#Orthopterida is a superorder including orthoptera and phasmida
pred.data$Order2[pred.data$Order2 == "Orthoptera" | pred.data$Order2 == "Phasmida"] <- "Orthopterida"

#Psocodea includes the suborder Phthiraptera 
pred.data$Order2[pred.data$Order2 == "Phthiraptera"] <- "Psocodea"

# grouping others that don't fit elsewhere and groups with less than 250 sites into "Other"
pred.data$Order2[pred.data$Order2 == "Siphonaptera" | pred.data$Order2 == "Zoraptera" | pred.data$Order2 == "Zygentoma"] <- "Other"
pred.data$Order2[pred.data$Order2 == "Embioptera" | pred.data$Order2 == "Mecoptera"] <- "Other"

# trying out grouping things with a freshwater life stage
pred.data$Order2[pred.data$Order2 == "Trichoptera" | pred.data$Order2 == "Ephemeroptera" | pred.data$Order2 == "Odonata"] <- "Freshwater"

# using code from Kyra:
# Here is the function I used to split the PREDICTS data set by insect order
# establish the naming convention for the data frames it will produce
# I wanted the data frames to be named by insect order and nothing else, so I left the first part empty ("")
# if I had wanted the names to be"Insect_Order", I would have put # OrderName <- paste0("Insect_",predicts$Order)
OrderName <- paste0("",pred.data$Order2) # 19 orders including "Other"

# create a list of data frames
by_Order <- split(pred.data,OrderName)

# check
by_Order

# extract data frames from list into global environment to work with them individually
#list2env(by_Order,globalenv())

# putting it all back together
#sites <- rbind(Blattodea,Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Neuroptera,Orthoptera,Other,Thysanoptera,Trichoptera)

# apply sampling effort correction to each subset
by_Order_CSE <- lapply(by_Order, CorrectSamplingEffort)

# apply merge sites to each subset
by_Order_MS <- lapply(by_Order_CSE, MergeSites, merge.extra = c("Order2"))

# apply site metrics to each subset
by_Order_SM <- lapply(by_Order_MS, SiteMetrics, extra.cols = c("Predominant_land_use", "SSB", "SSBS", "Order2"))

# bring back together into one dataset
sites.sub <- do.call("rbind", by_Order_SM)

table(sites.sub$Order2)

# # number of sites for each order
# # those with few sites will need to be combined somehow. 
# Archaeognatha    Coleoptera    Dermaptera   Dictyoptera       Diptera    Embioptera Ephemeroptera     Hemiptera   Hymenoptera 
# 323          2638           363           437           877            97            85          1006          4379 
# Isoptera   Lepidoptera     Mecoptera    Neuroptera       Odonata  Orthopterida         Other      Psocodea  Thysanoptera 
# 281          2079            77           538           103           664           145           488           556 
# Trichoptera 
# 526 

# Archaeognatha    Coleoptera    Dermaptera   Dictyoptera       Diptera    Freshwater     Hemiptera   Hymenoptera      Isoptera 
# 323          2638           363           437           877           554          1006          4379           281 
# Lepidoptera    Neuroptera  Orthopterida         Other      Psocodea  Thysanoptera 
# 2079           538           664           204           488           556 




# which land uses are we interested in, I have kept plantation in here too
sites.sub <- sites.sub[!sites.sub$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide"), ]

# 11294 rows (different now to number of sites as smae sites will have summaries of dif orders)


# remove sites with NA in lat/long columns
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ] # 11288 rows


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

# not sure if using this at the moment but adding in

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
# 6565      3997 



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
nrow(sites.sub) # 10562



##%######################################################%##
#                                                          #
####      checking correlations between variables       ####
#                                                          #
##%######################################################%##



# # function to create plot
# panel.cor <- function(x, y, ...)
# {
#   par(usr = c(0, 1, 0, 1))
#   txt <- as.character(format(cor(x, y), digits=2))
#   text(0.5, 0.5, txt, cex = 2)
# }
# 
# # save figure
# pdf(file = paste0(outdir, "/Correlations_all_variables.pdf"), width =14, height = 9)
# 
# # correlations, including all nlandcovers buffers
# pairs(sites.sub[ , c(23:27)], 
#       upper.panel=panel.cor, 
#       diag.panel = panel.hist, 
#       main = "",
#       cex = 2)
# 
# dev.off()
# 
# # only high ones are the high and low pesticide (as expected) and NH with field size (-0.42).
# 

sites.sub <- droplevels(sites.sub)

# save the untransformed datasets including all variables
save(sites.sub, file = paste0(outdir, "/PREDICTS_dataset_inc_variables_ORDERS.rdata"))
#load(file = paste0(outdir, "/PREDICTS_dataset_inc_variables_ORDERS.rdata"))


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
write.csv(values, paste0(outdir, "/Scaling_values_ORDERS.csv"), row.names = FALSE)


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


# combine plantation and cropland sites
final.data.trans$LU <- final.data.trans$Predominant_land_use

# set cropland and plantation to "Agriculture"
final.data.trans$LU[final.data.trans$LU %in% c("Cropland", "Plantation")] <- "Agriculture"

# set factor levels of predominant land use
final.data.trans$LU <- factor(final.data.trans$LU,
                              levels=c("Primary vegetation","Secondary vegetation", "Agriculture"))

# nsites per land use
table(final.data.trans$LU)

# Primary vegetation Secondary vegetation          Agriculture 
#               3266                 3676                 2711


# pdf(file = paste0(outdir, "/Correlations_final_variables.pdf"), width =9, height = 9)
# # correlations of final set of variables - transformed
# pairs(final.data.trans[ , c(26:31)], 
#       upper.panel=panel.cor, 
#       diag.panel = panel.hist, 
#       main = "",
#       cex.labels = 1.3)
# 
# dev.off()


# edit field size names
final.data.trans$fields <- as.character(final.data.trans$fields)
final.data.trans$fields[final.data.trans$fields == "0"] <- "No fields"
final.data.trans$fields[final.data.trans$fields == "3502"] <- "Very large"
final.data.trans$fields[final.data.trans$fields == "3503"] <- "Large"
final.data.trans$fields[final.data.trans$fields == "3504"] <- "Medium"
final.data.trans$fields[final.data.trans$fields == "3505"] <- "Small"
final.data.trans$fields[final.data.trans$fields == "3506"] <- "Very small"

table(final.data.trans$fields)

final.data.trans$fields <- factor(final.data.trans$fields, levels = c("No fields", "Very small", "Small", "Medium", "Large"))



# transform abundance
final.data.trans$logAbun <- log(final.data.trans$Total_abundance + 1)


# set Order2 as a factor level
final.data.trans$Order2 <- as.factor(final.data.trans$Order2)



# save transformed dataset
save(final.data.trans, file = paste0(outdir, "/PREDICTS_dataset_TRANS_ORDERS.rdata"))




##%######################################################%##
#                                                          #
####                 Get data summaries                 ####
#                                                          #
##%######################################################%##

# get the number of studies, sites, etc for the paper

length(unique(final.data.trans$SS)) # 248 studies
length(unique(final.data.trans$SSBS)) # 4987 sites

table(final.data.trans$Order2, final.data.trans$LU)


# 
#               Primary vegetation Secondary vegetation Agriculture
# Archaeognatha                109                  120           0
# Coleoptera                   898                  502         275
# Dermaptera                   110                   80          57
# Dictyoptera                  124                  133          52
# Diptera                       95                  250         204
# Freshwater                   135                  152          69
# Hemiptera                    160                  211         117
# Hymenoptera                  564                  741        1411
# Isoptera                      45                  106          53
# Lepidoptera                  424                  685         183
# Neuroptera                   141                  145          67
# Orthopterida                 129                  168         120
# Other                         43                   89          33
# Psocodea                     144                  147           2
# Thysanoptera                 145                  147          68
