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








# From Kyra:

# Here is the function I used to split the PREDICTS data set by insect order



# establish the naming convention for the data frames it will produce
# I wanted the data frames to be named by insect order and nothing else, so I left the first part empty ("")
# if I had wanted the names to be"Insect_Order", I would have put # OrderName <- paste0("Insect_",predicts$Order)
OrderName <- paste0("",predicts$Order)

# create a list of data frames
by_Order <- split(predicts,OrderName)

# check
by_Order

# extract data frames from list into global environment to work with them individually
list2env(by_Order,globalenv())

# putting it all back together
sites <- rbind(Blattodea,Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Neuroptera,Orthoptera,Other,Thysanoptera,Trichoptera)