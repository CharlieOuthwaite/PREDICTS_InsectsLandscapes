##%######################################################%##
#                                                          #
####            Plotting modelled responses             ####
#                                                          #
##%######################################################%##

# plot results from each model and prepare figures for manuscript

rm(list = ls())

# load libraries
library(StatisticalModels)
library(ggplot2)
library(cowplot)
library(raster)

# set directories
datadir <- "2_MODEL_SELECTION"
datadir2 <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "3_FIGURES"
if(!dir.exists(outdir)) dir.create(outdir)

#source('functions/sort_data.r')
source('functions/rescale.r')
source('functions/unscale.r')



##%######################################################%##
#                                                          #
####                1. All insect plots                 ####
#                                                          #
##%######################################################%##


# load the models
load(paste0(datadir, "/RICHNESS_Model_Set_INSECTS.rdata")) # sr1
load(paste0(datadir, "/ABUNDANCE_Model_Set_INSECTS.rdata")) # sr1

# load the dataset
load(paste0(datadir2, "/PREDICTS_dataset_TRANS_INSECTS.rdata")) # final.data.trans

# load the scalers
scalers_ins <- read.csv(paste0(datadir2, "/Scaling_values_INSECTS.csv"))


##%######################################################%##
#                                                          #
####     1. Compare global and predicts data ranges     ####
#                                                          #
##%######################################################%##

range(final.data.trans$ncrop) # 0.0000 123
range(final.data.trans$pest_L) # 0.0000 139.7747
range(final.data.trans$pest_H) # 0.0000 166.4538
range(final.data.trans$percNH) # 0.992 99.988
table(final.data.trans$fields) # no fields to large fields (no very large fields)

# No fields Very small      Small     Medium      Large 
#      2645        298        663       1085        296 

# read in global maps
ncrops <- raster("Data/Earthstat_NCrops_newmethod.tif")
percNH <- raster("Data/PercentNatural.tif")
percNH <- percNH/10
pest_H <- raster("Data/Pesticide_totalAPR_High.tif")
pest_L <- raster("Data/Pesticide_totalAPR_Low.tif")
fields <- raster("Data/Global Field Sizes/dominant_field_size_categories.tif")

# take a look at range of values
ncrops # 1 127
percNH # 0 100
pest_H # 0, 328.0646
pest_L # 0, 231.9041
fields # very small, small, medium, large, very large




































##### ncrops model SR results #####


#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans$ncrop),
                       percNH = median(final.data.trans$percNH),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod_IN$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_IN$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_IN$data$Forest_biome) 
levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod_IN$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



## SR

# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),40),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))



# list to save plots in
p <- list()

p[[1]] <- recordPlot()


### percNH ###

# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans$ncrop),
                       percNH = median(final.data.trans$percNH),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)

levels(pred_tab$Predominant_land_use) <- levels(srmod_IN$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_IN$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_IN$data$Forest_biome) 
levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# save this basic table for later use
base_pred_tab <- pred_tab

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_IN$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))

# SR plot = full range
p[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,25)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("All insects")

p[[2]] <- ggplotGrob(p[[2]])



### ncrop ###

pred_tab <- base_pred_tab

from = 0
to = 123
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_IN$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans$ncrop, scale = scalers_ins[6, 2], centre = scalers_ins[6, 3], log = F))

# SR plot = full range
p[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD6600")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD6600"), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 130)) +
  xlab("Number of crops") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("All insects")

p[[3]] <- ggplotGrob(p[[3]])




### Land use : percNH ###

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_IN$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))


p[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("All insects")

p[[4]] <- ggplotGrob(p[[4]])



### tropical : percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Tropical'
n <- 2

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans[, fac])[2]


# predict the result
result <- PredictGLMER(model = srmod_IN$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))

percNH <- as.data.frame(final.data.trans[, c("percNH", "Tropical")])
percNH$percNH <- unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F)

p[[5]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = percNH, col = Tropical), size = 0.1) +
  ylim(c(0,65)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#009ACD", "#9932CC"))+
  scale_fill_manual(values = c("#009ACD", "#9932CC")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #legend.position = 'none', 
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

p[[5]] <- ggplotGrob(p[[5]])


# organise plots
plots <- plot_grid(plotlist = p, ncol = 2)

save_plot(filename = paste0(outdir, "/INSECTS_all_plots_ncrop_SR.pdf"), plot = plots, base_height = 12, base_width = 8)


#####################################################################################################################################

### pest_H model SR results ###

#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans$ncrop),
                       percNH = median(final.data.trans$percNH),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod_IN$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_IN$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_IN$data$Forest_biome) 
levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod_IN2$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),40),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))




# list to save plots in
q <- list()

q[[1]] <- recordPlot()


### percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_IN2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))

# SR plot = full range
q[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,25)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("All insects")

q[[2]] <- ggplotGrob(q[[2]])



### Land use : percNH ###

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_IN2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))


q[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("All insects")

q[[3]] <- ggplotGrob(q[[3]])



### tropical : percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Tropical'
n <- 2

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans[, fac])[2]


# predict the result
result <- PredictGLMER(model = srmod_IN2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))

percNH <- as.data.frame(final.data.trans[, c("percNH", "Tropical")])
percNH$percNH <- unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F)

q[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = percNH, col = Tropical), size = 0.1) +
  ylim(c(0,65)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#009ACD", "#9932CC"))+
  scale_fill_manual(values = c("#009ACD", "#9932CC")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #legend.position = 'none', 
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

q[[4]] <- ggplotGrob(q[[4]])


# organise plots
plots <- plot_grid(plotlist = q, ncol = 2)

save_plot(filename = paste0(outdir, "/INSECTS_all_plots_pestH_SR.pdf"), plot = plots, base_height = 8, base_width = 8)


################################################################################################################################

#### All insects, Abundance model results ####

# load the models
load(paste0(datadir, "/ABMOD_output_INSECTS_ncrop.rdata")) # abmod_IN

#### land use - use intensity interactions ####

pred_tab <- do.call("rbind", replicate(9, base_pred_tab, simplify = FALSE))

# add a column for pest_H
pred_tab$pest_H_log <- median(final.data.trans$pest_H_log)

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = abmod_IN$model, data = pred_tab)

# transform the results
result <- exp(result)-1

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


### plot ###
# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),80),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Abundance (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))




# list to save plots in
r <- list()

r[[1]] <- recordPlot()



#### pest_H ####

from = 0
to = 170
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'pest_H'
logval = T
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)

# add in a pesticide value to the base table
base_pred_tab$pest_H_log <- median(final.data.trans$pest_H_log)

# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,'pest_H_log'] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_IN$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##
#pest_H_log <- as.data.frame(unscale(final.data.trans$pest_H_log, scale = scalers_ins[7, 2], centre = scalers_ins[7, 3], log = T))

# SR plot = full range
r[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = pest_H), size = 0.1) +
  ylim(c(0,150)) +
  xlim(c(0, 200)) +
  xlab("Pesticide application") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("All insects")

r[[2]] <- ggplotGrob(r[[2]])




### Use intensity : pest_H ###

from = 0
to = 170
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'pest_H'
logval = T
fac <- 'Use_intensity'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_ins[scalers_ins$variable == variable, 'centre'], 
                     scale = scalers_ins[scalers_ins$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, 'pest_H_log'] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_IN$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
#percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))


r[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = pest_H), size = 0.1) +
  ylim(c(0,200)) +
  xlim(c(0, 200)) +
  xlab("Pesticide application") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#6E8B3D", "#FFD700", "#CD3700"))+
  scale_fill_manual(values = c("#6E8B3D", "#FFD700", "#CD3700")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("All insects")

r[[3]] <- ggplotGrob(r[[3]])


# organise plots
plots <- plot_grid(plotlist = r, ncol = 2)

save_plot(filename = paste0(outdir, "/INSECTS_all_plots_AB.pdf"), plot = plots, base_height = 12, base_width = 8)






##%######################################################%##
#                                                          #
####                2. Pollinators plots                ####
#                                                          #
##%######################################################%##




# load the models
load(paste0(datadir, "/SRMOD_output_POLLINATORS.rdata")) # srmod_PO2
load(paste0(datadir, "/SRMOD_output_POLLINATORS_ncrop.rdata")) # srmod_PO

# load the dataset
load(paste0(datadir2, "/PREDICTS_dataset_TRANS_POLS.rdata"))

# load the scalers
scalers_pols <- read.csv(paste0(datadir2, "/Scaling_values_POLS.csv"))


##### ncrops model SR results #####


#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pols$ncrop),
                       percNH = median(final.data.trans.pols$percNH),
                       #Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod_PO$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PO$data$Use_intensity) 
#levels(pred_tab$Forest_biome) <- levels(srmod_IN$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod_PO$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



## SR

# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),80),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))



# list to save plots in
p <- list()

p[[1]] <- recordPlot()


### percNH ###

# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pols$ncrop),
                       percNH = median(final.data.trans.pols$percNH),
                       #Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)

levels(pred_tab$Predominant_land_use) <- levels(srmod_PO$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PO$data$Use_intensity) 
#levels(pred_tab$Forest_biome) <- levels(srmod_PO$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_PO$data$Tropical) 

# save this basic table for later use
base_pred_tab <- pred_tab

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_pols[scalers_pols$variable == variable, 'centre'], 
                     scale = scalers_pols[scalers_pols$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PO$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pols$percNH, scale = scalers_pols[5, 2], centre = scalers_pols[5, 3], log = F))

# SR plot = full range
p[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,15)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pollinators")

p[[2]] <- ggplotGrob(p[[2]])



### ncrop ###

pred_tab <- base_pred_tab

from = 0
to = 123
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_pols[scalers_pols$variable == variable, 'centre'], 
                     scale = scalers_pols[scalers_pols$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PO$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans.pols$ncrop, scale = scalers_pols[6, 2], centre = scalers_pols[6, 3], log = F))

# SR plot = full range
p[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD6600")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD6600"), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,15)) +
  xlim(c(0, 130)) +
  xlab("Number of crops") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pollinators")

p[[3]] <- ggplotGrob(p[[3]])


# organise plots
plots <- plot_grid(plotlist = p, ncol = 2)

save_plot(filename = paste0(outdir, "/POLLINATORS_all_plots_ncrop_SR.pdf"), plot = plots, base_height = 8, base_width = 8)



#####################################################################################################################################

### pest_H model SR results ###

#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(
  #ncrop = median(final.data.trans.pols$ncrop),
  percNH = median(final.data.trans.pols$percNH),
  #Forest_biome = "Temperate Broadleaf & Mixed Forests",
  Use_intensity = "Minimal use",
  Predominant_land_use = "Primary vegetation",
  #Tropical = "Temperate",
  Species_richness = 0,
  logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod_PO2$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PO2$data$Use_intensity) 
#levels(pred_tab$Forest_biome) <- levels(srmod_IN$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod_PO2$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),80),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))




# list to save plots in
q <- list()

q[[1]] <- recordPlot()


### percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_pols[scalers_pols$variable == variable, 'centre'], 
                     scale = scalers_pols[scalers_pols$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PO2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pols$percNH, scale = scalers_ins[5, 2], centre = scalers_ins[5, 3], log = F))

# SR plot = full range
q[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,15)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pollinators")

q[[2]] <- ggplotGrob(q[[2]])


# organise plots
plots <- plot_grid(plotlist = q, ncol = 2)

save_plot(filename = paste0(outdir, "/POLLINATORS_all_plots_pestH_SR.pdf"), plot = plots, base_height = 4, base_width = 8)



################################################################################################################################

#### Pollinators, Abundance model results ####

# load the models
load(paste0(datadir, "/ABMOD_output_POLLINATORS.rdata")) # abmod_PO

#### land use - use intensity interactions ####

pred_tab <- do.call("rbind", replicate(9, base_pred_tab, simplify = FALSE))

# add a column for pest_H
pred_tab$pest_H_log <- median(final.data.trans.pols$pest_H_log)

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = abmod_PO$model, data = pred_tab)

# transform the results
result <- exp(result)-1

result <- sweep(x = result, MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


pdf(file = paste0(outdir, "/POLLINATORS_all_plots_AB.pdf"))

### plot ###
# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),150),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Abundance (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))


dev.off()




##%######################################################%##
#                                                          #
####              3. Pest controller plots              ####
#                                                          #
##%######################################################%##



# load the models
load(paste0(datadir, "/SRMOD_output_PESTC.rdata")) # srmod_PC2
load(paste0(datadir, "/SRMOD_output_PESTC_ncrop.rdata")) # srmod_PC

# load the dataset
load(paste0(datadir2, "/PREDICTS_dataset_TRANS_PC.rdata"))

# load the scalers
scalers_PC <- read.csv(paste0(datadir2, "/Scaling_values_PC.csv"))


##### ncrops model SR results #####


#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod_PC$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PC$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_PC$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod_PC$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



## SR

# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),40),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))

title("interaction not selected though")

# list to save plots in
p <- list()

p[[1]] <- recordPlot()


### percNH ###

# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)

levels(pred_tab$Predominant_land_use) <- levels(srmod_PC$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PC$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_PC$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# save this basic table for later use
base_pred_tab <- pred_tab

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))

# SR plot = full range
p[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

p[[2]] <- ggplotGrob(p[[2]])



### ncrop ###

pred_tab <- base_pred_tab

from = 0
to = 123
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans.pc$ncrop, scale = scalers_PC[6, 2], centre = scalers_PC[6, 3], log = F))

# SR plot = full range
p[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD6600")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD6600"), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 130)) +
  xlab("Number of crops") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

p[[3]] <- ggplotGrob(p[[3]])



### nlandcovers ###

pred_tab <- base_pred_tab

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))

# SR plot = full range
p[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 10)) +
  xlab("Number of landcovers") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

p[[4]] <- ggplotGrob(p[[4]])



### Land use : ncrop ###

from = 0
to = 120
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans.pc$ncrop, scale = scalers_PC[6, 2], centre = scalers_PC[6, 3], log = F))


p[[5]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,25)) +
  xlim(c(0, 120)) +
  xlab("Number of crops") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

p[[5]] <- ggplotGrob(p[[5]])




### Land use : landcovers ###

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))


p[[6]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0, 20)) +
  xlim(c(0, 10)) +
  xlab("Number of landcovers") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

p[[6]] <- ggplotGrob(p[[6]])



### Use intensity : ncrop ###

from = 0
to = 120
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- 'Use_intensity'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, variable ] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans.pc$ncrop, scale = scalers_PC[6, 2], centre = scalers_PC[6, 3], log = F))


p[[7]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 120)) +
  xlab("Number of crops") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#6E8B3D", "#FFD700", "#CD3700"))+
  scale_fill_manual(values = c("#6E8B3D", "#FFD700", "#CD3700")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

p[[7]] <- ggplotGrob(p[[7]])



### Use intensity : percNH ###

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Use_intensity'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, variable ] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))


p[[8]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#6E8B3D", "#FFD700", "#CD3700"))+
  scale_fill_manual(values = c("#6E8B3D", "#FFD700", "#CD3700")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

p[[8]] <- ggplotGrob(p[[8]])


# organise plots
plots <- plot_grid(plotlist = p, ncol = 3, nrow = 3)

save_plot(filename = paste0(outdir, "/PESTCS_all_plots_SR_ncrop.pdf"), plot = plots, base_height = 12, base_width = 12)




##%######################################################%##
#                                                          #
####              3. Pest controller plots              ####
#                                                          #
##%######################################################%##



# load the models
load(paste0(datadir, "/SRMOD_output_PESTC.rdata")) # srmod_PC22
load(paste0(datadir, "/SRMOD_output_PESTC_ncrop.rdata")) # srmod_PC2

# load the dataset
load(paste0(datadir2, "/PREDICTS_dataset_TRANS_PC.rdata"))

# load the scalers
scalers_PC <- read.csv(paste0(datadir2, "/Scaling_values_PC.csv"))


##### pest_H model SR results #####


#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod_PC2$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PC2$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_PC2$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod_PC2$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



## SR

# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),40),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))

title("interaction not selected though")

# list to save plots in
q <- list()

q[[1]] <- recordPlot()


### percNH ###

# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)

levels(pred_tab$Predominant_land_use) <- levels(srmod_PC2$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_PC2$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_PC2$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# save this basic table for later use
base_pred_tab <- pred_tab

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))

# SR plot = full range
q[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

q[[2]] <- ggplotGrob(q[[2]])



### pest_H ###

pred_tab <- base_pred_tab

from = 0
to = 170
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'pest_H'
logval = T
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,'pest_H_log'] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# SR plot = full range
q[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = final.data.trans.pc, aes(x = pest_H), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 170)) +
  xlab("Pesticide application") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

q[[3]] <- ggplotGrob(q[[3]])



### nlandcovers ###

pred_tab <- base_pred_tab

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = srmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))

# SR plot = full range
q[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 10)) +
  xlab("Number of landcovers") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

q[[4]] <- ggplotGrob(q[[4]])



### Land use : pest_H ###

from = 0
to = 170
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'pest_H'
logval = T
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##


q[[5]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans.pc, aes(x = pest_H), size = 0.1) +
  ylim(c(0,25)) +
  xlim(c(0, 170)) +
  xlab("Pesticide application") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

q[[5]] <- ggplotGrob(q[[5]])





### Use intensity : percNH ###

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Use_intensity'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, variable ] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = srmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))


q[[6]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,20)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#6E8B3D", "#FFD700", "#CD3700"))+
  scale_fill_manual(values = c("#6E8B3D", "#FFD700", "#CD3700")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

q[[6]] <- ggplotGrob(q[[6]])


# organise plots
plots <- plot_grid(plotlist = q, ncol = 2, nrow = 3)

save_plot(filename = paste0(outdir, "/PESTCS_all_plots_SR_pestH.pdf"), plot = plots, base_height = 12, base_width = 8)



###############################################################################################################

#### Abundance models, ncrop ####


# load the models
load(paste0(datadir, "/ABMOD_output_PESTC.rdata")) # abmod_PC22
load(paste0(datadir, "/ABMOD_output_PESTC_ncrop.rdata")) # abmod_PC2



##%######################################################%##
#                                                          #
####              3. Pest controller plots              ####
#                                                          #
##%######################################################%##



# load the models
load(paste0(datadir, "/ABMOD_output_PESTC.rdata")) # abmod_PC2
load(paste0(datadir, "/ABMOD_output_PESTC_ncrop.rdata")) # abmod_PC

# load the dataset
load(paste0(datadir2, "/PREDICTS_dataset_TRANS_PC.rdata"))

# load the scalers
scalers_PC <- read.csv(paste0(datadir2, "/Scaling_values_PC.csv"))




#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(abmod_PC$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_PC$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_PC$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### ab predictions ####


# predict the result
result <- PredictGLMERRandIter(model = abmod_PC$model, data = pred_tab)

# transform the results
result <- exp(result)-1

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100




# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),250),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Total Abundance (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))


# list to save plots in
r <- list()

r[[1]] <- recordPlot()


### percNH ###

# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)

levels(pred_tab$Predominant_land_use) <- levels(abmod_PC$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_PC$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_PC$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(srmod_IN$data$Tropical) 

# save this basic table for later use
base_pred_tab <- pred_tab

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))

# SR plot = full range
r[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,100)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

r[[2]] <- ggplotGrob(r[[2]])



### ncrop ###

pred_tab <- base_pred_tab

from = 0
to = 120
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##

# SR plot = full range
r[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD6600")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD6600"), alpha = 0.3) +
  geom_rug(data = final.data.trans.pc, aes(x = pest_H), size = 0.1) +
  ylim(c(0,100)) +
  xlim(c(0, 170)) +
  xlab("Number of crops") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

r[[3]] <- ggplotGrob(r[[3]])



### nlandcovers ###

pred_tab <- base_pred_tab

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))

# SR plot = full range
r[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1), size = 0.1) +
  ylim(c(0,100)) +
  xlim(c(0, 10)) +
  xlab("Number of landcovers") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

r[[4]] <- ggplotGrob(r[[4]])



### Land use : ncrop ###

from = 0
to = 120
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans.pc$ncrop, scale = scalers_PC[6, 2], centre = scalers_PC[6, 3], log = F))


r[[5]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,170)) +
  xlim(c(0, 120)) +
  xlab("Number of crops") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

r[[5]] <- ggplotGrob(r[[5]])



### Land use : landcovers ###

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))


r[[6]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = landcovers.5k, aes(x = V1), size = 0.1) +
  ylim(c(0, 200)) +
  xlim(c(0, 10)) +
  xlab("Number of landcovers") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

r[[6]] <- ggplotGrob(r[[6]])

### Land use : percNH ###

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))


r[[7]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0, 200)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

r[[7]] <- ggplotGrob(r[[7]])

### Use intensity : ncrop ###

from = 0
to = 120
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'ncrop'
logval = F
fac <- 'Use_intensity'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, variable ] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
ncrop <- as.data.frame(unscale(final.data.trans.pc$ncrop, scale = scalers_PC[6, 2], centre = scalers_PC[6, 3], log = F))


r[[8]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = ncrop, aes(x = V1), size = 0.1) +
  ylim(c(0,120)) +
  xlim(c(0, 120)) +
  xlab("Number of crops") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#6E8B3D", "#FFD700", "#CD3700"))+
  scale_fill_manual(values = c("#6E8B3D", "#FFD700", "#CD3700")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

r[[8]] <- ggplotGrob(r[[8]])



# organise plots
plots <- plot_grid(plotlist = r, ncol = 3, nrow = 3)

save_plot(filename = paste0(outdir, "/PESTCS_all_plots_AB_ncrop.pdf"), plot = plots, base_height = 12, base_width = 12)


##########################################################################################################################

#### Abundance model, pest_H



#### land use - use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(abmod_PC2$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_PC2$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_PC2$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(abmod_IN$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))

pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = abmod_PC2$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



## SR

# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),200),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Total Abundance (%)",xlab="",bty="l")


axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))


# list to save plots in
q <- list()

q[[1]] <- recordPlot()


### percNH ###

# basic table of median values and reference factors
pred_tab <- data.frame(ncrop = median(final.data.trans.pc$ncrop),
                       percNH = median(final.data.trans.pc$percNH),
                       pest_H_log = median(final.data.trans.pc$pest_H_log),
                       landcovers.5k = median(final.data.trans.pc$landcovers.5k),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)

levels(pred_tab$Predominant_land_use) <- levels(abmod_PC2$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_PC2$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_PC2$data$Forest_biome) 
#levels(pred_tab$Tropical) <- levels(abmod_IN$data$Tropical) 

# save this basic table for later use
base_pred_tab <- pred_tab

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))

# SR plot = full range
q[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,100)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

q[[2]] <- ggplotGrob(q[[2]])



### pest_H ###

pred_tab <- base_pred_tab

from = 0
to = 170
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'pest_H'
logval = T
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,'pest_H_log'] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##

# SR plot = full range
q[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = final.data.trans.pc, aes(x = pest_H), size = 0.1) +
  ylim(c(0,100)) +
  xlim(c(0, 170)) +
  xlab("Pesticide application") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

q[[3]] <- ggplotGrob(p[[3]])



### nlandcovers ###

pred_tab <- base_pred_tab

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- NULL
n <- NULL

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = abmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))

# SR plot = full range
q[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1), size = 0.1) +
  ylim(c(0,100)) +
  xlim(c(0, 10)) +
  xlab("Number of landcovers") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Pest Controllers")

q[[4]] <- ggplotGrob(q[[4]])



### UI : pest_H ###

from = 0
to = 170
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'pest_H'
logval = T
fac <- 'Use_intensity'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[,variable] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##


q[[5]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans.pc, aes(x = pest_H), size = 0.1) +
  ylim(c(0,250)) +
  xlim(c(0, 170)) +
  xlab("Pesticide application") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#6E8B3D", "#FFD700", "#CD3700"))+
  scale_fill_manual(values = c("#6E8B3D", "#FFD700", "#CD3700")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

q[[5]] <- ggplotGrob(q[[5]])





### LU : percNH ###

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, variable ] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans.pc$percNH, scale = scalers_PC[5, 2], centre = scalers_PC[5, 3], log = F))


q[[6]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,250)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

q[[6]] <- ggplotGrob(q[[6]])



### LU : landcovers ###

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
logval = F
fac <- 'Predominant_land_use'
n <- 3

# transform the values
vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers_PC[scalers_PC$variable == variable, 'centre'], 
                     scale = scalers_PC[scalers_PC$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab <- do.call("rbind", replicate(length(vals), base_pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab[, variable ] <- vals_trans


# add in sections to pred_data for land uses
pred_tab <- do.call("rbind", replicate(n, pred_tab, simplify = FALSE))
pred_tab[, fac][(nrow(pred_tab)/n + 1):(nrow(pred_tab)/n + 1000)] <- levels(final.data.trans.pc[, fac])[2]
pred_tab[, fac][(nrow(pred_tab)/n + 1001):(nrow(pred_tab)/n + 2000)] <- levels(final.data.trans.pc[, fac])[3]


# predict the result
result <- PredictGLMER(model = abmod_PC2$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


## organise data for plotting ##
landcovers.5k <- as.data.frame(unscale(final.data.trans.pc$landcovers.5k, scale = scalers_PC[3, 2], centre = scalers_PC[3, 3], log = F))


q[[7]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,250)) +
  xlim(c(0, 10)) +
  xlab("Number of landcover") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("Pest Controllers")

q[[7]] <- ggplotGrob(q[[7]])

# organise plots
plots <- plot_grid(plotlist = q, ncol = 3, nrow = 3)

save_plot(filename = paste0(outdir, "/PESTCS_all_plots_AB_pestH.pdf"), plot = plots, base_height = 12, base_width = 8)


