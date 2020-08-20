##%######################################################%##
#                                                          #
####                 2. Model selection                 ####
#                                                          #
##%######################################################%##

# This script carries out the model selection process using the dataset generated
# in the previous script. Models are run for 3 biodiversity metrics. 

rm(list = ls())


# load libraries
library(devtools)
install_github(repo = "timnewbold/StatisticalModels")
library(StatisticalModels)
library(roquefort)
library(cowplot)
library(gridGraphics)


# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "2_MODEL_SELECTION"
dir.create(outdir)

# read in the PREDICTS datasets with landscape variables

load(paste0(datadir, "/PREDICTS_dataset_TRANS_INSECTS.rdata")) # final.data.trans
load(paste0(datadir, "/PREDICTS_dataset_TRANS_POLS.rdata")) # final.data.trans.pols 
load(paste0(datadir, "/PREDICTS_dataset_TRANS_PC.rdata")) # final.data.trans.pc



#### Run model selection process including landscape variables that are not correlated ####


##%######################################################%##
#                                                          #
####                1. Species richness                 ####
#                                                          #
##%######################################################%##


#### 1. Species richness models ####



### A. ALL INSECTS ###

# run the model selection process

system.time({sr1 <- GLMERSelect(modelData = final.data.trans, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(pest_H_log = 1, landcovers.5k = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_log", 
                                                      #"Predominant_land_use:ncrop",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:pest_H_log", 
                                                      #"Use_intensity:ncrop",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:percNH"), verbose = F)}) 
# save the output
save(sr1, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection_INSECTS.rdata"))

# take a look at the model output
summary(sr1$model)

# extract the stats produced as part of the model selection process
sr1stats <- as.data.frame(sr1$stats)

# save these
write.csv(sr1stats, file = paste0(outdir, "/sr_stats_INSECTS.csv"), row.names = F)



### B. POLLINATORS ###

system.time({sr2 <- GLMERSelect(modelData = final.data.trans.pols, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(pest_H_log = 1, landcovers.5k = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_log", 
                                                      #"Predominant_land_use:ncrop",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:pest_H_log", 
                                                      #"Use_intensity:ncrop",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:percNH"), verbose = F)}) 
# save the output
save(sr2, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection_POLLINATORS.rdata"))

# take a look at the model output
summary(sr2$model)

# extract the stats produced as part of the model selection process
sr2stats <- as.data.frame(sr2$stats)

# save these
write.csv(sr2stats, file = paste0(outdir, "/sr_stats_POLLINATORS.csv"), row.names = F)





### C. PEST CONTROLLERS ###


system.time({sr3 <- GLMERSelect(modelData = final.data.trans.pc, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(pest_H_log = 1, landcovers.5k = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_log", 
                                                      #"Predominant_land_use:ncrop",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:pest_H_log", 
                                                      #"Use_intensity:ncrop",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:percNH"), verbose = F)}) 

save(sr3, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection_PESTC.rdata"))

# take a look at the model output
summary(sr3$model)

# extract the stats produced as part of the model selection process
sr3stats <- as.data.frame(sr3$stats)

# save these
write.csv(sr3stats, file = paste0(outdir, "/sr_stats_PESTC.csv"), row.names = F)




##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##


final.data.trans$logAbun <- log(final.data.trans$Total_abundance + 1)
final.data.trans.pols$logAbun <- log(final.data.trans.pols$Total_abundance + 1)
final.data.trans.pc$logAbun <- log(final.data.trans.pc$Total_abundance + 1)



### A. ALL INSECTS ###


system.time({ab1 <- GLMERSelect(modelData = final.data.trans, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(pest_H_log = 1, landcovers.5k = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_log", 
                                                      #"Predominant_land_use:ncrop",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:pest_H_log", 
                                                      #"Use_intensity:ncrop",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:percNH"), verbose = F)})

# take a look at the model output
summary(ab1$model)

# save the output
save(ab1, file = paste0(outdir, "/ABUNDANCE_Model_Selection_INSECTS.rdata"))

# extract the stats produced as part of the model selection process
ab1stats <- as.data.frame(ab1$stats)

# save these
write.csv(ab1stats, file = paste0(outdir, "/ab_stats_INSECTS.csv"), row.names = F)



### B. POLLINATORS ###


system.time({ab2 <- GLMERSelect(modelData = final.data.trans.pols, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(pest_H_log = 1, landcovers.5k = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_log", 
                                                      #"Predominant_land_use:ncrop",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:pest_H_log", 
                                                      #"Use_intensity:ncrop",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:percNH"), verbose = F)})

# take a look at the model output
summary(ab2$model)

# save the output
save(ab2, file = paste0(outdir, "/ABUNDANCE_Model_Selection_POLLINATORS.rdata"))

# extract the stats produced as part of the model selection process
ab2stats <- as.data.frame(ab2$stats)

# save these
write.csv(ab2stats, file = paste0(outdir, "/ab_stats_POLLINATORS.csv"), row.names = F)




### C. PEST CONTROLLERS ###


system.time({ab3 <- GLMERSelect(modelData = final.data.trans.pc, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(pest_H_log = 1, landcovers.5k = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_log", 
                                                      #"Predominant_land_use:ncrop",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:pest_H_log", 
                                                      #"Use_intensity:ncrop",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:percNH"), verbose = F)})

# take a look at the model output
summary(ab3$model)

# save the output
save(ab3, file = paste0(outdir, "/ABUNDANCE_Model_Selection_PESTC.rdata"))

# extract the stats produced as part of the model selection process
ab3stats <- as.data.frame(ab3$stats)

# save these
write.csv(ab3stats, file = paste0(outdir, "/ab_stats_PESTC.csv"), row.names = F)




##%######################################################%##
#                                                          #
####        Run final selected model using REML         ####
#                                                          #
##%######################################################%##



### selected model for species richness ###

# one model with ncrop, one without #

### A. ALL INSECTS ###

# rerun the selected models, using REML

srmod_IN <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + percNH + ncrop + Predominant_land_use:percNH + Predominant_land_use:Use_intensity + Tropical:percNH",
               randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)

srmod_IN2 <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
                   fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + percNH + Predominant_land_use:percNH + Predominant_land_use:Use_intensity + Tropical:percNH",
                   randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)



# Warning messages: none



# take a look at the model output
summary(srmod_IN$model)
summary(srmod_IN2$model)

# extract the coefficents of the model
coefs <- fixef(srmod_IN$model)
coefs2 <- fixef(srmod_IN2$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs_INSECTS_ncrops.csv"), row.names = F)
write.csv(coefs2, file = paste0(outdir, "/SPECIESRICHNESS_coefs_INSECTS.csv"), row.names = F)

# save the model output
save(srmod_IN, file = paste0(outdir, "/SRMOD_output_INSECTS_ncrops.rdata"))
save(srmod_IN2, file = paste0(outdir, "/SRMOD_output_INSECTS.rdata"))





### B. POLLINATORS ###

# selected model for species richness:


# rerun the selected models, using REML

srmod_PO <- GLMER(modelData = final.data.trans.pols, responseVar = "Species_richness", fitFamily = "poisson",
                  fixedStruct = "Predominant_land_use + Use_intensity + percNH + ncrop + Predominant_land_use:Use_intensity",
                  randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)

srmod_PO2 <- GLMER(modelData = final.data.trans.pols, responseVar = "Species_richness", fitFamily = "poisson",
                   fixedStruct = "Predominant_land_use + Use_intensity + percNH + Predominant_land_use:Use_intensity",
                   randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)


# Warning messages: none


# take a look at the model output
summary(srmod_PO$model)
summary(srmod_PO2$model)

# extract the coefficents of the model
coefs <- fixef(srmod_PO$model)
coefs2 <- fixef(srmod_PO2$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs_POLLINATORS_ncrop.csv"), row.names = F)
write.csv(coefs2, file = paste0(outdir, "/SPECIESRICHNESS_coefs_POLLINATORS.csv"), row.names = F)

# save the model output
save(srmod_PO, file = paste0(outdir, "/SRMOD_output_POLLINATORS_ncrop.rdata"))
save(srmod_PO2, file = paste0(outdir, "/SRMOD_output_POLLINATORS.rdata"))




### C. PEST CONTROLLERS ###

# selected model for species richness:

# rerun the selected models, using REML

srmod_PC <- GLMER(modelData = final.data.trans.pc, responseVar = "Species_richness", fitFamily = "poisson",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + ncrop + percNH + landcovers.5k + Predominant_land_use:ncrop + Predominant_land_use:landcovers.5k +  Use_intensity:ncrop + Use_intensity:percNH",
                  randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)

srmod_PC2 <- GLMER(modelData = final.data.trans.pc, responseVar = "Species_richness", fitFamily = "poisson",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + pest_H_log + percNH + landcovers.5k + Predominant_land_use:pest_H_log + Use_intensity:percNH",
                  randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)

# Warning messages: none



# take a look at the model output
summary(srmod_PC$model)
summary(srmod_PC2$model)

# extract the coefficents of the model
coefs <- fixef(srmod_PC$model)
coefs2 <- fixef(srmod_PC2$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs_PESTC_ncrop.csv"), row.names = F)
write.csv(coefs2, file = paste0(outdir, "/SPECIESRICHNESS_coefs_PESTC.csv"), row.names = F)

# save the model output
save(srmod_PC, file = paste0(outdir, "/SRMOD_output_PESTC_ncrop.rdata"))
save(srmod_PC2, file = paste0(outdir, "/SRMOD_output_PESTC.rdata"))



### Selected model for abundance ###


### A. ALL INSECTS ###

# selected model:


# run selected model with REML

abmod_IN <- GLMER(modelData = final.data.trans, responseVar = "logAbun", fitFamily = "gaussian",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + pest_H_log + Use_intensity:pest_H_log + Predominant_land_use:Use_intensity",
               randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

abmod_IN2 <- GLMER(modelData = final.data.trans, responseVar = "logAbun", fitFamily = "gaussian",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + pest_H_log + Use_intensity:pest_H_log + Predominant_land_use:Use_intensity",
                  randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# warnings: none


# take a look at the model output
summary(abmod_IN$model)
summary(abmod_IN2$model)

# extract the coefficents of the model
coefs <- fixef(abmod_IN$model)
coefs2 <- fixef(abmod_IN2$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/ABUNDANCE_coefs_INSECTS_ncrop.csv"), row.names = F)
write.csv(coefs2, file = paste0(outdir, "/ABUNDANCE_coefs_INSECTS.csv"), row.names = F)

# save the model output
save(abmod_IN, file = paste0(outdir, "/ABMOD_output_INSECTS_ncrop.rdata"))
save(abmod_IN, file = paste0(outdir, "/ABMOD_output_INSECTS.rdata"))




### B. POLLINATORS ###

# selected model:



# run selected model with REML

abmod_PO <- GLMER(modelData = final.data.trans.pols, responseVar = "logAbun", fitFamily = "gaussian",
                  fixedStruct = "Predominant_land_use + Use_intensity + Predominant_land_use:Use_intensity",
                  randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# both models were the same with and without ncrop variable

# warnings: none



# take a look at the model output
summary(abmod_PO$model)

# extract the coefficents of the model
coefs <- fixef(abmod_PO$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/ABUNDANCE_coefs_POLLINATORS.csv"), row.names = F)

# save the model output
save(abmod_PO, file = paste0(outdir, "/ABMOD_output_POLLINATORS.rdata"))



### C. PEST CONTROLLERS ###

# selected model:




# run selected model with REML

abmod_PC <- GLMER(modelData = final.data.trans.pc, responseVar = "logAbun", fitFamily = "gaussian",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + percNH + ncrop + landcovers.5k + Predominant_land_use:ncrop + Predominant_land_use:landcovers.5k + Predominant_land_use:percNH + Use_intensity:ncrop + Predominant_land_use:Use_intensity",
                  randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

abmod_PC2 <- GLMER(modelData = final.data.trans.pc, responseVar = "logAbun", fitFamily = "gaussian",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + percNH + pest_H_log + landcovers.5k + Predominant_land_use:percNH + Predominant_land_use:landcovers.5k + Use_intensity:pest_H_log + Predominant_land_use:Use_intensity",
                  randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# warnings:



# take a look at the model output
summary(abmod_PC$model)
summary(abmod_PC2$model)

# extract the coefficents of the model
coefs <- fixef(abmod_PC$model)
coefs2 <- fixef(abmod_PC2$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/ABUNDANCE_coefs_PESTC_ncrop.csv"), row.names = F)
write.csv(coefs2, file = paste0(outdir, "/ABUNDANCE_coefs_PESTC.csv"), row.names = F)

# save the model output
save(abmod_PC, file = paste0(outdir, "/ABMOD_output_PESTC_ncrop.rdata"))
save(abmod_PC2, file = paste0(outdir, "/ABMOD_output_PESTC.rdata"))

