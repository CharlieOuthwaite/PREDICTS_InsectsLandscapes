##%######################################################%##
#                                                          #
####                 2. Model selection                 ####
#                                                          #
##%######################################################%##

# This script carries out the model selection process using the dataset generated
# in the previous script. Models are run for 2 biodiversity metrics. 

rm(list = ls())


# load libraries
library(devtools)
#install_github(repo = "timnewbold/StatisticalModels")
library(StatisticalModels)
library(cowplot)
library(gridGraphics)


# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES/"
outdir <- "2_MODEL_SELECTION/"
if(!dir.exists(outdir)) dir.create(outdir)

# read in the PREDICTS datasets with landscape variables

load(paste0(datadir, "/PREDICTS_dataset_TRANS_INSECTS.rdata")) # final.data.trans


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
                                fixedFactors = c("Predominant_land_use", "Use_intensity", "fields"),
                                fixedTerms = list(pest_H_RS = 1, ncrop_RS = 1, percNH_RS = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_RS", 
                                                      "Predominant_land_use:ncrop_RS",
                                                      "Predominant_land_use:percNH_RS",
                                                      "Predominant_land_use:fields",
                                                      "Use_intensity:pest_H_RS", 
                                                      "Use_intensity:ncrop_RS",
                                                      "Use_intensity:percNH_RS",
                                                      "Use_intensity:fields",
                                                      "Predominant_land_use:Use_intensity"), verbose = F)}) 
# save the output
save(sr1, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection_INSECTS.rdata"))

# take a look at the model output
summary(sr1$model)

# extract the stats produced as part of the model selection process
sr1stats <- as.data.frame(sr1$stats)

# save these
write.csv(sr1stats, file = paste0(outdir, "/sr_stats_INSECTS.csv"), row.names = F)


# Non selection version to assess all variables

sr2 <- GLMER(modelData = final.data.trans,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "Predominant_land_use +  Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + Predominant_land_use:pest_H_RS + Predominant_land_use:ncrop_RS + Predominant_land_use:percNH_RS + Predominant_land_use:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS+ Use_intensity:fields + Predominant_land_use:Use_intensity",
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS")

summary(sr2$model)

# save the output
save(sr2, file = paste0(outdir, "/RICHNESS_Model_Set_INSECTS.rdata"))


##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##

model_data <- final.data.trans[!is.na(final.data.trans$logAbun), ] # 4681 rows

length(unique(model_data$SS)) # 229
length(unique(model_data$SSBS)) # 4681

### A. ALL INSECTS ######SSBS A. ALL INSECTS ###


system.time({ab1 <- GLMERSelect(modelData = final.data.trans, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Use_intensity", "fields"),
                                fixedTerms = list(pest_H_RS = 1, ncrop_RS = 1, percNH_RS = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:pest_H_RS", 
                                                      "Predominant_land_use:ncrop_RS",
                                                      "Predominant_land_use:percNH_RS",
                                                      "Predominant_land_use:fields",
                                                      "Use_intensity:pest_H_RS", 
                                                      "Use_intensity:ncrop_RS",
                                                      "Use_intensity:percNH_RS",
                                                      "Use_intensity:fields",
                                                      "Predominant_land_use:Use_intensity"), verbose = F)}) 
# take a look at the model output
summary(ab1$model)

# save the output
save(ab1, file = paste0(outdir, "/ABUNDANCE_Model_Selection_INSECTS.rdata"))

# extract the stats produced as part of the model selection process
ab1stats <- as.data.frame(ab1$stats)

# save these
write.csv(ab1stats, file = paste0(outdir, "/ab_stats_INSECTS.csv"), row.names = F)


# Non selection version to assess all variables

ab2 <- GLMER(modelData = final.data.trans,responseVar = "logAbun",fitFamily = "gaussian",
      fixedStruct = "Predominant_land_use +  Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + Predominant_land_use:pest_H_RS + Predominant_land_use:ncrop_RS + Predominant_land_use:percNH_RS + Predominant_land_use:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS+ Use_intensity:fields + Predominant_land_use:Use_intensity",
      randomStruct = "(1|SS)+(1|SSB)")

summary(ab2$model)

# save the output
save(ab2, file = paste0(outdir, "/ABUNDANCE_Model_Set_INSECTS.rdata"))




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

