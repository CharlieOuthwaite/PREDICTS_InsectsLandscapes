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

load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS_INSECTS.rdata")) # final.data.trans
load(paste0(datadir, "/PREDICTS_abun_subset_INSECTS.rdata")) # final.data.abun
load(paste0(datadir, "/PREDICTS_rcar_subset_INSECTS.rdata")) # final.data.rcar



#### Run the model selection process including all landscape variables ####


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
                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                      "Predominant_land_use:Hansen_mindist_log",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:homogen", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:fert.total_log", 
                                                      "Use_intensity:Hansen_mindist_log",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:homogen",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:Hansen_mindist_log"), verbose = F)}) 
# save the output
save(sr1, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection_INSECTS.rdata"))

# take a look at the model output
summary(sr1$model)

# extract the stats produced as part of the model selection process
sr1stats <- as.data.frame(sr1$stats)

# save these
write.csv(sr1stats, file = paste0(outdir, "/sr_stats_INSECTS.csv"), row.names = F)



### B. POLLINATORS ###

system.time({sr2 <- GLMERSelect(modelData = final.data.trans, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                      "Predominant_land_use:Hansen_mindist_log",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:homogen", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:fert.total_log", 
                                                      "Use_intensity:Hansen_mindist_log",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:homogen",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:Hansen_mindist_log"), verbose = F)}) 
# save the output
save(sr2, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection_POLLINATORS.rdata"))

# take a look at the model output
summary(sr2$model)

# extract the stats produced as part of the model selection process
sr2stats <- as.data.frame(sr2$stats)

# save these
write.csv(sr2stats, file = paste0(outdir, "/sr_stats_POLLINATORS.csv"), row.names = F)





### C. PEST CONTROLLERS ###


system.time({sr3 <- GLMERSelect(modelData = final.data.trans, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                      "Predominant_land_use:Hansen_mindist_log",
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:homogen", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:fert.total_log", 
                                                      "Use_intensity:Hansen_mindist_log",
                                                      "Use_intensity:landcovers.5k",
                                                      "Use_intensity:homogen",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity",
                                                      "Tropical:Hansen_mindist_log"), verbose = F)}) 
# save the output
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


### A. ALL INSECTS ###


system.time({ab1 <- GLMERSelect(modelData = final.data.abun, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                      "Predominant_land_use:Hansen_mindist_log", 
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:homogen", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:fert.total_log", 
                                                      "Use_intensity:Hansen_mindist_log", 
                                                      "Use_intensity:landcovers.5k", 
                                                      "Use_intensity:homogen",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity", 
                                                      "Tropical:Hansen_mindist_log"), verbose = F)})

# take a look at the model output
summary(ab1$model)

# save the output
save(ab1, file = paste0(outdir, "/ABUNDANCE_Model_Selection_INSECTS.rdata"))

# extract the stats produced as part of the model selection process
ab1stats <- as.data.frame(ab1$stats)

# save these
write.csv(ab1stats, file = paste0(outdir, "/ab_stats_INSECTS.csv"), row.names = F)



### B. POLLINATORS ###


system.time({ab2 <- GLMERSelect(modelData = final.data.abun, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                      "Predominant_land_use:Hansen_mindist_log", 
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:homogen", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:fert.total_log", 
                                                      "Use_intensity:Hansen_mindist_log", 
                                                      "Use_intensity:landcovers.5k", 
                                                      "Use_intensity:homogen",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity", 
                                                      "Tropical:Hansen_mindist_log"), verbose = F)})

# take a look at the model output
summary(ab2$model)

# save the output
save(ab2, file = paste0(outdir, "/ABUNDANCE_Model_Selection_POLLINATORS.rdata"))

# extract the stats produced as part of the model selection process
ab2stats <- as.data.frame(ab2$stats)

# save these
write.csv(ab2stats, file = paste0(outdir, "/ab_stats_POLLINATORS.csv"), row.names = F)




### C. PEST CONTROLLERS ###


system.time({ab3 <- GLMERSelect(modelData = final.data.abun, 
                                responseVar = "logAbun",
                                fitFamily = "gaussian", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                randomStruct = "(1|SS)+(1|SSB)", 
                                fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                      "Predominant_land_use:Hansen_mindist_log", 
                                                      "Predominant_land_use:landcovers.5k", 
                                                      "Predominant_land_use:homogen", 
                                                      "Predominant_land_use:percNH",
                                                      "Use_intensity:fert.total_log", 
                                                      "Use_intensity:Hansen_mindist_log", 
                                                      "Use_intensity:landcovers.5k", 
                                                      "Use_intensity:homogen",
                                                      "Use_intensity:percNH",
                                                      "Predominant_land_use:Use_intensity", 
                                                      "Tropical:Hansen_mindist_log"), verbose = F)})

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
####                      3. RCAR                       ####
#                                                          #
##%######################################################%##


### A. ALL INSECTS ###


system.time({rcar1 <- GLMERSelect(modelData = final.data.rcar, responseVar = "RCAR_110km",
                                  fitFamily = "gaussian", 
                                  fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                  fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                  randomStruct = "(1|SS)+(1|SSB)", 
                                  fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                        "Predominant_land_use:Hansen_mindist_log", 
                                                        "Predominant_land_use:landcovers.5k", 
                                                        "Predominant_land_use:homogen", 
                                                        "Predominant_land_use:percNH",
                                                        "Use_intensity:fert.total_log", 
                                                        "Use_intensity:Hansen_mindist_log", 
                                                        "Use_intensity:landcovers.5k",
                                                        "Use_intensity:homogen",
                                                        "Use_intensity:percNH",
                                                        "Predominant_land_use:Use_intensity",
                                                        "Tropical:Hansen_mindist_log"), verbose = F)}
)

# take a look at the model output
summary(rcar1$model)

# save the model output
save(rcar1, file = paste0(outdir, "/RCAR_Model_Selection_INSECTS.rdata"))

# extract the stats produced as part of the model selection process
rcar1stats <- as.data.frame(rcar1$stats)

# save these
write.csv(rcar1stats, file = paste0(outdir, "/rcar_stats_INSECTS.csv"), row.names = F)



### B. POLLINATORS ###


system.time({rcar2 <- GLMERSelect(modelData = final.data.rcar, responseVar = "RCAR_110km",
                                  fitFamily = "gaussian", 
                                  fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                  fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                  randomStruct = "(1|SS)+(1|SSB)", 
                                  fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                        "Predominant_land_use:Hansen_mindist_log", 
                                                        "Predominant_land_use:landcovers.5k", 
                                                        "Predominant_land_use:homogen", 
                                                        "Predominant_land_use:percNH",
                                                        "Use_intensity:fert.total_log", 
                                                        "Use_intensity:Hansen_mindist_log", 
                                                        "Use_intensity:landcovers.5k",
                                                        "Use_intensity:homogen",
                                                        "Use_intensity:percNH",
                                                        "Predominant_land_use:Use_intensity",
                                                        "Tropical:Hansen_mindist_log"), verbose = F)}
)

# take a look at the model output
summary(rcar2$model)

# save the model output
save(rcar2, file = paste0(outdir, "/RCAR_Model_Selection_POLLINATORS.rdata"))

# extract the stats produced as part of the model selection process
rcar2stats <- as.data.frame(rcar2$stats)

# save these
write.csv(rcar2stats, file = paste0(outdir, "/rcar_stats_POLLINATORS.csv"), row.names = F)




### C. PEST CONTROLLERS ###


system.time({rcar3 <- GLMERSelect(modelData = final.data.rcar, responseVar = "RCAR_110km",
                                  fitFamily = "gaussian", 
                                  fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                  fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                  randomStruct = "(1|SS)+(1|SSB)", 
                                  fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                        "Predominant_land_use:Hansen_mindist_log", 
                                                        "Predominant_land_use:landcovers.5k", 
                                                        "Predominant_land_use:homogen", 
                                                        "Predominant_land_use:percNH",
                                                        "Use_intensity:fert.total_log", 
                                                        "Use_intensity:Hansen_mindist_log", 
                                                        "Use_intensity:landcovers.5k",
                                                        "Use_intensity:homogen",
                                                        "Use_intensity:percNH",
                                                        "Predominant_land_use:Use_intensity",
                                                        "Tropical:Hansen_mindist_log"), verbose = F)}
)

# take a look at the model output
summary(rcar3$model)

# save the model output
save(rcar3, file = paste0(outdir, "/RCAR_Model_Selection_POLLINATORS.rdata"))

# extract the stats produced as part of the model selection process
rcar3stats <- as.data.frame(rcar3$stats)

# save these
write.csv(rcar3stats, file = paste0(outdir, "/rcar_stats_POLLINATORS.csv"), row.names = F)





##%######################################################%##
#                                                          #
####        Run final selected model using REML         ####
#                                                          #
##%######################################################%##



### selected model for species richness ###

### A. ALL INSECTS ###

# rerun the selected models, using REML

srmod_IN <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:landcovers.5k + Predominant_land_use:homogen +  Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
               randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)


# Warning messages:




# take a look at the model output
summary(srmod_IN$model)

# extract the coefficents of the model
coefs <- fixef(srmod_IN$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs_INSECTS.csv"), row.names = F)

# save the model output
save(srmod_IN, file = paste0(outdir, "/SRMOD_output_INSECTS.rdata"))





### B. POLLINATORS ###

# selected model for species richness:



# rerun the selected models, using REML

srmod_PO <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:landcovers.5k + Predominant_land_use:homogen +  Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
                  randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)


# Warning messages:




# take a look at the model output
summary(srmod_PO$model)

# extract the coefficents of the model
coefs <- fixef(srmod_PO$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs_POLLINATORS.csv"), row.names = F)

# save the model output
save(srmod_PO, file = paste0(outdir, "/SRMOD_output_POLLINATORS.rdata"))






### C. PEST CONTROLLERS ###

# selected model for species richness:



# rerun the selected models, using REML

srmod_PC <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:landcovers.5k + Predominant_land_use:homogen +  Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
                  randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)


# Warning messages:




# take a look at the model output
summary(srmod_PC$model)

# extract the coefficents of the model
coefs <- fixef(srmod_PC$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs_PESTC.csv"), row.names = F)

# save the model output
save(srmod_PC, file = paste0(outdir, "/SRMOD_output_PESTC.rdata"))



### Selected model for abundance ###


### A. ALL INSECTS ###

# selected model:




# run selected model with REML

abmod_IN <- GLMER(modelData = final.data.abun, responseVar = "logAbun", fitFamily = "gaussian",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:percNH + Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
               randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)



# warnings:



# take a look at the model output
summary(abmod_IN$model)

# extract the coefficents of the model
coefs <- fixef(abmod_IN$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/ABUNDANCE_coefs_INSECTS.csv"), row.names = F)

# save the model output
save(abmod_IN, file = paste0(outdir, "/ABMOD_output_INSECTS.rdata"))




### B. POLLINATORS ###

# selected model:




# run selected model with REML

abmod_PO <- GLMER(modelData = final.data.abun, responseVar = "logAbun", fitFamily = "gaussian",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:percNH + Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
                  randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)



# warnings:



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

abmod_PC <- GLMER(modelData = final.data.abun, responseVar = "logAbun", fitFamily = "gaussian",
                  fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:percNH + Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
                  randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)



# warnings:



# take a look at the model output
summary(abmod_PC$model)

# extract the coefficents of the model
coefs <- fixef(abmod_PC$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/ABUNDANCE_coefs_PESTC.csv"), row.names = F)

# save the model output
save(abmod_PC, file = paste0(outdir, "/ABMOD_output_PESTC.rdata"))

