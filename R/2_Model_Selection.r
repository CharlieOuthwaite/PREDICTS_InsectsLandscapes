##%######################################################%##
#                                                          #
####                   2. Run Models                    ####
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


####################### ALL INSECTS ########################

##%######################################################%##
#                                                          #
####                1. Species richness                 ####
#                                                          #
##%######################################################%##


# Not using model selection as we are interested in each variable

model_struc <- "Predominant_land_use + Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + Predominant_land_use:pest_H_RS + Predominant_land_use:ncrop_RS + Predominant_land_use:percNH_RS + Predominant_land_use:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + Predominant_land_use:Use_intensity"

sr1 <- GLMER(modelData = final.data.trans,
             responseVar = "Species_richness",
             fitFamily = "poisson",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
             optimizer = "Nelder_Mead", 
             maxIters = 50000)
# trying an alternative optimizer to remove warning

# Warning message (with optimizer = bobyqa):
#   In commonArgs(par, fn, control, environment()) :
#   maxfun < 10 * length(par)^2 is not recommended.

# Warning message (with optimizer = Nelder_Mead):
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00802261 (tol = 0.002, component 1)

summary(sr1$model)

# save the output
save(sr1, file = paste0(outdir, "/RICHNESS_Model_Set_INSECTS.rdata"))




##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##

model_data <- final.data.trans[!is.na(final.data.trans$logAbun), ] # 4681 rows

model_data <- droplevels(model_data)

length(unique(model_data$SS)) # 229
length(unique(model_data$SSBS)) # 4681



model_struc <-"Predominant_land_use +  Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + Predominant_land_use:pest_H_RS + Predominant_land_use:ncrop_RS + Predominant_land_use:percNH_RS + Predominant_land_use:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + Predominant_land_use:Use_intensity"
  
# Non selection version to assess all variables

ab1 <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
      fixedStruct = model_struc,
      randomStruct = "(1|SS)+(1|SSB)")

summary(ab1$model)

# save the output
save(ab1, file = paste0(outdir, "/ABUNDANCE_Model_Set_INSECTS.rdata"))





####################### BY ORDER ####################### 

#  load in the by order dataset

load(paste0(datadir, "/PREDICTS_dataset_TRANS_ORDERS.rdata")) # 10562 rows



##%######################################################%##
#                                                          #
####                1. Species richness                 ####
#                                                          #
##%######################################################%##


model_struc <- "Predominant_land_use + Use_intensity + fields + Order2 + pest_H_RS + ncrop_RS + percNH_RS +Predominant_land_use:Order2 + Use_intensity:Order2 + fields:Order2 + pest_H_RS:Order2 + ncrop_RS:Order2 + percNH_RS:Order2 + Predominant_land_use:pest_H_RS:Order2 + Predominant_land_use:ncrop_RS:Order2 + Predominant_land_use:percNH_RS:Order2 + Predominant_land_use:fields:Order2 + Use_intensity:pest_H_RS:Order2 + Use_intensity:ncrop_RS:Order2 + Use_intensity:percNH_RS:Order2 + Use_intensity:fields:Order2 + Predominant_land_use:Use_intensity:Order2"

sr1.or <- GLMER(modelData = final.data.trans,
             responseVar = "Species_richness",
             fitFamily = "poisson",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)"
             #, 
             #optimizer = "Nelder_Mead", 
             #maxIters = 50000
             )

# Error
# fixed-effect model matrix is rank deficient so dropping 183 columns / coefficients
# Error in pwrssUpdate(pp, resp, tol = tolPwrss, GQmat = GHrule(0L), compDev = compDev,  : 
#                        Downdated VtV is not positive definite

summary(sr1.or$model)

# save the output
save(sr1.or, file = paste0(outdir, "/RICHNESS_Model_Set_ORDERS.rdata"))







##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##

model_data <- final.data.trans[!is.na(final.data.trans$logAbun), ] # 4681 rows

model_data <- droplevels(model_data)

length(unique(model_data$SS)) # 229
length(unique(model_data$SSBS)) # 4681


model_struc <- "Predominant_land_use + Use_intensity + fields + Order2 + pest_H_RS + ncrop_RS + percNH_RS +Predominant_land_use:Order2 + Use_intensity:Order2 + fields:Order2 + pest_H_RS:Order2 + ncrop_RS:Order2 + percNH_RS:Order2 + Predominant_land_use:pest_H_RS:Order2 + Predominant_land_use:ncrop_RS:Order2 + Predominant_land_use:percNH_RS:Order2 + Predominant_land_use:fields:Order2 + Use_intensity:pest_H_RS:Order2 + Use_intensity:ncrop_RS:Order2 + Use_intensity:percNH_RS:Order2 + Use_intensity:fields:Order2 + Predominant_land_use:Use_intensity:Order2"

# Non selection version to assess all variables

ab1.or <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)")

# fixed-effect model matrix is rank deficient so dropping 235 columns / coefficients

summary(ab1.or$model)

# save the output
save(ab1, file = paste0(outdir, "/ABUNDANCE_Model_Set_INSECTS.rdata"))






#################### SERVICE PROVIDERS ##################### 

##%######################################################%##
#                                                          #
####                1. Species richness                 ####
#                                                          #
##%######################################################%##








##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##





