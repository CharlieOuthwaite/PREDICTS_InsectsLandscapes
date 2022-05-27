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

model_struc <- "LU + Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + LU:pest_H_RS + LU:ncrop_RS + LU:percNH_RS + LU:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + LU:Use_intensity"

sr1 <- GLMER(modelData = final.data.trans,
             responseVar = "Species_richness",
             fitFamily = "poisson",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
             #optimizer = "Nelder_Mead", 
             maxIters = 50000)


# no warnings!

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



model_struc <- "LU + Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + LU:pest_H_RS + LU:ncrop_RS + LU:percNH_RS + LU:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + LU:Use_intensity"
  
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


model_struc <- "LU + Use_intensity + fields + Order2 + pest_H_RS + ncrop_RS + percNH_RS +LU:Order2 + Use_intensity:Order2 + fields:Order2 + pest_H_RS:Order2 + ncrop_RS:Order2 + percNH_RS:Order2 + LU:pest_H_RS:Order2 + LU:ncrop_RS:Order2 + LU:percNH_RS:Order2 + LU:fields:Order2 + Use_intensity:pest_H_RS:Order2 + Use_intensity:ncrop_RS:Order2 + Use_intensity:percNH_RS:Order2 + Use_intensity:fields:Order2 + LU:Use_intensity:Order2"

# including fewer interactions
model_struc2 <- "LU + Use_intensity + fields + Order2 + pest_H_RS + ncrop_RS + percNH_RS +LU:Order2 + Use_intensity:Order2 + fields:Order2 + pest_H_RS:Order2 + ncrop_RS:Order2 + percNH_RS:Order2 + LU:pest_H_RS:Order2 + LU:ncrop_RS:Order2 + LU:percNH_RS:Order2 + LU:fields:Order2 + Use_intensity:pest_H_RS:Order2 + Use_intensity:ncrop_RS:Order2 + Use_intensity:percNH_RS:Order2 + Use_intensity:fields:Order2 + LU:Use_intensity:Order2"

sr1.or <- GLMER(modelData = final.data.trans,
             responseVar = "Species_richness",
             fitFamily = "poisson",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)"
             , 
             #optimizer = "Nelder_Mead", 
             maxIters = 50000
             )

# Warning message:
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00818785 (tol = 0.002, component 1)

summary(sr1.or$model)

# save the output
save(sr1.or, file = paste0(outdir, "/RICHNESS_Model_Set_ORDERS.rdata"))







##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##

model_data <- final.data.trans[!is.na(final.data.trans$logAbun), ] # 10016 rows

model_data <- droplevels(model_data)

length(unique(model_data$SS)) # 229
length(unique(model_data$SSBS)) # 4681


model_struc <- "LU + Use_intensity + fields + Order2 + pest_H_RS + ncrop_RS + percNH_RS +LU:Order2 + Use_intensity:Order2 + fields:Order2 + pest_H_RS:Order2 + ncrop_RS:Order2 + percNH_RS:Order2 + LU:pest_H_RS:Order2 + LU:ncrop_RS:Order2 + LU:percNH_RS:Order2 + LU:fields:Order2 + Use_intensity:pest_H_RS:Order2 + Use_intensity:ncrop_RS:Order2 + Use_intensity:percNH_RS:Order2 + Use_intensity:fields:Order2 + LU:Use_intensity:Order2"

# Non selection version to assess all variables

ab1.or <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)")

# fixed-effect model matrix is rank deficient so dropping 194 columns / coefficients

summary(ab1.or$model)

# save the output
save(ab1.or, file = paste0(outdir, "/ABUNDANCE_Model_Set_ORDERS.rdata"))






#################### SERVICE PROVIDERS ##################### 

# load datasets
load(paste0(datadir, "/PREDICTS_dataset_TRANS_POLS.rdata")) # final.data.trans
load(paste0(datadir, "/PREDICTS_dataset_TRANS_PCS.rdata")) # final.data.trans


##%######################################################%##
#                                                          #
####                1. Species richness                 ####
#                                                          #
##%######################################################%##


## pollinators


model_struc <- "LU + Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + LU:pest_H_RS + LU:ncrop_RS + LU:percNH_RS + LU:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + LU:Use_intensity"

sr1.pols <- GLMER(modelData = final.data.trans.pols,
             responseVar = "Species_richness",
             fitFamily = "poisson",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
             #optimizer = "Nelder_Mead", 
             maxIters = 50000)

summary(sr1.pols$model)

# no warnings!

# save the output
save(sr1.pols, file = paste0(outdir, "/RICHNESS_Model_Set_POLS.rdata"))


## pest controllers


model_struc <- "LU + Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + LU:pest_H_RS + LU:ncrop_RS + LU:percNH_RS + LU:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + LU:Use_intensity"

sr1.pc <- GLMER(modelData = final.data.trans.pc,
             responseVar = "Species_richness",
             fitFamily = "poisson",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
             #optimizer = "Nelder_Mead", 
             maxIters = 50000)

summary(sr1.pc$model)

# no warnings!

# save the output
save(sr1.pc, file = paste0(outdir, "/RICHNESS_Model_Set_PCS.rdata"))







##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##


## pollinators


model_data <- final.data.trans.pols[!is.na(final.data.trans.pols$logAbun), ] # 3401 rows

model_data <- droplevels(model_data)

length(unique(model_data$SS)) # 139
length(unique(model_data$SSBS)) # 3401



model_struc <- "LU +  Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + LU:pest_H_RS + LU:ncrop_RS + LU:percNH_RS + LU:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + LU:Use_intensity"

# Non selection version to assess all variables

ab1.pols <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)")

summary(ab1.pols$model)

# no warnings!

# save the output
save(ab1.pols, file = paste0(outdir, "/ABUNDANCE_Model_Set_POLS.rdata"))


## pest controllers


model_data <- final.data.trans.pc[!is.na(final.data.trans.pc$logAbun), ] # 2158 rows

model_data <- droplevels(model_data)

length(unique(model_data$SS)) # 116
length(unique(model_data$SSBS)) # 2158



model_struc <- "LU +  Use_intensity + fields + pest_H_RS + ncrop_RS + percNH_RS + LU:pest_H_RS + LU:ncrop_RS + LU:percNH_RS + LU:fields + Use_intensity:pest_H_RS + Use_intensity:ncrop_RS + Use_intensity:percNH_RS + Use_intensity:fields + LU:Use_intensity"

# Non selection version to assess all variables

ab1.pc <- GLMER(modelData = model_data,responseVar = "logAbun",fitFamily = "gaussian",
             fixedStruct = model_struc,
             randomStruct = "(1|SS)+(1|SSB)")

summary(ab1.pc$model)

#  no warnings!

# save the output
save(ab1.pc, file = paste0(outdir, "/ABUNDANCE_Model_Set_PCS.rdata"))





