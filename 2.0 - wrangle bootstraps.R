
# 2.0 Combine all bootstrapped outputs into single dataframes

############
# Packages #
############
library(dplyr)
library(data.table)

# Useful vectors (color will be everything other than these) #
pattern <- c("band","stripe","blotch","stippling","spot")
ecology <- c("monophagous", "oligophagous", "polyphagous", "woody.shrub.tree","forb","grass", "leaf", "reproductive.tissue", "interior", "not.live.plants")

#######################
# Load all bootstraps #
#######################

##### Pattern ~ color associations
setwd("~/MSdata/bootstraps/patternColor-outputs") # wd should match output from bootstrap code (script 1.1)
# beta coefficients
patternColor.r.list <- list.files(pattern='patternColor_phyloGLM.r-*'); length(patternColor.r.list) 
patternColor.r.dataFiles = lapply(patternColor.r.list, read.csv); length(patternColor.r.dataFiles)
patternColor.r <- data.frame(rbindlist(patternColor.r.dataFiles, use.names=TRUE, fill = TRUE))
# p-values (raw)
patternColor.p.list <- list.files(pattern='patternColor_phyloGLM.p-*'); length(patternColor.p.list) 
patternColor.p.dataFiles = lapply(patternColor.p.list, read.csv); length(patternColor.p.dataFiles)
patternColor.p <- data.frame(rbindlist(patternColor.p.dataFiles, use.names=TRUE, fill = TRUE))

##### Ecology ~ coloration associations # wd should match output from bootstrap code (script 1.2)
setwd("~/MSdata/bootstraps/ecologyColoration-outputs")
# Regression coefficients / odds ratios
ecologyColoration.r.list <- list.files(pattern='ecologyColoration_phyloGLM.r-*'); length(ecologyColoration.r.list)
ecologyColoration.r.dataFiles = lapply(ecologyColoration.r.list, read.csv); length(ecologyColoration.r.dataFiles) 
ecologyColoration.r <- data.frame(rbindlist(ecologyColoration.r.dataFiles, use.names=TRUE, fill = TRUE))
# p-values (raw)
ecologyColoration.p.list <- list.files(pattern='ecologyColoration_phyloGLM.p-*'); length(ecologyColoration.p.list) 
ecologyColoration.p.dataFiles = lapply(ecologyColoration.p.list, read.csv); length(ecologyColoration.p.dataFiles) 
ecologyColoration.p <- data.frame(rbindlist(ecologyColoration.p.dataFiles, use.names=TRUE, fill = TRUE))

###### Also pull the bootstrapped data files 
setwd("~/MSdata/bootstraps/patternColor-outputs")
patternColor.list <- list.files(pattern='patternColor_dataBootstraps-*')
patternColor.dataFiles = lapply(patternColor.list, read.csv) 
patternColorBootstraps <- data.frame(rbindlist(patternColor.dataFiles, use.names=TRUE, fill = TRUE))

setwd("~/MSdata/bootstraps/ecologyColoration-outputs")
ecologyColoration.list <- list.files(pattern='ecologyColoration_dataBootstraps-*')
ecologyColoration.dataFiles = lapply(ecologyColoration.list, read.csv)
ecologyColorationBootstraps <- data.frame(rbindlist(ecologyColoration.dataFiles, use.names=TRUE, fill = TRUE))

# Write .csv files
setwd("~/MSdata/bootstraps-wrangled") # Create and set this new destination folder *NOTE: the current contents of this file are the wrangled bootstraps from our iterations. If you re-generate the bootstrapped results, this will/should overwrite the files already in the folder*
write.csv(patternColor.r, "patternColor_betas.csv")
write.csv(patternColor.p, "patternColor_p.csv")

write.csv(ecologyColoration.r, "ecologyPatternColor_betas.csv")
write.csv(ecologyColoration.p, "ecologyPatternColor_p.csv")

write.csv(patternColorBootstraps, "patternColor_bootstraps.csv")
write.csv(ecologyColorationBootstraps, "patternColor_bootstraps.csv")
