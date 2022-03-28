
# 3.0 Summarize stats across bootstraps, penalize for multiple comparisons, generate results tables
# This scripts generates Tables S1, S4, and S5

############
# Packages #
############
library(dplyr)
library(data.table)

########
# Data # 
########
setwd("~/MSdata/bootstraps-wrangled") # Set to destination folder written to in script 2.0

# Beta coefficients
patCol_beta0 <- read.csv("patternColor_betas.csv")
ecoPatCol_beta0 <- read.csv("ecologyPatternColor_betas.csv")
# P-values 
patCol_p0 <- read.csv("patternColor_p.csv")
ecoPatCol_p0 <- read.csv("ecologyPatternColor_p.csv")
# Bootstrapped raw data
patCol_boot0 <- read.csv("patternColor_bootstraps.csv")
ecoPatCol_boot0 <- read.csv("ecologyPatternColor_bootstraps.csv")

# list of a-priori hypotheses
aPriori0 <- read.csv("~/MSdata/data/aPrioris.csv")
aPriori <- subset(aPriori0, response != ""); nrow(aPriori)
head(aPriori)
# Make sure no typos in any of the a-prioris
aPrioriCheck <- unique(c(as.character(unique(aPriori$response)), unique(as.character(aPriori$predictor))))
allTraits <- unique(c(as.character(unique(ecoPatCol_p0$predictor)), unique(as.character(ecoPatCol_p0$response))))
setdiff(aPrioriCheck, allTraits) # should be 0. If not, sometimes wrote order of colors incorrectly in a-prioris (e.g. "green.brown" instead of "brown.green"); fix manually in a-prioris

# Useful vectors (color will be everything other than these) #
pattern <- c("band","stripe","blotch","stippling","spot")
ecology <- c("monophagous", "oligophagous", "polyphagous", "woody.shrub.tree","forb","grass", "leaf", "reproductive.tissue", "interior", "not.live.plants")

##############################################
# Summarize counts/frequencies of each trait # # this makes Table S1
##############################################
# Average frequency of each trait across all the bootstraps
# pattern ~ color
boot <- unique(patCol_boot0$bootstrap)
nboot <- length(boot)

tally_patCol0 <- list()

for(i in 1:nboot){
  #  temp <- subset(patCol_boot0, bootstrap == 5)
  temp <- subset(patCol_boot0, bootstrap == boot[i])
  temp$X <- NULL; temp$bootstrap <- NULL; temp$genus_species <- NULL; temp$numCol <- NULL
  # column sums, ignoring NAs
  tally <- data.frame(colSums(temp, na.rm=TRUE))
  tally$trait <- rownames(tally)
  names(tally) <- c("count", "trait")
  tally$bootstrap <- i
  # save output
  tally_patCol0 [[i]] <- tally
}

patCol_count0 <- do.call(rbind, tally_patCol0); nrow(patCol_count0); head(patCol_count0)
patCol_count <- patCol_count0[!grepl('X',patCol_count0$trait),]; head(patCol_count) # the "X" column gets carried over sometimes
patCol_count_summ <- do.call(data.frame, aggregate(count ~ trait, data=patCol_count, function(x) mean(x))); head(patCol_count_summ); nrow(patCol_count_summ)

# ecology ~ pattern & color
boot <- unique(ecoPatCol_boot0$bootstrap)
nboot <- length(boot)

tally_ecoPatCol0 <- list()

for(i in 1:nboot){
  #  temp <- subset(ecoPatCol_boot0, bootstrap == 5)
  temp <- subset(ecoPatCol_boot0, bootstrap == boot[i])
  temp$X <- NULL; temp$bootstrap <- NULL; temp$genus_species <- NULL; temp$numCol <- NULL
  # column sums, ignoring NAs
  tally <- data.frame(colSums(temp, na.rm=TRUE))
  tally$trait <- rownames(tally)
  names(tally) <- c("count", "trait")
  tally$bootstrap <- i
  # save output
  tally_ecoPatCol0 [[i]] <- tally
}

ecoPatCol_count0 <- do.call(rbind, tally_ecoPatCol0); nrow(ecoPatCol_count0); head(ecoPatCol_count0)
ecoPatCol_count <- ecoPatCol_count0[!grepl('X',ecoPatCol_count0$trait),]; head(ecoPatCol_count) # the "X" column gets carried over sometimes
ecoPatCol_count_summ <- do.call(data.frame, aggregate(count ~ trait, data=ecoPatCol_count, function(x) mean(x))); head(ecoPatCol_count_summ); nrow(ecoPatCol_count_summ)

# rbind & take the average
both <- rbind(ecoPatCol_count, patCol_count); head(both)
both_summ0 <- do.call(data.frame, aggregate(count ~ trait, data=both, function(x) mean(x))); head(both_summ0); nrow(both_summ0) # 140 traits total. 

# Identify low-frequency traits (present in 5 or fewer species, on average, across bootstraps)
both_summ <- both_summ0[order(both_summ0$count),]
nrow(subset(both_summ, count <= 5.4)) # 9 traits fit in this category 
min5 <- unique(subset(both_summ, count <= 5.4)$trait) # remove these traits
min5_df <- data.frame(min5); head(min5_df)
min5_df$numCol <- str_count(min5_df$min5, "\\.")
min5_df2 <- subset(min5_df, numCol != 2)

# Calculate proportions of caterpillars with different traits; export
# 1) Round counts to nearest whole number (integer)
# 2) Remove traits with < 5 instances
# 3) Calculate proportion, out of 1808 species
# 4) After export, remember to remove all 3-color combos
traitTally <- subset(both_summ, !trait %in% min5); nrow(traitTally)
traitTally$count_int <- round(traitTally$count, digits = 0); head(traitTally)
traitTally$propOfSpecies <- traitTally$count_int/1808
traitTally$percentOfSpecies <- round((traitTally$propOfSpecies*100), digits=0)
traitTally$numCols0 <- str_count(traitTally$trait, "\\."); head(traitTally)
traitTally$numCols <- traitTally$numCols0+1
traitTally_export <- subset(traitTally, select=c("trait", "count_int", "percentOfSpecies"))
names(traitTally_export) <- c("trait", "containedIn_nSpecies", "percentOfSpecies")

setwd("~/MSdata/outputs")
write.csv(traitTally_export, "traitFrequencies - TableS1.csv")

################################################
################################################
# Summarize beta coefficients & calculate odds # 
################################################
################################################
head(patCol_beta0)

# Pattern ~ Color
patCol_beta_summ0 <- do.call(data.frame, aggregate(beta ~ response + predictor, data=patCol_beta0, function(x) c(nBeta = length(x), mean = mean(x), sd = sd(x), min = min(x), max=max(x)))); head(patCol_beta_summ0) 
patCol_beta_summ0$odds <- (exp(patCol_beta_summ0$beta.mean)-1)*100
names(patCol_beta_summ0) <- c("response","predictor","nBeta","meanBeta","sdBeta", "minBeta", "maxBeta", "oddsRatio")
patCol_beta_summ <- subset(patCol_beta_summ0, !predictor %in% min5); nrow(patCol_beta_summ0); nrow(patCol_beta_summ) # N = 245 pattern ~ color models

# Ecology ~ Pattern / Color
ecoPatCol_beta_summ0 <- do.call(data.frame, aggregate(beta ~ response + predictor, data=ecoPatCol_beta0, function(x) c(nBeta = length(x), mean = mean(x), sd = sd(x), min = min(x), max=max(x)))); head(ecoPatCol_beta_summ0); nrow(ecoPatCol_beta_summ0)
ecoPatCol_beta_summ0$odds <- (exp(ecoPatCol_beta_summ0$beta.mean)-1)*100
names(ecoPatCol_beta_summ0) <- c("response","predictor","nBeta","meanBeta","sdBeta", "minBeta", "maxBeta", "oddsRatio")
ecoPatCol_beta_summ <- subset(ecoPatCol_beta_summ0, !predictor %in% min5); nrow(ecoPatCol_beta_summ0) # N = 540 ecology ~ paattern/color models

################
################
# Significance #
################
################

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse between a-Priori and exploratory #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Do this for each of the two questions (pattern ~ color, m1; ecology ~ color/pattern, m3)
head(aPriori); nrow(aPriori)
aPriori$hypothesis <- paste(aPriori$response, "~", aPriori$predictor, sep=" "); head(aPriori)

# Pattern ~ color
patCol_p <- subset(patCol_p0, !predictor %in% min5) # remove rare traits
patCol_p$hypothesis <- paste(patCol_p$response, "~", patCol_p$predictor, sep=" "); head(patCol_p)
allHypotheses_patCol <- unique(patCol_p$hypothesis)
aPriori_patCol <- subset(allHypotheses_patCol,  allHypotheses_patCol %in% aPriori$hypothesis) # N = 45 a-prior 
explore_patCol <- subset(allHypotheses_patCol,  !allHypotheses_patCol %in% aPriori$hypothesis) # N = 200 exploratory 

# Ecology ~ pattern / color
ecoPatCol_p <- subset(ecoPatCol_p, !predictor %in% min5) # remove rare traits
ecoPatCol_p$hypothesis <- paste(ecoPatCol_p$response, "~", ecoPatCol_p$predictor, sep=" ")
allHypotheses_ecoPatCol <- unique(ecoPatCol_p$hypothesis)
aPriori_ecoPatCol <- subset(allHypotheses_ecoPatCol,  allHypotheses_ecoPatCol %in% aPriori$hypothesis) # 84 a-priori 
explore_ecoPatCol <- subset(allHypotheses_ecoPatCol,  !allHypotheses_ecoPatCol %in% aPriori$hypothesis) # 456 exploratory 

# CHECK
ifelse(length(aPriori_patCol) + length(aPriori_ecoPatCol) == length(aPriori$hypothesis), "ONWARD HO", "STOP")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#      Pattern ~ color      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# For each set of models (pattern ~ color, ecology ~ pattern/color)...
# - Loop through each bootstrap; subset just to a-priori hypotheses & just exploratory hypotheses; penalize p-vals within those subsets

#-------------------------#
#------- penalize --------#
#-------------------------#
boot <- unique(patCol_p0$bootstrap)
nboot <- length(boot) # 1000 bootstraps

aPriori_patCol0 <- list()
explore_patCol0 <- list()

for(i in 1:nboot){
  temp <- subset(patCol_p, bootstrap == boot[i])
  # subset each bootstrap to a-prioris & penalize
  temp_aPriori <- subset(temp, temp$hypothesis %in% aPriori_patCol)
  temp_aPriori$p_fdr <- p.adjust(temp_aPriori$p_raw, method="fdr", n=length(temp_aPriori$p_raw))
  temp_aPriori$numComp <- length(temp_aPriori$p_raw)
  # subset each bootstrap to exploratories & penalize
  temp_explore <- subset(temp, temp$hypothesis %in% explore_patCol)
  temp_explore$p_fdr <- p.adjust(temp_explore$p_raw, method="fdr", n=length(temp_explore$p_raw))
  temp_explore$numComp <- length(temp_explore$p_raw)
  # save outputs
  aPriori_patCol0[[i]] <- temp_aPriori
  explore_patCol0[[i]] <- temp_explore
}

aPriori_patCol_p_all <- do.call(rbind, aPriori_patCol0); aPriori_patCol_p_all$hypothesisType <- "aPriori"
explore_patCol_p_all <- do.call(rbind, explore_patCol0); head(explore_patCol_p_all); explore_patCol_p_all$hypothesisType <- "exploratory"

# Merge them together
patCol_p_all <- rbind(aPriori_patCol_p_all, explore_patCol_p_all); patCol_p_all$question <- "pattern.color"; head(patCol_p_all); nrow(patCol_p_all)

#--------------------------------------------------------------------------#
#------- Calculate of % of p-values in different significance bins --------#
#--------------------------------------------------------------------------#
# Raw p-values
patCol_p_all$p_raw_05 <- 0
patCol_p_all$p_raw_05[which(patCol_p_all$p_raw <= .05)] <- 1
patCol_p_all$p_raw_01 <- 0
patCol_p_all$p_raw_01[which(patCol_p_all$p_raw <= .01)] <- 1
patCol_p_all$p_raw_001 <- 0
patCol_p_all$p_raw_001[which(patCol_p_all$p_raw <= .001)] <- 1

head(patCol_p_all)
patCol_p05_summ0 <- do.call(data.frame, aggregate(p_raw_05 ~ question + hypothesis + hypothesisType + response + predictor, data=patCol_p_all, function(x) propSig05_raw = sum(x)/length(x))) 
patCol_p01_summ0 <- do.call(data.frame, aggregate(p_raw_01 ~ question + hypothesis + hypothesisType + response + predictor, data=patCol_p_all, function(x) propSig05_raw = sum(x)/length(x))) 
patCol_p001_summ0 <- do.call(data.frame, aggregate(p_raw_001 ~ question + hypothesis + hypothesisType + response + predictor, data=patCol_p_all, function(x) propSig05_raw = sum(x)/length(x)))
patCol_summ_raw0 <- merge(patCol_p05_summ0, patCol_p01_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
patCol_summ_raw <- merge(patCol_summ_raw0, patCol_p001_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
head(patCol_summ_raw); nrow(patCol_summ_raw) # 245 total pattern ~ color comparisons 

# Adjusted p-values
patCol_p_all$p_fdr_05 <- 0
patCol_p_all$p_fdr_05[which(patCol_p_all$p_fdr <= .05)] <- 1
patCol_p_all$p_fdr_01 <- 0
patCol_p_all$p_fdr_01[which(patCol_p_all$p_fdr <= .01)] <- 1
patCol_p_all$p_fdr_001 <- 0
patCol_p_all$p_fdr_001[which(patCol_p_all$p_fdr <= .001)] <- 1

patCol_p05_summ0 <- do.call(data.frame, aggregate(p_fdr_05 ~ question + hypothesis + hypothesisType + response + predictor, data=patCol_p_all, function(x) propSig05_fdr = sum(x)/length(x)))
patCol_p01_summ0 <- do.call(data.frame, aggregate(p_fdr_01 ~ question + hypothesis + hypothesisType + response + predictor, data=patCol_p_all, function(x) propSig05_fdr = sum(x)/length(x)))
patCol_p001_summ0 <- do.call(data.frame, aggregate(p_fdr_001 ~ question + hypothesis + hypothesisType + response + predictor, data=patCol_p_all, function(x) propSig05_fdr = sum(x)/length(x)))
patCol_summ_fdr0 <- merge(patCol_p05_summ0, patCol_p01_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
patCol_summ_fdr <- merge(patCol_summ_fdr0, patCol_p001_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
head(patCol_summ_fdr); nrow(patCol_summ_fdr) # 245 total pattern ~ color comparisons

# Add asterisks. All significance is based on the fdr-corrected values now
patCol_summ_p0 <- merge(patCol_summ_raw, patCol_summ_fdr); head(patCol_summ_p0)
patCol_summ_p0$sig <- ""
patCol_summ_p0$sig[which(patCol_summ_p0$p_fdr_05 >= .9)] <- "*"
patCol_summ_p0$sig[which(patCol_summ_p0$p_fdr_01 >= .9)] <- "**"
patCol_summ_p0$sig[which(patCol_summ_p0$p_fdr_001 >= .9)] <- "***"

patCol_summ_p <- patCol_summ_p0[order(-patCol_summ_p0$p_fdr_05),]
head(patCol_summ_p); nrow(patCol_summ_p) # 245 total (40 a-priori + 205 exploratory)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#      Ecology ~ Pattern / color      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-------------------------#
#------- penalize --------#
#-------------------------#
boot <- unique(ecoPatCol_p0$bootstrap)
nboot <- length(boot) # 1000 bootstraps

aPriori_ecoPatCol0 <- list()
explore_ecoPatCol0 <- list()

for(i in 1:nboot){
  temp <- subset(ecoPatCol_p, bootstrap == boot[i])
  # subset each bootstrap to a-prioris & penalize
  temp_aPriori <- subset(temp, temp$hypothesis %in% aPriori_ecoPatCol)
  temp_aPriori$p_fdr <- p.adjust(temp_aPriori$p_raw, method="fdr", n=length(temp_aPriori$p_raw))
  temp_aPriori$numComp <- length(temp_aPriori$p_raw)
  # subset each bootstrap to exploratories & penalize
  temp_explore <- subset(temp, temp$hypothesis %in% explore_ecoPatCol)
  temp_explore$p_fdr <- p.adjust(temp_explore$p_raw, method="fdr", n=length(temp_explore$p_raw))
  temp_explore$numComp <- length(temp_explore$p_raw)
  # save outputs
  aPriori_ecoPatCol0[[i]] <- temp_aPriori
  explore_ecoPatCol0[[i]] <- temp_explore
}

aPriori_ecoPatCol_p_all <- do.call(rbind, aPriori_ecoPatCol0); aPriori_ecoPatCol_p_all$hypothesisType <- "aPriori"
explore_ecoPatCol_p_all <- do.call(rbind, explore_ecoPatCol0); explore_ecoPatCol_p_all$hypothesisType <- "exploratory"

# Merge them together
ecoPatCol_p_all <- rbind(aPriori_ecoPatCol_p_all, explore_ecoPatCol_p_all); ecoPatCol_p_all$question <- "ecology.pattern.color"; head(ecoPatCol_p_all); nrow(ecoPatCol_p_all)

#--------------------------------------------------------------------------#
#------- Calculate of % of p-values in different significance bins --------#
#--------------------------------------------------------------------------#
# Raw p-values
ecoPatCol_p_all$p_raw_05 <- 0
ecoPatCol_p_all$p_raw_05[which(ecoPatCol_p_all$p_raw <= .05)] <- 1
ecoPatCol_p_all$p_raw_01 <- 0
ecoPatCol_p_all$p_raw_01[which(ecoPatCol_p_all$p_raw <= .01)] <- 1
ecoPatCol_p_all$p_raw_001 <- 0
ecoPatCol_p_all$p_raw_001[which(ecoPatCol_p_all$p_raw <= .001)] <- 1

head(ecoPatCol_p_all)
ecoPatCol_p05_summ0 <- do.call(data.frame, aggregate(p_raw_05 ~ question + hypothesis + hypothesisType + response + predictor, data=ecoPatCol_p_all, function(x) propSig05_raw = sum(x)/length(x)))
ecoPatCol_p01_summ0 <- do.call(data.frame, aggregate(p_raw_01 ~ question + hypothesis + hypothesisType + response + predictor, data=ecoPatCol_p_all, function(x) propSig05_raw = sum(x)/length(x)))
ecoPatCol_p001_summ0 <- do.call(data.frame, aggregate(p_raw_001 ~ question + hypothesis + hypothesisType + response + predictor, data=ecoPatCol_p_all, function(x) propSig05_raw = sum(x)/length(x)))
ecoPatCol_summ_raw0 <- merge(ecoPatCol_p05_summ0, ecoPatCol_p01_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
ecoPatCol_summ_raw <- merge(ecoPatCol_summ_raw0, ecoPatCol_p001_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
head(ecoPatCol_summ_raw); nrow(ecoPatCol_summ_raw) # 540 total ecology ~ pattern & color comparisons 

# Adjusted p-values
ecoPatCol_p_all$p_fdr_05 <- 0
ecoPatCol_p_all$p_fdr_05[which(ecoPatCol_p_all$p_fdr <= .05)] <- 1
ecoPatCol_p_all$p_fdr_01 <- 0
ecoPatCol_p_all$p_fdr_01[which(ecoPatCol_p_all$p_fdr <= .01)] <- 1
ecoPatCol_p_all$p_fdr_001 <- 0
ecoPatCol_p_all$p_fdr_001[which(ecoPatCol_p_all$p_fdr <= .001)] <- 1

ecoPatCol_p05_summ0 <- do.call(data.frame, aggregate(p_fdr_05 ~ question + hypothesis + hypothesisType + response + predictor, data=ecoPatCol_p_all, function(x) propSig05_fdr = sum(x)/length(x)))
ecoPatCol_p01_summ0 <- do.call(data.frame, aggregate(p_fdr_01 ~ question + hypothesis + hypothesisType + response + predictor, data=ecoPatCol_p_all, function(x) propSig05_fdr = sum(x)/length(x)))
ecoPatCol_p001_summ0 <- do.call(data.frame, aggregate(p_fdr_001 ~ question + hypothesis + hypothesisType + response + predictor, data=ecoPatCol_p_all, function(x) propSig05_fdr = sum(x)/length(x))) 
ecoPatCol_summ_fdr0 <- merge(ecoPatCol_p05_summ0, ecoPatCol_p01_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
ecoPatCol_summ_fdr <- merge(ecoPatCol_summ_fdr0, ecoPatCol_p001_summ0, by=c("question", "hypothesis", "hypothesisType", "response", "predictor"))
head(ecoPatCol_summ_fdr); nrow(ecoPatCol_summ_fdr) # N=540 total ecology ~ pattern/color comparisons (83 aPriori + 457 exploratory)

# Add asterisks. All significance is based on the fdr-corrected values now, corrected at two different levels
ecoPatCol_summ_p0 <- merge(ecoPatCol_summ_raw, ecoPatCol_summ_fdr); head(ecoPatCol_summ_p0)
ecoPatCol_summ_p0$sig <- ""
ecoPatCol_summ_p0$sig[which(ecoPatCol_summ_p0$p_fdr_05 >= .9)] <- "*"
ecoPatCol_summ_p0$sig[which(ecoPatCol_summ_p0$p_fdr_01 >= .9)] <- "**"
ecoPatCol_summ_p0$sig[which(ecoPatCol_summ_p0$p_fdr_001 >= .9)] <- "***"

ecoPatCol_summ_p <- ecoPatCol_summ_p0[order(-ecoPatCol_summ_p0$p_fdr_05),]
head(ecoPatCol_summ_p); nrow(ecoPatCol_summ_p) # 540 total (84 a-priori + 456 exploratory)

################################
################################
# Merge effects + significance # Generates Table S4 & S5
################################
################################

##### pattern ~ color #####
head(patCol_summ_p); nrow(patCol_summ_p)
head(patCol_beta_summ); nrow(patCol_beta_summ)
patCol_beta_sig0 <- merge(patCol_summ_p, patCol_beta_summ, by=c("response", "predictor")); nrow(patCol_beta_sig0); head(patCol_beta_sig0)
patCol_beta_sig0$question <- NULL
patCol_beta_sig1 <- patCol_beta_sig0[c("hypothesis", "hypothesisType", "response", "predictor", "nBeta", "meanBeta", "sdBeta", "minBeta", "maxBeta", "oddsRatio", "p_raw_05", "p_fdr_05", "p_raw_01", "p_fdr_01", "p_raw_001", "p_fdr_001", "sig")]; head(patCol_beta_sig1)
patCol_beta_sig <- patCol_beta_sig1[order(-patCol_beta_sig1$p_fdr_05),]; head(patCol_beta_sig)

# Prettier version of the table for supplement
head(patCol_beta_sig)
patCol_beta_sig$meanBeta_SD <- paste(round(patCol_beta_sig$meanBeta, digits=2), " ", "(", round(patCol_beta_sig$sdBeta, digits=2), ")" , sep="")
patCol_beta_sig$betaRange <- paste(round(patCol_beta_sig$minBeta, digits=2), "/", round(patCol_beta_sig$maxBeta, digits=2), sep="")
patCol_beta_sig$oddsRatio_rounded <- round(patCol_beta_sig$oddsRatio, digits=2)
patCol_beta_sig$p_raw_05_perc <- paste(round(patCol_beta_sig$p_raw_05*100), "%", sep="")
patCol_beta_sig$p_raw_01_perc <- paste(round(patCol_beta_sig$p_raw_01*100), "%", sep="")
patCol_beta_sig$p_raw_001_perc <- paste(round(patCol_beta_sig$p_raw_001*100), "%", sep="")
patCol_beta_sig$p_fdr_05_perc <- paste(round(patCol_beta_sig$p_fdr_05*100), "%", sep="")
patCol_beta_sig$p_fdr_01_perc <- paste(round(patCol_beta_sig$p_fdr_01*100), "%", sep="")
patCol_beta_sig$p_fdr_001_perc <- paste(round(patCol_beta_sig$p_fdr_001*100), "%", sep="")
patCol_beta_sig$hypothesisTypeAbbr <- ifelse(patCol_beta_sig$hypothesisType == "aPriori", "aP", "Ex")
patCol_beta_sig$negPos <- ifelse(patCol_beta_sig$meanBeta < 0, "-", "+")
patCol_beta_sig$negPos <- factor(patCol_beta_sig$negPos, levels=c("+", "-"))
patCol_beta_sig_print0 <- subset(patCol_beta_sig, select=c("response", "predictor", "hypothesisTypeAbbr", "meanBeta_SD", "negPos", "betaRange", "oddsRatio_rounded", "p_raw_05_perc", "p_raw_01_perc", "p_raw_001_perc", "p_fdr_05_perc", "p_fdr_01_perc", "p_fdr_001_perc", "sig")); head(patCol_beta_sig_print0)
patCol_beta_sig_print0$response <- factor(patCol_beta_sig_print0$response, levels=c("band", "spot", "stripe", "blotch", "stippling"))
patCol_beta_sig_print0$sig <- factor(patCol_beta_sig_print0$sig, levels=c("***", "**", "*", ""))
patCol_beta_sig_print1 <- patCol_beta_sig_print0[order(patCol_beta_sig_print0$response, patCol_beta_sig_print0$hypothesisTypeAbbr,patCol_beta_sig_print0$negPos, patCol_beta_sig_print0$sig),] # sort first by response, and then into the a-priori and exploratory bins; within that, separate positive and negative, and finally sort by significance. I kind of like this; can work through each response, one by one.
patCol_beta_sig_print <- patCol_beta_sig_print1 %>%
  mutate(across(everything(), as.character)); head(patCol_beta_sig_print)
patCol_beta_sig_print$negPos <- NULL # remove this; was useful for ordering only
# See tutorial here for Word export process (https://sejdemyr.github.io/r-tutorials/basics/tables-in-r/)
# How many of each type of hypothesis?
nrow(subset(patCol_beta_sig_print, hypothesisTypeAbbr == "aP")) 
nrow(subset(patCol_beta_sig_print, hypothesisTypeAbbr == "Ex")) 

setwd("~/MSdata/outputs")
write.table(patCol_beta_sig_print,"patternColor-resultsSorted.txt", sep="," , quote = FALSE, row.names = F)


###### ecology ~ pattern & color ######
ecoPatCol_beta_sig0 <- merge(ecoPatCol_summ_p, ecoPatCol_beta_summ, by=c("response", "predictor"))
ecoPatCol_beta_sig0$question <- NULL
ecoPatCol_beta_sig1 <- ecoPatCol_beta_sig0[c("hypothesis", "hypothesisType", "response", "predictor", "nBeta", "meanBeta", "sdBeta", "minBeta", "maxBeta", "oddsRatio", "p_raw_05", "p_fdr_05", "p_raw_01", "p_fdr_01", "p_raw_001", "p_fdr_001", "sig")]
ecoPatCol_beta_sig <- ecoPatCol_beta_sig1[order(-ecoPatCol_beta_sig1$p_fdr_05),]

# Prettier version of the table for supplement
ecoPatCol_beta_sig$meanBeta_SD <- paste(round(ecoPatCol_beta_sig$meanBeta, digits=2), " ", "(", round(ecoPatCol_beta_sig$sdBeta, digits=2), ")" , sep="")
ecoPatCol_beta_sig$betaRange <- paste(round(ecoPatCol_beta_sig$minBeta, digits=2), "/", round(ecoPatCol_beta_sig$maxBeta, digits=2), sep="")
ecoPatCol_beta_sig$oddsRatio_rounded <- round(ecoPatCol_beta_sig$oddsRatio, digits=2)
ecoPatCol_beta_sig$p_raw_05_perc <- paste(round(ecoPatCol_beta_sig$p_raw_05*100), "%", sep="")
ecoPatCol_beta_sig$p_raw_01_perc <- paste(round(ecoPatCol_beta_sig$p_raw_01*100), "%", sep="")
ecoPatCol_beta_sig$p_raw_001_perc <- paste(round(ecoPatCol_beta_sig$p_raw_001*100), "%", sep="")
ecoPatCol_beta_sig$p_fdr_05_perc <- paste(round(ecoPatCol_beta_sig$p_fdr_05*100), "%", sep="")
ecoPatCol_beta_sig$p_fdr_01_perc <- paste(round(ecoPatCol_beta_sig$p_fdr_01*100), "%", sep="")
ecoPatCol_beta_sig$p_fdr_001_perc <- paste(round(ecoPatCol_beta_sig$p_fdr_001*100), "%", sep="")
ecoPatCol_beta_sig$hypothesisTypeAbbr <- ifelse(ecoPatCol_beta_sig$hypothesisType == "aPriori", "aP", "Ex")
ecoPatCol_beta_sig$negPos <- ifelse(ecoPatCol_beta_sig$meanBeta < 0, "-", "+")
ecoPatCol_beta_sig$negPos <- factor(ecoPatCol_beta_sig$negPos, levels=c("+", "-"))
ecoPatCol_beta_sig_print0 <- subset(ecoPatCol_beta_sig, select=c("response", "predictor", "hypothesisTypeAbbr", "meanBeta_SD", "negPos", "betaRange", "oddsRatio_rounded", "p_raw_05_perc", "p_raw_01_perc", "p_raw_001_perc", "p_fdr_05_perc", "p_fdr_01_perc", "p_fdr_001_perc", "sig")); head(ecoPatCol_beta_sig_print0)
ecoPatCol_beta_sig_print0$response <- factor(ecoPatCol_beta_sig_print0$response, levels=c("woody.shrub.tree", "forb", "grass","not.live.plants", "monophagous", "oligophagous", "polyphagous", "leaf", "reproductive.tissue", "interior"))
ecoPatCol_beta_sig_print0$sig <- factor(ecoPatCol_beta_sig_print0$sig, levels=c("***", "**", "*", ""))
ecoPatCol_beta_sig_print1 <- ecoPatCol_beta_sig_print0[order(ecoPatCol_beta_sig_print0$response, ecoPatCol_beta_sig_print0$hypothesisTypeAbbr,ecoPatCol_beta_sig_print0$negPos, ecoPatCol_beta_sig_print0$sig),] # sort first by response, and then into the a-priori and exploratory bins; within that, separate positive and negative, and finally sort by significance. I kind of like this; can work through each response, one by one.
ecoPatCol_beta_sig_print <- ecoPatCol_beta_sig_print1 %>%
  mutate(across(everything(), as.character)); head(ecoPatCol_beta_sig_print)
ecoPatCol_beta_sig_print$negPos <- NULL # remove this; was useful for ordering only
# See tutorial here for Word export process (https://sejdemyr.github.io/r-tutorials/basics/tables-in-r/)
# How many of each type of hypothesis?
nrow(subset(ecoPatCol_beta_sig_print, hypothesisTypeAbbr == "aP"))
nrow(subset(ecoPatCol_beta_sig_print, hypothesisTypeAbbr == "Ex"))

setwd("~/MSdata/outputs")
write.table(ecoPatCol_beta_sig_print,"ecologyPatternColor-resultsSorted.txt", sep="," , quote = FALSE, row.names = F)

