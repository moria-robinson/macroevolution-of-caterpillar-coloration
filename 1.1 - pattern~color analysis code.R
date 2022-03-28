
# 1.1 - Phylogenetic logistic regression: Associations between pattern (response) and colors (single colors and combinations)
 
# Moria L. Robinson

#----------#
# Packages #
#----------#
rm(list=ls()) # Clear everything else in the background
require(ape)
require(geiger)
require(splitstackshape)
require(phylolm)
require(reshape2)

#------- data --------#
data <- read.csv("~/MSdata/data.csv") # two rows for each species (two different observations of color, by two independent observers. Same (non-subjective) diet breadth data)
tree <- read.tree("~/MSdata/data/raxml.tree.analysis.labels_Feb2020.tree")

# Slice out diet breadth data, for use later in script
diet0 <- subset(data, select=c("genus_species", "monophagous", "oligophagous", "polyphagous", "not.live.plants", "grass", "forb", "shrub", "tree", "woody.shrub.tree","interior", "leaf", "flower", "fruit", "seed", "reproductive.tissue"))
diet <- diet0[!duplicated(diet0),] # nrow(diet) # 1808 spp

#################################################################
# TREE: Trim out a few subspecies & synonyms that were included #
#################################################################
# in all cases, removed synonym/subsp. tips were sister to full species; e.g. removal does not affect topology
all.species <- tree$tip.label
remove <- c("Catocala_badia_coelebs", "Oligia_illocata", "Cyclargus_ammon")
keep.species <- subset(all.species, !all.species %in% remove)
pruned.tree <- drop.tip(tree, tree$tip.label[-match(keep.species, tree$tip.label)])
str(pruned.tree) # 1808 spp

# RENAME three tips
pruned.tree$tip.label[pruned.tree$tip.label=="Drasteria_graphica_atlantica"] <- "Drasteria_graphica"
pruned.tree$tip.label[pruned.tree$tip.label=="Nymphalis_l-album"] <- "Nymphalis_I-album"
pruned.tree$tip.label[pruned.tree$tip.label=="Chlorostrymon_maesites"] <- "Chlorostrymon_telea" # all trait data also swapped for this congener

# Confirm no duplicate tips
all.species[duplicated(all.species)] 

# Use "pruned.tree" in analyses
# Remove temp files (will build up when iterating loop)
rm(tree, all.species, keep.species, remove, diet0)

###############################
# 1) Pattern ~ color analysis #
###############################
nrow(data) # 3616; two sets of 1808 species, each observed by 2 independent observers

# What this loop does: 
# - makes a bootstrapped dataset, selecting color data for each species from one of the two observers
# - Runs a phylogenetic logistic regression using the bootstrapped data
# - saves the output

# this is computationally intensive; recommend running in parallel

names(data)
# Subset to columns we'll be using - easier to check in on dataframe along the way 
temp0 <- subset(data, select=c("genus_species", "observer", "scan.color.1", "scan.color.1_alt", "scan.color.2", "scan.color.2_alt", "scan.color.3", "scan.color.3_alt", "stippling", "band", "stripe", "spot", "blotch"))

#### BOOTSTRAPPING DATA STARTS HERE

corr.p <- list()
corr.r <- list()
data <- list()

nreps <- 1000 # number of bootstrapped datasets and analyses to run

for (x in 1:nreps){
  setwd("~/MSdata/bootstraps/patternColor-outputs") # establish a destination folder for bootstrap outputs
  cat(file = "patternColor-runs.txt","i am on", x, append = T)
  tStart<-Sys.time()

  # Randomize across observers: randomly select one row per species; this effectively randomizes across the two observers 
  temp1 <- data.frame(stratified(temp0, "genus_species", 1)) 
  row.names(temp1) <- temp1$genus_species

  # Randomize across alternate colors
  # This is just for observations from MF, who was the only observer who did this (this was for colors that were between two categories - like yellowish-white. This observer made an alternative color, rather than choosing one of the two; we want to average over this)
  temp1_keep0 <- subset(temp1, observer != "MF") # other observers don't need this randomization
  temp1_keep1 <- subset(temp1_keep0, select=c("genus_species", "observer", "scan.color.1","scan.color.2", "scan.color.3", "stippling", "band", "stripe", "spot", "blotch")) 
  temp1_ran0 <- subset(temp1, observer == "MF") 
  temp1_ran0$index <- sample(1:2,replace=T, nrow(temp1_ran0)) # generate 1s and 2s randomly for each row
  ### Randomly choose between the main or alternate color, provided that there is an (not NA) alternate given:
  # scan color 1
  temp1_ran0$scan.color.1_final <- ifelse(temp1_ran0$index == 1 & !is.na(temp1_ran0$scan.color.1_alt), as.character(temp1_ran0$scan.color.1_alt), as.character(temp1_ran0$scan.color.1)) 
  # scan color 2
  temp1_ran0$index <- sample(1:2,replace=T, nrow(temp1_ran0)) # re-randomize the index
  temp1_ran0$scan.color.2_final <- ifelse(temp1_ran0$index == 1 & !is.na(temp1_ran0$scan.color.2_alt), as.character(temp1_ran0$scan.color.2_alt), as.character(temp1_ran0$scan.color.2)) 
  # scan color 3
  temp1_ran0$index <- sample(1:2,replace=T, nrow(temp1_ran0)) # re-randomize the index
  temp1_ran0$scan.color.3_final <- ifelse(temp1_ran0$index == 1 & !is.na(temp1_ran0$scan.color.3_alt), as.character(temp1_ran0$scan.color.3_alt), as.character(temp1_ran0$scan.color.3)) 
  # Check that it worked, re-order column to better compare
  temp1_ran1 <- temp1_ran0[,c("genus_species", "observer", "scan.color.1", "scan.color.1_alt", "scan.color.1_final", "scan.color.2", "scan.color.2_alt", "scan.color.2_final", "scan.color.3", "scan.color.3_alt", "scan.color.3_final", "stippling", "band", "stripe", "spot", "blotch")]# View(temp1_ran1)
  
  # Lastly - because of the randomization across alternative colors, there is a unique potential for the same color to appear twice. E.g. could have randomly selected an alternate that is the same as an already-entered scan color. 
  # Correct for this...
  temp1_ran1$scan.color.2_final2 <- ifelse(temp1_ran1$scan.color.2_final == temp1_ran1$scan.color.1_final, NA, temp1_ran1$scan.color.2_final)
  temp1_ran1$scan.color.3_final2 <- ifelse(temp1_ran1$scan.color.3_final == temp1_ran1$scan.color.1_final | temp1_ran1$scan.color.3_final == temp1_ran1$scan.color.2_final, NA, temp1_ran1$scan.color.3_final)

  # Subset and rename columns
  temp1_ran2 <- subset(temp1_ran1, select = c("genus_species", "observer","scan.color.1_final", "scan.color.2_final2", "scan.color.3_final2", "stippling", "band", "stripe", "spot", "blotch"))
  names(temp1_ran2) <- c("genus_species", "observer","scan.color.1", "scan.color.2", "scan.color.3", "stippling", "band", "stripe", "spot", "blotch")
  
  # Recombine the new randomized dataset with the non-randomized dataset (the other observer). 
  temp2 <- rbind(temp1_keep1, temp1_ran2) # This is a dataset that has been 1) randomized across observers and 2) randomized across alternative colors. 
  
  # Remove temp files - will build up when iterating loop
  rm(temp0, temp1, temp1_keep0, temp1_keep1, temp1_ran0, temp1_ran1, temp1_ran2)
  
  #----------------------------------------------------------------------------#
  # 4) wrangle all colors and color combinations into single-column 0/1 traits # 
  #----------------------------------------------------------------------------#
  # This is far from the most elegent approach, but it works

  # SINGLE COLORS: Make columns for each of the 11 colors... 
  temp2$black <- 0; temp2$gray <- 0; temp2$brown <- 0; temp2$white <- 0; temp2$red <- 0; temp2$pink <- 0; temp2$orange <- 0; temp2$yellow <- 0; temp2$green <- 0; temp2$blue <- 0; temp2$purple <- 0     
  # and now fill with "1" if the color is in any of the 3 scan color columns.
  temp2$black[which(temp2$scan.color.1 == "black" | temp2$scan.color.2 == "black" | temp2$scan.color.3 == "black")] <- 1
  temp2$gray[which(temp2$scan.color.1 == "gray" | temp2$scan.color.2 == "gray" | temp2$scan.color.3 == "gray")] <- 1
  temp2$brown[which(temp2$scan.color.1 == "brown" | temp2$scan.color.2 == "brown" | temp2$scan.color.3 == "brown")] <- 1
  temp2$white[which(temp2$scan.color.1 == "white" | temp2$scan.color.2 == "white" | temp2$scan.color.3 == "white")] <- 1
  temp2$red[which(temp2$scan.color.1 == "red" | temp2$scan.color.2 == "red" | temp2$scan.color.3 == "red")] <- 1
  temp2$pink[which(temp2$scan.color.1 == "pink" | temp2$scan.color.2 == "pink" | temp2$scan.color.3 == "pink")] <- 1
  temp2$orange[which(temp2$scan.color.1 == "orange" | temp2$scan.color.2 == "orange" | temp2$scan.color.3 == "orange")] <- 1
  temp2$yellow[which(temp2$scan.color.1 == "yellow" | temp2$scan.color.2 == "yellow" | temp2$scan.color.3 == "yellow")] <- 1
  temp2$green[which(temp2$scan.color.1 == "green" | temp2$scan.color.2 == "green" | temp2$scan.color.3 == "green")] <- 1
  temp2$blue[which(temp2$scan.color.1 == "blue" | temp2$scan.color.2 == "blue" | temp2$scan.color.3 == "blue")] <- 1
  temp2$purple[which(temp2$scan.color.1 == "purple" | temp2$scan.color.2 == "purple" | temp2$scan.color.3 == "purple")] <- 1
  
  # TWO-COLOR COMBINATIONS - clunky and long code, but it works
  # add ALL two- color combos PRESENT in each caterpillar body. So, for spp with three colors, this will tally each of the 3 2-color combos also present in the body. These are all possible 2-color combos.
  
  # Brown + others
  temp2$brown.green <- ifelse(temp2$brown == 1 & temp2$green == 1, 1, 0)
  temp2$brown.white <- ifelse(temp2$brown == 1 & temp2$white == 1, 1, 0)
  temp2$brown.red <- ifelse(temp2$brown == 1 & temp2$red == 1, 1, 0)
  temp2$brown.orange <- ifelse(temp2$brown == 1 & temp2$orange == 1, 1, 0)
  temp2$brown.yellow <- ifelse(temp2$brown == 1 & temp2$yellow == 1, 1, 0)
  temp2$brown.black <- ifelse(temp2$brown == 1 & temp2$black == 1, 1, 0)
  temp2$brown.gray <- ifelse(temp2$brown == 1 & temp2$gray == 1, 1, 0)
  temp2$brown.pink <- ifelse(temp2$brown == 1 & temp2$pink == 1, 1, 0)
  temp2$brown.blue <- ifelse(temp2$brown == 1 & temp2$blue == 1, 1, 0)
  temp2$brown.purple <- ifelse(temp2$brown == 1 & temp2$purple == 1, 1, 0)
  # Green + others
  temp2$green.white <- ifelse(temp2$green == 1 & temp2$white == 1, 1, 0)
  temp2$green.red <- ifelse(temp2$green == 1 & temp2$red == 1, 1, 0)
  temp2$green.orange <- ifelse(temp2$green == 1 & temp2$orange == 1, 1, 0)
  temp2$green.yellow <- ifelse(temp2$green == 1 & temp2$yellow == 1, 1, 0)
  temp2$green.black <- ifelse(temp2$green == 1 & temp2$black == 1, 1, 0)
  temp2$green.gray <- ifelse(temp2$green == 1 & temp2$gray == 1, 1, 0)
  temp2$green.pink <- ifelse(temp2$green == 1 & temp2$pink == 1, 1, 0)
  temp2$green.blue <- ifelse(temp2$green == 1 & temp2$blue == 1, 1, 0)
  temp2$green.purple <- ifelse(temp2$green == 1 & temp2$purple == 1, 1, 0)
  # White + others
  temp2$white.red <- ifelse(temp2$white == 1 & temp2$red == 1, 1, 0)
  temp2$white.orange <- ifelse(temp2$white == 1 & temp2$orange == 1, 1, 0)
  temp2$white.yellow <- ifelse(temp2$white == 1 & temp2$yellow == 1, 1, 0)
  temp2$white.black <- ifelse(temp2$white == 1 & temp2$black == 1, 1, 0)
  temp2$white.gray <- ifelse(temp2$white == 1 & temp2$gray == 1, 1, 0)
  temp2$white.pink <- ifelse(temp2$white == 1 & temp2$pink == 1, 1, 0)
  temp2$white.blue <- ifelse(temp2$white == 1 & temp2$blue == 1, 1, 0)
  temp2$white.purple <- ifelse(temp2$white == 1 & temp2$purple == 1, 1, 0)
  # Red + others
  temp2$red.orange <- ifelse(temp2$red == 1 & temp2$orange == 1, 1, 0)
  temp2$red.yellow <- ifelse(temp2$red == 1 & temp2$yellow == 1, 1, 0)
  temp2$red.black <- ifelse(temp2$red == 1 & temp2$black == 1, 1, 0)
  temp2$red.gray <- ifelse(temp2$red == 1 & temp2$gray == 1, 1, 0)
  temp2$red.pink <- ifelse(temp2$red == 1 & temp2$pink == 1, 1, 0)
  temp2$red.blue <- ifelse(temp2$red == 1 & temp2$blue == 1, 1, 0)
  temp2$red.purple <- ifelse(temp2$red == 1 & temp2$purple == 1, 1, 0)
  # Orange + others
  temp2$orange.yellow <- ifelse(temp2$orange == 1 & temp2$yellow == 1, 1, 0)
  temp2$orange.black <- ifelse(temp2$orange == 1 & temp2$black == 1, 1, 0)
  temp2$orange.gray <- ifelse(temp2$orange == 1 & temp2$gray == 1, 1, 0)
  temp2$orange.pink <- ifelse(temp2$orange == 1 & temp2$pink == 1, 1, 0)
  temp2$orange.blue <- ifelse(temp2$orange == 1 & temp2$blue == 1, 1, 0)
  temp2$orange.purple <- ifelse(temp2$orange == 1 & temp2$purple == 1, 1, 0)
  # Yellow + others
  temp2$yellow.black <- ifelse(temp2$yellow == 1 & temp2$black == 1, 1, 0)
  temp2$yellow.gray <- ifelse(temp2$yellow == 1 & temp2$gray == 1, 1, 0)
  temp2$yellow.pink <- ifelse(temp2$yellow == 1 & temp2$pink == 1, 1, 0)
  temp2$yellow.blue <- ifelse(temp2$yellow == 1 & temp2$blue == 1, 1, 0)
  temp2$yellow.purple <- ifelse(temp2$yellow == 1 & temp2$purple == 1, 1, 0)
  # Black + others
  temp2$black.gray <- ifelse(temp2$black == 1 & temp2$gray == 1, 1, 0)
  temp2$black.pink <- ifelse(temp2$black == 1 & temp2$pink == 1, 1, 0)
  temp2$black.blue <- ifelse(temp2$black == 1 & temp2$blue == 1, 1, 0)
  temp2$black.purple <- ifelse(temp2$black == 1 & temp2$purple == 1, 1, 0)
  # Gray + others
  temp2$gray.pink <- ifelse(temp2$gray == 1 & temp2$pink == 1, 1, 0)
  temp2$gray.blue <- ifelse(temp2$gray == 1 & temp2$blue == 1, 1, 0)
  temp2$gray.purple <- ifelse(temp2$gray == 1 & temp2$purple == 1, 1, 0)
  # Pink + others
  temp2$pink.blue <- ifelse(temp2$pink == 1 & temp2$blue == 1, 1, 0)
  temp2$pink.purple <- ifelse(temp2$pink == 1 & temp2$purple == 1, 1, 0)
  # Blue + others
  temp2$blue.purple <- ifelse(temp2$blue == 1 & temp2$purple == 1, 1, 0)

  # Add a column for the # of colors on each caterpillar (numCol)
  temp2$numCol<- temp2$black +  temp2$gray + temp2$brown + temp2$white + temp2$red + temp2$pink + temp2$orange + temp2$yellow + temp2$green + temp2$blue + temp2$purple
  
  #-------------------------#
  # 5) Merge with diet data #
  #-------------------------#
  temp3 <- merge(temp2, diet, by="genus_species") 
  
  #----------------------#
  # 6) Prep for phyloGLM #
  #----------------------#
  # Make sure levels & structure are OK; all traits of interest in numeric form; subset only to relevant columns
  # remove traits present in < 3 spp 
  # Assign species as row names to the data
  
  # Remove columns that will not be in analysis
  temp4 <- temp3[ , -which(names(temp3) %in% c("X", "observer", "scan.color.1", "scan.color.2","scan.color.3", "tree", "shrub", "flower", "fruit", "seed"))]
  
  # Make all columns (except genus_species) numeric
  temp5 <- temp4
  temp5[,-which(names(temp5) %in% c("genus_species"))] <- sapply(temp5[, -which(names(temp5) %in% c("genus_species"))], as.numeric) # str(temp5)
  
  # Drop traits that occur 3 or fewer (makes the loop crash)
  frequencies <- data.frame(colSums(temp5[,-which(names(temp5) %in% c("genus_species"))])); frequencies$trait <- rownames(frequencies); names(frequencies) <- c("count", "trait"); frequencies2 <- frequencies[order(frequencies$count),] # frequencies2
  keep <- subset(frequencies2, count >= 3)$trait # keep
  temp6 <- temp5[names(temp5)[names(temp5) %in% c("genus_species", keep)]] 
  
  # Establish which columns are the (1) predictors (color); (2) response (pattern); and (3) covariate/factor (numCol). 
  pattern <- c("band","stripe","blotch","stippling","spot")
  ecology <- c("monophagous", "oligophagous", "polyphagous", "woody.shrub.tree","forb","grass", "leaf", "reproductive.tissue", "interior", "not.live.plants")
  
  resp0 <- temp6[,which(names(temp6) %in% c("genus_species", pattern))] 
  pred0 <- temp6[,-which(names(temp6) %in% c(ecology, pattern, "numCol"))] 
  colorCount0 <- temp6[,which(names(temp6) %in% c("genus_species", "numCol"))]
  
  # row.names based on genus_species for each 
  row.names(resp0) <- resp0$genus_species
  row.names(pred0) <- pred0$genus_species
  row.names(colorCount0) <- colorCount0$genus_species
  
  # Remove genus_species
  resp <- resp0[,-which(names(resp0) %in% c("genus_species"))] 
  pred <- pred0[,-which(names(pred0) %in% c("genus_species"))] 
  colorCount <- colorCount0; colorCount$genus_species <- NULL 
  
  # Make sure all will map to tree
  #name.check(pruned.tree, resp) 
  #name.check(pruned.tree, pred) 
  #name.check(pruned.tree, colorCount) 
  
  # Remove temporary files; will build up thru loop iterations
  rm(colorCount0, frequencies, frequencies2,  pred0, resp0, temp2, temp3, temp4, temp5, keep, ecology, pattern)
  
  #---------------#
  # 7) RUN MODEL  #
  #---------------#
  predTraits <-colnames(pred)
  respTraits <-colnames(resp) 
  # REMEMBER:column = predictor(j, trait2) and row = response (i, trait 1)!
  phyloGLM.r<-matrix(NA,nrow=length(respTraits),ncol=length(predTraits)) 
  row.names(phyloGLM.r)<-respTraits
  colnames(phyloGLM.r)<-predTraits # Makes an empty matrix to fill with correl coefs
  phyloGLM.p<-phyloGLM.r # Replicates this empty matrix to fill with p-values
  
  for (i in 1:length(respTraits)){
    for (j in 1:length(predTraits)){
      if(respTraits[i]==predTraits[j]){phyloGLM.r[i,j]<-1}
      else{
        trait1<-resp[,i] # for a particular column (i).... treat (i) as trait 1 (will be the response)
        trait2<-pred[,j] # and a particular column (j).... treat (j) as trait 2 (will be the predictor)
        numCol <- colorCount$numCol
        names(trait1)<-row.names(resp)
        names(trait2)<-row.names(pred)
        names(numCol)<-row.names(colorCount)
        fit.temp<-phyloglm(trait1 ~ trait2 + as.numeric(numCol), phy=pruned.tree, method="logistic_IG10")
        summary.temp <- summary(fit.temp)
        phyloGLM.r[i,j]<-summary.temp$coefficients[2]
        phyloGLM.p[i,j]<-summary.temp$coefficients[11] 
        rm(fit.temp,summary.temp,trait1,trait2)
      }
    }
  }
  
  
  # convert to long-form dataframes; save the bootstrap iteration for each
  # correlation coefficients
  phyloGLM.r2 <- data.frame(phyloGLM.r)
  phyloGLM.r2$response <- rownames(phyloGLM.r2)
  phyloGLM.r3 <- reshape2::melt(phyloGLM.r2, id.vars="response", variable.name = "predictor", value.name="r")
  # correlation coefficients
  phyloGLM.p2 <- data.frame(phyloGLM.p)
  phyloGLM.p2$response <- rownames(phyloGLM.p2)
  phyloGLM.p3 <- reshape2::melt(phyloGLM.p2, id.vars="response", variable.name = "predictor", value.name="p_raw")
  # Add the bootstrap identifier as a column
  phyloGLM.r3$bootstrap <- x
  phyloGLM.p3$bootstrap <- x
  temp6$bootstrap <- x
  
  # Save each bootstrap (destination established at beginning of loop)
  write.csv(phyloGLM.p3, file = paste0("patternColor_phyloGLM.p-",x,".csv"))
  write.csv(phyloGLM.r3, file = paste0("patternColor_phyloGLM.r-",x,".csv"))
  write.csv(temp9, file = paste0("patternColor_dataBootstraps-",x,".csv"))
  
  # Remove all the objects
  rm(temp9,phyloGLM.r,phyloGLM.r1,phyloGLM.r2,phyloGLM.r3,phyloGLM.p,phyloGLM.p1,phyloGLM.p2,phyloGLM.p3)
  
  # Add time elapsed to the .txt file
  cat(file = "patternColor-runs.txt", ", hey that took ",print(Sys.time()-tStart), "mins" ,"\r", append = T)
}

# Next scripts: summarize results and statistics (mean, max, min of r) across bootstraps. For significance, penalize for multiple tests within a-priori and exploratory subsets; calculate % sig Ps for each association, across bootstraps