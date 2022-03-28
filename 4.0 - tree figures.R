####################
# Make tree figure #
####################

# The tree figure in the paper is a composite: it is a radial tree without any symbols, and then 'call outs' are clades subset manually (photoshop) from a vertical tree with symbols mapped to the tips. This script generates both trees. 

# Add pie charts, pattern symbols onto the tree
# https://guangchuangyu.github.io/software/ggtree/vignettes/ggtree-ggimage.html
# https://yulab-smu.top/treedata-book/chapter8.html
# plant silhouettes from http://phylopic.org/

##### Packages ######
rm(list=ls())
library(ape)
library(ggtree)
library(ggplot2)
library(ggimage)
library(stringr)
library(dplyr)
library(tidyr)

##### Data #####
tree <- read.tree("~/MSdata/data/raxml.tree.analysis.labels_Feb2020.tree")
traits <- read.csv("~/MSdata/figures/traits.csv") # this is one of the bootstrapped datasets, randomly-chosen (first one)

# APPROACH
# For each trait category, loop through the species and list the traits that are present / we want to see
# Merge them together with spacers in-between (so that, for each species, we'll see their colors, then a spacer and then the patterns, then a spacer and then the growth forms, etc)

# COLORS
sp <- unique(traits$genus_species)
nSp <- length(unique(traits$genus_species)) # 1808 spp
colors.out <- list()

for(i in 1:nSp){
  temp <- subset(traits, genus_species == sp[i])
#  temp <- subset(traits, genus_species == "Papilio_rutulus")
  temp.colors <- subset(temp, select=c("green", "yellow", "white", "brown", "gray", "black", "pink", "red", "orange", "blue", "purple"))
  colors <- c("green", "yellow", "white", "brown", "gray", "black", "pink", "red", "orange", "blue", "purple")
  colors.present <- colors[colSums(temp.colors[colors]) != 0]
  colors.present_ordered0 <- sort(colors.present)
  colors.present_ordered <- c(colors.present_ordered0, "spacer")
  numColors <- length(colors.present)
  colors_fileNames0 <- paste0(colors.present_ordered, ".png")
  colors_fileNames <- paste(colors_fileNames0, collapse=',')
  #colors_fileNames <- paste0('"', paste(colors_fileNames0, collapse='", "'), '"')
  genus_species <- temp$genus_species[1]
  colors.summary <- data.frame(genus_species, numColors, colors_fileNames)
  colors.out[[i]] <- colors.summary
  rm(temp, temp.colors, colors.present, colors, colors.present_ordered0, colors.present_ordered, numColors, colors_fileNames0, colors_fileNames, genus_species, colors.summary)
}

colors.summary <- do.call(rbind, colors.out); head(colors.summary)
# Make any spp without any colors (shouldn't be any cases like this) NA in this column
colors.summary$colors_fileNames[which(colors.summary$colors_fileNames == "spacer.png")] <- NA; head(colors.summary)

# PATTERNS
sp <- unique(traits$genus_species)
nSp <- length(unique(traits$genus_species))
patterns.out <- list()

for(i in 1:nSp){
  temp <- subset(traits, genus_species == sp[i])
#  temp <- subset(traits, genus_species == "Papilio_rutulus")
  temp.patterns <- subset(temp, select=c("stripe", "blotch", "stippling", "band", "spot"))
  patterns <- c("stripe", "blotch", "stippling", "band", "spot")
  patterns.present <- patterns[colSums(temp.patterns[patterns]) != 0]
  patterns.present_ordered0 <- as.character(factor(patterns.present, levels=patterns))
  patterns.present_ordered <- c(patterns.present_ordered0, "spacer")
  numPatterns <- length(patterns.present)
  patterns_fileNames0 <- paste0(patterns.present_ordered, ".png")
  patterns_fileNames <- paste(patterns_fileNames0, collapse=',')
  #patterns_fileNames <- paste0('"', paste(patterns_fileNames0, collapse='", "'), '"')
  genus_species <- temp$genus_species[1]
  patterns.summary <- data.frame(genus_species, numPatterns, patterns_fileNames)
  patterns.out[[i]] <- patterns.summary
  rm(temp, temp.patterns, patterns.present, patterns, patterns.present_ordered0, patterns.present_ordered, numPatterns, patterns_fileNames0, patterns_fileNames, genus_species, patterns.summary)
}

patterns.summary <- do.call(rbind, patterns.out); head(patterns.summary)
# Make any spp without any patterns NA in this column
patterns.summary$patterns_fileNames[which(patterns.summary$patterns_fileNames == "spacer.png")] <- NA; head(patterns.summary)


# DIET BREADTH
sp <- unique(traits$genus_species)
nSp <- length(unique(traits$genus_species))
dietBreadth.out <- list()

for(i in 1:nSp){
  temp <- subset(traits, genus_species == sp[i])
#  temp <- subset(traits, genus_species == "Papilio_rutulus")
  temp.diet <- subset(temp, select=c("monophagous", "oligophagous", "polyphagous"))
  diets <- c("monophagous", "oligophagous", "polyphagous")
  diet.present0 <- diets[colSums(temp.diet[diets]) != 0]
  diet.present <- c(diet.present0, "spacer")
  numDiets <- length(diet.present)
  diet_fileNames0 <- as.character(paste0(diet.present, ".png"))
  diet_fileNames <- paste(diet_fileNames0, collapse=',')
  #diet_fileNames <- paste0('"', paste(diet_fileNames0, collapse='", "'), '"')
  genus_species <- temp$genus_species[1]
  diet.summary <- data.frame(genus_species, numDiets, diet_fileNames)
  dietBreadth.out[[i]] <- diet.summary
  rm(temp, temp.diet, diets, diet.present, numDiets, diet_fileNames, genus_species, diet.summary)
}

dietBreadth.summary <- do.call(rbind, dietBreadth.out); head(dietBreadth.summary)
dietBreadth.summary$diet_fileNames[which(dietBreadth.summary$diet_fileNames == "spacer.png")] <- NA; head(dietBreadth.summary)

# PLANT GROWTH FORM
sp <- unique(traits$genus_species)
nSp <- length(unique(traits$genus_species))
growthForms.out <- list()

for(i in 1:nSp){
  temp <- subset(traits, genus_species == sp[i])
#  temp <- subset(traits, genus_species == "Papilio_rutulus")
  temp.growthForms <- subset(temp, select=c("grass", "forb", "woody.shrub.tree", "not.live.plants"))
  growthForms <- c("grass", "forb", "woody.shrub.tree", "not.live.plants")
  growthForms.present <- growthForms[colSums(temp.growthForms[growthForms]) != 0]
  growthForms.present_ordered0 <- as.character(factor(growthForms.present, levels=growthForms))
  growthForms.present_ordered <- c(growthForms.present_ordered0, "spacer")
  numGrowthForms <- length(growthForms.present)
  growthForms_fileNames0 <- paste0(growthForms.present_ordered, ".png")
  growthForms_fileNames <- paste(growthForms_fileNames0, collapse=',')
  #growthForms_fileNames <- paste0('"', paste(growthForms_fileNames0, collapse='", "'), '"')
  genus_species <- temp$genus_species[1]
  growthForms.summary <- data.frame(genus_species, numGrowthForms, growthForms_fileNames)
  growthForms.out[[i]] <- growthForms.summary
  rm(temp, temp.growthForms, growthForms.present, growthForms, growthForms.present_ordered0, growthForms.present_ordered, numGrowthForms, growthForms_fileNames0, growthForms_fileNames, genus_species, growthForms.summary)
}

growthForm.summary <- do.call(rbind, growthForms.out); head(growthForm.summary)
growthForm.summary$growthForms_fileNames[which(growthForm.summary$growthForms_fileNames == "spacer.png")] <- NA; head(growthForm.summary)


# PRIMARY TISSUE CONSUMED
sp <- unique(traits$genus_species)
nSp <- length(unique(traits$genus_species))
tissue.out <- list()

for(i in 1:nSp){
  temp <- subset(traits, genus_species == sp[i])
#  temp <- subset(traits, genus_species == "Papilio_rutulus")
  temp.tissue <- subset(temp, select=c("leaf", "reproductive.tissue", "interior"))
  tissues <- c("leaf", "reproductive.tissue", "interior")
  tissue.present0 <- tissues[colSums(temp.tissue[tissues]) != 0]
  tissue.present <- c(tissue.present0, "spacer")
  numTissues <- length(tissue.present)
  tissue_fileNames0 <- paste0(tissue.present, ".png")
  tissue_fileNames <- paste(tissue_fileNames0, collapse=',')
  #tissue_fileNames <- paste0('"', paste(tissue_fileNames0, collapse='", "'), '"')
  genus_species <- temp$genus_species[1]
  tissue.summary <- data.frame(genus_species, numTissues, tissue_fileNames)
  tissue.out[[i]] <- tissue.summary
  rm(temp, temp.tissue, tissues, tissue.present,tissue.present0, numTissues, tissue_fileNames, genus_species, tissue.summary)
}

tissue.summary <- do.call(rbind, tissue.out); head(tissue.summary)
tissue.summary$tissue_fileNames[which(tissue.summary$tissue_fileNames == "spacer.png")] <- NA; head(tissue.summary)


##############
allTraits1 <- merge(colors.summary, patterns.summary, by="genus_species"); 
nrow(colors.summary) == nrow(patterns.summary); nrow(patterns.summary) == nrow(allTraits1)
allTraits2 <- merge(allTraits1, dietBreadth.summary, by="genus_species"); 
nrow(dietBreadth.summary) == nrow(allTraits2)
allTraits3 <- merge(allTraits2, growthForm.summary, by="genus_species"); 
nrow(growthForm.summary) == nrow(allTraits3)
allTraits4 <- merge(allTraits3, tissue.summary, by="genus_species"); 
nrow(tissue.summary) == nrow(allTraits4)
head(allTraits4)

# Now, I want to paste all the columns together (ignoring NA columns) to make a single vector of all the traits, in order, with spacers. Can tweak order in which trait categories appear here

allTraits4$pastedTogether1 <- paste(allTraits4$colors_fileNames, allTraits4$patterns_fileNames, allTraits4$diet_fileNames, allTraits4$growthForms_fileNames, allTraits4$tissue_fileNames, sep=",")
#View(subset(allTraits4, select=c("pastedTogether1")))
allTraits4$pastedTogetherFinal <- gsub(",NA","", allTraits4$pastedTogether1)
#View(subset(allTraits4, select=c("pastedTogetherFinal")))

### OK. Now, need to split the contents of this column up into new columns...
allTraits5 <- subset(allTraits4, select=c("genus_species", "pastedTogetherFinal"))
#View(allTraits5)

# Maximum number of columns is 3 colors + 5 patterns + 1 diet breadth + 4 growth forms + 1 tissue + 5 spacers; 16 symbols
allTraits6 <- separate(data = allTraits5, col = pastedTogetherFinal, into=c("symbol1", "symbol2", "symbol3", "symbol4", "symbol5", "symbol6", "symbol7", "symbol8", "symbol9", "symbol10", "symbol11", "symbol12", "symbol13", "symbol14", "symbol15", "symbol16"), sep = ",")
head(allTraits6) # sweet!

# Now, make sure row order matches tip label order of tree...
spOrder <- tree$tip.label
allTraits6$genus_species <- factor(allTraits6$genus_species, levels=spOrder)
allTraits7 <- allTraits6[order(allTraits6$genus_species),]

# OK - now, the final step will be to concatenate down the columns, in order, for each symbol position. Here, we will KEEP the NAs in any species rows/ within the columns, since we want them to map on as blanks on those tips. Then, we need to remember to  add the required number of NAs at the end of the trait vector (1 less than the number of tips)
NAs <- rep(NA, nrow(allTraits7)-1)

symbol1 <- c(allTraits7$symbol1, NAs); length(symbol1) # 3615
symbol2 <- c(allTraits7$symbol2, NAs); length(symbol2) # 3615
symbol3 <- c(allTraits7$symbol3, NAs); length(symbol3) # 3615
symbol4 <- c(allTraits7$symbol4, NAs); length(symbol4) # 3615
symbol5 <- c(allTraits7$symbol5, NAs); length(symbol5) # 3615
symbol6 <- c(allTraits7$symbol6, NAs); length(symbol6) # 3615
symbol7 <- c(allTraits7$symbol7, NAs); length(symbol7) # 3615
symbol8 <- c(allTraits7$symbol8, NAs); length(symbol8) # 3615
symbol9 <- c(allTraits7$symbol9, NAs); length(symbol9) # 3615
symbol10 <- c(allTraits7$symbol10, NAs); length(symbol10) # 3615
symbol11 <- c(allTraits7$symbol11, NAs); length(symbol11) # 3615
symbol12 <- c(allTraits7$symbol12, NAs); length(symbol12) # 3615
symbol13 <- c(allTraits7$symbol13, NAs); length(symbol13) # 3615
symbol14 <- c(allTraits7$symbol14, NAs); length(symbol14) # 3615
symbol15 <- c(allTraits7$symbol15, NAs); length(symbol15) # 3615
symbol16 <- c(allTraits7$symbol16, NAs); length(symbol16) # 3615

setwd("~/MSdata/figures/symbols")

# Make the ginormous tree!
ggtree(tree, size=.25) + 
  geom_image(image=symbol1, nudge_x=.0165, size=.0004) +
  geom_image(image=symbol2, nudge_x=.0175, size=.0004) +
  geom_image(image=symbol3, nudge_x=.0185, size=.0004) +
  geom_image(image=symbol4, nudge_x=.0195, size=.0004) +
  geom_image(image=symbol5, nudge_x=.0205, size=.0004) +
  geom_image(image=symbol6, nudge_x=.0215, size=.0004) +
  geom_image(image=symbol7, nudge_x=.0225, size=.0004) +
  geom_image(image=symbol8, nudge_x=.0235, size=.0004) +
  geom_image(image=symbol9, nudge_x=.0245, size=.0004) +
  geom_image(image=symbol10, nudge_x=.0255, size=.0004) +
  geom_image(image=symbol11, nudge_x=.0265, size=.0004) +
  geom_image(image=symbol12, nudge_x=.0275, size=.0004) +
  geom_image(image=symbol13, nudge_x=.0285, size=.0004) +
  geom_image(image=symbol14, nudge_x=.0295, size=.0004) +
  geom_image(image=symbol15, nudge_x=.0305, size=.0004) +
  geom_image(image=symbol16, nudge_x=.0315, size=.0004) +
  geom_tiplab(size=1, offset=-.004, hjust=0, color="black", align=FALSE) +
  ggplot2::xlim(0,2)
# PDF 100 (height) x 100 (width)

#### Another way to set up the symbols would be to space them out, so that there was some visual separation between color, pattern, diet breadth *THIS IS THE FIGURE IN THE PAPER

# I think to do this we need to add blanks for each NA, plus a separator/spacer between each type of trait
# Color traits will occupy symbol positions 1,2,3,4 (4=separator)
colors.simp0 <- colors.summary; head(colors.simp0)
colors.simp0$colors_fileNames2 <- str_replace(colors.simp0$colors_fileNames, ",spacer.png", ""); head(colors.simp0)
colors.simp0$color.blanks <- ifelse(colors.simp0$numColors == 1, "blank.png,blank.png", "blank.png")
colors.simp0$color.blanks[which(colors.simp0$numColors == 3)] <- ""; colors.simp0[1:30,]
colors.simp0$colorsAndBlanks <- paste(colors.simp0$colors_fileNames2, colors.simp0$color.blanks, sep=","); colors.simp0[1:30,]
colors.simp0$colorsAndBlanks2 <- gsub(",$","",colors.simp0$colorsAndBlanks); colors.simp0[1:30,]
colors.simp2 <- separate(data = colors.simp0, col = colorsAndBlanks2, into=c("symbol1", "symbol2", "symbol3"), sep = ","); colors.simp2[1:30,]
colors.simp <- subset(colors.simp2, select=c("genus_species", "symbol1", "symbol2", "symbol3")); colors.simp[1:30,]
colors.simp$symbol4 <- "separator.png"; colors.simp[1:30,]

# Next up is pattern... pattern traits will occupy symbol positions 5,6,7,8,9,10 (10=separator)
head(patterns.summary); range(patterns.summary$numPatterns) # some spp have all the patterns.
patterns.simp0 <- patterns.summary; head(patterns.simp0)
patterns.simp0$patterns_fileNames2 <- str_replace(patterns.simp0$patterns_fileNames, ",spacer.png", ""); head(patterns.simp0)
patterns.simp0$patterns.blanks <- ""
patterns.simp0$patterns.blanks[which(patterns.simp0$numPatterns == 0)] <- "blank.png,blank.png,blank.png,blank.png,blank.png"
patterns.simp0$patterns.blanks[which(patterns.simp0$numPatterns == 1)] <- "blank.png,blank.png,blank.png,blank.png"
patterns.simp0$patterns.blanks[which(patterns.simp0$numPatterns == 2)] <- "blank.png,blank.png,blank.png"
patterns.simp0$patterns.blanks[which(patterns.simp0$numPatterns == 3)] <- "blank.png,blank.png"
patterns.simp0$patterns.blanks[which(patterns.simp0$numPatterns == 4)] <- "blank.png"
patterns.simp0$patterns.blanks[which(patterns.simp0$numPatterns == 5)] <- ""; patterns.simp0[1:30,]
patterns.simp0$patternsAndBlanks <- paste(patterns.simp0$patterns_fileNames2, patterns.simp0$patterns.blanks, sep=","); patterns.simp0[1:30,]
patterns.simp0$patternsAndBlanks2 <- gsub(",$","",patterns.simp0$patternsAndBlanks); patterns.simp0[1:30,]
patterns.simp0$patternsAndBlanks2[which(patterns.simp0$numPatterns == 0)] <- "blank.png,blank.png,blank.png,blank.png,blank.png"
patterns.simp2 <- separate(data = patterns.simp0, col = patternsAndBlanks2, into=c("symbol5", "symbol6", "symbol7", "symbol8", "symbol9"), sep = ","); patterns.simp2[1:30,]
patterns.simp <- subset(patterns.simp2, select=c("genus_species", "symbol5", "symbol6", "symbol7", "symbol8", "symbol9")); patterns.simp[1:30,]
patterns.simp$symbol10 <- "separator.png"; patterns.simp[1:30,]

# Next up is diet breadth - occupy symbol positions 11, 12 (12=separator)
head(dietBreadth.summary)
dietBreadth.simp0 <- dietBreadth.summary; head(dietBreadth.simp0)
dietBreadth.simp0$dietBreadth_fileNames2 <- str_replace(dietBreadth.simp0$diet_fileNames, ",spacer.png", ""); head(dietBreadth.simp0)
dietBreadth.simp0$dietBreadth_fileNames3 <- ifelse(dietBreadth.simp0$numDiets == 1, "blank.png", dietBreadth.simp0$dietBreadth_fileNames2); head(dietBreadth.simp0)
dietBreadth.simp <- subset(dietBreadth.simp0, select=c("genus_species", "dietBreadth_fileNames3")); names(dietBreadth.simp) <- c("genus_species", "symbol11"); dietBreadth.simp[1:30,]
dietBreadth.simp$symbol12 <- "separator.png"; dietBreadth.simp[1:30,]

# Next up is growth form traits. Occupy symbol positions 13,14,15,16 (16=separator)
head(growthForm.summary); range(growthForm.summary$numGrowthForms) # ranges from 1 to 3. 
growthForm.simp0 <- growthForm.summary; head(growthForm.simp0)
growthForm.simp0$growthForms_fileNames2 <- str_replace(growthForm.simp0$growthForms_fileNames, ",spacer.png", ""); growthForm.simp0[1:30,]
growthForm.simp0$growthForms.blanks <- ""
growthForm.simp0$growthForms.blanks[which(growthForm.simp0$numGrowthForms == 1)] <- "blank.png,blank.png"
growthForm.simp0$growthForms.blanks[which(growthForm.simp0$numGrowthForms == 2)] <- "blank.png"
growthForm.simp0$growthForms.blanks[which(growthForm.simp0$numGrowthForms == 3)] <- "" ; growthForm.simp0[1:30,]
growthForm.simp0$patternsAndBlanks <- paste(growthForm.simp0$growthForms_fileNames2, growthForm.simp0$growthForms.blanks, sep=","); growthForm.simp0[1:30,]
growthForm.simp0$patternsAndBlanks2 <- gsub(",$","",growthForm.simp0$patternsAndBlanks); growthForm.simp0[1:30,]
growthForm.simp2 <- separate(data = growthForm.simp0, col = patternsAndBlanks2, into=c("symbol13", "symbol14", "symbol15"), sep = ","); growthForm.simp2[1:30,]
growthForm.simp <- subset(growthForm.simp2, select=c("genus_species", "symbol13", "symbol14", "symbol15")); growthForm.simp[1:30,]
growthForm.simp$symbol16 <- "separator.png"; growthForm.simp[1:30,]

# Last is tissue consumed. Occupy symbol positions 17,18,19 (19=separator)
head(tissue.summary); range(tissue.summary$numTissues)
tissue.simp0 <- tissue.summary; head(tissue.simp0)
tissue.simp0$tissues_fileNames2 <- str_replace(tissue.simp0$tissue_fileNames, ",spacer.png", ""); tissue.simp0[1:30,]
tissue.simp0$tissues.blanks <- ""
tissue.simp0$tissues.blanks[which(tissue.simp0$numTissues == 1)] <- "blank.png,blank.png"
tissue.simp0$tissues.blanks[which(tissue.simp0$numTissues == 2)] <- "blank.png"
tissue.simp0$tissues.blanks[which(tissue.simp0$numTissues == 3)] <- ""
tissue.simp0$patternsAndBlanks <- paste(tissue.simp0$tissues_fileNames2, tissue.simp0$tissues.blanks, sep=","); tissue.simp0[1:30,]
tissue.simp0$patternsAndBlanks[which(tissue.simp0$numTissues == 1)] <- "blank.png,blank.png"; tissue.simp0[1:30,]
tissue.simp0$patternsAndBlanks2 <- gsub(",$","",tissue.simp0$patternsAndBlanks); tissue.simp0[1:30,]
tissue.simp2 <- separate(data = tissue.simp0, col = patternsAndBlanks2, into=c("symbol17", "symbol18"), sep = ","); tissue.simp2[1:30,]
tissue.simp <- subset(tissue.simp2, select=c("genus_species", "symbol17", "symbol18")); tissue.simp[1:30,]
tissue.simp$symbol19 <- "separator.png"; tissue.simp[1:30,]

###### Merge all of these together into a single dataframe...
colors.patterns <- merge(colors.simp, patterns.simp, by="genus_species"); nrow(colors.patterns)
colors.patterns.diet <- merge(colors.patterns, dietBreadth.simp, by="genus_species"); nrow(colors.patterns.diet)
colors.patterns.diet.growthForm <- merge(colors.patterns.diet, growthForm.simp, by="genus_species"); nrow(colors.patterns.diet.growthForm)
allTraits8 <- merge(colors.patterns.diet.growthForm, tissue.simp, by="genus_species"); nrow(allTraits8)

# Now, make sure row order matches tip label order of tree...
spOrder <- tree$tip.label
allTraits8$genus_species <- factor(allTraits8$genus_species, levels=spOrder)
allTraits9 <- allTraits8[order(allTraits8$genus_species),]

# OK - now, the final step will be to concatenate down the columns, in order, for each symbol position. Here, we will KEEP the NAs in any species rows/ within the columns, since we want them to map on as blanks on those tips. Then, we need to remember to  add the required number of NAs at the end of the trait vector (1 less than the number of tips)
NAs <- rep(NA, nrow(allTraits9)-1)

symbol1 <- c(allTraits9$symbol1, NAs); length(symbol1); unique(allTraits9$symbol1) # 3615
symbol2 <- c(allTraits9$symbol2, NAs); length(symbol2); unique(allTraits9$symbol2) # 3615
symbol3 <- c(allTraits9$symbol3, NAs); length(symbol3); unique(allTraits9$symbol3) # 3615
symbol4 <- c(allTraits9$symbol4, NAs); length(symbol4); unique(allTraits9$symbol4) # 3615
symbol5 <- c(allTraits9$symbol5, NAs); length(symbol5); unique(allTraits9$symbol5) # 3615
symbol6 <- c(allTraits9$symbol6, NAs); length(symbol6); unique(allTraits9$symbol6) # 3615
symbol7 <- c(allTraits9$symbol7, NAs); length(symbol7); unique(allTraits9$symbol7) # 3615
symbol8 <- c(allTraits9$symbol8, NAs); length(symbol8); unique(allTraits9$symbol8) # 3615
symbol9 <- c(allTraits9$symbol9, NAs); length(symbol9); unique(allTraits9$symbol9) # 3615
symbol10 <- c(allTraits9$symbol10, NAs); length(symbol10); unique(allTraits9$symbol10) # 3615
symbol11 <- c(allTraits9$symbol11, NAs); length(symbol11); unique(allTraits9$symbol11) # 3615
symbol12 <- c(allTraits9$symbol12, NAs); length(symbol12); unique(allTraits9$symbol12) # 3615
symbol13 <- c(allTraits9$symbol13, NAs); length(symbol13); unique(allTraits9$symbol13) # 3615
symbol14 <- c(allTraits9$symbol14, NAs); length(symbol14); unique(allTraits9$symbol14) # 3615
symbol15 <- c(allTraits9$symbol15, NAs); length(symbol15); unique(allTraits9$symbol15) # 3615
symbol16 <- c(allTraits9$symbol16, NAs); length(symbol16); unique(allTraits9$symbol16) # 3615
symbol17 <- c(allTraits9$symbol17, NAs); length(symbol17); unique(allTraits9$symbol17) # 3615
symbol18 <- c(allTraits9$symbol18, NAs); length(symbol18); unique(allTraits9$symbol18) # 3615
symbol19 <- c(allTraits9$symbol19, NAs); length(symbol19); unique(allTraits9$symbol19) # 3615

# Make the ginormous tree!
tree_vert <- ggtree(tree, size=.25) + 
  geom_image(image=symbol1, nudge_x=.0165, size=.0004) +
  geom_image(image=symbol2, nudge_x=.0175, size=.0004) +
  geom_image(image=symbol3, nudge_x=.0185, size=.0004) +
  geom_image(image=symbol4, nudge_x=.0195, size=.0004) + # separator
  geom_image(image=symbol5, nudge_x=.0205, size=.0004) +
  geom_image(image=symbol6, nudge_x=.0215, size=.0004) +
  geom_image(image=symbol7, nudge_x=.0225, size=.0004) +
  geom_image(image=symbol8, nudge_x=.0235, size=.0004) +
  geom_image(image=symbol9, nudge_x=.0245, size=.0004) +
  geom_image(image=symbol10, nudge_x=.0255, size=.0004) + # separator
  geom_image(image=symbol11, nudge_x=.0265, size=.0004) +
  geom_image(image=symbol12, nudge_x=.0275, size=.0004) + # separator
  geom_image(image=symbol13, nudge_x=.0285, size=.0004) +
  geom_image(image=symbol14, nudge_x=.0295, size=.0004) +
  geom_image(image=symbol15, nudge_x=.0305, size=.0004) +
  geom_image(image=symbol16, nudge_x=.0315, size=.0004) + #separator
  geom_image(image=symbol17, nudge_x=.0325, size=.0004) +
  geom_image(image=symbol18, nudge_x=.0335, size=.0004) +
  geom_image(image=symbol19, nudge_x=.0345, size=.0004) + # separator
#  geom_vline(xintercept=c(.75, 1, 1.25, 1.5, 1.75), linetype="dotted") +
  geom_tiplab(size=1, offset=-.004, hjust=0, color="black", align=FALSE) +
  ggplot2::xlim(0,2)
# PDF 100 (height) x 100 (width)


# Circular tree without symbols
tree_rad <- ggtree(tree, layout="circular") +
  geom_tiplab(size=1.5, offset=0.01, color="black", align=FALSE) +
  ggplot2::xlim(-1.02, 1.02)
# 60 x 60 PDF