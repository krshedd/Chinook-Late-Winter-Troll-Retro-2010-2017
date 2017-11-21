#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### SEAK Chinook Late Winter Troll 2010-2017 ####
# Kyle Shedd Tue Nov 21 11:00:08 2017
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to revisit Chinook salmon mixtures from the SEAK
# commercial late winter troll harvests from 2010-2017 using the GAPS3.0
# baseline containing 357 populations in 26 reporting groups characterized by 
# 13 uSATs. All mixtures are to be analyzed with the program BAYES.
# This re-analysis is similar to what was done with Spring Troll, this is to
# aid in the development of action plans for the new stocks of concern in 2018
# BoF.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Specific Objectives ####
# This script will:
# 1) Import mixture data
# 2) Add attribute data
# 3) Define spatio-temporal strata
# 4) Perform a data QC on mixtures
# 5) Prepare BAYES input files
# 6) Summarize BAYES results
# 7) Generate plots and tables of results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Late Winter Troll 2010-2017")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
username <- "krshedd"
password <- "********"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
LWinter2010Mixtures <- c("KWINT10LC", "KWINT10LJ", "KWINT10LK", "KWINT10LPG", "KWINT10LS", "KWINT10LY")
LWinter2011Mixtures <- c("KWINT11LC", "KWINT11LJ", "KWINT11LK", "KWINT11LPG", "KWINT11LS", "KWINT11LY")
LWinter2012Mixtures <- c("KWINT12LJ", "KWINT12LK", "KWINT12LPG", "KWINT12LS", "KWINT12LY")  # "KWINT12LC", 
LWinter2013Mixtures <- c("KWINT13LC", "KWINT13LJ", "KWINT13LK", "KWINT13LPG", "KWINT13LS", "KWINT13LY")
LWinter2014Mixtures <- c("KWINT14LC", "KWINT14LJ", "KWINT14LK", "KWINT14LPG", "KWINT14LS")  # , "KWINT14LY"
LWinter2015Mixtures <- c("KWINT15LC", "KWINT15LJ", "KWINT15LK", "KWINT15LPG", "KWINT15LS")
LWinter2016Mixtures <- c("KWINT16LC", "KWINT16LJ", "KWINT16LK", "KWINT16LPG", "KWINT16LS")
LWinter2017Mixtures <- c("KTROL17LW")

## Pull genotypes
LOKI2R_GAPS.GCL(sillyvec = unlist(sapply(objects(pattern = "LWinter"), get)), username = username, password = password); rm(password)


## Save unaltered .gcls
# dir.create("Raw genotypes")
# dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(unlist(sapply(objects(pattern = "LWinter"), get)), function(silly) {dput(x = get(paste0(silly, ".gcl")), file = paste0("Raw genotypes/OriginalCollections/" , silly, ".txt"))} )); beep(8)

# dir.create("Objects")
dput(x = LocusControl, file = "Objects/LocusControl.txt")
invisible(sapply(objects(pattern = "Mixtures"), function(mix) {dput(x = get(mix), file = paste0("Objects/", mix, ".txt"))}))
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GAPSLoci_reordered.txt", to = "Objects")
GAPSLoci_reordered <- dget(file = "Objects/GAPSLoci_reordered.txt")

dimnames(KTROL17LW.gcl$counts)[[2]]
GAPSLoci_reordered

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pool into a single silly per year
PoolCollections.GCL(collections = LWinter2010Mixtures, loci = GAPSLoci_reordered, newname = "KTROL10LW")
PoolCollections.GCL(collections = LWinter2011Mixtures, loci = GAPSLoci_reordered, newname = "KTROL11LW")
PoolCollections.GCL(collections = LWinter2012Mixtures, loci = GAPSLoci_reordered, newname = "KTROL12LW")
PoolCollections.GCL(collections = LWinter2013Mixtures, loci = GAPSLoci_reordered, newname = "KTROL13LW")
PoolCollections.GCL(collections = LWinter2014Mixtures, loci = GAPSLoci_reordered, newname = "KTROL14LW")
PoolCollections.GCL(collections = LWinter2015Mixtures, loci = GAPSLoci_reordered, newname = "KTROL15LW")
PoolCollections.GCL(collections = LWinter2016Mixtures, loci = GAPSLoci_reordered, newname = "KTROL16LW")
PoolCollections.GCL(collections = LWinter2017Mixtures, loci = GAPSLoci_reordered, newname = "KTROL17LW")

dimnames(KTROL16LW.gcl$counts)[[2]]

sapply(paste0("KTROL", 10:17, "LW"), function(silly) {get(paste0(silly, ".gcl"))$n} )
# KTROL10LW KTROL11LW KTROL12LW KTROL13LW KTROL14LW KTROL15LW KTROL16LW KTROL17LW 
# 522       517       491       520       489       567       572       499 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Change FK_FISH_ID to the back end of SillySource
str(KTROL10LW.gcl$attributes$FK_FISH_ID)
str(KTROL17LW.gcl$attributes$FK_FISH_ID)

KTROL10LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL10LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL11LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL11LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL12LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL12LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL13LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL13LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL14LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL14LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL15LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL15LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL16LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL16LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )
KTROL17LW.gcl$attributes$FK_FISH_ID <- sapply(as.character(KTROL17LW.gcl$attributes$SillySource), function(ind) {as.numeric(unlist(strsplit(x = ind, split = "_"))[2])} )


# dir.create("Raw genotypes/PooledCollections")
invisible(sapply(paste0("KTROL", 10:17, "LW"), function(silly) {dput(x = get(paste0(silly, ".gcl")), file = paste0("Raw genotypes/PooledCollections/" , silly, ".txt"))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Late Winter Troll 2010-2017/")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
SEAKobjects <- list.files(path = "Objects", recursive = FALSE)
# SEAKobjects <- SEAKobjects[-which(SEAKobjects == "Vials" | SEAKobjects == "OLD_BAD_LOCUSCONTROL")]
SEAKobjects

invisible(sapply(SEAKobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

## Get un-altered mixtures
invisible(sapply(paste0("KTROL", 10:17, "LW"), function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pair with metadata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

latewinter_troll.df <- read.csv(file = "Late Winter Troll ASL 2010-2017.csv", header = TRUE)
str(latewinter_troll.df)
addmargins(table(latewinter_troll.df$Year, latewinter_troll.df$Stat.Week))  # all samples


ids <- sapply(paste0("KTROL", 10:17, "LW"), function(silly) {get(paste0(silly, ".gcl"))$attributes$FK_FISH_ID} )
str(ids)

# Are we missing metadata for fish we have genotyped?
table(unlist(ids) %in% latewinter_troll.df$Dna.Specimen.No)
# FALSE  TRUE 
#     5  4172
# Got this donw to 1 fish

# Which years are missing metadata
table(sapply(names(unlist(ids)[!unlist(ids) %in% latewinter_troll.df$Dna.Specimen.No]), function(id) {
  unlist(strsplit(x = unlist(strsplit(x = id, split = "KTROL"))[2], split = "LW"))[1]
} ))
# 13 14 17 
#  1  1  3
# Fixed all but the 2014 fish

# Paste the missing fish into ASL .csv to see what project they are from
writeClipboard(as.character(unlist(ids)[!unlist(ids) %in% latewinter_troll.df$Dna.Specimen.No]))
# "171164" "215537" "421201" "735501" "785701"

# Match data up by year and look at district breakdowns
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "LW.gcl"))
  match.yr <- match(my.gcl$attributes$FK_FISH_ID, latewinter_troll.df$Dna.Specimen.No)
  # table(latewinter_troll.df$Year[match.yr])
  my.gcl$attributes$StatWeek <- latewinter_troll.df$Stat.Week[match.yr]
  my.gcl$attributes$Port <- latewinter_troll.df$Port.Code[match.yr]
  my.gcl$attributes$Quadrant <- latewinter_troll.df$District[match.yr]
  my.gcl$attributes$Age <- latewinter_troll.df$Age.European[match.yr]
  my.gcl$attributes$LengthType <- latewinter_troll.df$Length.Type[match.yr]
  my.gcl$attributes$Length <- latewinter_troll.df$Average.Length.mm[match.yr]
  
  assign(x = paste0("match.20", yr), value = match.yr)
  assign(x = paste0("KTROL", yr, "LW.gcl"), value = my.gcl)
}

sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "LW.gcl"))
  addmargins(table(my.gcl$attributes$Quadrant, useNA = "always"))
} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

LWint10_17_Strata <- paste0("KTROL", 10:17, "LW")
dput(x = LWint10_17_Strata, file = "Objects/LWint10_17_Strata.txt")

LWint10_17_Strata_SampleSizes <- matrix(data = NA, nrow = length(LWint10_17_Strata), ncol = 4, 
                                         dimnames = list(LWint10_17_Strata, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_LWint10_17_Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = LWint10_17_Strata, loci = GAPSLoci_reordered)
min(Original_LWint10_17_Strata_SampleSizebyLocus)  ## 477
apply(Original_LWint10_17_Strata_SampleSizebyLocus, 1, min) / apply(Original_LWint10_17_Strata_SampleSizebyLocus, 1, max)  ## Good, 0.952

Original_LWint10_17_Strata_PercentbyLocus <- apply(Original_LWint10_17_Strata_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_LWint10_17_Strata_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_LWint10_17_Strata_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares

#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_LWint10_17_Strata_ColSize <- sapply(paste0(LWint10_17_Strata, ".gcl"), function(x) get(x)$n)
LWint10_17_Strata_SampleSizes[, "Genotyped"] <- Original_LWint10_17_Strata_ColSize

### Missing
## Remove individuals with >20% missing data
LWint10_17_Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = LWint10_17_Strata, proportion = 0.8)
dput(x = LWint10_17_Strata_MissLoci, file = "Objects/LWint10_17_Strata_MissLoci.txt")

## Get number of individuals per silly after removing missing loci individuals
ColSize_LWint10_17_Strata_PostMissLoci <- sapply(paste0(LWint10_17_Strata, ".gcl"), function(x) get(x)$n)
LWint10_17_Strata_SampleSizes[, "Missing"] <- Original_LWint10_17_Strata_ColSize - ColSize_LWint10_17_Strata_PostMissLoci

### Duplicate
## Check within collections for duplicate individuals.
LWint10_17_Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = LWint10_17_Strata, loci = GAPSLoci_reordered, quantile = NULL, minproportion = 0.95)
LWint10_17_Strata_DuplicateCheckReportSummary <- sapply(LWint10_17_Strata, function(x) LWint10_17_Strata_DuplicateCheck95MinProportion[[x]]$report)
LWint10_17_Strata_DuplicateCheckReportSummary
dput(x = LWint10_17_Strata_DuplicateCheckReportSummary, file = "Objects/LWint10_17_Strata_DuplicateCheckReportSummary.txt")

## Remove duplicate individuals
LWint10_17_Strata_RemovedDups <- RemoveDups.GCL(LWint10_17_Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_LWint10_17_Strata_PostDuplicate <- sapply(paste0(LWint10_17_Strata, ".gcl"), function(x) get(x)$n)
LWint10_17_Strata_SampleSizes[, "Duplicate"] <- ColSize_LWint10_17_Strata_PostMissLoci-ColSize_LWint10_17_Strata_PostDuplicate

### Final
LWint10_17_Strata_SampleSizes[, "Final"] <- ColSize_LWint10_17_Strata_PostDuplicate
LWint10_17_Strata_SampleSizes

dput(x = LWint10_17_Strata_SampleSizes, file = "Objects/LWint10_17_Strata_SampleSizes.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC .gcl's as back-up:
# dir.create("Raw genotypes/PooledCollections_PostQC")
invisible(sapply(LWint10_17_Strata, function(silly) {
  dput(x = get(paste(silly, ".gcl", sep = '')), file = paste0("Raw genotypes/PooledCollections_PostQC/" , silly, ".txt"))
} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Add Genotyped TRUE/FALSE to ASL.df
ids.genotpyed <- sapply(paste0("KTROL", 10:17, "LW"), function(silly) {get(paste0(silly, ".gcl"))$attributes$FK_FISH_ID} )
dput(x = ids.genotpyed, file = "Objects/ids.genotyped.txt")


latewinter_troll.df$Genotyped <- latewinter_troll.df$Dna.Specimen.No %in% unlist(ids.genotpyed)
write.csv(x = latewinter_troll.df, file = "Late Winter Troll ASL 2010-2017_new.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/Spring Troll 2010-2017/")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
SEAKobjects <- list.files(path = "Objects", recursive = FALSE)
# SEAKobjects <- SEAKobjects[-which(SEAKobjects == "Vials" | SEAKobjects == "OLD_BAD_LOCUSCONTROL")]
SEAKobjects

invisible(sapply(SEAKobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

## Get un-altered mixtures
invisible(sapply(paste0("KTROL", 10:17, "SP"), function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections_PostQC/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")
