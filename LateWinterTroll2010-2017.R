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

# Add Month variable
latewinter_troll.df$Date <- as.Date(latewinter_troll.df$Sample.Date, format = "%m/%d/%Y")
latewinter_troll.df$Month <- format(x = latewinter_troll.df$Date, format = "%m")
latewinter_troll.df$Month.abb <- factor(x = month.abb[as.numeric(latewinter_troll.df$Month)], levels = month.abb)

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
  my.gcl$attributes$Date <- latewinter_troll.df$Date[match.yr]
  my.gcl$attributes$Month <- latewinter_troll.df$Month[match.yr]
  my.gcl$attributes$Month.abb <- latewinter_troll.df$Month.abb[match.yr]
  
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
# dir.create("Raw genotypes/PooledCollections_PostQC_Metadata")
invisible(sapply(LWint10_17_Strata, function(silly) {
  dput(x = get(paste(silly, ".gcl", sep = '')), file = paste0("Raw genotypes/PooledCollections_PostQC_Metadata/" , silly, ".txt"))
} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Add Genotyped TRUE/FALSE to ASL.df
ids.genotyped <- sapply(paste0("KTROL", 10:17, "LW"), function(silly) {get(paste0(silly, ".gcl"))$attributes$FK_FISH_ID} )
dput(x = ids.genotyped, file = "Objects/ids.genotyped.txt")


latewinter_troll.df$Genotyped <- latewinter_troll.df$Dna.Specimen.No %in% unlist(ids.genotpyed)
write.csv(x = latewinter_troll.df, file = "Late Winter Troll ASL 2010-2017_new.csv", row.names = FALSE)

# Pivot of # genotyped by year by month
require(reshape)
cast(data = latewinter_troll.df, Year ~ Month.abb, value = "Genotyped", subset = Genotyped == TRUE, margins = TRUE)
dput(x = latewinter_troll.df, file = "Objects/latewinter_troll.df.txt")


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

## Get the postQC_metadata mixtures
invisible(sapply(LWint10_17_Strata, function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0("Raw genotypes/PooledCollections_PostQC_Metadata/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pivot of years by mixtures with NA
addmargins(t(sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "LW.gcl"))
  table(my.gcl$attributes$Month.abb, useNA = "always")
} )))


# Add new factor with mixtures
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "LW.gcl"))
  
  my.gcl$attributes$Mixture <- NA
  my.gcl$attributes$Mixture[my.gcl$attributes$Month.abb %in% "Jan"] <- "Jan"
  my.gcl$attributes$Mixture[my.gcl$attributes$Month.abb %in% "Feb"] <- "Feb"
  my.gcl$attributes$Mixture[my.gcl$attributes$Month.abb %in% "Mar"] <- "Mar"
  my.gcl$attributes$Mixture[my.gcl$attributes$Month.abb %in% c("Apr", "May")] <- "Apr"

  my.gcl$attributes$Mixture <- factor(x = my.gcl$attributes$Mixture, levels = month.abb[1:4])
  
  assign(x = paste0("KTROL", yr, "LW.gcl"), value = my.gcl)
}


# Pivot of years by mixtures with NA

mixture.n.mat <- t(sapply(as.character(10:17), function(yr) {
  my.gcl <- get(paste0("KTROL", yr, "LW.gcl"))
  table(my.gcl$attributes$Mixture, useNA = "always")
} ))

addmargins(mixture.n.mat)


mixtures.yr <- apply(mixture.n.mat, 1, function(yr) {which(yr > 100)} )
dput(x = mixtures.yr, file = "Objects/mixtures.yr.txt")

mixtures <- levels(KTROL10LW.gcl$attributes$Mixture)
dput(x = mixtures, file = "Objects/mixtures.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES")
# dir.create("BAYES/Mixture")

# Loop over years and mixtures
for(yr in 10:17){
  my.gcl <- get(paste0("KTROL", yr, "LW.gcl"))
  
  for(mix in mixtures.yr[[as.character(yr)]]) {
    IDs <- list("my" = na.omit(AttributesToIDs.GCL(silly = "my", attribute = "Mixture", matching = mixtures[mix])))
    if(length(IDs[["my"]])) {
      invisible(CreateMixture.GCL(sillys = "my", loci = GAPSLoci_reordered, IDs = IDs, 
                                  mixname = paste0("AllQuad", mixtures[mix], "Troll_20", yr), 
                                  dir = "BAYES/Mixture/", type = "BAYES", PT = FALSE))
    }  # if IDS
  }  # mixture within year
  
}  # year

all.mixtures <- unlist(sapply(names(mixtures.yr), function(yr) {paste0("AllQuad", names(mixtures.yr[[yr]]), "Troll_20", yr)} ))
all.mixtures <- unname(all.mixtures)
dput(x = all.mixtures, file = "Objects/all.mixtures.txt")

# Double check mixture files

sapply(list.files(path = "BAYES/Mixture/", full.names = TRUE), function(fle) {nrow(read.table(file = fle, header = FALSE))} )
mixture.n.mat


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES/Baseline")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/BAYES/Baseline/GAPS357pops13loci.bse", to = "BAYES/Baseline/")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupVec26RG_357.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupNames26.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/SEAKPops357.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GAPS357PopsInits.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/WASSIPSockeyeSeeds.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/mixfortran.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/bayesfortran_357.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/Spring Troll 2010-2017/Objects/GAPS357PopFlatPrior.txt", to = "Objects")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES/Control")

SEAKobjects <- list.files(path = "Objects", recursive = FALSE)
invisible(sapply(SEAKobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); rm(SEAKobjects); beep(2)

sapply(all.mixtures, function(Mix) {
  invisible(CreateControlFile.GCL(sillyvec = SEAKPops357, loci = GAPSLoci_reordered, mixname = Mix, basename = "GAPS357pops13loci", suffix = "", nreps = 40000, nchains = 5,
                                  groupvec = GroupVec26RG_357, priorvec = GAPS357PopFlatPrior, initmat = GAPS357PopsInits, dir = "BAYES/Control",
                                  seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = mixfortran, basefortran = bayesfortran_357, switches = "F T F T T T F"))
} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BAYES/Output")
sapply(all.mixtures, function(Mix) {
  invisible(dir.create(path = paste0("BAYES/Output/", Mix)))
})




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LWint10_17_26RG_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:26, groupnames = GroupNames26, maindir = "BAYES/Output", mixvec = all.mixtures,
                               prior = '', ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)
# dir.create("Estimates objects")
dput(x = LWint10_17_26RG_Estimates, file = "Estimates objects/LWint10_17_26RG_Estimates.txt")
sapply(LWint10_17_26RG_Estimates, function(mix) {table(mix[, "GR"] > 1.2)})  # March2012 and March 2017 had 1) RG each over 1.2
sapply(c("AllQuadMarTroll_2012", "AllQuadMarTroll_2017"), function(mix) {LWint10_17_26RG_Estimates[[mix]][, "GR"]})


file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupNames4.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupNames4Pub.txt", to = "Objects")
file.copy(from = "V:/Analysis/1_SEAK/Chinook/Mixture/SEAK17/Objects/GroupVec4.txt", to = "Objects")

LWint10_17_4RG_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec4, groupnames = GroupNames4Pub, maindir = "BAYES/Output", mixvec = all.mixtures,
                               prior = '', ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)
dput(x = LWint10_17_4RG_Estimates, file = "Estimates objects/LWint10_17_4RG_Estimates.txt")
LWint10_17_4RG_Estimates <- dget(file = "Estimates objects/LWint10_17_4RG_Estimates.txt")

sapply(LWint10_17_4RG_Estimates, function(mix) {table(mix[, "GR"] > 1.2)})
sapply(LWint10_17_4RG_Estimates, function(mix) {mix[c("Alaska", "TBR"), "mean"]})

sapply(LWint10_17_4RG_Estimates, function(mix) {mix[, "mean"]} )
sapply(LWint10_17_4RG_Estimates, function(mix) {mix[, "sd"] / mix[, "mean"]} )  # CV
sapply(LWint10_17_4RG_Estimates, function(mix) {sum(mix[, "sd"] / mix[, "mean"] < 0.20)} )  # how many RGs with CV < 20%
sapply(LWint10_17_4RG_Estimates, function(mix) {mix[, "95%"] - mix[, "5%"]} )  # 90% CI range


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Read in Harvest and Sample Size Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(xlsx)
harvest.df <- read.xlsx(file = "Late Winter Troll ASL 2010-2017.xlsx", sheetName = "CE000678", startRow = 23, header = TRUE)
str(harvest.df)

harvest.df$Mixture <- NA
harvest.df$Mixture[harvest.df$District %in% 101:102] <- "101/102"
harvest.df$Mixture[harvest.df$District %in% 103] <- "103"
harvest.df$Mixture[harvest.df$District %in% 106:108 & harvest.df$Area.Value != 10643] <- "106/107/108"
harvest.df$Mixture[harvest.df$District %in% c(109:110, 112) & harvest.df$Area.Value != 11265] <- "109/110/112"
harvest.df$Mixture[harvest.df$District %in% 113] <- "113"
harvest.df$Mixture[harvest.df$District %in% 114 | harvest.df$Area.Value %in% c(11265, 11395, 11397)] <- "114"
harvest.df$Mixture[harvest.df$District %in% 183] <- "183"

harvest.df$Mixture <- factor(x = harvest.df$Mixture, levels = c("101/102", "103", "106/107/108", "109/110/112", "113", "114", "183"))

dput(x = harvest.df, file = "Objects/harvest.df.txt")

require(reshape)
harvest_mix.df <- aggregate(N.Catch ~ Year + Mixture, data = harvest.df, sum)
str(harvest_mix.df)
cast(harvest_mix.df, Year ~ Mixture, value = "N.Catch")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample sizes
all.mixtures.samplesize <- sapply(all.mixtures, function(mix) {dim(read.table(file = paste0("BAYES/Mixture/", mix, ".mix")))[1]} )
dput(x = all.mixtures.samplesize, file = "Objects/all.mixtures.samplesize.txt")

mixtures.names
mixtures.names2 <- names(mixtures.names)
names(mixtures.names2) <- mixtures.names

mixtures.df <- as.data.frame(t(sapply(all.mixtures, function(mix) {unlist(strsplit(x = mix, split = "_"))} )), stringsAsFactors = FALSE)
names(mixtures.df) <- c("Mixname", "Year")
mixtures.df$Year <- as.numeric(mixtures.df$Year)
mixtures.df$Full.Mixname <- all.mixtures
mixtures.df$Mixture <- factor(x = mixtures.names2[mixtures.df$Mixname], levels = levels(harvest_mix.df$Mixture))
mixtures.df$n <- all.mixtures.samplesize

dput(x = mixtures.df, file = "Objects/mixtures.df.txt")


str(mixtures.df)

all.mixtures.n100 <- names(which(all.mixtures.samplesize >= 100))
dput(x = all.mixtures.n100, file = "Objects/all.mixtures.n100.txt")

# Subset data for n >= 100
Spring10_17_4RG_Estimates_n100 <- Spring10_17_4RG_Estimates[all.mixtures.n100]
dput(x = Spring10_17_4RG_Estimates_n100, file = "Estimates objects/Spring10_17_4RG_Estimates_n100.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge harvest and sample size
spring.estimates.df <- merge(x = harvest_mix.df, y = mixtures.df, by = c("Mixture", "Year"), all = FALSE)
spring.estimates.df$Alaska.mean.p <- sapply(Spring10_17_4RG_Estimates, function(mix) {mix["Alaska", "mean"]})
spring.estimates.df$TBR.mean.p <- sapply(Spring10_17_4RG_Estimates, function(mix) {mix["TBR", "mean"]})
spring.estimates.df$Alaska.mean.C <- spring.estimates.df$Alaska.mean.p * spring.estimates.df$N.Catch
spring.estimates.df$TBR.mean.C <- spring.estimates.df$TBR.mean.p * spring.estimates.df$N.Catch

round(cast(spring.estimates.df, Year ~ Mixture, value = "Alaska.mean.C"))
round(cast(spring.estimates.df, Year ~ Mixture, value = "TBR.mean.C"))
round(cast(spring.estimates.df, Year ~ Mixture, value = "n"))

dput(x = spring.estimates.df, file = "Objects/spring.estimates.df.txt")

harvest <- setNames(object = spring.estimates.df$N.Catch, nm = spring.estimates.df$Full.Mixname)
dput(x = harvest, file = "Objects/harvest.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Subset to only mixtures with >= 100 fish
spring.estimates.n100.df <- spring.estimates.df
spring.estimates.n100.df[spring.estimates.n100.df$n < 100, c("Alaska.mean.p", "Alaska.mean.C", "TBR.mean.p", "TBR.mean.C")] <- NA

# Heatmap of total Catch
require(lattice)
new.colors <- colorRampPalette(c("white", "darkgreen"))
data.mat <- as.matrix(round(cast(spring.estimates.n100.df, Year ~ Mixture, value = "N.Catch")))
# data.mat[is.na(data.mat)] <- 0
levelplot(data.mat, 
          col.regions = new.colors, 
          at = seq(from = 0, to = max(data.mat, na.rm = TRUE), length.out = 100), 
          main = "Total Catch", xlab = "Year", ylab = "District Area", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill",
          panel = function(...) {
            panel.fill("black")
            panel.levelplot(...)}
)  # aspect = "iso" will make squares


# Heatmap of mean Alaska
require(lattice)
new.colors <- colorRampPalette(c("white", "darkblue"))
data.mat <- as.matrix(cast(spring.estimates.n100.df, Year ~ Mixture, value = "Alaska.mean.p")) * 100
levelplot(data.mat, 
          col.regions = new.colors, 
          at = seq(from = 0, to = 100, length.out = 100), 
          main = "Mean Alaska %", xlab = "Year", ylab = "District Area", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill", 
          panel = function(...) {
            panel.fill("black")
            panel.levelplot(...)}
)  # aspect = "iso" will make squares

# Heatmap of mean TBR %
require(lattice)
new.colors <- colorRampPalette(c("white", "darkblue"))
data.mat <- as.matrix(cast(spring.estimates.n100.df, Year ~ Mixture, value = "TBR.mean.p")) * 100
levelplot(data.mat, 
          col.regions = new.colors, 
          at = seq(from = 0, to = 100, length.out = 100), 
          main = "Mean TBR %", xlab = "Year", ylab = "District Area", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill", 
          panel = function(...) {
            panel.fill("black")
            panel.levelplot(...)}
)  # aspect = "iso" will make squares


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create 4RG Summary Tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dir.create("Estimates tables")
require(xlsx)

EstimatesStats <- Spring10_17_4RG_Estimates_n100
SampSizes <- all.mixtures.samplesize
HarvestVec <- harvest
PubNames <- setNames(object = paste("Spring Troll", spring.estimates.df$Year, "District(s)", spring.estimates.df$Mixture),
                     nm = spring.estimates.df$Full.Mixname)

for(mix in all.mixtures.n100) {
  
  TableX <- matrix(data = "", nrow = 7, ncol = 7)
  TableX[1, 1] <- paste0(PubNames[mix], " (n=", SampSizes[mix], ", catch=", formatC(x = HarvestVec[mix], digits = 0, big.mark = ",", format = "f"), ")")
  TableX[2, 6] <- "90% CI"
  TableX[3, 2:7] <- c("Reporting Group", "Mean", "SD", "Median", "5%", "95%")
  TableX[4:7, 1] <- 1:4
  TableX[4:7, 2] <- rownames(EstimatesStats[[mix]])
  TableX[4:7, 3:7] <- formatC(x = EstimatesStats[[mix]][, c("mean", "sd", "median", "5%", "95%")], digits = 3, format = "f")
  
  write.xlsx(x = TableX, file = "Estimates tables/SpringTroll2017_4RG_Estimates.xlsx",
             col.names = FALSE, row.names = FALSE, sheetName = paste(mix, " 4RG"), append = TRUE)
  
}


# save.image("SpringTroll2010-2017.RData")
