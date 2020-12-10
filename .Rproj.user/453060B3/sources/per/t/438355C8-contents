# clear R workspace
rm(list=ls())

############## load libraries #################
library(mvabund) # for model fitting
library(cluster) # for gower similarity and pam
library(plyr) # for rbind.fill
library(dplyr) # for data cleaning & piping
library(Cairo) # to save plot image
library(tidyr) # for gather()
library(fmsb) # for radar plot

# list data types
datas<-c("Australia","Spain")

##############################################
# read and clean raw data
source("raw.data/clean.raw.data.2019.R")

# run rare neighbor species analysis
source("Species.grouping/Rare.grouping/rare.species.GROUPED.2019.kfold.by.focal.R")
source("Species.grouping/Rare.grouping/rare.species.REMOVED.2019.kfold.by.focal.R")

# create clusters
source("Species.grouping/Functional.trait.grouping/create.clusters.2019.R")

# create data list with various neighbor groupings
source("Neigh.ID/making.data.list.groupings.2019.R")

# run neighbor ID analysis
source("Neigh.ID/run.neigh.ID.whole.model.2019.direct.kfolds.by.focal.R")
source("Neigh.ID/run.neigh.ID.whole.model.2019.HOI.kfolds.by.focal.R")
