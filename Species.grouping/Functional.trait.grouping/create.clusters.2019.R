## read in trait data
traits<-read.csv("Species.grouping/Functional.trait.grouping/species.traits.Australia.Spain.2019.csv")
# find clusters and make plots (in supps)
source("Species.grouping/Functional.trait.grouping/cluster.functions.2019.R")

# list the columns to remove from the cluster analysis
#   all = all quantiative traits and qualitative traits (exotic/naitve, life history, leaf type, grass/forb)
cols.rm.all<-c("Data","SP.CODE","Species","family","NOTES","TRY.ID","NOTES.TRY","species.dataset.code","NOTES.1","exotic")

#returns the traits dataframe with the cluster values attached
# remove should be the respective column remove vector as the type
traits.all<-find.clusters(traits=traits,remove=cols.rm.all,type="all")

# write out final large dataset
write.csv(traits.all,"Species.grouping/Functional.trait.grouping/output/traits.cluster.all.results.2019.csv")

# write out dataset with just the cluster and ID that will run in the model
model.cols<-c("Data","SP.CODE","species.dataset.code","cluster.all")
traits.model<-traits.all[,model.cols]
write.csv(traits.model,"Species.grouping/Functional.trait.grouping/output/clusters.for.model.2019.csv")
print("Creating trait clustering group done")
