# identify leading columns that hold non-species information
X.M<-sum(colSums(Australia.data[5:49]))
# identify neighbor species that have less than 1% of the total number of individuals across all plots
group.nm.M<-names(which(colSums(Australia.data[5:49])<=round(X.M*0.01)))
# summ across rows for rare species 
rares.M<-rowSums(Australia.data[,which(names(Australia.data) %in% group.nm.M)])
# remove rare species from dataframe
Australia.data.sub<-Australia.data[,-which(names(Australia.data) %in% group.nm.M)]
# rename dataframe with rares removed
Australia.data.sub.rm<-Australia.data.sub


# identify leading columns that hold non-species information
X.G<-sum(colSums(Spain.data[4:20]))
# identify neighbor species that have less than 1% of the total number of individuals across all plots
group.nm.G<-names(which(colSums(Spain.data[4:20])<=round(X.G*0.01)))
# summ across rows for rare species 
rares.G<-rowSums(Spain.data[,which(names(Spain.data) %in% group.nm.G)])
# remove rare species from dataframe
Spain.data.sub<-Spain.data[,-which(names(Spain.data) %in% group.nm.G)]
# rename dataframe with rares removed
Spain.data.sub.rm<-Spain.data.sub

# make data list
data.list.W<-list()
data.list.W[["Spain"]]<-Spain.data.sub.rm
data.list.W[["Australia"]]<-Australia.data.sub.rm

########## whole model removed rare species
# make empty summary dataframe
sum<-data.frame(Model=rep(100,5000),beta=rep("F",5000), 
                coef=rep(0,5000),AIC=rep(0,5000),RMSE=rep(0,5000))
sum$beta<-as.character(sum$beta)
b<-1

# beta flag
beta<-T

# source model could
source("Species.grouping/Rare.grouping/fecundity.model.NEW.AB.focal.mvabund.R")

for (d in datas) { # for each data type
  
  fecundity.data<-data.list.W[[d]] # read in data
  
  # run direct-only model
  alpha.models.AB <- fecundity.model(fecundity.data, fit.alphas=TRUE, fit.betas=FALSE)
  
  # run HOI inclusive model
  if (beta==T) {
    alpha.beta.models.AB <- fecundity.model(fecundity.data, fit.alphas=TRUE, fit.betas=TRUE)
  }
  
  #### fill in summary dataframe
  sum$Model[b]<-d
  sum$beta[b]<-paste0("F")
  sum$coef[b]<-length(alpha.models.AB$coefficients)
  sum$AIC[b]<-AIC(alpha.models.AB)
  sum$RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
  sum$likeli[b]<-logLik(alpha.models.AB)
  
  b<-b+1
  
  if (beta==T){
    sum$Model[b]<-d
    sum$beta[b]<-paste0("T")
    sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
    sum$AIC[b]<-AIC(alpha.beta.models.AB)
    sum$RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
    sum$likeli[b]<-logLik(alpha.beta.models.AB)
    
    b<-b+1
  }
  rm(alpha.beta.models.AB)
  rm(alpha.models.AB)
}

# clean summary dataframe
sum<-unique(sum)
sum<-sum[-dim(sum)[1],]
write.csv(sum,paste0("Species.grouping/Rare.grouping/output/rare.removed.summary.2019.csv"),row.names = FALSE)


################# by focal #################
# make empty summary dataframe
sum<-data.frame(Model=rep(100,5000),beta=rep("F",5000),focal=as.character(rep("O",5000)), 
                coef=rep(0,5000),AIC=rep(0,5000),RMSE=rep(0,5000))
sum$focal<-as.character(sum$focal)
sum$beta<-as.character(sum$beta)
b<-1

# beta flag
beta<-T

# read in model function for running model by species
source("Species.grouping/Rare.grouping/fecundity.model.NEW.AB.focal.mvabund.by.species.R")

for (d in datas) { # for each data type
  print(d)
  fecundity.data.W<-data.list[[d]] # read in non-pre-removed rare species data
  
  # get a list of focal species
  focals<-unique(fecundity.data.W$focal)
  
  for(f in focals) { # for each focal species
    
    fecundity.data<-fecundity.data.W[which(fecundity.data.W$focal==f),] # get data for just focal species
    
    # remove rare species for focal specific dataset
    X.M<-sum(colSums(fecundity.data[5:dim(fecundity.data)[2]]))
    group.nm.M<-names(which(colSums(fecundity.data[5:dim(fecundity.data)[2]])<=round(X.M*0.01)))
    rares.M<-rowSums(fecundity.data[,which(names(fecundity.data) %in% group.nm.M)])
    fecundity.data.sub<-fecundity.data[,-which(names(fecundity.data) %in% group.nm.M)]
    
    # run direct-only models  
    alpha.models.AB <- fecundity.model(fecundity.data.sub, fit.alphas=TRUE, fit.betas=FALSE)
    
    
    # run HOI models
    if (beta==T) {
      alpha.beta.models.AB <- fecundity.model(fecundity.data.sub, fit.alphas=TRUE, fit.betas=TRUE)
    }
    
    
    #### fill in summary dataframe
    sum$Model[b]<-d
    sum$beta[b]<-paste0("F")
    sum$focal[b]<-paste0(f)
    sum$coef[b]<-length(alpha.models.AB$coefficients)
    sum$AIC[b]<-AIC(alpha.models.AB)
    sum$RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
    sum$likeli[b]<-logLik(alpha.models.AB)
    b<-b+1
    
    if (beta==T){
      sum$Model[b]<-d
      sum$beta[b]<-paste0("T")
      sum$focal[b]<-paste0(f)
      sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
      sum$AIC[b]<-AIC(alpha.beta.models.AB)
      sum$RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
      sum$likeli[b]<-logLik(alpha.beta.models.AB)
      
      b<-b+1
    }
    rm(alpha.beta.models.AB)
    rm(alpha.models.AB)
  }
}

# clean and save summary dataframe
sum<-unique(sum)
sum<-sum[-dim(sum)[1],]
write.csv(sum,paste0("Species.grouping/Rare.grouping/output/rare.removed.by.species.2019.csv"),row.names = FALSE)

print("Rare species removed done")