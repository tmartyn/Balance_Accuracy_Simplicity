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
datas<-c("Spain","Australia")
########## whole model removed rare species
# make empty summary dataframe
sum<-data.frame(Model=rep("O",5000),beta=rep("F",5000), 
                coef=rep(0,5000),AIC=rep(0,5000),train.RMSE=rep(0,5000),
                predict.RMSE=rep(0,5000),
                BIC=rep(0,5000))
sum$beta<-as.character(sum$beta)
sum$Model<-as.character(sum$Model)
b<-1

# beta flag
beta<-T

# source model could
source("Species.grouping/Rare.grouping/fecundity.model.NEW.AB.focal.mvabund.R")
source("predict.manyglm.R")

for (d in datas) { # for each data type
  
  fecundity.data<-data.list.W[[d]] # read in data
  data.length<-dim(fecundity.data)[1]
  require(caret)
  flds <- createFolds(seq(1:data.length), k = 5, list = TRUE, returnTrain = FALSE)
  
  for (k in 1:5){
    
    fecundity.data2<-fecundity.data[-flds[[k]],]
    
    # run direct-only model
    alpha.models.AB <- fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=FALSE)
    newdata<-fecundity.data[flds[[k]],]
    newdata$fecundity<-newdata$focal
    predict.alpha<-predict.manyglm2(alpha.models.AB,newdata=newdata)
    predict.alpha.2<-predict.alpha
    observed.alpha<-log(fecundity.data$seeds[flds[[k]]]+1)
    predict.RMSE.alpha<-sqrt(mean((observed.alpha-predict.alpha.2)^2))
    
    # run HOI inclusive model
    #if (beta==T) {
    alpha.beta.models.AB <- fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=TRUE)
    newdata<-fecundity.data[flds[[k]],]
    newdata$fecundity<-newdata$focal
    predict.beta<-predict.manyglm2(alpha.beta.models.AB,newdata=newdata)
    predict.beta.2<-predict.beta
    observed.beta<-log(fecundity.data$seeds[flds[[k]]]+1)
    predict.RMSE.beta<-sqrt(mean((observed.beta-predict.beta.2)^2))
    #}
    
    sum$Model[b]<-paste0(d)
    sum$beta[b]<-paste0("F")
    sum$fold[b]<-k
    sum$train.RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
    sum$predict.RMSE[b]<-predict.RMSE.alpha
    sum$coef[b]<-length(alpha.models.AB$coefficients)
    sum$AIC[b]<-AIC(alpha.models.AB)
    sum$BIC[b]<-BIC(alpha.models.AB)
    sum$likeli[b]<-logLik(alpha.models.AB)
    ## increment summary DF
    b<-b+1
    
    ## print model statistics to summary DF for beta models
    #if (beta==T){
    sum$Model[b]<-paste0(d)
    sum$beta[b]<-paste0("T")
    sum$fold[b]<-k
    sum$train.RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
    sum$predict.RMSE[b]<-predict.RMSE.beta
    sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
    sum$AIC[b]<-AIC(alpha.beta.models.AB)
    sum$BIC[b]<-BIC(alpha.beta.models.AB)
    sum$likeli[b]<-logLik(alpha.beta.models.AB)
    b<-b+1
    
    
    #}
    rm(alpha.beta.models.AB)
    rm(alpha.models.AB)
  } # end for each fold
}
# clean summary dataframe
sum<-unique(sum)
sum<-sum[-dim(sum)[1],]
write.csv(sum,paste0("Species.grouping/Rare.grouping/output/rare.removed.summary.kfolds.2019.csv"),row.names = FALSE)


################# by focal #################
# make empty summary dataframe
sum<-data.frame(Model=rep(100,5000),beta=rep("F",5000),focal=as.character(rep("O",5000)), 
                coef=rep(0,5000),AIC=rep(0,5000),train.RMSE=rep(0,5000),predict.RMSE=rep(0,5000),
                BIC=rep(0,5000))
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
    
    require(caret)
    flds <- createFolds(seq(1:data.length), k = 5, list = TRUE, returnTrain = FALSE)
    
    for (k in 1:5){
      
      
      
      # remove rare species for focal specific dataset
      X.M<-sum(colSums(fecundity.data[5:dim(fecundity.data)[2]]))
      group.nm.M<-names(which(colSums(fecundity.data[5:dim(fecundity.data)[2]])<=round(X.M*0.01)))
      rares.M<-rowSums(fecundity.data[,which(names(fecundity.data) %in% group.nm.M)])
      fecundity.data.sub<-fecundity.data[,-which(names(fecundity.data) %in% group.nm.M)]
      
      fecundity.data2<-fecundity.data.sub[-flds[[k]],]
      
      # run direct-only model
      alpha.models.AB <- fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=FALSE)
      newdata<-na.omit(fecundity.data[flds[[k]],])
      newdata$fecundity<-newdata$focal
      predict.alpha<-predict.manyglm2(alpha.models.AB,newdata=newdata)
      predict.alpha.2<-predict.alpha
      observed.alpha<-log(na.omit(fecundity.data$seeds[flds[[k]]])+1)
      predict.RMSE.alpha<-sqrt(mean((observed.alpha-predict.alpha.2)^2))
      
      #### fill in summary dataframe
      
      sum$Model[b]<-paste0(d)
      sum$beta[b]<-paste0("F")
      sum$fold[b]<-k
      sum$focal[b]<-paste0(f)
      sum$train.RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
      sum$predict.RMSE[b]<-predict.RMSE.alpha
      sum$coef[b]<-length(alpha.models.AB$coefficients)
      sum$AIC[b]<-AIC(alpha.models.AB)
      sum$BIC[b]<-BIC(alpha.models.AB)
      sum$likeli[b]<-logLik(alpha.models.AB)
      ## increment summary DF
      b<-b+1
      
      # run HOI inclusive model
      #if (beta==T) {
      try(alpha.beta.models.AB <- fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=TRUE),TRUE)
      if ("alpha.beta.models.AB" %in% ls()) {
        newdata<-na.omit(fecundity.data[flds[[k]],])
        newdata$fecundity<-newdata$focal
        predict.beta<-predict.manyglm2(alpha.beta.models.AB,newdata=newdata)
        predict.beta.2<-predict.beta
        observed.beta<-log(na.omit(fecundity.data$seeds[flds[[k]]])+1)
        predict.RMSE.beta<-sqrt(mean((observed.beta-predict.beta.2)^2))
        #}
        
        ## print model statistics to summary DF for beta models
        #if (beta==T){
        sum$Model[b]<-paste0(d)
        sum$beta[b]<-paste0("T")
        sum$fold[b]<-k
        sum$focal[b]<-paste0(f)
        sum$train.RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
        sum$predict.RMSE[b]<-predict.RMSE.beta
        sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
        sum$AIC[b]<-AIC(alpha.beta.models.AB)
        sum$BIC[b]<-BIC(alpha.beta.models.AB)
        sum$likeli[b]<-logLik(alpha.beta.models.AB)
        b<-b+1
      } else {
        sum$Model[b]<-paste0(d)
        sum$beta[b]<-paste0("T")
        sum$fold[b]<-k
        sum$focal[b]<-paste0(f)
        sum$train.RMSE[b]<-NA
        sum$predict.RMSE[b]<-NA
        sum$coef[b]<-NA
        sum$AIC[b]<-NA
        sum$BIC[b]<-NA
        sum$likeli[b]<-NA
        b<-b+1
      }
      
      
      #}
      #}
      rm(alpha.beta.models.AB)
      rm(alpha.models.AB)
      
      
    } # end for each fold
  }
}
# clean and save summary dataframe
sum<-unique(sum)
sum<-sum[-dim(sum)[1],]
write.csv(sum,paste0("Species.grouping/Rare.grouping/output/rare.removed.by.species.kfold.2019.csv"),row.names = FALSE)

print("Rare species removed done")