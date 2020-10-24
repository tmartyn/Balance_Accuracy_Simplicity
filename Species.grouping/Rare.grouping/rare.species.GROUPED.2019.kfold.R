
# identify leading columns that hold non-species information
X.M<-sum(colSums(Australia.data[5:49]))
# identify neighbor species that have less than 1% of the total number of individuals across all plots
group.nm.M<-names(which(colSums(Australia.data[5:49])<=round(X.M*0.01)))
# summ across rows for rare species 
rares.M<-rowSums(Australia.data[,which(names(Australia.data) %in% group.nm.M)])
# remove rare species from dataframe
Australia.data.sub<-Australia.data[,-which(names(Australia.data) %in% group.nm.M)]
# add rare species column to dataframe
Australia.data.sub$Rare.sp<-rares.M

# identify leading columns that hold non-species information
X.G<-sum(colSums(Spain.data[4:20]))
# identify neighbor species that have less than 1% of the total number of individuals across all plots
group.nm.G<-names(which(colSums(Spain.data[4:20])<=round(X.G*0.01)))
# summ across rows for rare species 
rares.G<-rowSums(Spain.data[,which(names(Spain.data) %in% group.nm.G)])
# remove rare species from dataframe
Spain.data.sub<-Spain.data[,-which(names(Spain.data) %in% group.nm.G)]
# add rare species column to dataframe
Spain.data.sub$Rare.sp<-rares.G

# make data list
data.list.W<-list()
data.list.W[["Spain"]]<-Spain.data.sub
data.list.W[["Australia"]]<-Australia.data.sub
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
write.csv(sum,paste0("Species.grouping/Rare.grouping/output/rare.grouped.summary.kfolds.2019.csv"),row.names = FALSE)


print("Rare species grouped done")