source("predict.manyglm.R")
# list focal identity types
type<-c("NA","NA.part","A","AB")

for (d in datas) {
  # read the dataframe
  print(d)
  fecundity.data<-data.list[[d]]
  
  
  data.length<-dim(fecundity.data)[1]
  
  require(caret)
  flds <- createFolds(seq(1:data.length), k = 5, list = TRUE, returnTrain = FALSE)
  
  # create an empty dataframe to catch output
  sum<-data.frame(Model=rep(100,5000),beta=rep("F",5000),fold=rep(0,5000),
                  coef=rep(0,5000),
                  AIC=rep(0,5000),train.RMSE=rep(0,5000),
                  predict.RMSE=rep(0,5000),BIC=rep(0,5000))
  #sum$focal<-as.character(sum$focal)
  sum$beta<-as.character(sum$beta)
  b<-1
  
  for (k in 1:5) {
    print(k)
    fecundity.data2<-fecundity.data[-flds[[k]],]
    
    for (t in type) {
      print(t)

      # source model (this is the best fitting of the HOI models so far...)
      source(paste0("Focal.ID/model/fecundity.model.NEW.",t,".focal.mvabund.R"))
      
      # if the model type is the ones below - run the direct-only ("alpha") model
      if (t %in% c("A","NA","NA.part")) {
        #   print(paste0("Running model alpha ",t," for data ",d))
        try(alpha.models.AB<-fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=FALSE),TRUE)
        if ("alpha.models.AB" %in% ls()) {
          newdata<-fecundity.data[flds[[k]],]
          newdata$fecundity<-newdata$focal
          predict.alpha<-predict.manyglm2(alpha.models.AB,newdata=newdata)
          predict.alpha.2<-predict.alpha
          observed.alpha<-log(fecundity.data$seeds[flds[[k]]]+1)
          predict.RMSE.alpha<-sqrt(mean((observed.alpha-predict.alpha.2)^2))
          sum$Model[b]<-paste0(t)
          sum$coef[b]<-length(alpha.models.AB$coefficients)
          sum$AIC[b]<-AIC(alpha.models.AB)
          sum$fold[b]<-k
          sum$BIC[b]<-BIC(alpha.models.AB)
          sum$train.RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
          sum$predict.RMSE[b]<-predict.RMSE.alpha
          sum$likeli[b]<-logLik(alpha.models.AB)
          sum$beta[b]<-"FALSE"
          b<-b+1
          
          saveRDS(alpha.models.AB,
                  paste0("Focal.ID/output/",d,".mvabund.",t,".alpha.model.fold.",k,".24.10.2020.rds"))
          
          rm(alpha.models.AB)
          
        } else {     
          print("Direct-only model did not fit this fold")
          sum$Model[b]<-paste0(t)
          sum$coef[b]<-NA
          sum$AIC[b]<-NA
          sum$fold[b]<-k
          sum$BIC[b]<-NA
          sum$train.RMSE[b]<-NA
          sum$predict.RMSE[b]<-NA
          sum$likeli[b]<-NA
          sum$beta[b]<-"FALSE"
          b<-b+1
          rm(alpha.models.AB)
        }
      } 
      
      # if the model includes HOIs - run the HOI-inclusive ("alpha.beta") model
      if (t %in% c("A","AB","B","NA","NA.part")) {
        # print(paste0("Running model alpha beta ",t," for data ",d))
        try(alpha.beta.models.AB <- fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=TRUE),TRUE)
        if ("alpha.beta.models.AB" %in% ls()) {
          #predict.manyglm2(alpha.beta.models.AB,new.data=fecundity.data[flds[[1]],])
          newdata<-fecundity.data[flds[[k]],]
          newdata$fecundity<-newdata$focal
          try(predict.beta<-predict.manyglm2(alpha.beta.models.AB,newdata=newdata),TRUE)
          if ("predict.beta" %in% ls()) {
          predict.beta.2<-predict.beta
          observed.beta<-log(fecundity.data$seeds[flds[[k]]]+1)
          predict.RMSE.beta<-sqrt(mean((observed.beta-predict.beta.2)^2))
          rm(predict.beta)
          } else {
            predict.RMSE.beta<-0
          }
          sum$Model[b]<-paste0(t)
          sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
          sum$AIC[b]<-AIC(alpha.beta.models.AB)
          sum$fold[b]<-k
          sum$BIC[b]<-BIC(alpha.beta.models.AB)
          sum$train.RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
          sum$predict.RMSE[b]<-predict.RMSE.beta
          sum$likeli[b]<-logLik(alpha.beta.models.AB)
          sum$beta[b]<-"TRUE"
          b<-b+1
          
          saveRDS(alpha.beta.models.AB,
                  paste0("Focal.ID/output/",d,".mvabund.",t,".alpha.beta.model.fold.",k,".24.10.2020.rds"))
          
          rm(alpha.beta.models.AB)
          
        } else { 
          print("HOI-inclusive model did not fit this fold")
          sum$Model[b]<-paste0(t)
          sum$coef[b]<-NA
          sum$AIC[b]<-NA
          sum$fold[b]<-k
          sum$BIC[b]<-NA
          sum$train.RMSE[b]<-NA
          sum$predict.RMSE[b]<-NA
          sum$likeli[b]<-NA
          sum$beta[b]<-"TRUE"
          b<-b+1
          rm(alpha.beta.models.AB)
        }
      }
      

    } # end for each type
    
    
  }# end for each kfold
  
  # clean up the filled-in dataframe
  sum<-unique(sum)
  sum<-sum[-dim(sum)[1],]
  
  # write dataframe as output csv
  write.csv(sum,paste0("Focal.ID/output/",d,".focal.id.model.summary.kfolds.2020.csv"))
  
  
}# end for each data
print("Focal ID k-folds analysis done")
