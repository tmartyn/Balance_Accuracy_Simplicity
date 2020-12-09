
source("Neigh.ID/model/fecundity.model.NEW.AB.focal.mvabund.by.focal.R")
source("predict.manyglm.R")
## beta flag
beta<-F

## start model

for (D in datas){
  
  print(D)
  
  ## read in data
  data.listt<-readRDS(paste0("Neigh.ID/data/",D,".data.neighbor.grouping.list.2019.rds"))
  
  ## list different types
  mode.types<-names(data.listt)
  
  ## create empty lists
  alpha.model.list<-list()
  alpha.beta.model.list<-list()
  
  ## make an empty summary dataframe
  sum<-data.frame(Model=rep(100,5000),beta=rep("F",5000),fold=rep(0,5000),
                  focal=as.character(rep("O",5000)), coef=rep(0,5000),
                  AIC=rep(0,5000),train.RMSE=rep(0,5000),
                  predict.RMSE=rep(0,5000),BIC=rep(0,5000))
  sum$focal<-as.character(sum$focal)
  sum$beta<-as.character(sum$beta)
  
  
  ## DF filler start
  b<-1
  
  ## identify collumns to ignore for some quantificaitons
  if (D=="Australia") {
    cols.ignore<-c("focal","seeds","quadrat","site")  
  } else if (D=="Spain") {
    cols.ignore<-c("focal","site","seeds")
  }
  
  focals<-unique(data.listt[[1]]$focal)
  
  ## for each focal species
  for (f in focals) {
    
    for (t in mode.types) {
      
      DF<-data.listt[[t]][which(data.listt[[t]]$focal==f),]
      
      ## special case for rm.rares dataset
      if (t=="rm.rares") {
        DF<-DF[,-which(names(DF)=="rare.count")]
      }
      DF.sp<-DF[,-which(names(DF) %in% cols.ignore)]
      
      ## remove columns that have all 0's
      if (is.null(dim(data.frame(DF.sp)))) {
        
        rm.names<-names(DF.sp[which(colSums(DF.sp)==0)])
        
        if (length(rm.names)>0) {
          DF<-DF[,-which(names(DF)%in%rm.names)]
        } else {DF<-DF}
        
      }
      
      ## identify fecundity.data dataset
      fecundity.data<-DF
      fecundity.data[is.na(fecundity.data)]<-0
      
      data.length<-dim(fecundity.data)[1]
      
      require(caret)
      flds <- createFolds(seq(1:data.length), k = 5, list = TRUE, returnTrain = FALSE)
      
      for (k in 1:5){
        
        fecundity.data2<-fecundity.data[-flds[[k]],]
        
        if (beta==F) {
        ## run the alpha models
        alpha.models.AB <-fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=FALSE)
        newdata<-fecundity.data[flds[[k]],]
        newdata$fecundity<-newdata$focal
        predict.alpha<-predict.manyglm2(alpha.models.AB,newdata=newdata)
        predict.alpha.2<-predict.alpha
        observed.alpha<-log(fecundity.data$seeds[flds[[k]]]+1)
        predict.RMSE.alpha<-sqrt(mean((observed.alpha-predict.alpha.2)^2))
        
        saveRDS(alpha.models.AB,paste0("Neigh.ID/output/by.focal/",D,".alpha.model.list.direct.",t,"fold.",k,"focal.",f,".2019.rds"))
        
        ## print model statistics to summary DF for alpha models
        sum$Model[b]<-paste0(t)
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
        }
        ## add alpha models to list
        
        ## if running beta moels
        if (beta==T) {
          
          ## run beta models
          alpha.beta.models.AB <-fecundity.model(fecundity.data2, fit.alphas=TRUE, fit.betas=TRUE)
          newdata<-fecundity.data[flds[[k]],]
          newdata$fecundity<-newdata$focal
          predict.beta<-predict.manyglm2(alpha.beta.models.AB,newdata=newdata)
          predict.beta.2<-predict.beta
          observed.beta<-log(fecundity.data$seeds[flds[[k]]]+1)
          predict.RMSE.beta<-sqrt(mean((observed.beta-predict.beta.2)^2))
          
          
          # end run beta models 
          
          saveRDS(alpha.beta.models.AB,paste0("Neigh.ID/output/by.focal/",D,".alpha.beta.model.list.direct.",t,"fold.",k,"focal.",f,".2019.rds"))  
          ## add beta models to list
          
          
          ## print model statistics to summary DF for beta models
            sum$Model[b]<-paste0(t)
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
        }
        
 
        
 
          
          
        
        ##
      }# end for each kfold
    } # end for each type t
    
  } # end for each focal
  
  ###################
  ## clean summary DF
  sum<-unique(sum)
  sum<-sum[-dim(sum)[1],]
  
  ## write summary csv
  write.csv(sum,paste0("Neigh.ID/output/by.focal/",D,".neighbor.ID.summary.2019.direct.kfolds.by.focal.csv"))
  
} # close for each data D

print("Neighbor ID analysis done - direct-only")