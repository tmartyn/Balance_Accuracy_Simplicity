source("Neigh.ID/model/fecundity.model.NEW.AB.focal.mvabund.by.species.R")

## beta flag
beta<-T

## start model
for (D in datas){
  
  print(D)
  
  ## read in data
  data.listt<-readRDS(paste0("Neigh.ID/data/",D,".data.neighbor.grouping.list.2019.rds"))
  
  ## list different types
  mode.types<-names(data.listt)
  
  ## list focals for this dataset
  focal<-unique(data.listt[[1]]$focal)
  
  ## create empty lists
  alpha.model.list<-list()
  alpha.beta.model.list<-list()
  
  ## make an empty summary dataframe
  sum<-data.frame(Model=rep("A",500),beta=rep("0",500),focal=rep("O",500), coef=rep(0,500),AIC=rep(0,500),likeli=rep(0,500), RMSE=rep(0,500))
  ## make sure these columns are character
  sum$Model<-as.character(sum$Model)
  sum$beta<-as.character(sum$beta)
  sum$focal<-as.character(sum$focal)
  
  ## DF filler start
  b<-1
  
  ## identify collumns to ignore for some quantificaitons
  if (D=="Australia") {
    cols.ignore<-c("focal","seeds","quadrat","site")  
  } else if (D=="Spain") {
    cols.ignore<-c("focal","site","seeds")
  }
  
  ## for each focal species
  for (f in focal) {
    
    ## for each data type
    for (t in mode.types) {
      
      DF<-data.listt[[t]][which(data.listt[[t]]$focal==f),]
      
      ## remove columns that have all NA
      cols.rm<-which(apply(DF,2, function(x) any(is.na(x)))==TRUE)
      
      if (length(cols.rm)>0) {
        DF<-DF[,-cols.rm]
      }
      
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
      
      ## run the alpha models
      
      alpha.models.AB <- fecundity.model(fecundity.data, fit.alphas=TRUE, fit.betas=FALSE)
      
      ## add alpha models to list
      alpha.model.list[[paste0(t,".",f,".","alpha.model")]]<-alpha.models.AB
      
      ## if running beta moels
      if (beta==T) {
        
        alpha.beta.models.AB <- fecundity.model(fecundity.data, fit.alphas=TRUE, fit.betas=TRUE)
        
      }
      
      ## add beta models to list
      alpha.beta.model.list[[paste0(t,".",f,".","beta.model")]]<-alpha.beta.models.AB
      
      ## print model statistics to summary DF for alpha models
      sum$Model[b]<-paste0(t)
      sum$beta[b]<-paste0("F")
      sum$focal[b]<-paste0(f)
      sum$coef[b]<-length(alpha.models.AB$coefficients)
      sum$RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
      sum$AIC[b]<-AIC(alpha.models.AB)
      sum$likeli[b]<-logLik(alpha.models.AB)
      ## increment summary DF
      b<-b+1
      
      ## print model statistics to summary DF for beta models
      if (beta==T){
        sum$Model[b]<-paste0(t)
        sum$beta[b]<-paste0("T")
        sum$focal[b]<-paste0(f)
        sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
        sum$RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
        sum$AIC[b]<-AIC(alpha.beta.models.AB)
        sum$likeli[b]<-logLik(alpha.beta.models.AB)
        ## increment summary DF
        b<-b+1
        
        rm(alpha.beta.models.AB)
        rm(alpha.models.AB)
        
      }
      ##
      
    } # end for each type t
  } # end for each focal f
  
  ###################
  ## clean summary DF
  sum<-unique(sum)
  sum<-sum[-dim(sum)[1],]
  
  ## write summary csv
  write.csv(sum,paste0("Neigh.ID/output/",D,".neighbor.ID.summary.2019.by.species.csv"))
  
  ## save model lists
  saveRDS(alpha.model.list,paste0("Neigh.ID/output/",D,".alpha.model.list.by.species.rds"))
  saveRDS(alpha.beta.model.list,paste0("Neigh.ID/output/",D,".alpha.beta.model.list.by.species.rds"))
} # close for each data D

print("Neighbor ID analysis done - by species")