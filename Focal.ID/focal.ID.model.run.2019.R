# list focal identity types
type<-c("NA","NA.part","A","AB")

for (d in datas) {
  # read the dataframe
  print(d)
  fecundity.data<-data.list[[d]]
  
  # create an empty dataframe to catch output
  sum<-data.frame(Model=rep(100,5000),beta=rep("F",5000),focal=as.character(rep("O",5000)), coef=rep(0,5000),AIC=rep(0,5000),RMSE=rep(0,5000))
  sum$focal<-as.character(sum$focal)
  sum$beta<-as.character(sum$beta)
  b<-1
  
  for (t in type) {
    
    # make empty model lists
    alpha.models.AB<-list()
    alpha.beta.models.AB<-list()
    
    #print(paste0("Testing Focal: ",t))
    
    # source model (this is the best fitting of the HOI models so far...)
    source(paste0("Focal.ID/model/fecundity.model.NEW.",t,".focal.mvabund.R"))
    
    # if the model type is the ones below - run the direct-only ("alpha") model
    if (t %in% c("A","NA","NA.part")) {
      #   print(paste0("Running model alpha ",t," for data ",d))
      alpha.models.AB<-fecundity.model(fecundity.data, fit.alphas=TRUE, fit.betas=FALSE)
    } 
    
    # if the model includes HOIs - run the HOI-inclusive ("alpha.beta") model
    if (t %in% c("A","AB","B","NA","NA.part")) {
      # print(paste0("Running model alpha beta ",t," for data ",d))
      alpha.beta.models.AB <- fecundity.model(fecundity.data, fit.alphas=TRUE, fit.betas=TRUE)
    }
    
    #print("Saving data")
    
    # save model lists as output
    if (t %in% c("A","NA","NA.part")) {
      saveRDS(alpha.models.AB,paste0("Focal.ID/output/",d,".mvabund.",t,".alpha.model.09.09.2018.rds"))
    }
    
    if (t %in% c("A","AB","B","NA","NA.part")) {
      saveRDS(alpha.beta.models.AB,paste0("Focal.ID/output/",d,".mvabund.",t,".alpha.beta.model.09.09.2018.rds"))
    }
    
    # fill in empty dataframe created above
    if (t %in% c("A","NA","NA.part")) {
      sum$Model[b]<-paste0(t)
      sum$coef[b]<-length(alpha.models.AB$coefficients)
      sum$AIC[b]<-AIC(alpha.models.AB)
      sum$RMSE[b]<-sqrt(mean(alpha.models.AB$residuals^2))
      sum$likeli[b]<-logLik(alpha.models.AB)
      sum$beta[b]<-"FALSE"
      b<-b+1}
    
    if (t %in% c("A","AB","B","NA","NA.part")) {
      sum$Model[b]<-paste0(t)
      sum$coef[b]<-length(alpha.beta.models.AB$coefficients)
      sum$AIC[b]<-AIC(alpha.beta.models.AB)
      sum$RMSE[b]<-sqrt(mean(alpha.beta.models.AB$residuals^2))
      sum$likeli[b]<-logLik(alpha.beta.models.AB)
      sum$beta[b]<-"TRUE"
      b<-b+1}
    
    
    rm(alpha.beta.models.AB)
    rm(alpha.models.AB)
  }
  
  # clean up the filled-in dataframe
  sum<-unique(sum)
  sum<-sum[-dim(sum)[1],]
  
  # write dataframe as output csv
  write.csv(sum,paste0("Focal.ID/output/",d,".focal.id.model.summary.2019.csv"))
  
}

print("Focal ID analysis done")