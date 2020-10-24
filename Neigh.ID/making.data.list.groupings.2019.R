# read in trait and cluster groupings
TR<-read.csv("Species.grouping/Functional.trait.grouping/output/traits.cluster.all.results.2019.csv")

################### make large dataset list ###########
# for each dataset
for (d in datas) {
  print(d)
  
  fecundity.data<-data.list[[d]]
  # identitfy unique focal species 
  focals.list<-unique(fecundity.data$focal)
  
  #####################################################################
  # create large list for dataframes
  data.list.groups<-list()
  data.list.groups[["full"]]<-fecundity.data
  
  ############################################################
  ### count number of competitors and add to dataframe ####
  if (d=="Australia") {
    cols.ignore<-c("focal","seeds","quadrat","site")
    
  } else if (d=="Spain") {
    
    cols.ignore<-c("focal","site","seeds")
  }

  
  # pull out the non-species columns in the dataframe
  DF.starter<-fecundity.data[,cols.ignore]
  # pull out the species-only columns in the dataframe
  species.df<-fecundity.data[,-which(names(fecundity.data) %in% cols.ignore)]
  
  
  DF.starter$total.num.competitor<-rowSums(species.df)
  # for each row of the dataframe
  
  data.list.groups[["number.indivs"]]<-DF.starter[,c(cols.ignore,"total.num.competitor")]
  
  ################################################################
  
  #################################################
  fecundity.data2<-fecundity.data
  
  
  #######################################################################
  ### exoctics vs natives #############################
  
  # get list of species' identity of native vs exotic
  exotic<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="exotic")]
  native<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="native")]
  
  DF.starter<-fecundity.data[,cols.ignore]
  DF.starter$exotic<-0
  DF.starter$native<-0
  
  # assign and summ exotic or native individuals
  for (b in 1:dim(fecundity.data2)[1]) {
    
    if (length(exotic)>0) {
      ex.sum<-sum(fecundity.data2[b,which(names(fecundity.data2) %in% exotic)])
    } else {ex.sum<-0}
    if (length(native)>0) {
      nt.sum<-sum(fecundity.data2[b,which(names(fecundity.data2) %in% native)])
    } else {nt.sum<-0}
    
    DF.starter$exotic[b]<-ex.sum
    DF.starter$native[b]<-nt.sum
  }
  
  data.list.groups[["native.exotic"]]<-DF.starter
  
  ######################################################################################
  
  ###############################################################
  ########## grass vs forb groups #################################
  
  # get list of species assigned grass or forb
  grass<-TR$species.dataset.code[which(TR$Data==d & TR$grass_forb=="grass")]
  forb<-TR$species.dataset.code[which(TR$Data==d & TR$grass_forb=="forb")]
  
  DF.starter<-fecundity.data[,cols.ignore]
  DF.starter$grass<-0
  DF.starter$forb<-0
  
  # assign and sum grass or forb individuals
  for (b in 1:dim(fecundity.data2)[1]) {
    
    gr.sum<-sum(fecundity.data2[b,which(names(fecundity.data2) %in% grass)])
    fb.sum<-sum(fecundity.data2[b,which(names(fecundity.data2) %in% forb)])
    
    DF.starter$grass[b]<-gr.sum
    DF.starter$forb[b]<-fb.sum
  }
  
  data.list.groups[["grass.forb"]]<-DF.starter
  
  ######################################################################################
  
  
  ###############################################################
  ########## functoinal groups X exotic native #################
  func.group.list<-list()
  
  # assign species to functional groups
  func.group.list[["native.p.forb"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="native" & TR$annual_perennial=="perennial" & TR$grass_forb=="forb")]
  func.group.list[["native.p.grass"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="native" & TR$annual_perennial=="perennial" & TR$grass_forb=="grass")]
  func.group.list[["exotic.a.forb"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="exotic" & TR$annual_perennial=="annual" & TR$grass_forb=="forb")]
  func.group.list[["exotic.a.grass"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="exotic" & TR$annual_perennial=="annual" & TR$grass_forb=="grass")]
  func.group.list[["native.a.forb"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="native" & TR$annual_perennial=="annual" & TR$grass_forb=="forb")]
  func.group.list[["native.a.grass"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="native" & TR$annual_perennial=="annual" & TR$grass_forb=="grass")]
  func.group.list[["exotic.p.forb"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="exotic" & TR$annual_perennial=="perennial" & TR$grass_forb=="forb")]
  func.group.list[["exotic.p.grass"]]<-TR$species.dataset.code[which(TR$Data==d & TR$native_extoic=="exotic" & TR$annual_perennial=="perennial" & TR$grass_forb=="grass")]
  
  # make empty vector for groups that have no species
  ele.rm<-vector()
  
  # remove groups with no species
  for (ii in 1:length(func.group.list)) {
    if (length(func.group.list[[ii]])==0){
      ele.rm<-append(ele.rm, ii)
    } else {next}
  }
  func.group.list<-func.group.list[-ele.rm]
  
  # list types of functional groups
  type.func<-names(func.group.list)
  
  # remove type of group that occurs exactly once in all the data (too rare to get an accurate result for this particular functional group)
  if(d=="Australia") {type.func<-type.func[-which(type.func=="native.p.forb")]}
  
  # create starter dataframe to full
  DF.starter<-fecundity.data[,cols.ignore]
  for (tf in type.func) {
    DF.starter[,paste0(tf)]<-0
  }
  
  
  # fill in the dataframe summing for each functional group
  for (tfp in type.func) {
    
    names.f<-func.group.list[[tfp]]
    for (b in 1:dim(fecundity.data2)[1]) {
      
      sum<-sum(fecundity.data2[b,which(names(fecundity.data2) %in% names.f)])
      
      DF.starter[[tfp]][b]<-sum
    }
  }
  
  data.list.groups[["functional"]]<-DF.starter
  
  ###############################################3
  
  
  ###############################################################
  ########## families #################
  
  fam.list<-list()
  families<-unique(TR$family[which(TR$Data==d)])
  
  # remove specific families of neighbor species that have been since cleaned from the datasets
  if (d=="Spain") {families<-families[-which(families=="Brassicaeae")]}
  if (d=="Australia") {families<-families[-which(families=="Carophyllaceae")]}
  
  # create dataframe to fill
  DF.starter<-fecundity.data[,cols.ignore]
  for (fm in families) {
    fam.list[[fm]]<-TR$species.dataset.code[which(TR$family==fm & TR$Data==d)]
    DF.starter[,paste0(fm)]<-0
  }
  
  temp.df<-fecundity.data2[,-which(names(fecundity.data2)%in%cols.ignore)]
  
  # summ and fill dataframe
  for (b in 1:dim(temp.df)[1]) {
    
    for (fm in families) {
      
      fm.vec<-fam.list[[fm]]
      
      count.col<-vector()
      
      fam.count.list<-list()
      
      fm.count<-sum(temp.df[b,which(names(temp.df) %in% fm.vec)])
      
      DF.starter[b,which(names(DF.starter)==fm)]<-fm.count
    }
  }
  
  data.list.groups[["families"]]<-DF.starter
  
  #####################################
  
  
  #####################################
  ### funcational groups ##############
  list.type<-list()
  
  # read in trait cluster
  nm.all<-unique(TR$cluster.all[which(TR$Data==d)])
  
  list.type[["all"]]<-nm.all
  
  type<-names(list.type)
  
  # sum and fill in number of individuals for each cluster
  for (tp in type) {
    
    clust<-list.type[[tp]]
    col.name<-paste0("cluster.",tp)
    cols<-c(col.name,"species.dataset.code")
    DF.starter<-fecundity.data[,cols.ignore]
    
    for (cl in clust) {
      FF<-TR[which(TR$Data==d),cols]
      DF.starter[[paste0(tp,".",cl)]]<-0
      
      for (b in 1:dim(fecundity.data2)[1]) {
        tp.spec<-FF$species.dataset.code[which(FF[[col.name]]==cl)]
        tp.count<-sum(fecundity.data2[b,which(names(fecundity.data2) %in% tp.spec)])
        DF.starter[[paste0(tp,".",cl)]][b]<-tp.count
      }
    }
  }
  
  data.list.groups[["cluster"]]<-DF.starter
  ###############################################################
  
  # save data list to use in neighbor species abundance
  saveRDS(data.list.groups,paste0("Neigh.ID/data/",d,".data.neighbor.grouping.list.2019.rds"))
  
}
print("Neighbor grouping data collating done")
