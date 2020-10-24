# type will help identify the column of cluster (choices are "quant","leaves","all") for quantitiative only, qunatiative with leaves, and all traits
find.clusters<-function(traits,remove,type) {
  
  ## Spain trait data
  data.Spain<-traits[which(traits$Data=="Spain"),]
  
  ## remove data in which there are NA's 
  data.Spain.2 <- data.Spain[,colSums(is.na(data.Spain))<nrow(data.Spain)]
  col.rm.G<-c(remove)
  data.Spain.3 <- data.Spain.2[,-which(names(data.Spain.2)%in% col.rm.G)]
  data.Spain.3<-data.Spain.3[-which(names(data.Spain.3)=="X")]
  
  ## Australia trait data
  data.Australia<-traits[which(traits$Data=="Australia"),]
  
  ## remove data in which there are NA's 
  data.Australia.2 <- data.Australia[,colSums(is.na(data.Australia))<nrow(data.Australia)]
  col.rm.May<-c("leaf.area.for.sla","leaf.dry.mass.for.sla")
  col.rm.M<-c(remove,col.rm.May)
  data.Australia.3 <- data.Australia.2[,-which(names(data.Australia.2)%in% col.rm.M)]
  data.Australia.3<-na.omit(data.Australia.3)
  data.Australia.3<-data.Australia.3[-which(names(data.Australia.3)=="X")]
  
  #######################################################################################
  #### cluster analysis Spain.data #####
  
  ## Calculate gower distances
  gower_dist.G <- daisy(data.Spain.3,
                        metric = "gower",
                        type = list(logratio = 3))
  
  summary(gower_dist.G)
  
  
  ##################################
  # change distance lower triangler to matrix
  gower_mat <- as.matrix(gower_dist.G)
  
  # Output most similar pair
  data.Spain.3[
    which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
          arr.ind = TRUE)[1, ], ]
  
  ############################
  # Calculate silhouette width for many clusters using PAM
  sil_width <- vector()
  num_clusters<-2:15 
  
  # for the number of clusters
  for(i in num_clusters){
    
    pam_fit <- pam(gower_dist.G,
                   diss = TRUE,
                   k = i)
    
    sil_width[i] <- pam_fit$silinfo$avg.width
    
  }
  
  sil_width<-na.omit(sil_width)
  sils<-data.frame(cbind(num_clusters,sil_width))
  
  # Plot sihouette width (higher is better) 
  
  plot(sils[,1], sils[,2],
       xlab = "Number of clusters",
       ylab = "Silhouette Width", type="b", col="blue",ylim=c(0,.70))
  
  max_clust.G<-sils$num_clusters[which(sils[,2]==max(sils[,2]))]
  
  
  ###################################
  pam_fit.G <- pam(gower_dist.G, diss = TRUE, k = max_clust.G)
  
  pam_results.G <- data.Spain.3 %>%
    mutate(cluster = pam_fit.G$clustering) %>%
    group_by(cluster) %>%
    do(the_summary = summary(.))
  
  clust.names<-paste0("G",pam_fit.G$clustering)
  data.clust.Spain<-cbind(data.Spain,clust.names)
  colnames(data.clust.Spain) [colnames(data.clust.Spain)=="clust.names"] <-paste0("cluster.",type)
  ##########################
  
  
  ####################################################################################
  #### cluster analysis Australia.data #####
  
  ## Calculate gower distances
  gower_dist <- daisy(data.Australia.3,
                      metric = "gower",
                      type = list(logratio = 3))
  
  summary(gower_dist)
  
  
  ##################################
  gower_mat <- as.matrix(gower_dist)
  
  # Output most similar pair
  data.Australia.3[
    which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
          arr.ind = TRUE)[1, ], ]
  
  ############################
  # Calculate silhouette width for many k using PAM
  
  sil_width <- vector()
  num_clusters<-2:15 
  
  for(i in num_clusters){
    
    pam_fit <- pam(gower_dist,
                   diss = TRUE,
                   k = i)
    
    sil_width[i] <- pam_fit$silinfo$avg.width
    
  }
  
  sil_width<-na.omit(sil_width)
  sils<-data.frame(cbind(num_clusters,sil_width))
  
  # Plot sihouette width (higher is better) 
  points(sils[,1], sils[,2])
  lines(sils[,1], sils[,2])
  
  legend("topleft",col=c("black","blue"),legend=c("Australia","Spain"), pch=c(1,1))
  
  
  # find what cluster value is the highest
  max_clust<-sils$num_clusters[which(sils[,2]==max(sils[,2]))]
  
  ###################################
  pam_fit <- pam(gower_dist, diss = TRUE, k = max_clust)
  pam_results <- data.Australia.3 %>%
    mutate(cluster = pam_fit$clustering) %>%
    group_by(cluster) %>%
    do(the_summary = summary(.))
  
  clust.names<-paste0("M",pam_fit$clustering)
  data.clust.Australia<-cbind(data.Australia,clust.names)
  
  colnames(data.clust.Australia) [colnames(data.clust.Australia)=="clust.names"] <-paste0("cluster.",type)
  
  # make sliter plot for each dataset
  clusplot(pam_fit.G, main="Spain dataset", lines=0, color=T)
  clusplot(pam_fit, main="Australia dataset", lines=0, color=T)
  
  # make grouped dataframe
  traits.clustered<-rbind.fill(data.clust.Australia,data.clust.Spain)
  
  # return dataframe
  return(traits.clustered)
}
#########################
