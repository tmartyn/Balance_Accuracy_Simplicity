# calculate Pielou eveness

AUS<-readRDS("Neigh.ID/data/Australia.data.neighbor.grouping.list.2019.rds")
SPN<-readRDS("Neigh.ID/data/Spain.data.neighbor.grouping.list.2019.rds")

AUS.full<-AUS[["full"]]
SPN.full<-SPN[["full"]]

cols.rm<-c("focal","seeds","quadrat","site")

### Australia evenness
AUS.species.num<-length(names(AUS.full)[-which(names(AUS.full)%in%cols.rm)])
AUS.H.max<-log(AUS.species.num)
AUS.temp<-AUS.full[-which(names(AUS.full)%in% cols.rm)]
species.count<-colSums(AUS.temp)
sum.sp.ct<-sum(species.count)
species.prop<-species.count/sum.sp.ct
AUS.H<--1*sum(species.prop*log(species.prop))
AUS.H/AUS.H.max #more even than SPN

## by site
# site K (Kunjin)
AUS.K<-AUS.full[AUS.full$site=="K",]
cols.K.rm<-names(AUS.K[5:49])[which(colSums(AUS.K[,5:49])==0)]
AUS.K<-AUS.K[,-which(names(AUS.K)%in%cols.K.rm)]



AUS.K.species.num<-length(names(AUS.K)[-which(names(AUS.K)%in%cols.rm)])
AUS.K.H.max<-log(AUS.K.species.num)
AUS.K.temp<-AUS.K[-which(names(AUS.K)%in% cols.rm)]
species.count.K<-colSums(AUS.K.temp)
sum.sp.ct.K<-sum(species.count.K)
species.prop.K<-species.count.K/sum.sp.ct.K
AUS.H.K<--1*sum(species.prop.K*log(species.prop.K))
AUS.H.K/AUS.K.H.max #more even than SPN

# site B (Bendering)
AUS.B<-AUS.full[AUS.full$site=="B",]
cols.B.rm<-names(AUS.B[5:49])[which(colSums(AUS.B[,5:49])==0)]
AUS.B<-AUS.B[,-which(names(AUS.B)%in%cols.B.rm)]


AUS.B.species.num<-length(names(AUS.B)[-which(names(AUS.B)%in%cols.rm)])
AUS.B.H.max<-log(AUS.B.species.num)
AUS.B.temp<-AUS.B[-which(names(AUS.B)%in% cols.rm)]
species.count.B<-colSums(AUS.B.temp)
sum.sp.ct.B<-sum(species.count.B)
species.prop.B<-species.count.B/sum.sp.ct.B
AUS.H.B<--1*sum(species.prop.B*log(species.prop.B))
AUS.H.B/AUS.B.H.max #more even than SPN


########### spain eveness
SPN.species.num<-length(names(SPN.full)[-which(names(SPN.full)%in%cols.rm)])
SPN.H.max<-log(SPN.species.num)
SPN.temp<-SPN.full[-which(names(SPN.full)%in% cols.rm)]
species.count<-colSums(SPN.temp)
sum.sp.ct<-sum(species.count)
species.prop<-species.count/sum.sp.ct
SPN.H<--1*sum(species.prop*log(species.prop))
SPN.H/SPN.H.max #less even than AUS

# driven by ~15000 Hordeum individuals (at least 10X more individuals than the next most abundant)