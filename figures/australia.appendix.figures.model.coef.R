
################# Aussie Native exotic ###############################
alpha.model.list<-list()
alpha.beta.model.list<-list()

rm(df.all)

for (i in 1:5){
  alpha.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.model.list.native.exoticfold.",i,".2019.rds"))
  alpha.beta.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.beta.model.list.native.exoticfold.",i,".2019.rds"))
  alpha.coef<-coef(alpha.model.list[[i]])
  alpha.beta.coef<-coef(alpha.beta.model.list[[i]])[1:(length(alpha.coef))]
  temp.df<-cbind(alpha.coef,alpha.beta.coef)
  temp.df<-data.frame(temp.df[-grep("fecundity",rownames(temp.df)),])
  temp.df$interaction<-rownames(temp.df)
  if ("df.all" %in% ls()) {
    df.all<-rbind(df.all,temp.df)
  } else {df.all<-temp.df}
}


std <- function(x) sd(x)/sqrt(length(x))

df.mean<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),mean))
df.std<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),std))
    
df.mean$alpha.beta.std<-df.std$alpha.beta.coef 
df.mean$alpha.std<-df.std$seeds

names(df.mean)<-c("interaction","alpha.beta.mean","alpha.mean","alpha.beta.std","alpha.std")

library(ggplot2)
library(Cairo)
CairoPDF("figures/Australia_Orgin_Status.pdf",height=10,width=10)
colors<-rep(c("red","blue"),6)
ggplot(data=df.mean,aes(x=alpha.mean,y=alpha.beta.mean,col=interaction))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1,size=1.5)+
  ylim(-0.3,0.3)+xlim(-0.3,0.3)+
  geom_errorbarh(aes(xmin=alpha.mean-2*alpha.std,xmax=alpha.mean+2*alpha.std))+
  geom_errorbar(aes(ymin=alpha.beta.mean-2*alpha.beta.std,ymax=alpha.beta.mean+2*alpha.beta.std)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia Origin Status")+
  geom_vline(xintercept=0, col="grey60")+
  geom_hline(yintercept=0, col="grey60")+
  scale_color_manual(values=colors)+
  guides(col=F)+
  theme_minimal(base_size = 20)+theme(panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank())
dev.off()
################# Aussie grass forb #####################
alpha.model.list<-list()
alpha.beta.model.list<-list()

rm(df.all)

for (i in 1:5){
  alpha.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.model.list.grass.forbfold.",i,".2019.rds"))
  alpha.beta.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.beta.model.list.grass.forbfold.",i,".2019.rds"))
  alpha.coef<-coef(alpha.model.list[[i]])
  alpha.beta.coef<-coef(alpha.beta.model.list[[i]])[1:(length(alpha.coef))]
  temp.df<-cbind(alpha.coef,alpha.beta.coef)
  temp.df<-data.frame(temp.df[-grep("fecundity",rownames(temp.df)),])
  temp.df$interaction<-rownames(temp.df)
  if ("df.all" %in% ls()) {
    df.all<-rbind(df.all,temp.df)
  } else {df.all<-temp.df}
}

std <- function(x) sd(x)/sqrt(length(x))

df.mean<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),mean))
df.std<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),std))

df.mean$alpha.beta.std<-df.std$alpha.beta.coef 
df.mean$alpha.std<-df.std$seeds

names(df.mean)<-c("interaction","alpha.beta.mean","alpha.mean","alpha.beta.std","alpha.std")


CairoPDF("figures/Australia_Life_Form.pdf",height=10,width=10)
colors<-rep(c("red","blue"),6)
ggplot(data=df.mean,aes(x=alpha.mean,y=alpha.beta.mean,col=interaction))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1,size=1.5)+
  #facet_wrap(vars(Metric.name), scales="free")+
  # geom_blank(data=dummy)+
  ylim(-0.3,0.3)+xlim(-0.3,0.3)+
  geom_errorbarh(aes(xmin=alpha.mean-2*alpha.std,xmax=alpha.mean+2*alpha.std))+
  geom_errorbar(aes(ymin=alpha.beta.mean-2*alpha.beta.std,ymax=alpha.beta.mean+2*alpha.beta.std)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia Life Form")+
  geom_vline(xintercept=0, col="grey60")+
  geom_hline(yintercept=0, col="grey60")+
  scale_color_manual(values=colors)+
  guides(col=F)+
  theme_minimal(base_size = 20)+theme(panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank())
dev.off()


################# Aussie
alpha.model.list<-list()
alpha.beta.model.list<-list()

rm(df.all)

for (i in 1:5){
  alpha.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.model.list.number.indivsfold.",i,".2019.rds"))
  alpha.beta.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.beta.model.list.number.indivsfold.",i,".2019.rds"))
  alpha.coef<-coef(alpha.model.list[[i]])
  alpha.beta.coef<-coef(alpha.beta.model.list[[i]])[1:(length(alpha.coef))]
  temp.df<-cbind(alpha.coef,alpha.beta.coef)
  temp.df<-data.frame(temp.df[-grep("fecundity",rownames(temp.df)),])
  temp.df$interaction<-rownames(temp.df)
  if ("df.all" %in% ls()) {
    df.all<-rbind(df.all,temp.df)
  } else {df.all<-temp.df}
}


std <- function(x) sd(x)/sqrt(length(x))

df.mean<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),mean))
df.std<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),std))

df.mean$alpha.beta.std<-df.std$alpha.beta.coef 
df.mean$alpha.std<-df.std$seeds

names(df.mean)<-c("interaction","alpha.beta.mean","alpha.mean","alpha.beta.std","alpha.std")

CairoPDF("figures/Australia_Abundance.pdf",height=10,width=10)
colors<-rep(c("red"),6)
ggplot(data=df.mean,aes(x=alpha.mean,y=alpha.beta.mean,col=interaction))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1,size=1.5)+
  ylim(-0.3,0.3)+xlim(-0.3,0.3)+
  geom_errorbarh(aes(xmin=alpha.mean-2*alpha.std,xmax=alpha.mean+2*alpha.std))+
  geom_errorbar(aes(ymin=alpha.beta.mean-2*alpha.beta.std,ymax=alpha.beta.mean+2*alpha.beta.std)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia Abundance")+
  geom_vline(xintercept=0, col="grey60")+
  geom_hline(yintercept=0, col="grey60")+
  scale_color_manual(values=colors)+
  guides(col=F)+
  theme_minimal(base_size = 20)+theme(panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank())
dev.off()


################# Aussie fucntional type
alpha.model.list<-list()
alpha.beta.model.list<-list()

rm(df.all)

for (i in 1:5){
  alpha.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.model.list.clusterfold.",i,".2019.rds"))
  alpha.beta.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Australia.alpha.beta.model.list.clusterfold.",i,".2019.rds"))
  alpha.coef<-coef(alpha.model.list[[i]])
  alpha.beta.coef<-coef(alpha.beta.model.list[[i]])[1:(length(alpha.coef))]
  temp.df<-cbind(alpha.coef,alpha.beta.coef)
  temp.df<-data.frame(temp.df[-grep("fecundity",rownames(temp.df)),])
  temp.df$interaction<-rownames(temp.df)
  if ("df.all" %in% ls()) {
    df.all<-rbind(df.all,temp.df)
  } else {df.all<-temp.df}
}


std <- function(x) sd(x)/sqrt(length(x))

df.mean<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),mean))
df.std<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),std))

df.mean$alpha.beta.std<-df.std$alpha.beta.coef 
df.mean$alpha.std<-df.std$seeds

names(df.mean)<-c("interaction","alpha.beta.mean","alpha.mean","alpha.beta.std","alpha.std")

CairoPDF("figures/Australia_Trait_complex.pdf",height=10,width=10)
colors<-rep(c("red","green","orange","hotpink","blue","black","grey40"),6)
ggplot(data=df.mean,aes(x=alpha.mean,y=alpha.beta.mean,col=interaction))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1,size=1.5)+
  ylim(-3,3)+xlim(-3,3)+
  geom_errorbarh(aes(xmin=alpha.mean-2*alpha.std,xmax=alpha.mean+2*alpha.std))+
  geom_errorbar(aes(ymin=alpha.beta.mean-2*alpha.beta.std,ymax=alpha.beta.mean+2*alpha.beta.std)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia Trait Complex")+
  geom_vline(xintercept=0, col="grey60")+
  geom_hline(yintercept=0, col="grey60")+
  scale_color_manual(values=colors)+
  guides(col=F)+
  theme_minimal(base_size = 20)+theme(panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank())

dev.off()

rm(df.all)
