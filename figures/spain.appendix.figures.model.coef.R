################# Spain grass forb #####################
alpha.model.list<-list()
alpha.beta.model.list<-list()

rm(df.all)

for (i in 1:5){
  alpha.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Spain.alpha.model.list.grass.forbfold.",i,".2019.rds"))
  alpha.beta.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Spain.alpha.beta.model.list.grass.forbfold.",i,".2019.rds"))
  alpha.coef<-coef(alpha.model.list[[i]])
  alpha.beta.coef<-coef(alpha.beta.model.list[[i]])[1:(length(alpha.coef))]
  temp.df<-cbind(alpha.coef,alpha.beta.coef)
  temp.df<-data.frame(temp.df[-grep("fecundity",rownames(temp.df)),])
  temp.df$interaction<-rownames(temp.df)
  if ("df.all" %in% ls()) {
    df.all<-rbind(df.all,temp.df)
  } else {df.all<-temp.df}
}


with(df.all,plot(alpha.beta.coef,seeds))
abline(a=0,b=1)

std <- function(x) sd(x)/sqrt(length(x))

df.mean<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),mean))
df.std<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),std))

df.mean$alpha.beta.std<-df.std$alpha.beta.coef 
df.mean$alpha.std<-df.std$seeds

names(df.mean)<-c("interaction","alpha.beta.mean","alpha.mean","alpha.beta.std","alpha.std")

library(ggplot2)
library(Cairo)
CairoPDF("figures/Spain_Grass_Forb.pdf",height=10,width=10)
colors<-rep(c("red","blue"),12)
ggplot(data=df.mean,aes(x=alpha.mean,y=alpha.beta.mean,col=interaction))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1,size=1.5)+
  ylim(-0.8,1.5)+xlim(-0.8,1.5)+
  geom_errorbarh(aes(xmin=alpha.mean-2*alpha.std,xmax=alpha.mean+2*alpha.std))+
  geom_errorbar(aes(ymin=alpha.beta.mean-2*alpha.beta.std,ymax=alpha.beta.mean+2*alpha.beta.std)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Spain Life Form")+
  geom_vline(xintercept=0, col="grey60")+
  geom_hline(yintercept=0, col="grey60")+
  scale_color_manual(values=colors)+
  guides(col=F)+
  theme_minimal(base_size = 20)+theme(panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank())

dev.off()

################# Spain
alpha.model.list<-list()
alpha.beta.model.list<-list()

rm(df.all)

for (i in 1:5){
  alpha.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Spain.alpha.model.list.number.indivsfold.",i,".2019.rds"))
  alpha.beta.model.list[[i]]<-readRDS(paste0("Neigh.ID/output/Spain.alpha.beta.model.list.number.indivsfold.",i,".2019.rds"))
  alpha.coef<-coef(alpha.model.list[[i]])
  alpha.beta.coef<-coef(alpha.beta.model.list[[i]])[1:(length(alpha.coef))]
  temp.df<-cbind(alpha.coef,alpha.beta.coef)
  temp.df<-data.frame(temp.df[-grep("fecundity",rownames(temp.df)),])
  temp.df$interaction<-rownames(temp.df)
  if ("df.all" %in% ls()) {
    df.all<-rbind(df.all,temp.df)
  } else {df.all<-temp.df}
}

with(df.all,plot(alpha.beta.coef,seeds))
abline(a=0,b=1)

std <- function(x) sd(x)/sqrt(length(x))

df.mean<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),mean))
df.std<-with(df.all,aggregate(cbind(alpha.beta.coef,seeds),by=list(interaction),std))

df.mean$alpha.beta.std<-df.std$alpha.beta.coef 
df.mean$alpha.std<-df.std$seeds

names(df.mean)<-c("interaction","alpha.beta.mean","alpha.mean","alpha.beta.std","alpha.std")

CairoPDF("figures/Spain_Abundance.pdf",height=10,width=10)
colors<-rep(c("red"),12)
ggplot(data=df.mean,aes(x=alpha.mean,y=alpha.beta.mean,col=interaction))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1,size=1.5)+
  ylim(-0.5,1.0)+xlim(-0.5,1.0)+
  geom_errorbarh(aes(xmin=alpha.mean-2*alpha.std,xmax=alpha.mean+2*alpha.std))+
  geom_errorbar(aes(ymin=alpha.beta.mean-2*alpha.beta.std,ymax=alpha.beta.mean+2*alpha.beta.std)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Spain Abundance")+
  geom_vline(xintercept=0, col="grey60")+
  geom_hline(yintercept=0, col="grey60")+
  scale_color_manual(values=colors)+
  guides(col=F)+
  theme_minimal(base_size = 20)+theme(panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank())

dev.off()

