rm(list=ls())
library(ggplot2)
library(Cairo)
library(tidyverse)
rares.G<-read.csv("Species.grouping/Rare.grouping/output/rare.grouped.summary.kfolds.2019.csv")
rares.R<-read.csv("Species.grouping/Rare.grouping/output/rare.removed.summary.kfolds.2019.csv")
names(rares.G)<-c("Data","beta","coef","AIC","train.RMSE","predict.RMSE","BIC","fold","likeli")
names(rares.R)<-c("Data","beta","coef","AIC","train.RMSE","predict.RMSE","BIC","fold","likeli")
rares.G$group<-"rares.grouped"
rares.R$group<-"rares.removed"
rares.G<-rares.G[,sort(names(rares.G))]
rares.R<-rares.R[,sort(names(rares.R))]
spain.neigh<-read.csv("Neigh.ID/output/Spain.neighbor.ID.summary.2019.HOI.kfolds.csv")
spain.neigh2<-spain.neigh[,-1]
spain.neigh2<-spain.neigh2[,-which(names(spain.neigh2)=="focal")]
names(spain.neigh2)<-c("group","beta","fold","coef","AIC","train.RMSE","predict.RMSE","BIC","likeli")
spain.neigh2$Data<-"Spain"
spain.neigh3<-spain.neigh2[,sort((names(spain.neigh2)))]
aus.neigh<-read.csv("Neigh.ID/output/Australia.neighbor.ID.summary.2019.HOI.kfolds.csv")
aus.neigh2<-aus.neigh[,-1]
aus.neigh2<-aus.neigh2[,-which(names(aus.neigh2)=="focal")]
names(aus.neigh2)<-c("group","beta","fold","coef","AIC","train.RMSE","predict.RMSE","BIC","likeli")
aus.neigh2$Data<-"Australia"
aus.neigh3<-aus.neigh2[,sort((names(aus.neigh2)))]

###############################################
grouped.1<-rbind(aus.neigh3,spain.neigh3)
grouped.2<-rbind(rares.G,rares.R)
grouped<-rbind(grouped.1,grouped.2)

groups<-unique(grouped$group)
group.shapes<-c("SP","AB","OS","LF","FG","FAM","CL","RG","RR")

grouped.HOI<-grouped[which(grouped$beta==T),]
long.grouped.HOI<- grouped.HOI %>% 
  gather("AIC","BIC","likeli","predict.RMSE","train.RMSE",key="Metric.name",value="Metric")


grouped.direct<-grouped[which(grouped$beta==F),]
long.grouped.direct<- grouped.direct %>% 
  gather("AIC","BIC","likeli","predict.RMSE","train.RMSE",key="Metric.name",value="Metric")



m.grouped<-merge(long.grouped.HOI,long.grouped.direct,by=c("Data","fold","group","Metric.name"),suffixes=c(".HOI",".direct"))

m.grouped.mean<-with(m.grouped,aggregate(cbind(Metric.HOI,Metric.direct),by=list(Data,group,Metric.name),mean))
names(m.grouped.mean)[1:3]<-c("Data","group","Metric.name")
std <- function(x) sd(x)/sqrt(length(x))
m.grouped.sd<-with(m.grouped,aggregate(cbind(Metric.HOI,Metric.direct),by=list(Data,group,Metric.name),std))
names(m.grouped.sd)[1:3]<-c("Data","group","Metric.name")
m.grouped.coef<-with(m.grouped,aggregate(cbind(coef.HOI,coef.direct),by=list(Data,group,Metric.name),mean))
names(m.grouped.coef)[1:3]<-c("Data","group","Metric.name")

m.grouped.2<-merge(m.grouped.mean,m.grouped.sd,by=c("Data","group","Metric.name"),suffixes=c(".mean",".sd"))
m.grouped.2<-merge(m.grouped.2,m.grouped.coef,by=c("Data","group","Metric.name"))
#m.grouped.2<-m.grouped.2[-which(m.grouped.2$group=="functional"),]
m.grouped.2<-m.grouped.2[-which(m.grouped.2$Metric.name=="likeli"),]
gp.rm<-c("full","rares.grouped","rares.removed")
gp.rm.A<-c("full","rares.grouped","rares.removed","families")
m.grouped.3<-m.grouped.2[-which(m.grouped.2$Metric.name=="predict.RMSE"& m.grouped.2$group %in% gp.rm&m.grouped.2$Data=="Spain"),]
m.grouped.3<-m.grouped.3[-which(m.grouped.3$Metric.name=="predict.RMSE"& m.grouped.3$group %in% gp.rm.A&m.grouped.3$Data=="Australia"),]

unique(m.grouped.3$group)

m.grouped.3$group<-factor(m.grouped.3$group,
                          levels=c("full","rares.grouped","rares.removed",
                                   "number.indivs","families","cluster","functional",
                                   "grass.forb","native.exotic"),
                          labels=c("Species Identity","Rare Species Grouped","Rare Species Removed",
                                   "Abundance","Plant Family","Trait Complex","Functional Type",
                                   "Life Form","Origin Status"))

m.grouped.3s<-m.grouped.3

m.grouped.3$Metric.name<-factor(m.grouped.3$Metric.name,
                                levels=c("AIC","BIC","train.RMSE","predict.RMSE"),
                                labels=c("(A) AIC","(B) BIC","(C) RMSE training","(D) RMSE testing"))


colors<-c("black","grey40","grey65",
          "red","darkgoldenrod1","chocolate2","darkmagenta",
          "seagreen","royalblue3")
#################### AUstralia ###################
CairoPDF("figures/Australia_neigh_ID_Aug_2020.pdf",height=10,width=10)
dummy <- data.frame(Metric.direct.mean = c(7900, 11100, 7800, 8300,0.94,1.2,0.7,2.5),
                    Metric.HOI.mean = c(7900, 11100, 7800, 8300,0.94,1.2,0.7,2.5),
                    Metric.name = c("(B) BIC", "(B) BIC", "(A) AIC", "(A) AIC","(C) RMSE training","(C) RMSE training","(D) RMSE testing","(D) RMSE testing"),
                    group=rep("Life Form",8))
p.aus<-ggplot(data=m.grouped.3[which(m.grouped.3$Data=="Australia"),],aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1)+
  facet_wrap(vars(Metric.name), scales="free")+
  geom_blank(data=dummy)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia")+
  scale_color_manual(values=colors)+
  theme_minimal(base_size = 20)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      strip.text.x = element_text(hjust = -.01),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=guide_legend(nrow=3))#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")
p.aus
dev.off()

CairoPDF("figures/Australia_neigh_ID_Aug_2020.pdf",height=10,width=10)
dummy <- data.frame(Metric.direct.mean = c(7900, 11100, 7800, 8300,0.94,1.2,0.7,2.5),
                    Metric.HOI.mean = c(7900, 11100, 7800, 8300,0.94,1.2,0.7,2.5),
                    Metric.name = c("(B) BIC", "(B) BIC", "(A) AIC", "(A) AIC","(C) RMSE training","(C) RMSE training","(D) RMSE testing","(D) RMSE testing"),
                    group=rep("Life Form",8))
p.aus<-ggplot(data=m.grouped.3[which(m.grouped.3$Data=="Australia"),],aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1)+
  facet_wrap(vars(Metric.name), scales="free")+
  geom_blank(data=dummy)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia")+
  scale_color_manual(values=colors)+
  theme_minimal(base_size = 22)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=20),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      strip.text.x = element_text(hjust = -.01),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=guide_legend(nrow=3))#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")
p.aus
dev.off()

################### Spain #################################

S.colors<-c("black","grey40","grey65",
            "red","darkgoldenrod1","chocolate2","darkmagenta",
            "seagreen")

m.grouped.3.S<-m.grouped.3[-which(m.grouped.3$Data=="Spain"& m.grouped.3$group=="Origin Status"),]


CairoPDF("figures/Spain_neigh_ID_Aug_2020.pdf",height=10,width=10)
dummy <- data.frame(Metric.direct.mean = c(13400,16500, 13200,14000,1.05,1.5,0.5,1.6),
                    Metric.HOI.mean = c(13400,16500, 13200,14000,1.05,1.5,0.5,1.6),
                    Metric.name = c("(B) BIC", "(B) BIC", "(A) AIC", "(A) AIC","(C) RMSE training","(C) RMSE training","(D) RMSE testing","(D) RMSE testing"),
                    group=rep("Life Form",8))
p.spain<-ggplot(data=m.grouped.3.S[which(m.grouped.3.S$Data=="Spain"),],aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1)+
  facet_wrap(vars(Metric.name), scales="free")+
  geom_blank(data=dummy)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Spain")+
  scale_color_manual(values=S.colors)+
  theme_minimal(base_size = 20)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      strip.text.x = element_text(hjust = -.01),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=guide_legend(nrow=3))#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")
p.aus
dev.off()

m.grouped.3.S$Metric.name<-factor(m.grouped.3.S$Metric.name,
                                levels=c("(A) AIC","(B) BIC","(C) RMSE training","(D) RMSE testing"),
                                labels=c("(E) AIC","(F) BIC","(G) RMSE training","(H) RMSE testing"))

CairoPDF("figures/Spain_neigh_ID_Aug_2020.pdf",height=10,width=10)
dummy <- data.frame(Metric.direct.mean = c(13400,16500, 13200,14000,1.05,1.5,0.5,1.6),
                    Metric.HOI.mean = c(13400,16500, 13200,14000,1.05,1.5,0.5,1.6),
                    Metric.name = c("(F) BIC", "(F) BIC", "(E) AIC", "(E) AIC","(G) RMSE training","(G) RMSE training","(H) RMSE testing","(H) RMSE testing"),
                    group=rep("Life Form",8))
p.spain<-ggplot(data=m.grouped.3.S[which(m.grouped.3.S$Data=="Spain"),],aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1)+
  facet_wrap(vars(Metric.name), scales="free")+
  geom_blank(data=dummy)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Spain")+
  scale_color_manual(values=S.colors)+
  theme_minimal(base_size = 22)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=20),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      strip.text.x = element_text(hjust = -.01),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=guide_legend(nrow=3))#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")
p.spain
dev.off()

########### predict.RMSE plot ###################
m.grouped.2$group<-factor(m.grouped.2$group,
                          levels=c("full","rares.grouped","rares.removed",
                                   "number.indivs","families","cluster","functional",
                                   "grass.forb","native.exotic"),
                          labels=c("Species Identity","Rare Species Grouped","Rare Species Removed",
                                   "Abundance","Plant Family","Trait Complex","Functional Type",
                                   "Life Form","Origin Status"))
m.grouped.2$Metric.name<-factor(m.grouped.2$Metric.name,
                                levels=c("AIC","BIC","train.RMSE","predict.RMSE"),
                                labels=c("AIC","BIC","RMSE training","RMSE testing"))


colors<-c("black","grey40","grey65",
          "red","darkgoldenrod1","chocolate2","darkmagenta",
          "seagreen","royalblue3")


CairoPDF("figures/Australia_RMSE_test_all.pdf",width=8,height=8)
ggplot(data=m.grouped.2[which(m.grouped.2$Data=="Australia" & m.grouped.2$Metric.name=="RMSE testing"),],
       aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+
  geom_abline(slope=1)+ylim(-150,650)+xlim(-150,650)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI model")+ xlab("Direct model")+
  ggtitle("Australia RMSE testing")+
  scale_color_manual(values=colors)+
  theme_minimal(base_size = 20)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=guide_legend(nrow=3))#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")

dev.off()


########### Metric plot ###################
S.colors<-c("black","grey40","grey65",
            "red","darkgoldenrod1","chocolate2","darkmagenta",
            "seagreen")

m.grouped.2.S<-m.grouped.2[-which(m.grouped.2$Data=="Spain"& m.grouped.2$group=="Origin Status"),]

CairoPDF("figures/Spain_RMSE_test_all.pdf",width=8,height=8)
ggplot(data=m.grouped.2.S[which(m.grouped.2.S$Data=="Spain" & m.grouped.2.S$Metric.name=="RMSE testing"),],
       aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+
  geom_abline(slope=1)+ylim(0,7)+xlim(0,7)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI model")+ xlab("Direct model")+
  ggtitle("Spain RMSE testing")+
  scale_color_manual(values=S.colors)+
  theme_minimal(base_size = 20)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=guide_legend(nrow=3))#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")

dev.off()

##################################################
source("multiplot.function.R")
#############################
dummy <- data.frame(Metric.direct.mean = c(7900, 11100, 7800, 8300,0.94,1.2,0.7,2.5),
                    Metric.HOI.mean = c(7900, 11100, 7800, 8300,0.94,1.2,0.7,2.5),
                    Metric.name = c("(B) BIC", "(B) BIC", "(A) AIC", "(A) AIC","(C) RMSE training","(C) RMSE training","(D) RMSE testing","(D) RMSE testing"),
                    group=rep("Life Form",8))
p.aus2<-ggplot(data=m.grouped.3[which(m.grouped.3$Data=="Australia"),],aes(x=Metric.direct.mean,y=Metric.HOI.mean,col=group))+
  geom_point(shape=21,size=5,stroke=1.5)+geom_abline(slope=1)+
  facet_wrap(vars(Metric.name), scales="free")+
  geom_blank(data=dummy)+
  geom_errorbarh(aes(xmin=Metric.direct.mean-2*Metric.direct.sd,xmax=Metric.direct.mean+2*Metric.direct.sd))+
  geom_errorbar(aes(ymin=Metric.HOI.mean-2*Metric.HOI.sd,ymax=Metric.HOI.mean+2*Metric.HOI.sd)) +
  ylab("HOI Model")+ xlab("Direct Model")+
  ggtitle("Australia")+
  scale_color_manual(values=colors)+
  theme_minimal(base_size = 20)+theme(legend.direction="horizontal",
                                      legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),   # Major grid lines
                                      panel.grid.minor = element_blank(),
                                      strip.text.x = element_text(hjust = -.01),
                                      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                      panel.margin = unit(2, "lines"))+
  guides(col=F)#+guides(col=F)#+theme(legend.direction="horizontal",legend.position = "top")
p.aus2

CairoPDF("figures/Figure_2.pdf",height=20,width=12)
multiplot(p.aus2,p.spain, cols=1)
dev.off()
  