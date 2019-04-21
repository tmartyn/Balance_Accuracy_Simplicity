# read in rare species grouping results (Australia)
rares<-read.csv("Species.grouping/Rare.grouping/output/rare.grouped.by.species.2019.csv")

# read in other species grouping results (Australia)
Australia.summary<-read.csv("Neigh.ID/output/Australia.neighbor.ID.summary.2019.by.species.csv")
Australia.summary<-Australia.summary[,-1]

# merge rare grouping with other groupings
rares.M<-rares[which(rares$Model=="Australia"),]
rares.M$Model<-"group.rares"
Australia.summary<-rbind(Australia.summary,rares.M)

# read in rare speies grouping results (Spain)
rares.G<-rares[which(rares$Model=="Spain"),]
rares.G$Model<-"group.rares"

# read in other grouping resluts (Spain)
Spain.summary<-read.csv("Neigh.ID/output/Spain.neighbor.ID.summary.2019.by.species.csv")
Spain.summary<-Spain.summary[,-1]

# merge results together (Spain)
Spain.summary<-rbind(Spain.summary,rares.G)


##########################
# catch plot as pdf
CairoPDF(width=17,height=8,"figures/BOTH_radar_updated.pdf", bg="transparent")

# transform dataframe to long
MM<-Australia.summary %>%
  gather(key=variable,value,-c(Model,focal,beta)) 
MM<-spread(MM,Model,value)

# list species names fully
M.nm<-c("A.caryophylla","H.glabra","P.gnaphalioides","T.ornata","U.anthemoides","W.acuminata")

# don't plot results for these groups as they are summarized elsewhere
groups<-unique(Australia.summary$Model)
groups<-groups[!groups%in%c("full","number.indivs")]

#subtract values from full value
subs<-MM[,which(names(MM) %in% groups)]-MM$full
# attach data together
JJ<-cbind(subs,focal=MM$focal,beta=MM$beta,variable=MM$variable)
# only pull out the data with includes AIC and HOIs
JJ<-JJ[which(JJ$beta=="TRUE"& JJ$variable=="AIC"),]
# only use the functional groups in groups
KK<-JJ[,which(names(JJ) %in% groups)]
# rownames are focal species names
rownames(KK)<-M.nm
# columns names are functional grouping names
names(KK)<-c("Trait Complex","Plant Family","Functional Type","Life Form","Native/Exotic","Group Rares")

# build dataframe specific to radar plot specificaton where the irst 3 rows need to be the max lim of plot,
#   the min lim of plots, and any lines you want to add to the plots (i.e. the dashed line a 0)
data<-data.frame(rbind(rep(20,dim(KK)[2]) , rep(-160,dim(KK)[2]) , rep(0,dim(KK)[2]),KK))

#calculate how many times each species X each functional group is less than 0 (will be center numbers)
HH<-apply(KK,2,function(x) length(which(x < 0)==T))

# assign colors for lines in plot
h<-c("firebrick2","darkorange1","gold","chartreuse2",
     "dodgerblue","purple")

# make plot
colors_border<-c("black",h)
line.wd<-c(2,2,2,2,2,2,2)
line.typ<-c(5,1,1,1,1,1,1)
point.type<-c(32,20,20,20,20,20,20)

par(mfrow=c(1,2))

radarchart( data  , axistype=1 , 
            #custom polygon
            pcol=colors_border , plwd=line.wd , plty=line.typ, pty=point.type,
            #custom the grid
            cglcol="grey70", cglty=1, axislabcol="grey40", seg=9, caxislabels=round(seq(-160,20,length.out=10)), cglwd=1.3,
            #custom labels
            vlcex=1.3, bg="grey"
)
legend("bottomleft", legend = rownames(data[-c(1:3),]), text.font=3,
       bty = "n", pch=20 , col=colors_border[2:10] , text.col = "grey", cex=1, pt.cex=3)

# make inner numbers
text.x<-c(0,-0.06,-.06,0,.06,.06)
text.y<-c(.07,.04,-.04,-.07,-.04,.04)
text.labels<-c(HH)
text(text.x,text.y,text.labels,cex=1)

######## Spain data ###########################
# make data long
GG<-Spain.summary %>%
  gather(key=variable,value,-c(Model,focal,beta)) 
GG<-spread(GG,Model,value)

# use only specific functional group results for plot
groups<-unique(Spain.summary$Model)
groups<-groups[!groups%in%c("full","number.indivs","native.exotic")]

# subtract values from full model
FF<-GG[,which(names(GG) %in% groups)]-GG$full
# bind data together
JJ<-cbind(FF,focal=GG$focal,beta=GG$beta,variable=GG$variable)
# pull out only AIC and HOI data
JJ<-JJ[which(JJ$beta=="TRUE"& JJ$variable=="AIC"),]
# use only the data in the selected functional groups
KK<-JJ[,which(names(JJ) %in% groups)]
colnames(KK)<-c("Trait Complex","Plant Family","Functional Type","Life Form","Group Rares")

# list species names full
G.nm<-c("C.tennuiflorum","C.fuscatum","C.mixtum","H.marinum","L.maroccanus","L.tribracteatum",
        "P.incurva","P.coronopus","P.maritimus","P.paludosa","S.soda","S.laciniata")
rownames(KK)<-G.nm # assign species names to rownames

# bind data together specifically for radar plot rad-in strucutre
data<-data.frame(rbind(rep(15,dim(KK)[2]) , rep(-60,dim(KK)[2]) , rep(0,dim(KK)[2]),KK))

# summarize number of focals X groupings under 0 (will be venter numbers)
HH<-apply(KK,2,function(x) length(which(x < 0)==T))

# make plot
h<-c("firebrick2","darkorange1","gold","chartreuse2","darkgreen","cyan2","dodgerblue","blue4","purple","magenta2","chocolate4","lightslategrey") 
colors_border<-c("black",h)
line.wd<-c(2,2,2,2,2,2,2,2,2,2,2,2,2)
line.typ<-c(5,1,1,1,1,1,1,1,1,1,1,1,1)
point.type<-c(32,20,20,20,20,20,20,20,20,20,20,20,20)

radarchart( data  , axistype=1 , 
            #custom polygon
            pcol=colors_border , plwd=line.wd , plty=line.typ, pty=point.type,
            #custom the grid
            cglcol="grey70", cglty=1, axislabcol="grey40", seg=5, caxislabels=round(seq(-60,15,length.out=6)), cglwd=1.3,
            #custom labels
            vlcex=1.2, bg="grey"
)

# add legend
legend("bottomleft", legend = rownames(data[-c(1:3),]), text.font=3,
       bty = "n", pch=20 , col=colors_border[2:13] , text.col = "grey", cex=1, pt.cex=3,ncol=1)

#make center numbers
text.x<-c(0,-0.085,-.065,.065,.085)
text.y<-c(.09,.04,-.07,-.07,.04)
text.labels<-c(HH)
text(text.x,text.y,text.labels,cex=1.0)

# end plot catch
dev.off()
print("Radar plot done")


