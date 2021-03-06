############### read in data #######################
# load the Australia data
load('raw.data/mayfield.Rdata')
# stitch the species back together into one large dataframe
Australia.data <- do.call(rbind, mayfield)
# remove neighbors with no individuals
rm<-names(which(colSums(Australia.data[5:51])==0))
Australia.data<-Australia.data[,-which(names(Australia.data)%in%rm)]


# load the Spain data
Spain.data<-read.csv("raw.data/Caracoles.competition.updated.14.05.2018.csv")
# remove columns not used in analysis
remove.col.G<-c("subplot","date","order","fruit")
Spain.data<-Spain.data[,-which(names(Spain.data) %in% remove.col.G)]
# rename a few columns 
colnames(Spain.data)[which(names(Spain.data)=="seed")]<-"seeds"
Spain.data$site<-"Caracoles"
# change any NA's to '0'
Spain.data[is.na(Spain.data)]<-0

table(Spain.data$focal)

# remove focals with less than 20 points
Spain.data<-Spain.data[-which(Spain.data$focal=="MEEL"),] # remove focal with only one point
Spain.data<-Spain.data[-which(Spain.data$focal=="MESU"),] # remove focal with only 14 points
Spain.data<-Spain.data[-which(Spain.data$focal=="POMO"),] # remove focal with only 12 points
Spain.data<-Spain.data[-which(Spain.data$focal=="SPRU"),] # remove focal with only 12 points
Spain.data<-Spain.data[-which(Spain.data$focal=="BEMA"),] # remove focal with only 18 points

colSums(Spain.data[4:21])

Spain.data<-Spain.data[,-which(names(Spain.data)=="MEEL")] # remove neighbor that now has no occurances after removing tha above focals

# create an empty data list and add the full data set for the Spain and Australia data to it
data.list<-list()
data.list[["Spain"]]<-Spain.data
data.list[["Australia"]]<-Australia.data

print("Reading and cleaning data done")