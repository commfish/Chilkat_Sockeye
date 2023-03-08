# this code gets sourced from this file, and creates the "dat" object, which is then used in
# the main script
# read in raw data----
rawdat<-as.data.frame(read.csv("2023_analysis/data/chilkat_sockeye.csv",header=T))
nyrs<-as.numeric(length(rawdat$year))
fyr<-min(rawdat$year)
lyr<-max(rawdat$year)
nages<-3
a.min<-4
a.max<-6
A<-3

# data clean----
year <- as.numeric(as.character(rawdat$year))
DS <- as.numeric(as.character(rawdat$DS))
DS.cv <- as.numeric(as.character(rawdat$DS.cv)) 
mr <- as.numeric(as.character(rawdat$mr))
mr.cv <- as.numeric(as.character(rawdat$mr.cv)) 
weir <- as.numeric(as.character(rawdat$weir)) 
weir.cv <- as.numeric(as.character(rawdat$weir.cv)) 
hbelow <- as.numeric(as.character(rawdat$hbelow))
hbelow.cv <- as.numeric(as.character(rawdat$hbelow.cv))
x<-as.matrix(rawdat[,substr(colnames(rawdat), 1,1)=="x"])#age comp count data matrix
colnames(x)<-NULL
n.a<-rowSums(x) # effective sample size


dat=list(Y = nyrs, A=nages, a.min=a.min, a.max=a.max,
         x=x, DS=DS, DS.cv=DS.cv, mr=mr, mr.cv=mr.cv, weir=weir,
         weir.cv=weir.cv, hbelow=hbelow, hbelow.cv=hbelow.cv,
         n.a=n.a)

