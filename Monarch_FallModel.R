#############################################################################
#############FALL MODEL using OW COLONY SIZES##############################
#############################################################################

#load jagsUI library
library(jagsUI)
library(reshape)

#set working directory

#read in colony sizes data
fallmon = read.csv("SurveyData_Fall.csv", header=TRUE,sep=',',na.strings=T)
fallmon[1:10,]
names(fallmon)
str(fallmon)

#Look at data (checking all years/sites present)
uyear = sort(unique(fallmon$Year))
usite = sort(unique(fallmon$SiteID))

#Reshape the survey data.
junk.melt=melt(fallmon,id.var=c("Year", "SiteID"), measure.var="Size")
fallmonarchs=cast(junk.melt, SiteID ~ Year, sum)
fallmonarchs=fallmonarchs[,2:13]
#table of 19 sites by 12 years

#read in SITE-YEAR data##################
siteyear = read.csv("SiteYear_Fall.csv", header=T,sep=',',na.strings=T)
str(siteyear)

#organize site by year
#report appropriate DENSE forest cover values according to siteID and year
forest_d = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)) {
  for (j in 1:length(usite)) {
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    forest_d[j,tt] = siteyear$Forest_Dense[a] 
  }} 

#then standardize
forestdmean=mean(forest_d,na.rm=T)
forestdsd=sd(as.vector(forest_d), na.rm=T)
forestd.st<-round((forest_d-forestdmean)/forestdsd,3)

#dim of forest.st: sites (19) by year (12)

#report appropriate DENSE+OPEN forest cover values according to siteID and year
forest_do = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)) {
  for (j in 1:length(usite)) {
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    forest_do[j,tt] = siteyear$Forest_D_O[a] 
  }} 

#then standardize
forestdomean=mean(forest_do,na.rm=T)
forestdosd=sd(as.vector(forest_do), na.rm=T)
forestdo.st<-round((forest_do-forestdomean)/forestdosd,3)

# include density dependence component (colony size in prev year)
size_1yr = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)) {
  for (j in 1:length(usite)) {
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    size_1yr[j,tt] = siteyear$Size_1yr[a] 
  }} 

#then standardize
ddmean=mean(size_1yr,na.rm=T)
ddsd=sd(as.vector(size_1yr), na.rm=T)
dd.st<-round((size_1yr-ddmean)/ddsd,3)

#convert to autoregressive presence/absence for hurdle part of model
autoreg <- matrix(NA, nrow = nrow(fallmonarchs), ncol = ncol(fallmonarchs))
autoreg <- ifelse(size_1yr > 0, 1, 0)

#read in YEAR data#####################
yeareffects=read.csv("Year_Fall.csv", header=T,sep=',',na.strings=T)
str(yeareffects)
names(yeareffects)

roosts=yeareffects$JNroosts
nectar1=yeareffects$NDVIR1
nectar2=yeareffects$NDVIR2
mintemp1=yeareffects$TMinR1
mintemp2=yeareffects$TMinR2
avgtemp1=yeareffects$TAvgR1
avgtemp2=yeareffects$TAvgR2
sumcount=yeareffects$SumCount
oeinf=yeareffects$OEInf
naba=yeareffects$NABA

#standardize all year covariates
#fall roosts
mroosts=mean(as.matrix(roosts))
sdroosts=sd(as.vector(roosts))
roosts.st=as.vector((roosts-mroosts)/sdroosts)

#nectar region 1 availability 
mnectar1=mean(as.matrix(nectar1))
sdnectar1=sd(as.vector(nectar1))
nectar.st1=as.vector((nectar1-mnectar1)/sdnectar1)

#nectar region 2 availability
mnectar2=mean(as.matrix(nectar2))
sdnectar2=sd(as.vector(nectar2))
nectar.st2=as.vector((nectar2-mnectar2)/sdnectar2)

#min temps in fall region 1
mmintemp1=mean(as.matrix(mintemp1))
sdmintemp1=sd(as.vector(mintemp1))
mintemp.st1=as.vector((mintemp1-mmintemp1)/sdmintemp1)

#min temps in fall region 2
mmintemp2=mean(as.matrix(mintemp2))
sdmintemp2=sd(as.vector(mintemp2))
mintemp.st2=as.vector((mintemp2-mmintemp2)/sdmintemp2)

#average temps in fall region 1
mavgtemp1=mean(as.matrix(avgtemp1))
sdavgtemp1=sd(as.vector(avgtemp1))
avgtemp.st1=as.vector((avgtemp1-mavgtemp1)/sdavgtemp1)

#average temps in fall region 2
mavgtemp2=mean(as.matrix(avgtemp2))
sdavgtemp2=sd(as.vector(avgtemp2))
avgtemp.st2=as.vector((avgtemp2-mavgtemp2)/sdavgtemp2)

#proportion of OE infection
#NOTE: 2015 value is imputed
moeinfect=mean(as.matrix(oeinf))
sdoeinfect=sd(as.vector(oeinf))
oeinfect.st=as.vector((oeinf-moeinfect)/sdoeinfect)

#NABA end of summer counts
mnaba=mean(as.matrix(naba))
sdnaba=sd(as.vector(naba))
naba.st=as.vector((naba-mnaba)/sdnaba)

#year (to test if any trend as FE)
myear=mean(uyear)
sdyear=sd(uyear)
year.st=(uyear-myear)/sdyear

#read in SITE data#######################
siteeffects=read.csv("Site_Fall.csv", header=T,sep=',',na.strings=T)
str(siteeffects)
dim(siteeffects) 

#in or out of reserve (1 is in reserve, 0 is out)
reserve=siteeffects$Reserve

#### hurdle model adjustments
# create separate matrix of 0s and 1s indicating whether colony present in each site-year
# convert to z in bugsdata below to include in hurdle model
non_zero <- matrix(NA, nrow = nrow(fallmonarchs), ncol = ncol(fallmonarchs))
non_zero <- ifelse(fallmonarchs > 0, 1, 0)

#create matrix of ones to feed in as data for the ones trick
ones <- matrix(NA, nrow = nrow(fallmonarchs), ncol = ncol(fallmonarchs))
ones[is.na(ones)] <- 1
######

####################################################################
#PREP JAGSUI DATA
bugsdata<-list(uyear=length(uyear), usite=length(usite),
               y=data.matrix(fallmonarchs), z=data.matrix(non_zero),ones=ones, roosts.st=roosts.st, nectar.st1=nectar.st1, nectar.st2=nectar.st2,
               mintemp.st1=mintemp.st1, mintemp.st2=mintemp.st2, avgtemp.st1=avgtemp.st1, avgtemp.st2=avgtemp.st2, oeinfect.st=oeinfect.st, 
               dd.st=dd.st, naba.st=naba.st, reserve=reserve, forestd.st=forestd.st, forestdo.st=forestdo.st, autoreg=autoreg, year.st=year.st) 

inits<-function(){
  list(a1=runif(1,1,10),a2=rnorm(1),a3=rnorm(1),a4=rnorm(1),a5=rnorm(1),a6=rnorm(1),g0=rnorm(1))
}

parameters<-c('a1', 'a2', 'a3', 'a4', 'a5', 'g0', 'sigma.a8', 'sigma.a9', 'w')

###################################################################
##RUN BUGS MODEL IN JAGSUI
###################################################################
#example run with top-supported covariates
output <-jags(data = bugsdata,
                     inits = inits,               
                     parameters.to.save = parameters,
                     model.file = 'top_model_covs.txt',
                     n.chains = 3,
                     n.adapt = 5000, 
                     n.iter = 100000, 
                     n.burnin = 50000, 
                     n.thin = 5, 
                     parallel = TRUE,
                     store.data = TRUE)
                    
print(output,digits=4)   
traceplot(output)
whiskerplot(output,parameters=c('a1', 'a2', 'a3', 'a4', 'a5', 'g0', 'sigma.a8', 'sigma.a9', 'w'))
