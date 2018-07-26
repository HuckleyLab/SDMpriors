#load libraries
library(dismo)  #see also zoon R package?
library(plyr)
library(rgbif)
library(GRaF)  #see methods paper here: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12523/pdf
library(pROC)
library(ROCR)
library(ncdf4)

#for cleaning data
library(biogeo) #https://cran.r-project.org/web/packages/biogeo/index.html
library(rgeospatialquality) #https://cran.r-project.org/web/packages/rgeospatialquality/

#tutorial here: https://rawgit.com/goldingn/intecol2013/master/tutorial/graf_workshop.htm

#Potential Resources 
#web: http://sdmdata.sdmserialsoftware.org

#--------------------------------
# load Sunday database
setwd("/Users/laurenbuckley/SDMpriors/")
dat= read.csv("Sundayetal_thermallimits.csv")
#start with reptiles and amphibians
dat= subset(dat, dat$phylum=="Chordata")
#start with species with Tmin and Tmax
dat= dat[!is.na(dat$tmax) & !is.na(dat$tmin),]
#subset to critical rather than lethal

#write out list
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/")
write.csv(dat,"SpeciesList.csv")

#species name to enable match
dat$spec = gsub("_"," ",dat$species)

#------------------------
# Query GBIF (R rgbif package) for specimen localities

#Write out localities
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/GBIF/")

# loop through species
for(spec.k in 29:43){ #nrow(dat)

#look up species
key <- name_suggest(q=dat$spec[spec.k], rank='species')$key[1]

occ <- occ_data(scientificName=dat$spec[spec.k], limit=1000)
occ <- occ$data

#write out
filename<-paste("GBIFloc_", dat$spec[spec.k],".csv", sep="")

write.csv(occ[,1:5],filename)

} #end looop species

#occ=occ_search(taxonKey=key, limit=2000, return="data") 
#fields=c('name','basisOfRecord','protocol')
#return: can get metadata, etc.

spec.k=55 #Sceloporus occidentalis
spec.k=44 #Uta
spec.k=56

#-------------------------
#map
#gbifmap(occ)

#library(ggmap)

#set up map
bbox <- ggmap::make_bbox(decimalLongitude, decimalLatitude, occ, f = 0.1)

map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
map1=ggmap(map_loc, margins=FALSE) #

map1 +geom_point(data=occ, aes(y=decimalLatitude, x=decimalLongitude) ) + coord_cartesian() 

#Using thet GBIF map web tile service, making a raster and visualizing it
#x <- map_fetch(search = "taxonKey", id = 3118771, year = 2010)
#library(raster)
#plot(x)

#======================================================
#Exploratory plots
#Plot prior and localities

#Performance Curve Function from Deutsch et al. 2008
TPC= function(T,Topt=33,CTmin=10.45, CTmax=42.62, bound="mean"){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #set negetative to zero
  F[F<0]<-0
  #set to 1 above or below min or max
  if(bound=="max") F[T<=Topt]=1
  if(bound=="min") F[T>=Topt]=1
  
  return(F)
}

#======================================================
# Use Worldclim bioclimatic variables (getData function in R raster library). 
BClim = getData("worldclim", var="bio", res=2.5)

#-----------------
#Try Kearney microcliamte data

#load microclim data
#currently use mean across hours

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/data/microclim/0_shade/")
#use july for max
temp= brick("TA1cm_soil_0_7.nc")
tmax_0= mean(temp) #or max
#use jan for min
temp= brick("TA1cm_soil_0_1.nc")
tmin_0= mean(temp)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/data/microclim/50_shade/")
#use july for max
temp= brick("TA1cm_soil_50_7.nc")
tmax_50= mean(temp)
#use jan for min
temp= brick("TA1cm_soil_50_1.nc")
tmin_50= mean(temp)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/data/microclim/100_shade/")
#use july for max
temp= brick("TA1cm_soil_100_7.nc")
tmax_100= mean(temp)
#use jan for min
temp= brick("TA1cm_soil_100_1.nc")
tmin_100= mean(temp)

#======================================================
#set up file
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/")
pdf("priors_microclim.pdf", height = 10, width = 12)

par(mfrow=c(4,4), cex=1.2, mar=c(3, 3, 1.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

for(spec.k in 1:58){ #nrow(dat)
  
#load localities
  setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/GBIF/")
  filename<-paste("GBIFloc_", dat$spec[spec.k],".csv", sep="")
  
  #check file size
  if( file.size(filename) > 3 ){
  locs<- read.csv(filename)
  
  #restrict to points with lat and lon
  locs<- locs[which(!is.na(locs$"decimalLongitude") & !is.na(locs$"decimalLatitude"))  ,]
  
#crop to limits
#crop to observed range 
ext = extent(rbind(range(locs$decimalLongitude), range(locs$decimalLatitude))) # define the extent
# extent
ext[1]= ext[1]-10; ext[2]= ext[2]+10; ext[3]=ext[3]-10; ext[4]=ext[4]+10

clim = crop(BClim, ext)
tmax0=  crop(tmax_0, ext)
tmin0=  crop(tmin_0, ext)
tmax50=  crop(tmax_50, ext)
tmin50=  crop(tmin_50, ext)
tmax100=  crop(tmax_100, ext)
tmin100=  crop(tmin_100, ext)
#rescale temp
clim=clim/10

#Set up prior
CTmin1= dat$tmin[spec.k]
CTmax1= dat$tmax[spec.k]
#approximate Topt, but fix based on data
Topt= CTmin1+ (CTmax1-CTmin1)*0.7
#-----
# sun to shade
# thermoregulation scenario

#max
tmax0.dif= abs(tmax0 - Topt) 
tmax50.dif= abs(tmax50 - Topt) 
tmax100.dif= abs(tmax100 - Topt) 
tmax.dif= stack(tmax0.dif, tmax50.dif, tmax100.dif)
tr.ind= which.min(tmax.dif)

tr<- tmax0
tr[]<-NA
tr[tr.ind==1]<- tmax0[tr.ind==1]
tr[tr.ind==2]<- tmax50[tr.ind==2]
tr[tr.ind==3]<- tmax100[tr.ind==3]
trmax=tr

#min
tmin0.dif= abs(tmin0 - Topt) 
tmin50.dif= abs(tmin50 - Topt) 
tmin100.dif= abs(tmin100 - Topt) 
tmin.dif= stack(tmin0.dif, tmin50.dif, tmin100.dif)
tr.ind= which.min(tmin.dif)

tr<- tmin0
tr[]<-NA
tr[tr.ind==1]<- tmin0[tr.ind==1]
tr[tr.ind==2]<- tmin50[tr.ind==2]
tr[tr.ind==3]<- tmin100[tr.ind==3]
trmin=tr

#----
#Calculate priors
#bioclim
mean.prior= calc(clim$bio1, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean")
max.prior= calc(clim$bio10, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="max") #bio10, bio5
min.prior= calc(clim$bio11, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="min") #bio11, bio6

#microclim
# mmean.prior0= calc(tmax0, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
# mmax.prior0= calc(tmax0, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="max") 
# mmin.prior0= calc(tmin0, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="min") 
mmean.prior50= calc(tmax50, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
mmax.prior50= calc(tmax50, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="max") 
mmin.prior50= calc(tmin50, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="min") 
# mmean.prior100= calc(tmax100, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
# mmax.prior100= calc(tmax100, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="max") 
# mmin.prior100= calc(tmin100, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="min") 

#thermoregulation functions
trmax.prior= calc(trmax, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
trmin.prior= calc(trmin, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 

#-----
#Plot 
# plot(clim$bio1, main=dat$spec[spec.k])
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
# 
# plot(mean.prior)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
# 
# plot(max.prior)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
# 
# plot(min.prior)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")

# plot(mmean.prior0)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
# plot(mmin.prior0)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")

plot(mmean.prior50)
points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
plot(mmin.prior50)
points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")

# plot(mmean.prior100)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
# plot(mmin.prior100)
# points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")

plot(trmax.prior)
points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
plot(trmin.prior)
points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")

} #check empty files
} #end loop speices

dev.off()

#====================
#run GRaF

mc_tr= stack(tmax50, trmin, trmax)
names(mc_tr)=c('tmax50', 'trmin', 'trmax')
#----------------------------
#generate pseudo absence

# define circles with a radius of 50 km around the subsampled points
x = circles(locs[,c("decimalLongitude","decimalLatitude")], d=50000, lonlat=TRUE)
# draw random points that must fall within the circles in object x
bg = spsample(x@polygons, 1000, type='random', iter=100)

#----
# extract environmental values
#occ_bc = extract(BClim, occ[,c("decimalLongitude","decimalLatitude")] ) # for the subsampled presence points
#bg_bc = extract(BClim, bg) # for the pseudo-absence points
occ_bc = extract(mc_tr, locs[,c("decimalLongitude","decimalLatitude")] ) # for the subsampled presence points
bg_bc = extract(mc_tr, bg) # for the pseudo-absence points

occ_bc = data.frame(lon=locs$decimalLongitude, lat=locs$decimalLatitude, occ_bc)
bgpoints = bg@coords
colnames(bgpoints) = c("lon","lat")
bg_bc = data.frame(cbind(bgpoints,bg_bc))

# Create dataframe from bioclim and presense/absance.
pres<-rep(1,dim(occ_bc)[1])
temp1<-data.frame(pres,occ_bc[,3:5])
pres<-rep(0,dim(bg_bc)[1])
temp2<-data.frame(pres,bg_bc[,3:5])
df<-rbind(temp1,temp2)
head(df,5)

#--------------------------------
# Implement Gaussian Random Fields

#covs <- df[, c("pres","bio1", "bio10","bio11")]#Pick variables #"bio5","bio6"
covs <- na.omit(df)

## 75% of the sample size
smp_size <- floor(0.75 * nrow(covs))
set.seed(123)
train_ind <- sample(seq_len(nrow(covs)), size = smp_size)
train <- covs[train_ind, ]
test <- covs[-train_ind, ]

pa_tr <- train$pres
pa_te <- test$pres
m1 <- graf(pa_tr, train[,2:ncol(train)])
pred_df<-data.frame(predict(m1,test[,2:ncol(train)]))

#plot
plot(m1)

#-------
#with priors

y= train[,1]
x= as.data.frame(train[,2:4])

# x in mean, max, min
e.prior= function(x, CTmin= CTmin1, CTmax= CTmax1){ 
  Topt= CTmin+ (CTmax-CTmin)*0.7
  P1= TPC(x[,1], Topt, CTmin, CTmax, bound="mean")
  P2= TPC(x[,2], Topt, CTmin, CTmax, bound="mean")
  P3= TPC(x[,3], Topt, CTmin, CTmax, bound="mean")
  
  return(cbind(P1,P2,P3))
}

#run model
m3 <- graf(y, x,opt.l = FALSE, prior = e.prior )
## fails for multiple variables

#-------------
#switch to one varialbe

y1= train[,1]
x1= as.data.frame(train[,2])

# x in mean, max, min
e.prior= function(x, CTmin= CTmin1, CTmax= CTmax1){ 
  Topt= CTmin+ (CTmax-CTmin)*0.7
  P1= TPC(x[,1], Topt, CTmin, CTmax, bound="mean")
  return(P1)
}

#run model
m3 <- graf(y1, x1,opt.l = FALSE, prior = e.prior )
plot(m3)

plot(x1[,1], e.prior(x1))

#---------------------
#Plot

#GP map
#https://github.com/goldingn/gp_sdm_paper/blob/master/figures/fig4.R

#predict currently not running
pred_me = predict(m3, mc_tr) # generate the predictions
# make a nice plot
plot(pred_me, 1, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main="Predicted presence of the species")
map("state", xlim=c(-119, -110), ylim=c(33.5, 38), fill=F, col="cornsilk", add=T)

# presence points
points(locs$lon, locs$lat, pch=20, cex=0.5, col="darkgreen")
# pseud-absence points
points(bg, cex=0.5, col="darkorange3")

# add axes
axis(1,las=1)
axis(2,las=1)

#============================================
#Maps with GRaF

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/")
pdf("priorsandGrAF.pdf", height = 10, width = 12)

par(mfrow=c(5,4), cex=1.2, mar=c(3, 3, 1.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

for(spec.k in 1:58){ #nrow(dat)
  
  #load localities
  setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/GBIF/")
  filename<-paste("GBIFloc_", dat$spec[spec.k],".csv", sep="")
  
  #check file size
  if( file.size(filename) > 3 ){
    locs<- read.csv(filename)
    
    #restrict to points with lat and lon
    locs<- locs[which(!is.na(locs$"decimalLongitude") & !is.na(locs$"decimalLatitude"))  ,]
    
    #crop to limits
    #crop to observed range 
    ext = extent(rbind(range(locs$decimalLongitude), range(locs$decimalLatitude))) # define the extent
    # extent
    ext[1]= ext[1]-10; ext[2]= ext[2]+10; ext[3]=ext[3]-10; ext[4]=ext[4]+10
    
    #clim = crop(BClim, ext)
    tmax0=  crop(tmax_0, ext)
    tmin0=  crop(tmin_0, ext)
    tmax50=  crop(tmax_50, ext)
    tmin50=  crop(tmin_50, ext)
    tmax100=  crop(tmax_100, ext)
    tmin100=  crop(tmin_100, ext)
    ##rescale temp
    #clim=clim/10
    
    #Set up prior
    CTmin1= dat$tmin[spec.k]
    CTmax1= dat$tmax[spec.k]
    #approximate Topt, but fix based on data
    Topt= CTmin1+ (CTmax1-CTmin1)*0.7
    #-----
    # sun to shade
    # thermoregulation scenario
    
    #max
    tmax0.dif= abs(tmax0 - Topt) 
    tmax50.dif= abs(tmax50 - Topt) 
    tmax100.dif= abs(tmax100 - Topt) 
    tmax.dif= stack(tmax0.dif, tmax50.dif, tmax100.dif)
    tr.ind= which.min(tmax.dif)
    
    tr<- tmax0
    tr[]<-NA
    tr[tr.ind==1]<- tmax0[tr.ind==1]
    tr[tr.ind==2]<- tmax50[tr.ind==2]
    tr[tr.ind==3]<- tmax100[tr.ind==3]
    trmax=tr
    
    #min
    tmin0.dif= abs(tmin0 - Topt) 
    tmin50.dif= abs(tmin50 - Topt) 
    tmin100.dif= abs(tmin100 - Topt) 
    tmin.dif= stack(tmin0.dif, tmin50.dif, tmin100.dif)
    tr.ind= which.min(tmin.dif)
    
    tr<- tmin0
    tr[]<-NA
    tr[tr.ind==1]<- tmin0[tr.ind==1]
    tr[tr.ind==2]<- tmin50[tr.ind==2]
    tr[tr.ind==3]<- tmin100[tr.ind==3]
    trmin=tr
    
    #----
    #Calculate priors
    
    #microclim
    mmean.prior50= calc(tmax50, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
   
    #thermoregulation functions
    trmax.prior= calc(trmax, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
    trmin.prior= calc(trmin, fun=TPC, CTmin= CTmin1, CTmax=CTmax1, bound="mean") 
    
    #-----
    #Plot 
    plot(mmean.prior50)
    points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
   
    plot(trmax.prior)
    points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
    #plot(trmin.prior)
    #points(locs$"decimalLongitude", locs$"decimalLatitude", pch=20, cex=0.5, col="darkgreen")
    
    #run GRaF
    mc_tr= stack(tmax50, trmin, trmax)
    names(mc_tr)=c('tmax50', 'trmin', 'trmax')
    #----------------------------
    #generate pseudo absence
    
    # define circles with a radius of 50 km around the subsampled points
    x = circles(locs[,c("decimalLongitude","decimalLatitude")], d=50000, lonlat=TRUE)
    # draw random points that must fall within the circles in object x
    bg = spsample(x@polygons, 1000, type='random', iter=100)
    
    #----
    # extract environmental values
    occ_bc = extract(mc_tr, locs[,c("decimalLongitude","decimalLatitude")] ) # for the subsampled presence points
    bg_bc = extract(mc_tr, bg) # for the pseudo-absence points
    
    occ_bc = data.frame(lon=locs$decimalLongitude, lat=locs$decimalLatitude, occ_bc)
    bgpoints = bg@coords
    colnames(bgpoints) = c("lon","lat")
    bg_bc = data.frame(cbind(bgpoints,bg_bc))
    
    # Create dataframe from bioclim and presense/absance.
    pres<-rep(1,dim(occ_bc)[1])
    temp1<-data.frame(pres,occ_bc[,3:5])
    pres<-rep(0,dim(bg_bc)[1])
    temp2<-data.frame(pres,bg_bc[,3:5])
    df<-rbind(temp1,temp2)
    head(df,5)
    
    #--------------------------------
    # Implement Gaussian Random Fields
    
    #covs <- df[, c("pres","bio1", "bio10","bio11")]#Pick variables #"bio5","bio6"
    covs <- na.omit(df)
    
    ## 75% of the sample size
    smp_size <- floor(0.75 * nrow(covs))
    set.seed(123)
    train_ind <- sample(seq_len(nrow(covs)), size = smp_size)
    train <- covs[train_ind, ]
    test <- covs[-train_ind, ]
    
    #-------------
    
    y1= train[,1]
    x1= as.data.frame(train[,2])
    
    # x in mean, max, min
    e.prior= function(x, CTmin= CTmin1, CTmax= CTmax1){ 
      Topt= CTmin+ (CTmax-CTmin)*0.7
      P1= TPC(x[,1], Topt, CTmin, CTmax, bound="mean")
      return(P1)
    }
    
    #no prior
    m1 <- graf(y1, x1, opt.l = TRUE)
    plot(m1)
    
    #with prior, one predictor
    m3 <- graf(y1, x1, prior = e.prior, opt.l = TRUE) #opt.l = TRUE ## adjust lengthscale l = 100,
    plot(m3, prior=TRUE)
    
  } #check empty files
} #end loop speices

dev.off()

#====================









