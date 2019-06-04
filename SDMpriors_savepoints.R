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
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/presabs/")
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

#------------
# Raster to xyz

micro.brick= brick(tmax_0, tmin_0, tmax_50, tmin_50, tmax_100, tmin_100)
micro.xyz= rasterToPoints(micro.brick)
#to data frame
micro.xyz = as.data.frame(micro.xyz)
#add names
names(micro.xyz)[3:8]=c('tmax_0', 'tmin_0', 'tmax_50', 'tmin_50', 'tmax_100', 'tmin_100')

#plot to check
library(ggplot2)
ggplot(micro.xyz, aes(micro.xyz[,1],micro.xyz[,2], color=micro.xyz[,3]) )+geom_tile()

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/presabs/")
write.csv(micro.xyz, "EnviDat.csv")

#============================================
#Write out presence / absence data

#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/presabs/")
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/SDMpriors/out/presabs/")

for(spec.k in 1:58){ #nrow(dat)
  
  #load localities
  setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/GBIF/")
  filename<-paste("GBIFloc_", dat$spec[spec.k],".csv", sep="")
  
  #check file size
  if( file.size(filename) > 200 ){
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
  
    #---------
    #Species tolerance data
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
    
    #predictors  
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
    temp1<-data.frame(pres,occ_bc[1:5])
    pres<-rep(0,dim(bg_bc)[1])
    temp2<-data.frame(pres,bg_bc[1:5])
    df<-rbind(temp1,temp2)
    head(df,5)
    
    covs <- na.omit(df)
  
    #-------------
    #write out
    setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/SDMpriors/out/presabs/")
    pa.filename<-paste("PresAbs_", dat$spec[spec.k],".csv", sep="")
    write.csv(covs,pa.filename)
      
  } #check empty files
} #end loop speices

#======================
#Check PA data

keep=NA

for(spec.k in 1:58){ #nrow(dat)
  
  #load presence absence
  setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/SDMpriors/out/presabs/")
  filename<-paste("PresAbs_", dat$spec[spec.k],".csv", sep="")
  
  #check file size
  if( !is.na(file.size(filename)) ) keep=c(keep, spec.k)
} #end loop species

#write out species with PA data
write.csv(dat[keep[2:length(keep)],], 'SpeciesList_PresAbs.csv')

