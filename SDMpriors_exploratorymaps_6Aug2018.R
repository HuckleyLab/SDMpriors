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
    plot(mmean.prior50, main=dat$spec[spec.k])
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
    m1 <- graf(y1, x1, l=100)
    plot(m1)
    
    #with prior, one predictor
    m3 <- graf(y1, x1, prior = e.prior, l=100) #opt.l = TRUE ## adjust lengthscale l = 100,
    plot(m3, prior=TRUE)
    
  } #check empty files
} #end loop speices

dev.off()

#====================









