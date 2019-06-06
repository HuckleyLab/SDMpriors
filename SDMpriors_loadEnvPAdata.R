#load libraries
library(ggplot2)

#--------------------------------
#setwd
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/SDMpriors/out/presabs/")

#load species data
#for species with presence absence data, currently 48 species
spec.dat= read.csv("SpeciesList_PresAbs.csv")
#tmin and tmax are critical thermal limits that we are thinking to use to define the prior
#prior functions we've explored are in SDMpriors_priors.R file

#load environmental data
env.dat= read.csv("EnviDat.csv")
#envi variables are as follows:
#data from microclim dataset, https://www.nature.com/articles/sdata20146
#tmax_0= annual daily maximum temperature with 0 shade, tmin_0= annual daily minimum temperature with 0 shade, 
#tmax_50= annual daily maximum temperature with 50% shade, tmin_50= annual daily minimum temperature with 50% shade, 
#tmax_100= annual daily maximum temperature with 100% shade, tmin_100= annual daily minimum temperature with 100% shade,  

#plot envi data to check
#ggplot(env.dat, aes(x,y, color=tmax_50) )+geom_tile()

#============================================
#Load presence / absence data for each species

for(spec.k in nrow(spec.dat)){ 
  
  #load presence absence
  setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/SDMpriors/out/presabs/")
  filename<-paste("PresAbs_", spec.dat$spec[spec.k],".csv", sep="")
  
  #load pa data
  pa<- read.csv(filename)
  #envi variables are as follows:
  #tmax50= annual daily maximum temperature with 50% shade, corresponds to tmax_50 in env.dat
  #trmin= annual daily minimum temperature with thermoregulation, selecting among 0, 50, and 100% shade to get as close to species' Topt as possible
  #trmax= annual daily maximum temperature with thermoregulation, selecting among 0, 50, and 100% shade to get as close to species' Topt as possible
  #can make species specific envi data file for each species if desired
  
  #Species tolerance data
  #For setting up prior
  CTmin1= spec.dat$tmin[spec.k]
  CTmax1= spec.dat$tmax[spec.k]
  #approximate Topt, but fix based on data
  Topt= CTmin1+ (CTmax1-CTmin1)*0.7
    
  #plot out to check
  pa$pres= as.factor(pa$pres)
  ggplot(pa, aes(lon,lat,color=tmax50,shape=pres) ) +geom_point()
  
} #end loop species




