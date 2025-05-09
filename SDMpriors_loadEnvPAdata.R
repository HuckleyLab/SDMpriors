#load libraries
library(ggplot2)

#--------------------------------
#setwd
desktop<- "y"
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")

#load species data
#for species with presence absence data, currently 48 species
spec.dat= read.csv("out/presabs/SpeciesList_PresAbs.csv")
#tmin and tmax are critical thermal limits that we are thinking to use to define the prior
#prior functions we've explored are in SDMpriors_priors.R file

#load environmental data
env.dat= read.csv("out/presabs/EnviDat.csv")
#envi variables are as follows:
#data from microclim dataset, https://www.nature.com/articles/sdata20146
#tmax0= annual daily maximum temperature with 0 shade, tmin_0= annual daily minimum temperature with 0 shade, 
#tmax50= annual daily maximum temperature with 50% shade, tmin_50= annual daily minimum temperature with 50% shade, 
#tmax100= annual daily maximum temperature with 100% shade, tmin_100= annual daily minimum temperature with 100% shade,  

#plot envi data to check
#ggplot(env.dat, aes(x,y, color=tmax50) )+geom_tile()

#============================================
#Load presence / absence data for each species

for(spec.k in nrow(spec.dat)){ 
  
  #load presence absence
  filename<-paste("out/presabs/PresAbs_", spec.dat$spec[spec.k],".csv", sep="")
  
  #load pa data
  pa<- read.csv(filename)
  #envi variables are as follows:
  #trmin= annual daily minimum temperature with thermoregulation, selecting among 0, 50, and 100% shade to get as close to species' Topt as possible
  #trmax= annual daily maximum temperature with thermoregulation, selecting among 0, 50, and 100% shade to get as close to species' Topt as possible
  #tmax0= annual daily maximum temperature with 0 shade, tmin0= annual daily minimum temperature with 0 shade, 
  #tmax50= annual daily maximum temperature with 50% shade, tmin50= annual daily minimum temperature with 50% shade, 
  #tmax100= annual daily maximum temperature with 100% shade, tmin100= annual daily minimum temperature with 100% shade,  
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




