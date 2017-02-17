#load libraries
library(dismo)
library(plyr)
library(rgbif)
library(GRaF)
library(pROC)

#for cleaning data
#https://cran.r-project.org/web/packages/biogeo/index.html; 
#https://cran.r-project.org/web/packages/rgeospatialquality/.
library(biogeo)
library(rgeospatialquality)

#--------------------------------
# load Sunday database
dat= read.csv("Sundayetal_thermallimits.csv")
#start with reptiles and amphibians
dat= subset(dat, dat$phylum=="Chordata")
#start with species with Tmin and Tmax
dat= dat[!is.na(dat$tmax) & !is.na(dat$tmin),]
#subset to critical rather than lethal

#species name to enable match
dat$spec = gsub("_"," ",dat$species)

# loop through species
#for(spec.k in 1:nrow(dat)){
spec.k=1

#------------------------
# Query GBIF (R rgbif package) for specimen localities

#look up species
key <- name_suggest(q=dat$spec[spec.k], rank='species')$key[1]
occ=occ_search(taxonKey=key, limit=2000, return="data") 
#fields=c('name','basisOfRecord','protocol')
#return: can get metadata, etc.

#map
gbifmap(occ)

#---------------------------
#clean up data

#http://onlinelibrary.wiley.com/doi/10.1111/ecog.02118/abstract
#errorcheck
#quickclean
#geo2envpca

#http://rpubs.com/jotegui/rgeospatialquality_rgbif

#---------------------------
# Use Worldclim bioclimatic variables (getData function in R raster library). They are widely used and easily available.

# Implement Gaussian Random Fields


# Compare SDM from above to SDM without physiological data and standard SDMS 

#===================================
#GENERATE PSEUDO-ABSENCE

require(raster)
# define circles with a radius of 50 km around the subsampled points
x = circles(subs[,c("lon","lat")], d=50000, lonlat=T)
# draw random points that must fall within the circles in object x
bg = spsample(x@polygons, 100, type='random', iter=100)
BClim = getData("worldclim", var="bio", res=2.5, path="data/")
YbrevRange = extent(-119.25,-112.75,33.25,38.25) # define the extent
BClim = crop(BClim, YbrevRange)
writeRaster(BClim, filename="./data/YbrevBC_2.5.grd", overwrite=T)
BClim = brick("data/YbrevBC_2.5.grd")
# this format plots 1 (of 19) variables stored in BClim
plot(BClim, 1, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main="Annual mean temperature (ºC x 10)")
map("state", xlim=c(-119, -110), ylim=c(33.5, 38), fill=F, col="cornsilk", add=T)
# state names
text(x=-117.5, y=35.5, "California", col=rgb(1,1,1,0.6), cex=3)
text(x=-116, y=37.5, "Nevada", col=rgb(1,1,1,0.6), cex=3)
text(x=-113, y=34.5, "Arizona", col=rgb(1,1,1,0.6), cex=3)
text(x=-113, y=37.75, "Utah", col=rgb(1,1,1,0.6), cex=3)
# plot the presence points
points(locs$lon, locs$lat, pch=20, cex=2, col="darkgreen")
# and the pseudo-absence points
points(bg, cex=0.5, col="darkorange3")
# add axes
axis(1,las=1)
axis(2,las=1)
box()

#---------------------------
#EXTRACT BIOCLIM VARIABLES FOR LOCATIONS

# pulling bioclim values
Ybrev_bc = extract(BClim, subs[,c("lon","lat")]) # for the subsampled presence points
bg_bc = extract(BClim, bg) # for the pseudo-absence points
Ybrev_bc = data.frame(lon=subs$lon, lat=subs$lat, Ybrev_bc)
bgpoints = bg@coords
colnames(bgpoints) = c("lon","lat")
bg_bc = data.frame(cbind(bgpoints,bg_bc))
length(which(is.na(bg_bc$bio1))) # double-check for missing data
## [1] 0
bg_bc = bg_bc[!is.na(bg_bc$bio1), ] # and pull out the missing lines

#---------------------------------

#SDM

#covs <- df[1:1037, c("pres","bio1", "bio12")]# Not sure these are the best variables.
covs <- df

## 75% of the sample size
smp_size <- floor(0.75 * nrow(covs))
set.seed(123)
train_ind <- sample(seq_len(nrow(covs)), size = smp_size)
train <- covs[train_ind, ]
test <- covs[-train_ind, ]

pa_tr <- train$pres
pa_te <- test$pres
m1 <- graf(pa_tr, train[,2:20])
pred_df<-data.frame(predict(m1,test[,2:20]))

#print(paste("Area under ROC with No knowledge of thermal niche : ",auc(pa_te, pred_df$posterior.mode)))

thresh <- function(x) ifelse(x$bio1 < 150 | x$bio1 > 350 ,0.3, 0.6)

# fit the model, optimising the lengthscale
# fit a linear model
m.lin <- glm(pa_tr ~ bio1, data=train, family = binomial)
# wrap the predict method up in a new function
lin <- function(temp) predict(m.lin, temp, type = "response")
m3 <- graf(pa_tr, train[, c(2,6,7), drop = FALSE],opt.l = TRUE, prior = lin)
pred_df<-data.frame(predict(m3,test[, c(2,6,7), drop = FALSE]))
print(paste("Area under ROC with prior knowledge of thermal niche : ",auc(pa_te, pred_df$posterior.mode)))


