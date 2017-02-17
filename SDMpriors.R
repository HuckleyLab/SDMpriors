#load libraries
library(dismo)
library(plyr)
library(rgbif)
library(GRaF)
library(pROC)
 
#for cleaning data
library(biogeo) #https://cran.r-project.org/web/packages/biogeo/index.html
library(rgeospatialquality) #https://cran.r-project.org/web/packages/rgeospatialquality/

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

occ <- occ_data(scientificName=dat$spec[spec.k], limit=1000, minimal=FALSE)
occ <- occ$data

#occ=occ_search(taxonKey=key, limit=2000, return="data") 
#fields=c('name','basisOfRecord','protocol')
#return: can get metadata, etc.

#map
gbifmap(occ)

#---------------------------
#clean up data

#restrict to points with lat and lon
occ<- occ[which(!is.na(occ$"decimalLongitude") & !is.na(occ$"decimalLatitude"))  ,]

#http://onlinelibrary.wiley.com/doi/10.1111/ecog.02118/abstract
#errorcheck(occ)
#quickclean
#geo2envpca

#USE rgeospatialquality_rgbif
#check names
"countryCode" %in% names(occ)
"scientificName" %in% names(occ)

##Add quality flags
#http://rpubs.com/jotegui/rgeospatialquality_rgbif
occ1 <- add_flags(occ)

#drop porblematic records
flags=occ1$flags
#drop several fields #REVISE
flags=flags[,-which(names(flags) %in% c("highPrecisionCoordinates","distanceToCountryInKm"))]

#keep records passing all quality checks
check= apply(flags, MARGIN=1, FUN=all, na.rm=TRUE)

occ= occ[check,]

#----------------------------
#generate pseudo absence

# define circles with a radius of 50 km around the subsampled points
x = circles(occ[,c("decimalLongitude","decimalLatitude")], d=50000, lonlat=T)
# draw random points that must fall within the circles in object x
bg = spsample(x@polygons, 100, type='random', iter=100)

#---------------------------
# Use Worldclim bioclimatic variables (getData function in R raster library). 
BClim = getData("worldclim", var="bio", res=2.5)

#crop to observed range #REVISIT
ext = extent(rbind(range(occ$decimalLongitude), range(occ$decimalLatitude))) # define the extent
BClim = crop(BClim, ext)

# pulling bioclim values
occ_bc = extract(BClim, occ[,c("decimalLongitude","decimalLatitude")] ) # for the subsampled presence points
bg_bc = extract(BClim, bg) # for the pseudo-absence points
occ_bc = data.frame(lon=occ$decimalLongitude, lat=occ$decimalLatitude, occ_bc)
bgpoints = bg@coords
colnames(bgpoints) = c("lon","lat")
bg_bc = data.frame(cbind(bgpoints,bg_bc))

# Create dataframe from bioclim and presense/absance.
pres<-rep(1,dim(occ_bc)[1])
temp1<-data.frame(pres,occ_bc[,3:21])
pres<-rep(0,dim(bg_bc)[1])
temp2<-data.frame(pres,bg_bc[,3:21])
df<-rbind(temp1,temp2)
head(df,5)

#--------------------------------
# Implement Gaussian Random Fields

covs <- df[, c("pres","bio1", "bio5","bio6")]#Pcik variables
#covs <- df

#divide var by 10
covs[,2:ncol(covs)]= covs[,2:ncol(covs)]/10

#remove NAs
covs= na.omit(covs)

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

#establish function incorporating priors
thresh <- function(x) ifelse(x$bio6 < dat$tmin[spec.k] | x$bio5 > dat$tmax[spec.k] ,0.1, 0.6)

# fit the model, optimising the lengthscale
# fit a linear model
m.lin <- glm(pa_tr ~ bio1, data=train, family = binomial)

# wrap the predict method up in a new function
lin <- function(temp) predict(m.lin, temp, type = "response")
m3 <- graf(pa_tr, train[,2:ncol(train), drop = FALSE],opt.l = TRUE, prior = lin)
pred_df<-data.frame(predict(m3,test[,2:ncol(train), drop = FALSE]))
print(paste("Area under ROC with prior knowledge of thermal niche : ",auc(pa_te, pred_df$posterior.mode)))

#--------------------------------
# Compare SDM from above to SDM without physiological data and standard SDMS 

