#load libraries
library(dismo)  #see also zoon R package?
library(plyr)
library(rgbif)
library(GRaF)  #see methods paper here: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12523/pdf
library(pROC)
library(ROCR)
 
#for cleaning data
library(biogeo) #https://cran.r-project.org/web/packages/biogeo/index.html
library(rgeospatialquality) #https://cran.r-project.org/web/packages/rgeospatialquality/

#tutorial here: https://rawgit.com/goldingn/intecol2013/master/tutorial/graf_workshop.htm
#--------------------------------
# load Sunday database
setwd("~/Dropbox/Projects/UW/SDMpriors/")
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
spec.k=55 #Sceloporus occidentalis

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
nrow(occ)
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
#ADD also run as presence only

# define circles with a radius of 50 km around the subsampled points
x = circles(occ[,c("decimalLongitude","decimalLatitude")], d=50000, lonlat=T)
x
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

covs <- df[, c("pres","bio1", "bio5","bio6")]#Pick variables
#covs <- df
head(covs, 5)
#divide var by 10
covs[,2:ncol(covs)]= covs[,2:ncol(covs)]/10
head(covs)

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

par(mfrow = c(1, 3))
plot(m1, prior=TRUE)
#---------------------------------------------
#establish function incorporating priors

thresh <- function(x) ifelse(x$bio1 > dat$tmin[spec.k] & x$bio1 < dat$tmax[spec.k] ,0.6, 0.1)

# fit the model, optimising the lengthscale
# fit a linear model
m.lin <- glm(pa_tr ~ bio1, data=train, family = binomial)

# wrap the predict method up in a new function
lin <- function(temp) predict(m.lin, temp, type = "response")

m3 <- graf(pa_tr, train[,2:ncol(train), drop = FALSE],opt.l = TRUE, prior = lin)
m4 <- graf(pa_tr, train[,2:ncol(train), drop = FALSE],opt.l = TRUE)
par(mfrow = c(1, 3))

plot(m3, prior=TRUE)


#NEW THRESHOLD FUNCTIONS

#threshold using bio1
thresh <- function(x) ifelse(x > dat$tmin[spec.k] & x$bio1 < dat$tmax[spec.k] ,0.6, 0.1)
m3 <- graf(pa_tr, train[,2, drop = FALSE],opt.l = TRUE, prior = thresh)
par(mfrow=c(1,3)) 
plot(m3, prior=T)

#normal curve using bio1
n.mean= (dat$tmin[spec.k]+dat$tmax[spec.k])/2
n.sd= (dat$tmax[spec.k] - dat$tmin[spec.k])/2/3 #set CTs as 3 sds

n.prior= function(x) dnorm(x[,1], mean=n.mean, sd=n.sd)* (1/dnorm(n.mean, mean=n.mean, sd=n.sd)) #normalized to peak at 1
#plot(1:60, n.prior(1:60))

m3 <- graf(pa_tr, train[,2, drop = FALSE],opt.l = TRUE, prior = n.prior) #drop=FALSE maintains matrix
### ERROR previously, but works with S. occidentalis

#plot prior
plot(m3, prior = TRUE)

#-------------------------------------------
#Threshold with multiple envi variables, needs fixing
#Threshold with bio5 and bio6
e.max<-function(x) ifelse(x<dat$tmax[spec.k], 0.8, 0.1) #max  
e.min<-function(x) ifelse(x<dat$tmin[spec.k], 0.1, 0.8) #min

#exponential based on bio5 and bio6
#start decline ten degrees above / below 
e.max<-function(x) ifelse(x<dat$tmax[spec.k]-10, 0.8, exp(-(x-dat$tmax[spec.k]+10)/5)) #max  
e.min<-function(x) ifelse(x<dat$tmin[spec.k]   , 0.1, 1- exp(-(x-(dat$tmax[spec.k])/10000) ) ) #min fix
                            

#bio1: mean, bio5:max, bio6:min (THIS IS A BAD PRIOR -- FOR TESTING)
e.prior = function(x) e.max(x[,2]) * e.min(x[,3])

m3 <- graf(pa_tr, train[,2:4, drop = FALSE],opt.l = TRUE, prior=e.prior)
plot(m3, prior=TRUE)


#---------------------------------------------------

pred_df<-data.frame(predict(m3,test[,2:ncol(train), drop = FALSE]))
print(paste("Area under ROC with prior knowledge of thermal niche : ",auc(pa_te, pred_df$posterior.mode)))

#Plot response curves
plot(m3, prior=TRUE)
m3
#plot3d(m3)

#--------------------------------
# Compare SDM from above to SDM without physiological data and standard SDMS 

prob <- pred_df$posterior.mode
pred <- prediction(prob, pa_te)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]

roc.data <- data.frame(fpr=unlist(perf@x.values),
                       tpr=unlist(perf@y.values),
                       model="GP")
plot(roc.data$fpr,roc.data$tpr,type="l",col="red",ylab="TPR",xlab="FPR",main="ROC for GP vs MaxEnt",lwd=3.5)

group_p = kfold(occ_bc, 5) # vector of group assignments splitting the Ybrev_bc into 5 groups
group_a = kfold(bg_bc, 5) # ditto for bg_bc

test = 3

train_p = occ_bc[group_p!=test, c("lon","lat")]
train_a = bg_bc[group_a!=test, c("lon","lat")]
test_p = occ_bc[group_p==test, c("lon","lat")]
test_a = bg_bc[group_a==test, c("lon","lat")]

me = maxent(BClim, p=train_p, a=train_a) #modify variables incorporated in maxent model
e = evaluate(test_p, test_a, me, BClim)

print(e)

#------------------------------------
#ROC plot

probs_me<-c(e@presence,e@absence)
class_me<-c(rep(1,length(e@presence)),rep(0,length(e@absence)))
pred_me <- prediction(probs_me, class_me)
perf_me <- performance(pred_me, measure = "tpr", x.measure = "fpr")
auc_me <- performance(pred_me, measure = "auc")
auc_me <- auc_me@y.values[[1]]

roc.data_me <- data.frame(fpr=unlist(perf_me@x.values),
                          tpr=unlist(perf_me@y.values),
                          model="ME")

lines(roc.data_me$fpr,roc.data_me$tpr,type="l",col="green",lwd=3.5)

legend(0.6,0.4, # places a legend at the appropriate place 
       c("GP","MaxEnt"), # puts text in the legend
       
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       
       lwd=c(2.5,2.5),col=c("red","green"))

#------------------------------------
#Maxent map

pred_me = predict(me, BClim) # generate the predictions
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

# restore the box around the mapaxi
box()

#-----------------------------------------
#GP map
#https://github.com/goldingn/gp_sdm_paper/blob/master/figures/fig4.R

#predict currently not running
pred_me = predict(m3, BClim) # generate the predictions
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



