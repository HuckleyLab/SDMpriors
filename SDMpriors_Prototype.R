# loads the dismo library
library(dismo)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
file
bradypus <- read.table(file,  header=TRUE,  sep=",")
> # inspect the values of the file
  > # first rows
head(bradypus)
bradypus <- bradypus[,2:3]
head(bradypus)
acaule = gbif("solanum", "acaule*", geo=FALSE)
data(acaule)
dim(acaule)
colnames(acaule)
acgeo <- subset(acaule, !is.na(lon) & !is.na(lat))
acgeo[1:4, c(1:5,7:10)]
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, 
     xlim=c(-80,70), 
     ylim=c(-60,60), 
     axes=TRUE, 
     col="light yellow")
box()
points(acgeo$lon, acgeo$lat, col='orange', pch=20, cex=0.75)
points(acgeo$lon, acgeo$lat, col='red', cex=0.75)


acaule[c(303,885),1:10]
lonzero = subset(acgeo, lon==0)
lonzero[, 1:13]

# which records are duplicates (only for the first 10 columns)?
dups <- duplicated(lonzero[, 1:10])
# remove duplicates
lonzero  <-  lonzero[dups, ]
lonzero[,1:13]

# differentiating by (sub) species
  > # dups2 <- duplicated(acgeo[, c('species', 'lon', 'lat')]) > # ignoring (sub) species and other naming variation
dups2 <- duplicated(acgeo[, c('lon', 'lat')])
> # number of duplicates
sum(dups2)
acg <- acgeo[!dups2, ]

i <- acg$lon > 0 & acg$lat > 0
acg$lon[i] <- -1 * acg$lon[i]
acg$lat[i] <- -1 * acg$lat[i]
acg <- acg[acg$lon < -50 & acg$lat > -50, ]
library(sp)
coordinates(acg) <- ~lon+lat
crs(acg) <- crs(wrld_simpl)
class(acg)
class(wrld_simpl)
ovr <- over(acg, wrld_simpl)
head(ovr)
cntr <- ovr$NAME
i <- which(is.na(cntr))
j <- which(cntr != acg$country)
cbind(cntr, acg$country)[j,]
plot(acg)
plot(wrld_simpl, add=T, border='blue', lwd=2)
points(acg[j, ], col='red', pch=20, cex=2)



georef <- subset(acaule, (is.na(lon) | is.na(lat)) & ! is.na(locality) )
georef$cloc[4]

b <- try(geocode(georef$cloc[4]))
b
acg
r <- raster(acg)
res(r) <- 1
r <- extend(r, extent(r)+1)
acsel <- gridSample(acg, r, n=1)
p <- rasterToPolygons(r)
plot(p, border='gray')
points(acg)
points(acsel, cex=1, col='red', pch='x')
file <- paste(system.file(package="dismo"), '/ex/acaule.csv', sep='')
acsel <- read.csv(file)


# get the file names
files <- list.files(path=paste(system.file(package="dismo"), '/ex',
                                   sep=''), pattern='grd', full.names=TRUE )
# we use the first file to create a RasterLayer
mask <- raster(files[1])
# select 500 random points
# set seed to assure that the examples will always
# have the same random sample.
set.seed(1963)
bg <- randomPoints(mask, 500 )

# set up the plotting area for two maps
par(mfrow=c(1,2))

plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)
# now we repeat the sampling, but limit
# the area of sampling using a spatial extent
e <- extent(-80, -53, -39, -22)
bg2 <- randomPoints(mask, 50, ext=e)
plot(!is.na(mask), legend=FALSE)
plot(e, add=TRUE, col='red')
points(bg2, cex=0.5)

file <- paste(system.file(package="dismo"), '/ex/acaule.csv', sep='')
ac <- read.csv(file)
coordinates(ac) <- ~lon+lat
projection(ac) <- CRS('+proj=longlat +datum=WGS84')

# circles with a radius of 50 km
x <- circles(ac, d=50000, lonlat=TRUE)
pol <- polygons(x)
# sample randomly from all circles
samp1 <- spsample(pol, 250, type='random', iter=25)
# get unique cells
cells <- cellFromXY(mask, samp1)
length(cells)
cells <- unique(cells)
length(cells)
xy <- xyFromCell(mask, cells)
#Plot to inspect the results:
plot(pol, axes=TRUE)
points(xy, cex=0.75, pch=20, col='blue')

spxy <- SpatialPoints(xy, proj4string=CRS('+proj=longlat +datum=WGS84'))
o <- over(spxy, geometry(x))
xyInside <- xy[!is.na(o), ]


# extract cell numbers for the circles
v <- extract(mask, x@polygons, cellnumbers=T)
# use rbind to combine the elements in list v
v <- do.call(rbind, v)
# get unique cell numbers from which you could sample
v <- unique(v[,1])
head(v)

# to display the results
m <- mask
m[] <- NA
m[v] <- 1
plot(m, ext=extent(x@polygons)+1)
plot(x@polygons, add=T)


files <- list.files(path=paste(system.file(package="dismo"),
                                 '/ex', sep=''), pattern='grd', full.names=TRUE )
# The above finds all the files with extension "grd" in the
# examples ("ex") directory of the dismo package. You do not
 # need such a complex statement to get your own files.
files
predictors <- stack(files)
names(predictors)
plot(predictors)


library(maptools)
data(wrld_simpl)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file, header=TRUE, sep=',')
# we do not need the first column
bradypus <- bradypus[,-1]
#And now plot:
# first layer of the RasterStack
plot(predictors, 1)
# note the "add=TRUE" argument with plot
plot(wrld_simpl, add=TRUE)
# with the points function, "add" is implicit
points(bradypus, col='blue')

presvals <- extract(predictors, bradypus)
# setting random seed to always create the same
# random set of points for this example
set.seed(0)
backgr <- randomPoints(predictors, 500)
absvals <- extract(predictors, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata[,'biome'] = as.factor(sdmdata[,'biome'])
head(sdmdata)

m1 <- glm(pb ~ bio1 + bio5 + bio12, data=sdmdata)
class(m1)
summary(m1)

m2 = glm(pb ~ ., data=sdmdata)
m2


bc <- bioclim(presvals[,c('bio1', 'bio5', 'bio12')])
class(bc)

pairs(bc)


bio1 = c(40, 150, 200)
bio5 = c(60, 115, 290)
bio12 = c(600, 1600, 1700)
pd = data.frame(cbind(bio1, bio5, bio12))
pd
predict(m1, pd)
predict(bc, pd)

names(predictors)
p <- predict(predictors, m1)
plot(p)


p <- rnorm(50, mean=0.7, sd=0.3)
a <- rnorm(50, mean=0.4, sd=0.4)
p <- rnorm(50, mean=0.7, sd=0.3)
a <- rnorm(50, mean=0.4, sd=0.4)
par(mfrow=c(1, 2))
plot(sort(p), col='red', pch=21)
points(sort(a), col='blue', pch=24)
legend(1, 0.95 * max(a,p), c('presence', 'absence'),
        pch=c(21,24), col=c('red', 'blue'))
comb = c(p,a)
group = c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb~group, col=c('blue', 'red'))



######



library(raster)
x <- raster()
x

x <- raster(ncol=36, nrow=18, xmn=-1000, xmx=1000, ymn=-100, ymx=900)
res(x)
res(x) <- 100
ncol(x)
projection(x) <- "+proj=utm +zone=48 +datum=WGS84"
x
r <- raster(ncol=10, nrow=10)
ncell(r)
hasValues(r)
values(r) <- 1:ncell(r)
set.seed(0)
values(r) <- runif(ncell(r))
hasValues(r)