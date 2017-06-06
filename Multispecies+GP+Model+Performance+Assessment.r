
library(dismo)  #see also zoon R package?
library(plyr)
library(rgbif)
library(GRaF)  #see methods paper here: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12523/pdf
library(pROC)
library(ROCR)
library(foreach)
library(doMC)
library(rJava)
library(optparse)
 
registerDoMC(cores=10)

phys = read.csv("Sundayetal_thermallimits.csv")
phys= phys[!is.na(phys$tmax) & !is.na(phys$tmin),]
phys$spec = gsub("_", " ", phys$species)
paste("We've got Tmin and TMax for", nrow(phys), "species.")

climVars = c("presence", "bio1", "bio5", "bio6")


runModels = function(physdata, speciesIdx){
    #
    # the purpose of this function is to run a series of 
    # species distribution models given physiological information
    # to inform them for a single species.
    #
    # It should really be 4 or 5 different functions. 
    #
    # each time it runs it produces a .csv file in the current
    # directory with {speciesname}-intermed.csv as the filename
    # which contains the results from the model runs
    # in order. 

    specName = phys$spec[speciesIdx]
    print(paste(specName, speciesIdx))
    ## get data from GBIF
    occs = occ_data(scientificName = specName, limit=1000, minimal=TRUE)$data
    if(is.null(occs)) {
        print(paste("No occurence information for ", specName, ', SKIPPING.'))
    	return(list(specName, NaN, NaN, NaN, NaN))
    }
    occs = occs[which(!is.na(occs$"decimalLongitude") & !is.na(occs$"decimalLatitude")),]
    
    
        
    ## generate presence + absence
    bufs = circles(occs[,c("decimalLongitude", "decimalLatitude")], d=50000, lonlat=TRUE)
    abs  = spsample(bufs@polygons, 100, type='random', iter=100)
    
    ## get climate data
    BClim = getData("worldclim", var='bio', res=2.5)
    specExt = extent(rbind(range(occs$decimalLongitude), range(occs$decimalLatitude)))
    BClim = crop(BClim, specExt)
    
    ## assign data to each presence ...
    clim_Pres = extract(BClim, occs[,c("decimalLongitude", "decimalLatitude")])
    if (all(is.na(clim_Pres))) {
        print(paste("No climate data for ", specName, ", SKIPPING"))
        return(list(specName, NaN, NaN, NaN, NaN))
    }

    clim_Pres = data.frame(lon=occs$decimalLongitude,
                           lat=occs$decimalLatitude,
                           clim_Pres)
    ## ..and absence point.
    clim_Abs  = extract(BClim, abs)
    clim_Abs  = data.frame(lon=abs@coords[,'x'], lat=abs@coords[,'y'], clim_Abs)
        
    ## assign binary presence (1) absence (0) flags.
    presence = rep(1,dim(clim_Pres)[1])
    presence_temp = data.frame(presence, clim_Pres[,3:ncol(clim_Pres)])
    presence = rep(0, dim(clim_Abs)[1])
    absence_temp = data.frame(presence, clim_Abs[,3:ncol(clim_Abs)])
    

    ## and combine them. 
    clim_PresAbs = rbind(presence_temp, absence_temp)
    
    ## extract and transform relevant climate information
    covs = clim_PresAbs[, climVars]
    covs[,2:ncol(covs)] = covs[,2:ncol(covs)]/10 ## (BClim needs to be divided by 10)
    # omit rows with NA
    covs = na.omit(covs)
    
    ## split data into train and test sets.
    train_size = floor(0.75*nrow(covs))
    trainPres = NULL
    testPres = NULL
    ## make sure that there are 2 classes in the test data
    numTries = 0 
    while (!(length(unique(testPres)) > 1 && length(unique(trainPres)) > 1) && numTries < 10){
        train_idxs = sample(seq_len(nrow(covs)), size=train_size)
        trainData = covs[train_idxs,]
        trainPres = trainData[,1]
        trainCovs = trainData[,2:ncol(trainData)]
        testData  = covs[-train_idxs,]
        testPres = testData[,1]
        testCovs = testData[,2:ncol(trainData)]
        numTries =+ 1
    }

    if (!(length(unique(testPres)) > 1 && length(unique(trainPres)) > 1)){
        print(paste("Bad train/test split for", specName, ", SKIPPING"))
        return(list(specName, NaN, NaN, NaN, NaN))

    }
    ## BUILD MODELS

    # simple model
    simple = graf(trainPres, trainCovs, opt.l=T)

    # threshold model 

    # define threshold functions and threshold prior.
    e.max<-function(x) ifelse(x<physdata$tmax[speciesIdx]-10, 0.9, exp(-(x-physdata$tmax[speciesIdx]+10)/5)) #max  
    e.min<-function(x) ifelse(x<physdata$tmin[speciesIdx]   , 0.1, 1- exp(-(x-(physdata$tmax[speciesIdx])/10000) ) ) #min fix
    e.prior = function(x) e.max(x[,2]) * e.min(x[,3])

    if (any(is.na(qnorm(e.prior(trainCovs))))){
        print(paste("Error in prior for species ", specName, ", SKIPPING"))
        return(list(specName, NaN, NaN, NaN, NaN))

#         results = append(results, 0)
        #return(0)
    }

    eModel = graf(trainPres, trainCovs, prior = e.prior, opt.l=T) 

    # normal prior 
    ct.mean = physdata$tmax[speciesIdx] - physdata$tmin[speciesIdx]
    ct.std  = ct.mean / 2
    ct.thresh <- function(x) ifelse(x$bio6 > physdata$tmin[speciesIdx] & x$bio5 < physdata$tmax[speciesIdx], 1, .00001)
    ct.prob = function(x) {
        dnorm(x$bio1, mean=ct.mean, sd = ct.std) * ct.thresh(x)
    }
    ct.cumprob = function(x){
        ct.thresh(x) * (pnorm(x$bio5, mean=ct.mean, sd=ct.std) - pnorm(x$bio6, mean=ct.mean, sd=ct.std))
    }

    ct.model = graf(trainPres, trainCovs, prior = ct.prob, opt.l=T) 

    # MaxEnt model
    me = maxent(trainCovs , p=trainPres)

    ## EVALUATE MODELS
    ## order of eval: Simple, e.prior, normal.prior, maxEnt
    intermed_results = c(specName)

    for (mod in list(simple, eModel, ct.model, me)){
	    prob = data.frame(predict(mod, testCovs))
	    if(class(mod) != "MaxEnt") prob = prob$posterior.mode
	    pred = prediction(prob, testPres)
	    auc  = performance(pred, measure='auc')
	    auc = auc@y.values[[1]]
	    intermed_results = c(intermed_results, auc)
	}
	write.table(t(matrix(unlist(intermed_results))), paste("./",gsub(" ", "-", specName),"-intermed.csv", sep=""), sep=",", col.names=FALSE)
	print(paste("finished index", speciesIdx))
	return(intermed_results)
}


# for most experiments this loop won't run all the way through. 
# Seems to get stuck at the end. Not sure why. We push through by
# segmenting input indices (1..30, 31..60, 61..90, 91..nrow(phys))

# this parallelized for loop runs the `runmodels` function for each species in the 
# physiological data file. 
results = foreach(speciesIdx=seq(1, nrow(phys)), .errorhandling = 'remove') %dopar% {
	return(runModels(phys, speciesIdx))
}



# # it's unlikely that the code reaches this point. 
# resultsDF = as.data.frame(t(matrix(unlist(results), nrow=length(unlist(results[1])))))
# colnames(resultsDF) = c("species", "simple", "eMinMax", "normdist", "maxent")
# print(resultsDF)
# write.csv(resultsDF, "./allSpecies_results.csv")


