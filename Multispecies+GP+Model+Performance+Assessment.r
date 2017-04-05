
library(dismo)  #see also zoon R package?
library(plyr)
library(rgbif)
library(GRaF)  #see methods paper here: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12523/pdf
library(pROC)
library(ROCR)
library(foreach)
library(doMC)
registerDoMC(cores=10)

phys = read.csv("Sundayetal_thermallimits.csv")
phys= phys[!is.na(phys$tmax) & !is.na(phys$tmin),]
phys$spec = gsub("_", " ", phys$species)
paste("We've got Tmin and TMax for", nrow(phys), "species.")

climVars = c("presence", "bio1", "bio5", "bio6")

results = foreach(speciesIdx=seq(1, 30), .errorhandling = 'pass') %dopar% {
    specName = phys$spec[speciesIdx]
    print(specName)
    ## get data from GBIF
    occs = occ_data(scientificName = specName, limit=1000, minimal=TRUE)$data
    if(is.null(occs)) {
        return(paste("No occurence information for ", specName, ', SKIPPING.'))
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
        return(paste("No climate data for ", specName, ", SKIPPING"))
#         append(results, 0)
        #return(0) 
    }

    clim_Pres = data.frame(lon=occs$decimalLongitude,
                           lat=occs$decimalLatitude,
                           clim_Pres)
    ## ..and absence point.
    clim_Abs  = extract(BClim, abs)
    clim_Abs  = data.frame(lon=abs@coords[,'x'], lat=abs@coords[,'y'], clim_Abs)
        
    presence = rep(1,dim(clim_Pres)[1])
    presence_temp = data.frame(presence, clim_Pres[,3:ncol(clim_Pres)])
    presence = rep(0, dim(clim_Abs)[1])
    absence_temp = data.frame(presence, clim_Abs[,3:ncol(clim_Abs)])
    

    ## and combine them. 
    clim_PresAbs = rbind(presence_temp, absence_temp)
    
    ## extract relevant information
    covs = clim_PresAbs[, climVars]
    covs[,2:ncol(covs)] = covs[,2:ncol(covs)]/10 ## (BClim needs to be divided by 10)
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
    if (numTries > 1) {print("Warning, many tries")}
    ## BUILD MODEL
    # define threshold functions
    e.max<-function(x) ifelse(x<phys$tmax[speciesIdx]-10, 0.9, exp(-(x-phys$tmax[speciesIdx]+10)/5)) #max  
    e.min<-function(x) ifelse(x<phys$tmin[speciesIdx]   , 0.1, 1- exp(-(x-(phys$tmax[speciesIdx])/10000) ) ) #min fix
    e.prior = function(x) e.max(x[,2]) * e.min(x[,3])

    if(any(is.na(qnorm(e.prior(trainCovs))))){
        return(paste("Error in prior for species ", specName, ", SKIPPING"))
#         results = append(results, 0)
        #return(0)
    }
        
    eModel = graf(trainPres, trainCovs, prior = e.prior)#, opt.l=T) 

    if (!(length(unique(testPres)) > 1 && length(unique(trainPres)) > 1)){
        return(paste("Bad train/test split for", specName, ", SKIPPING"))
        
        #results = append(results, 0)
        #return(0)
    }

    ## EVALUATE MODEL
    csimplePred = data.frame(predict(eModel, testCovs))
    prob = csimplePred$posterior.mode
    pred = prediction(prob, testPres)
    auc  = performance(pred, measure='auc')
    auc = auc@y.values[[1]]
    if (is.null(auc)){ return("NULL AUC")}
    return(auc)
}


print(results)


