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
source("gbifOccurrence.R")
source("climPresAbs.R")

runModels = function(physdata, speciesIdx, climVars, occs, plots=FALSE, output=TRUE, suffix=NULL){
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
    
	specName = physdata$spec[speciesIdx]

    print(paste(specName, speciesIdx))

    if (!is.null(suffix)){
        outputSuffix = paste(gsub(" ", "-", specName), "-", suffix, sep="")
    } else {
        outputSuffix = paste(gsub(" ", "-", specName), "-", "intermed", sep="")
    }

    ## get data from GBIF

    if(is.null(occs)) {
        print(paste("No occurence information for ", specName, ', SKIPPING.'))
    	return(list(specName, NaN, NaN, NaN, NaN))
    }
    
    
    ## generate presence + absence
    clim_PresAbs = gbifclimPresAbs(occs)
    if (is.null(clim_PresAbs)){
        print(paste("No climate data for ", specName, ", SKIPPING"))
        return(list(specName, NaN, NaN, NaN, NaN))
    }
    
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
    if(!is.null(plots)){
        pdf(paste(outputSuffix, "-eModel.pdf", sep=""))
        par(mfrow=c(1,length(climVars)-1)) 
        plot(eModel, prior=T)
        title(main=paste(outputSuffix, "eModel", sep=""))
        dev.off()
    }

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
    if(!is.null(plots)){
        pdf(paste(outputSuffix, "-ct.model.pdf", sep=""))
        par(mfrow=c(1,length(climVars)-1)) 
        plot(ct.model, prior=T)
        title(main=paste(outputSuffix, "ct.model", sep=""))
        dev.off()
    }

    # MaxEnt model
    me = maxent(trainCovs , p=trainPres)
    if(!is.null(plots)){
        pdf(paste(outputSuffix, "-MaxEnt.pdf", sep=""))
        par(mfrow=c(1,1)) 
        response(me)
        title(main=paste(outputSuffix, "MaxEnt", sep=""))
        dev.off()
    }

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

	write.table(t(matrix(unlist(intermed_results))), paste(outputSuffix, ".csv", sep=""), sep=",", col.names=FALSE)
	print(paste("finished index", speciesIdx))
	return(intermed_results)
}
