library(foreach)
library(doMC)
source("runModels.R")

 
registerDoMC(cores=10)

phys = read.csv("Sundayetal_thermallimits.csv")
phys= phys[!is.na(phys$tmax) & !is.na(phys$tmin),]
phys$spec = gsub("_", " ", phys$species)

climVars = c("presence", "bio1", "bio5", "bio6")



# for most experiments this loop won't run all the way through. 
# Seems to get stuck at the end. Not sure why. We push through by
# segmenting input indices (1..30, 31..60, 61..90, 91..nrow(phys))

# this parallelized for loop runs the `runmodels` function for each species in the 
# physiological data file. 
results = foreach(speciesIdx=seq(1, nrow(phys)), .errorhandling = 'remove') %dopar% {
	return(runModels(phys, speciesIdx, climVars))
}



# # it's unlikely that the code reaches this point. 
# resultsDF = as.data.frame(t(matrix(unlist(results), nrow=length(unlist(results[1])))))
# colnames(resultsDF) = c("species", "simple", "eMinMax", "normdist", "maxent")
# print(resultsDF)
# write.csv(resultsDF, "./allSpecies_results.csv")


