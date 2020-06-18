library(rgbif)

## get_centroids.R
##
# 	usage: rscript get_centroids.R <datafile> <outputfilename>
#
#   given a <datafile> with species names in the second column,
#	will use first 1000 GBIF records to compute centroid of
#	distribution and put results in a file in the same order 
#   as the input file. 
##

args = commandArgs(trailingOnly=TRUE)

res = read.csv(args[1], header=FALSE)

get.centroid = function(spec){
	occs = occ_data(scientificName = spec, limit=1000, minimal=TRUE)$data
	if (is.null(occs)){
		print(paste("No occurrence data found for ", spec))
	}
	occs = occs[which(!is.na(occs$"decimalLongitude") & !is.na(occs$"decimalLatitude")),]
	cent.lat = mean(occs$decimalLatitude)
	cent.lon = mean(occs$decimalLongitude)
	print(paste(spec, " finished", sep=""))
	return(c(cent.lat, cent.lon))
}

write.table(lapply(res[,2], get.centroid), args[2], sep=",")