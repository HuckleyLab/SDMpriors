library(rgbif)

res = read.csv("all_species.csv", header=FALSE)

get.centroid = function(spec){
	occs = occ_data(scientificName = spec, limit=1000, minimal=TRUE)$data
	if (is.null(occs)){
		print(paste("No occurrence data found for ", spec))
	}
	occs = occs[which(!is.na(occs$"decimalLongitude") & !is.na(occs$"decimalLatitude")),]
	cent.lat = mean(occs$decimalLatitude)
	cent.lon = mean(occs$decimalLongitude)
	return(c(cent.lat, cent.lon))
}

write.table(lapply(res[,2], get.centroid), "centroids.csv", sep=",", header=FALSE)