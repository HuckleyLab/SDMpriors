library(rgbif)

res = read.csv("all_species.csv", header=FALSE)

get.centroid = function(spec){
	occs = occ_data(scientificName = spec, limit=1000, minimal=TRUE)$data
	if (is.null(occs)){
		print(paste("No occurrence data found for ", spec))
	}
	cent.lat = mean(occs$decimalLatitude)
	cent.lon = mean(occs$decimalLongitude)
	return(c(cent.lat, cent.lon))
}

lapply(res[10:20,2], get.centroid)