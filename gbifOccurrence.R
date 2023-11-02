library(rgbif)

gbifOccurence <- function(specName){
	occs = occ_data(scientificName = specName, limit=1000, minimal=TRUE)$data
    occs = occs[which(!is.na(occs$"decimalLongitude") & !is.na(occs$"decimalLatitude")),]
    return(occs)
}