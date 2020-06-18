
gbifclimPresAbs <- function(occs){
    #
    # Given occurrence information, generates presence and absence
    # records with accompanying climate data from bioclim. 
    #
    # returns a dataframe containing a "presence" column, lat, lon, 
    # and all bioclim vars.

    bufs = circles(occs[,c("decimalLongitude", "decimalLatitude")], d=50000, lonlat=TRUE)
    abs  = spsample(bufs@polygons, 100, type='random', iter=100)
    
    ## get climate data
    BClim = getData("worldclim", var='bio', res=2.5)
    specExt = extent(rbind(range(occs$decimalLongitude), range(occs$decimalLatitude)))
    BClim = crop(BClim, specExt)
    
    ## assign data to each presence ...
    clim_Pres = extract(BClim, occs[,c("decimalLongitude", "decimalLatitude")])
    if (all(is.na(clim_Pres))) {
        return(NULL)
    }

    clim_Pres = data.frame(lon=occs$decimalLongitude,
                           lat=occs$decimalLatitude,
                           clim_Pres)
    ## ..and absence point.
    clim_Abs  = extract(BClim, abs)
    clim_Abs  = data.frame(lon=abs@coords[,'x'], lat=abs@coords[,'y'], clim_Abs)
        
    ## assign binary presence (1) absence (0) flags.
    presence = rep(1,dim(clim_Pres)[1])
    presence_temp = data.frame(presence, clim_Pres)
    presence = rep(0, dim(clim_Abs)[1])
    absence_temp = data.frame(presence, clim_Abs)

    

    ## and combine them. 
    clim_PresAbs = rbind(presence_temp, absence_temp)
    return(clim_PresAbs)
}