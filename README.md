### Disclaimer:
This is in-progress research software. We provide no warranty, guarantee, or description of functionality. All credit belongs to Tony Cannistra, Lauren Buckley, and the University of Washington. If you're interested in collaboration, contact [tonycan@uw.edu](mailto:tonycan@uw.edu)

# SDMpriors
Exploring the potential to incorporate physiological priors in species distribution models

Objectives:
We aim to test whether using laboratory thermal tolerance data (critical thermal minima and maxima) to inform priors in simple species distribution models can improve model performance in extrapolation to new time periods.

Realted Approaches:
Marine example: https://www.nature.com/articles/s41598-018-38416-3

multivariate GP spatio-temporal models (R package VAST and a list of references:  https://github.com/James-Thorson-NOAA/VAST#references). 

This Gaussian Random Field R package appears to be no longer maintained, but is a good fit to the application: https://github.com/goldingn/GRaF. Paper here: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12523.

Isaac Caruso thesis: https://github.com/icaruso21/bayes-gp-sdm-manuscript 

Source data:
* Species thermal tolerance data from Sunday et. al (2104, https://www.pnas.org/content/111/15/5610.short). Currently examining 48 species of amphibians and reptiles.
Species occurrence data from GBIF, pseudo absences were selected within a circle with a radius of 50km around the presence points
* Environmental data from microclim (https://www.nature.com/articles/sdata20146). Data are surface temperatures for three levels of shade. We also explored using Worldclim data, but surface temperatures appear to provide a better match to thermal tolerance data.
* Hindcasting data: We hope to test the models in extrapolation using resurvey data such as from the Grinnel Resurvey Project from Yosemite (http://mvz.berkeley.edu/Grinnell/).

1. Input data format: (SDMpriors_loadEnvPAdata.R loads and describes data)
Species data, "SpeciesList_PresAbs.csv": Master list of 47 species. Each row is a species with columns for priors (tmax and tmin: critical thermal minima and maxima (C), respectively) etc.

2. Observations:  a data frame for each species (e.g., "PresAbs_Bufo alvarius.csv") in the master list with presences and pseudoabsence data, where each presence/absence is a row, and with columns for Lat, Lon, Presence/Absence (1/0), and 8 columns for environmental data:
  trmin= annual daily minimum temperature (C) with thermoregulation, selecting among 0, 50, and 100% shade to get as close to species' Topt as possible
  trmax= annual daily maximum temperature with thermoregulation, selecting among 0, 50, and 100% shade to get as close to species' Topt as possible
tmax0= annual daily maximum temperature (C) with 0 shade,
tmin0= annual daily minimum temperature (C) with 0 shade, 
tmax50= annual daily maximum temperature (C) with 50% shade, 
tmin50= annual daily minimum temperature (C) with 50% shade, 
tmax100= annual daily maximum temperature (C) with 100% shade
tmin100= annual daily minimum temperature (C) with 100% shade,  

3.  Extrapolations, "EnviDat.csv":  a data frame encoding the raster that represents the spatial domain for interpreting results, where each row is the centroid of a raster grid cell, and columns for Grid ID, Lat (y), Lon (x), and then 6 environmental columns (as above without trmin and trmax since they are species specific due to modelling thermoregulation). 
I was thinking to start with only modelling “tmax50”, shared between observation and environmental data and a good fit to physiological data. I can make species specific environmental files if helpful.

Prior:
I believe we discussed trying out the model with a simple polynomial prior function (using CTmin and CTmax as  y-intercepts), but we’ve explored a number of prior functions with variable realism (https://github.com/HuckleyLab/SDMpriors/blob/master/SDMpriors_priors.R)

