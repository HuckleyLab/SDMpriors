
desktop<- "y"

library(ggplot2)
library(reshape)
library(viridis)
library(patchwork)

#ibis.iSDM package
# For Installation directly from github
install.packages("remotes")
remotes::install_github("IIASA/ibis.iSDM")

#https://iiasa.github.io/ibis.iSDM/articles/02_train_simple_model.html
#names(mod)
#"set_priors"

#https://iiasa.github.io/ibis.iSDM/articles/03_integrate_data.html
#Integration with priors
# Probabilistic priors with estimates placed on for example the mean (μ) and standard deviation (σ) or precision in the case of [engine_inla]

#https://iiasa.github.io/ibis.iSDM/articles/08_frequently-asked-questions.html
#How exactly do I add a prior to a model ?

# We have prior information that 'Forest' is important for a species
# In this case and for the INLA engine we define normal prior on the mean and precision
p <- INLAPrior(variable = "Forest",type = "normal",hyper = c(2, 10))
# This is then wrapped in a PriorList
pp <- priors(p)
print( pp )

# And can now added to the model
mod <- distribution(background, limits = zone) |>  
  add_biodiversity_poipo(species_data) |>  
  add_predictors(covariates) |>  
  add_priors(priors = pp)
engine_inlabru()

#====================
# Load the package
# Load the package
library(ibis.iSDM)
library(inlabru)
library(xgboost)
library(terra)
library(uuid)
library(assertthat)

# Don't print out as many messages
options("ibis.setupmessages" = FALSE)

# Background layer
background <- terra::rast(system.file("extdata/europegrid_50km.tif",package = "ibis.iSDM", mustWork = TRUE))
# Load virtual species points
virtual_species <- sf::st_read(system.file("extdata/input_data.gpkg",package = "ibis.iSDM", mustWork = TRUE), "points") 
#> Reading layer `points' from data source 
#>   `/home/runner/work/_temp/Library/ibis.iSDM/extdata/input_data.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 208 features and 5 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 4.109162 ymin: 48.7885 xmax: 24.47594 ymax: 64.69323
#> Geodetic CRS:  WGS 84
# Predictors
predictors <- terra::rast(list.files(system.file("extdata/predictors/", package = "ibis.iSDM", mustWork = TRUE), "*.tif",full.names = TRUE))
# Make use only of a few of them
predictors <- subset(predictors, c("bio01_mean_50km","bio03_mean_50km","bio19_mean_50km",
                                   "CLC3_112_mean_50km","CLC3_132_mean_50km",
                                   "CLC3_211_mean_50km","CLC3_312_mean_50km",
                                   "elevation_mean_50km"))

# First we define a distribution object using the background layer
mod <- distribution(background)

# Then lets add species data to it. 
# This data needs to be in sf format and key information is that
# the model knows where occurrence data is stored (e.g. how many observations per entry) as
# indicated by the field_occurrence field.
mod <- add_biodiversity_poipo(mod, virtual_species,
                              name = "Virtual test species",
                              field_occurrence = "Observed")

# Then lets add predictor information
# Here we are interested in basic transformations (scaling), but derivates (like quadratic)
# for now, but check options
mod <- add_predictors(mod, 
                      env = predictors,
                      transform = "scale", derivates = "none")

# Finally define the engine for the model
# This uses the default data currently backed in the model,
# !Note that any other data might require an adaptation of the default mesh parameters used by the engine!
mod <- engine_inlabru(mod)

# Print out the object to see the information that is now stored within
print(mod)
#> <Biodiversity distribution model>
#> Background extent: 
#>      xmin: -16.064, xmax: 34.95,
#>      ymin: 36.322, ymax: 71.535
#>    projection: +proj=longlat +datum=WGS84 +no_defs
#>  --------- 
#> Biodiversity data:
#>    Point - Presence only <208 records>
#>  --------- 
#>   predictors:     bio01_mean_50km, bio03_mean_50km, bio19_mean_50km, ... (8 predictors)
#>   priors:         <Default>
#>   latent:         None
#>   log:            <Console>
#>   engine:         <INLABRU>

print("Create model")

mod <- distribution(background) |> 
  add_biodiversity_poipo(virtual_species,
                         name = "Virtual test species",
                         field_occurrence = "Observed") |>  
  add_predictors(env = predictors, transform = "scale", derivates = "none") |> 
  engine_inlabru() 
