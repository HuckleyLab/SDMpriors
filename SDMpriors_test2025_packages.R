
desktop<- "y"

library(ggplot2)
library(reshape)
library(viridis)
library(patchwork)

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

#ibis.iSDM package
# For Installation directly from github
#install.packages("remotes")
#remotes::install_github("IIASA/ibis.iSDM")

#https://iiasa.github.io/ibis.iSDM/articles/02_train_simple_model.html
#names(mod)
#"set_priors"

#https://iiasa.github.io/ibis.iSDM/articles/03_integrate_data.html
#Integration with priors
# Probabilistic priors with estimates placed on for example the mean (μ) and standard deviation (σ) or precision in the case of [engine_inla]

#====================
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

#Create model
mod <- distribution(background) |> 
  add_biodiversity_poipo(virtual_species,
                         name = "Virtual test species",
                         field_occurrence = "Observed") |>  
  add_predictors(env = predictors, transform = "scale", derivates = "none") |>
  engine_inlabru() 

# Make visualization of the contained biodiversity data
plot(mod$biodiversity)

# Other options to explore
names(mod)

# Finally train
fit <- train(mod,
             runname =  "Test INLA run",
             aggregate_observations = FALSE, # Don't aggregate point counts per grid cell
             verbose = TRUE # Don't be chatty
)

# Plot the mean of the posterior predictions
plot(fit, "mean")

# Print out some summary statistics
summary(fit)

# Show the default effect plot from inlabru
effects(fit)

# To calculate a partial effect for a given variable
o <- partial(fit, x.var = "CLC3_312_mean_50km", plot = TRUE)

# The object o contains the data underlying this figure

# Similarly the partial effect can be visualized spatially as 'spartial'
s <- spartial(fit, x.var = "CLC3_312_mean_50km")
plot(s[[1]], col = rainbow(10), main = "Marginal effect of forest on the relative reporting rate")

# Calculate a threshold based on a 50% percentile criterion
fit <- threshold(fit, method = "percentile", value = 0.5)

# Notice that this is now indicated in the fit object
print(fit)

# There is also a convenient plotting function
fit$plot_threshold()

# It is also possible to use truncated thresholds, which removes non-suitable areas
# while retaining those that are suitable. These are then normalized to a range of [0-1]
fit <- threshold(fit, method = "percentile", value = 0.5, format = "normalize")
fit$plot_threshold()

#--------------------------
#xgboost algorithm

# We are going to fit two separate Poisson Process Models (PPMs) on presence-only data.

# Load the predictors again
predictors <- terra::rast(list.files(system.file("extdata/predictors/", package = "ibis.iSDM"), "*.tif",full.names = TRUE))
predictors <- subset(predictors, c("bio01_mean_50km","bio03_mean_50km","bio19_mean_50km",
                                   "CLC3_112_mean_50km","CLC3_132_mean_50km",
                                   "CLC3_211_mean_50km","CLC3_312_mean_50km",
                                   "elevation_mean_50km",
                                   "koeppen_50km"))
# One of them (Köppen) is a factor, we will now convert this to a true factor variable
predictors$koeppen_50km <- terra::as.factor(predictors$koeppen_50km)

# Create a distribution modelling pipeline
x <- distribution(background) |> 
  add_biodiversity_poipo(virtual_species, field_occurrence = 'Observed', name = 'Virtual points') |>
  add_predictors(predictors, transform = 'scale', derivates = "none") |>
  engine_xgboost(iter = 8000)

# Now train 2 models, one without and one with a spatial latent effect
mod_null <- train(x, runname = 'Normal PPM projection', only_linear = TRUE, verbose = FALSE)
# And with an added constrain
# Calculated as nearest neighbour distance (NND) between all input points
mod_dist <- train(x |> add_latent_spatial(method = "nnd"),
                  runname = 'PPM with NND constrain', only_linear = TRUE, verbose = FALSE)

# Compare both
plot(background, main = "Biodiversity data"); plot(virtual_species['Observed'], add = TRUE)

plot(mod_null)
plot(mod_dist)

# Create again a distribution object, but this time with limits (use the Köppen-geiger layer from above)
# The zones layer must be a factor layer (e.g. is.factor(layer) )

# Zone layers can be supplied directly to distribution(background, limits = zones)
# or through an extrapolation control as shown below.
x <- distribution(background) |> 
  add_biodiversity_poipo(virtual_species, field_occurrence = 'Observed', name = 'Virtual points') |>
  add_predictors(predictors, transform = 'scale', derivates = "none") |>
  # Since we are adding the koeppen layer as zonal layer, we disgard it from the predictors
  rm_predictors("koeppen_50km") |> 
  add_limits_extrapolation(layer = predictors$koeppen_50km, method = "zones") |> 
  engine_xgboost(iter = 3000, learning_rate = 0.01)

# Spatially limited prediction
mod_limited <- train(x, runname = 'Limited prediction background', only_linear = TRUE, verbose = FALSE)

# Compare the output
plot(mod_limited)

#==========================================
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


