#Gaussian Process Species Distribution Models with physiological priors, following the approach described in Golding & Purse (2016) but without relying on the GRaF package

# Core GP implementation with Laplace approximation
# Physiological prior functions
# Prediction functions
# Model evaluation
# Cross-validation
# Lengthscale optimization
# Visualization

#other packages
#Gaussian process regression: https://cran.r-project.org/web/packages/GauPro/vignettes/GauPro.html
#https://github.com/iiasa/ibis.iSDM
#https://link.springer.com/article/10.1007/s13253-023-00595-6, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13897

library(ggplot2)
library(reshape)
library(viridis)
library(patchwork)

#source priors
source("SDMpriors_priors.R")

# Install and load necessary packages
packages <- c("raster", "sp", "dplyr", "mgcv", "ROCR", "kernlab")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#' Gaussian Process Species Distribution Model with Physiological Priors
#' 
#' @param formula Formula specifying the model structure
#' @param data Data frame containing presence/absence and environmental variables
#' @param mean_function Custom mean function incorporating physiological priors (default: NULL)
#' @param kernel Covariance function (default: "rbf")
#' @param lengthscales Vector of lengthscales for each predictor
#' @param sigma Noise parameter
#' @param method Inference method ("Laplace" or "EP")
#' @return A fitted GP model
gp_sdm <- function(formula, data, mean_function = NULL, 
                   kernel = "rbf", lengthscales = NULL, 
                   sigma = 0.1, method = "Laplace") {
  
  # Extract response and predictor variables
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, mf)[, -1, drop = FALSE]  # Remove intercept
  
  # Set default lengthscales if not provided
  if (is.null(lengthscales)) {
    lengthscales <- rep(1, ncol(X))
    names(lengthscales) <- colnames(X)
  }
  
  # Set default mean function if not provided
  if (is.null(mean_function)) {
    mean_function <- function(X) rep(0.5, nrow(X)) #or set to 0, try constant 0.5
  }
  
  # Compute prior mean
  prior_mean <- mean_function(data)
  
  # Define kernel function
  if (kernel == "rbf") {
    kernel_fn <- function(X1, X2 = NULL) {
      if (is.null(X2)) X2 <- X1
      
      K <- matrix(0, nrow = nrow(X1), ncol = nrow(X2))
      
      for (i in 1:nrow(X1)) {
        for (j in 1:nrow(X2)) {
          # Compute squared distance with lengthscale weighting
          dist_sq <- sum(((X1[i,] - X2[j,])^2) / (lengthscales^2))
          K[i, j] <- exp(-0.5 * dist_sq)
        }
      }
      return(K)
    }
  } else {
    stop("Only RBF kernel is currently implemented")
  }
  
  # Compute kernel matrix
  K <- kernel_fn(X)
  
  # Add noise to diagonal for numerical stability
  K_y <- K + diag(sigma^2, nrow(K))
  
  # Fit GP model using specified inference method
  if (method == "Laplace") {
    model <- laplace_inference(y, X, K_y, prior_mean)
  } else if (method == "EP") {
    stop("EP inference not yet implemented")
  } else {
    stop("Method must be either 'Laplace' or 'EP'")
  }
  
  # Return model object
  result <- list(
    formula = formula,
    data = data,
    X = X,
    y = y,
    mean_function = mean_function,
    kernel_fn = kernel_fn,
    lengthscales = lengthscales,
    sigma = sigma,
    K = K,
    K_y = K_y,
    method = method,
    model = model,
    prior_mean = prior_mean
  )
  
  class(result) <- "gp_sdm"
  return(result)
}

#' Laplace Approximation for GP Inference
#' 
#' @param y Binary response variable (0/1)
#' @param X Predictor matrix
#' @param K Kernel matrix
#' @param prior_mean Prior mean vector
#' @return List containing posterior parameters
laplace_inference <- function(y, X, K, prior_mean) {
  n <- length(y)
  f <- prior_mean  # Initialize latent function at prior mean
  
  # Newton's method for finding posterior mode
  max_iter <- 100
  tol <- 1e-6
  
  for (iter in 1:max_iter) {
    # Compute first and second derivatives of log likelihood
    p <- 1 / (1 + exp(-f))
    W <- diag(p * (1 - p))
    
    # Gradient and Hessian of negative log posterior
    grad <- (y - p) - solve(K, f - prior_mean)
    hess <- -W - solve(K)
    
    # Newton update
    delta_f <- -solve(hess, grad)
    f_new <- f + delta_f
    
    # Check convergence
    if (max(abs(f_new - f)) < tol) {
      f <- f_new
      break
    }
    
    f <- f_new
  }
  
  # Posterior covariance (Laplace approximation)
  W_sqrt <- sqrt(W)
  B <- diag(n) + W_sqrt %*% K %*% W_sqrt
  L <- chol(B)
  V <- solve(L, W_sqrt %*% K)
  Sigma <- K - t(V) %*% V
  
  return(list(
    f_map = f,
    Sigma = Sigma,
    W = W
  ))
}

#' Predict Method for GP SDM
#' 
#' @param object Fitted GP SDM model
#' @param newdata New data for prediction
#' @param type Type of prediction ("link" or "response")
#' @return Vector of predictions
predict.gp_sdm <- function(object, newdata, type = "response") {
  # Extract model components
  X <- object$X
  y <- object$y
  K <- object$K
  f_map <- object$model$f_map
  Sigma <- object$model$Sigma
  
  # Prepare new data
  mf <- model.frame(object$formula, newdata, na.action = na.pass)
  X_new <- model.matrix(object$formula, newdata)[, -1, drop = FALSE]
  
  # Compute prior mean for new data
  prior_mean_new <- object$mean_function(newdata)
  
  # Compute cross-covariance
  K_s <- object$kernel_fn(X_new, X)
  
  # Predictive mean
  f_pred <- prior_mean_new + K_s %*% solve(object$K_y, f_map - object$prior_mean)
  
  # Predictive variance (optional, can be computationally expensive)
  # K_ss <- object$kernel_fn(X_new)
  # var_pred <- diag(K_ss - K_s %*% solve(object$K_y, t(K_s)))
  
  # Return predictions
  if (type == "link") {
    return(f_pred)
  } else if (type == "response") {
    return(1 / (1 + exp(-f_pred)))
  } else {
    stop("Type must be either 'link' or 'response'")
  }
}

#--------------------------------
#Physiological Prior Functions
#Now let's implement some physiological prior functions that can be used as mean functions in our GP model:
#see SDMpriors_test.R for combined priors

#' Temperature Response Function
#' 
#' @param data Data frame with environmental variables
#' @param temp_col Name of temperature column
#' @param ctmin critical thermal minima
#' @param ctmax critical thermal maxima
#' @param type prior type in "beta", "gaussian", "poly", "threshold", "sigmoid", "tpc"
#' @return Logit-transformed probabilities for use as GP prior mean
temperature_prior <- function(data, 
                              temp_col = "temperature", 
                              ctmin = 5, 
                              ctmax = 35,
                              type= "sigmoid") {
  
  # Extract environmental variables
  temp <- data[,temp_col]
  
  # Calculate individual responses
  if(type=="beta") temp_prob<- my.betaFun(temp, CTmin=ctmin, CTmax=ctmax, 0.2, 0.2)
  if(type=="gaussian") temp_prob<- my.custnorm(temp, CTmin=ctmin, CTmax=ctmax, prob=0.99)
  if(type=="poly") temp_prob<- poly.prior(temp, CTmin=ctmin, CTmax=ctmax, pmax=1)
  if(type=="threshold") temp_prob<- thresh.prior(temp, CTmin=ctmin, CTmax=ctmax)
  if(type=="sigmoid") temp_prob<- sigmoid.prior(temp, CTmin=ctmin, CTmax=ctmax, high_p=1)
  if(type=="sigmoid") temp_prob<- sigmoid.prior(temp, CTmin=ctmin, CTmax=ctmax, high_p=1)
  if(type=="tpc") temp_prob<- TPC.prior(temp, CTmin=ctmin, CTmax=ctmax) 
  
  # Convert to logit scale
  logit_prob <- log(temp_prob / (1 - temp_prob))
  
  return(logit_prob)
}

#-----------------------------------------------------------
# Evaluate model performance
evaluate_model <- function(predictions, observations) {
  # Calculate AUC
  pred_obj <- prediction(predictions, observations)
  auc_obj <- performance(pred_obj, "auc")
  auc <- auc_obj@y.values[[1]]
  
  # Calculate RMSE
  rmse <- sqrt(mean((predictions - observations)^2))
  
  return(list(AUC = auc, RMSE = rmse))
}

#Spatial Prediction with Raster Data
#' Predict to Raster
#' 
#' @param model Fitted GP SDM model
#' @param raster_stack RasterStack of environmental variables
#' @return Raster of predictions
predict_to_raster <- function(model, raster_stack) {
  # Extract variable names from model formula
  var_names <- all.vars(model$formula)[-1]  # Remove response variable
  
  # Check if all variables are in the raster stack
  if (!all(var_names %in% names(raster_stack))) {
    stop("Not all model variables found in raster stack")
  }
  
  # Extract raster data as a data frame
  raster_df <- as.data.frame(raster_stack, xy = TRUE)
  
  # Remove rows with NA
  raster_df_complete <- na.omit(raster_df)
  raster_df_complete$pres<-NA
  
  # Make predictions
  predictions <- predict(model, raster_df_complete, type = "response")
  
  # Add predictions to data frame
  raster_df_complete$prediction <- predictions
  
  # Convert back to raster
  pred_raster <- rasterFromXYZ(raster_df_complete[, c("x", "y", "prediction")])
  
  # Match extent and projection with input raster
  extent(pred_raster) <- extent(raster_stack)
  projection(pred_raster) <- projection(raster_stack)
  
  return(pred_raster)
}

# Example usage with raster data
# Assuming we have environmental rasters:
# env_rasters <- stack(temp_raster, precip_raster)
# pred_raster <- predict_to_raster(gp_prior, env_rasters)
# plot(pred_raster)

#--------------------
#Cross-Validation and Model Selection
#To perform k-fold cross-validation for model selection:

#' K-fold Cross-Validation for GP SDM
#' 
#' @param formula Model formula
#' @param data Full dataset
#' @param mean_functions List of mean functions to compare
#' @param k Number of folds
#' @param lengthscales Lengthscales for the GP kernel
#' @return Data frame of cross-validation results
cross_validate_gp <- function(formula, data, mean_functions, k = 5, 
                              lengthscales = NULL) {
  
  # Create folds
  set.seed(42)
  folds <- sample(1:k, nrow(data), replace = TRUE)
  
  # Initialize results
  results <- data.frame()
  
  # Loop through each fold
  for (i in 1:k) {
    # Split data
    train_data <- data[folds != i, ]
    test_data <- data[which(folds == i), ]
    
    # Loop through mean functions
    for (name in names(mean_functions)) {
      mean_fn <- mean_functions[[name]]
      
      # Fit model
      model <- gp_sdm(
        formula = formula,
        data = train_data,
        mean_function = mean_fn,
        lengthscales = lengthscales
      )
      
      # Make predictions
      predictions <- predict(model, test_data)
      
      # Evaluate
      eval_metrics <- evaluate_model(predictions, test_data$pres)
      
      # Store results
      fold_results <- data.frame(
        Fold = i,
        Model = name,
        AUC = eval_metrics$AUC,
        RMSE = eval_metrics$RMSE
      )
      
      results <- rbind(results, fold_results)
    }
  }
  
  # Summarize results
  summary <- results %>%
    group_by(Model) %>%
    summarize(
      Mean_AUC = mean(AUC),
      SD_AUC = sd(AUC),
      Mean_RMSE = mean(RMSE),
      SD_RMSE = sd(RMSE)
    )
  
  return(list(fold_results = results, summary = summary))
}

#---------------
#length scale optimization

#' Optimize GP Lengthscales
#' 
#' @param formula Model formula
#' @param data Training data
#' @param mean_function Mean function
#' @param init_lengthscales Initial lengthscales
#' @return Optimized lengthscales
optimize_lengthscales <- function(formula, data, mean_function = NULL, 
                                  init_lengthscales = NULL) {
  
  # Extract variable names
  var_names <- all.vars(formula)[-1]
  
  # Set initial lengthscales if not provided
  if (is.null(init_lengthscales)) {
    init_lengthscales <- rep(1, length(var_names))
    names(init_lengthscales) <- var_names
  }
  
  # Define objective function (negative log marginal likelihood)
  objective <- function(log_ls) {
    # Convert from log scale
    ls <- exp(log_ls)
    names(ls) <- var_names
    
    # Fit model with current lengthscales
    model <- try(gp_sdm(
      formula = formula,
      data = data,
      mean_function = mean_function,
      lengthscales = ls
    ), silent = TRUE)
    
    # Return high value if model fails
    if (inherits(model, "try-error")) {
      return(1e10)
    }
    
    # Calculate negative log marginal likelihood (approximation)
    f <- model$model$f_map
    K_inv <- solve(model$K_y)
    W <- model$model$W
    
    # Laplace approximation to log marginal likelihood
    p <- 1 / (1 + exp(-f))
    log_lik <- sum(data$presence * log(p) + (1 - data$presence) * log(1 - p))
    log_prior <- -0.5 * t(f - model$prior_mean) %*% K_inv %*% (f - model$prior_mean)
    log_det <- -0.5 * determinant(diag(nrow(data)) + W %*% model$K, logarithm = TRUE)$modulus
    
    lml <- as.numeric(log_lik + log_prior + log_det)
    
    return(-lml)  # Return negative for minimization
  }
  
  # Optimize lengthscales
  opt_result <- optim(
    par = log(init_lengthscales),
    fn = objective,
    method = "L-BFGS-B",
    control = list(maxit = 100)
  )
  
  # Convert optimized lengthscales back from log scale
  opt_lengthscales <- exp(opt_result$par)
  names(opt_lengthscales) <- var_names
  
  return(list(
    lengthscales = opt_lengthscales,
    convergence = opt_result$convergence,
    value = -opt_result$value  # Return positive log marginal likelihood
  ))
}

#===========================
#Example

# load physiological priors from Sunday database
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")

dat= read.csv("out/presabs/SpeciesList_PresAbs.csv")

#Load envi data
#DATA CAN BE DOWNLOADED TO SHARED DRIVE FROM HERE: https://figshare.com/collections/microclim_Global_estimates_of_hourly_microclimate_based_on_long_term_monthly_climate_averages/878253
#USED 1cm Air temperature

if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/")

#load all clim data
#use july for max
temp= brick("data/microclim/0_shade/TA1cm_soil_0_7.nc")
tmax_0= mean(temp) #or max
names(tmax_0)<-"tmax0"

#use july for max
temp= brick("data/microclim/50_shade/TA1cm_soil_50_7.nc")
tmax_50= mean(temp)
names(tmax_50)<-"tmax50"

#use july for max
temp= brick("data/microclim/100_shade/TA1cm_soil_100_7.nc")
tmax_100= mean(temp)
names(tmax_100)<-"tmax100"

#---------------
#set up data storage
#models<- array(NA, dim=c(nrow(dat), 3, 2), dimnames=list(dat$species, c("model","AUC","RMSE"),c("flat","prior") ))
models<- matrix(NA, nrow=nrow(dat), ncol=6)
rownames(models)<- dat$species
models<- as.data.frame(models)
colnames(models)<- c("mod1", "AUC1", "RMSE1", "mod2", "AUC2", "RMSE2")

#load presence absence
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")

for(spec.k in 1:nrow(dat)){

  print(spec.k)
#load presence absence
pa= read.csv(paste("out/presabs/PresAbs_",dat$spec[spec.k],".csv",sep=""))

#----
#plot localities and temperature
#crop to observed range 
ext = extent(rbind(range(pa$lon), range(pa$lat))) # define the extent
# extent
ext[1]= ext[1]-10; ext[2]= ext[2]+10; ext[3]=ext[3]-10; ext[4]=ext[4]+10
#crop
tmax0=  crop(tmax_0, ext)
tmax50=  crop(tmax_50, ext)
tmax100=  crop(tmax_100, ext)

#------
#model thermoregulation
#Set up prior
CTmin1= dat$tmin[spec.k]
CTmax1= dat$tmax[spec.k]
#approximate Topt, but fix based on data
Topt= CTmin1+ (CTmax1-CTmin1)*0.7

#-----
# sun to shade
# thermoregulation scenario

#max
tmax0.dif= abs(tmax0 - Topt) 
tmax50.dif= abs(tmax50 - Topt) 
tmax100.dif= abs(tmax100 - Topt) 
tmax.dif= stack(tmax0.dif, tmax50.dif, tmax100.dif)
tr.ind= which.min(tmax.dif)

tr<- tmax0
tr[]<-NA
tr[tr.ind==1]<- tmax0[tr.ind==1]
tr[tr.ind==2]<- tmax50[tr.ind==2]
tr[tr.ind==3]<- tmax100[tr.ind==3]
trmax=tr
names(trmax)<-"trmax"
#----
# Split into training and testing sets
train_idx <- sample(1:nrow(pa), 0.7 * nrow(pa))
train_data <- pa[train_idx, ]
test_data <- pa[-train_idx, ]

#------------------------------
# Define physiological prior
my_prior <- function(data, temp_col, ctmin=dat$tmin[spec.k], ctmax=dat$tmax[spec.k]) {
  temperature_prior(
    data, 
    temp_col="trmax",
    ctmin = ctmin, 
    ctmax = ctmax,
    type= "tpc"
  )
}

#-----------------
#make predictions
pred= temperature_prior(pa, temp_col="trmax", ctmin=dat$tmin[spec.k], ctmax=dat$tmax[spec.k], type="tpc")
#"beta", "gaussian", "poly", "threshold", "sigmoid", "tpc"

#convert back from logit scale
pa$pred= 1 / (1 + exp(-pred))

#plot prior
plot.prior<- ggplot() +
  geom_point(aes(trmax, pres), data=pa)+
  #add prior
  geom_line(aes(trmax, pred), data=pa)+
  xlim(0, 42)

#--------------------------------
# Optimize lengthscales

 # opt_result <- optimize_lengthscales(
 #   pres ~ trmax,
 #   data = train_data,
 #   mean_function = my_prior
 # )

# Fit models with optimized lengthscales
# Fit GP model without physiological prior
gp_flat <- gp_sdm(
  pres ~ trmax,
  data = train_data,
  lengthscales = 5 #opt_result$lengthscales #5
)

# Fit GP model with physiological prior
gp_prior <- gp_sdm(
  pres ~ trmax,
  data = train_data,
  mean_function = my_prior,
  lengthscales = 5 #opt_result$lengthscales
)

# Make predictions
pred_flat <- predict(gp_flat, test_data)
pred_prior <- predict(gp_prior, test_data)

# Compare models
eval_flat <- evaluate_model(pred_flat, test_data$pres)
eval_prior <- evaluate_model(pred_prior, test_data$pres)

results <- data.frame(
  Model = c("GP-Flat", "GP-Prior"),
  AUC = c(eval_flat$AUC, eval_prior$AUC),
  RMSE = c(eval_flat$RMSE, eval_prior$RMSE)
)

models[spec.k,1:3]<- c(results[1,])
models[spec.k,4:6]<- c(results[2,])

# # Cross validate
# mean_functions <- list(
#   "Flat" = function(data) rep(0, nrow(data)),
#   # "Temperature_Only" = function(data) {
#   #   temp_prob <- temp_response(data$trmax, optimal_temp = 22, tolerance = 6)
#   #   return(log(temp_prob / (1 - temp_prob)))
#   # },
#   "Full_Physiological" = my_prior
# )
# 
# # Run cross-validation
# cv_results <- cross_validate_gp(
#   pres ~ trmax,
#   data = train_data,
#   mean_functions = mean_functions,
#   k = 5,
#   lengthscales = opt_result$lengthscales
# )
# 
# # Print summary
# print(cv_results$summary)

#---------------------------
#point based approach
#extract values
pts= rasterToPoints(trmax)
colnames(pts)=c("lon","lat","trmax")
pts= as.data.frame(pts)
pts$pres<- 0

#predictions
pred_flat <- predict(gp_flat, newdata=pts, type="response")
pred_prior <- predict(gp_prior, newdata=pts, type="response")

#combine predictions
pts= cbind(pts, pred_flat[,1], pred_prior[,1])
colnames(pts)[(ncol(pts)-1):ncol(pts)]=c("occ_pp","occ_np")

#plot
#to long format
pts.l <- melt(pts[,which(names(pts)!="pres")], id=c("lon","lat","trmax"))
#presence points
pres= subset(pa, pa$pres=="1")
pres$variable="occ_pp"
#replicate pres for plotting
pres2<- pres
pres2$variable="occ_np"
pres.all= rbind(pres, pres2)

#trmax
tr.plot=ggplot(pts, aes(lon, lat, fill= trmax)) + 
  geom_tile()+scale_fill_viridis(na.value = 'grey')

#predictions
occ.plot= ggplot(pts.l, aes(lon, lat)) + 
  geom_tile(aes(fill= value))+scale_fill_viridis(na.value = 'grey') +facet_wrap(~variable, nrow=1)+
  ggtitle(dat$species[spec.k])
#add localities
occ.plot= occ.plot +geom_point(pres.all, mapping=aes(lon, lat, color="red"))

#combine and write out
design <- "ABBB"

#save figure 
pdfname <- paste("./figures/",dat$species[spec.k], ".pdf", sep="")
pdf(pdfname,height = 6, width = 10)
print(plot.prior + occ.plot + plot_layout(design = design))
dev.off()

} #end loop species

#-----------
#examine results 
plot(models$AUC1, models$AUC2)
abline(a=0, b=1)

plot(models$RMSE1, models$RMSE2)
abline(a=0, b=1)

#================
#alternative visualization of model predictions
# Visualize model predictions
#raster approach

# Add a dummy 'pres' layer (with all NA values)
dummy_pres <- trmax
values(dummy_pres) <- 0
names(dummy_pres) <- "pres"

# Stack them together
pred_stack <- stack(dummy_pres, trmax)

# Now predict
prediction_raster <- raster::predict(pred_stack, gp_flat)
plot(prediction_raster)

prediction_raster_prior <- raster::predict(pred_stack, gp_prior)
plot(prediction_raster_prior)

#add localities
#ppt= subset(pa, pa$pres=="0")
#points(ppt$lon, ppt$lat)
ppt= subset(pa, pa$pres=="1")
points(ppt$lon, ppt$lat, pch=16)

#-----------------

