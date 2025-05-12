#Gaussian Process Species Distribution Models with physiological priors, following the approach described in Golding & Purse (2016) but without relying on the GRaF package

# Core GP implementation with Laplace approximation
# Physiological prior functions
# Prediction functions
# Model evaluation
# Cross-validation
# Lengthscale optimization
# Visualization

library(ggplot2)

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
  X <- model.matrix(formula, data)[, -1, drop = FALSE]  # Remove intercept
  
  # Set default lengthscales if not provided
  if (is.null(lengthscales)) {
    lengthscales <- rep(1, ncol(X))
    names(lengthscales) <- colnames(X)
  }
  
  # Set default mean function if not provided
  if (is.null(mean_function)) {
    mean_function <- function(X) rep(0, nrow(X))
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

#' Temperature Response Function
#' 
#' @param temp Temperature values
#' @param optimal_temp Optimal temperature for the species
#' @param tolerance Temperature tolerance (width of response curve)
#' @return Probability values between 0 and 1
temp_response <- function(temp, optimal_temp = 20, tolerance = 5) {
  # A Gaussian response curve centered at optimal temperature
  prob <- exp(-((temp - optimal_temp)^2) / (2 * tolerance^2))
  return(prob)
}

#' Precipitation Response Function
#' 
#' @param precip Precipitation values
#' @param min_precip Minimum viable precipitation
#' @param max_precip Optimal maximum precipitation
#' @return Probability values between 0 and 1
precip_response <- function(precip, min_precip = 300, max_precip = 2000) {
  # Sigmoid response with decline after maximum
  prob <- 1 / (1 + exp(-0.01 * (precip - min_precip)))
  # Declining response after maximum
  high_idx <- which(precip > max_precip)
  if (length(high_idx) > 0) {
    decline_factor <- 1 - (precip[high_idx] - max_precip) / max_precip
    decline_factor <- pmax(0, pmin(1, decline_factor))  # Constrain to [0,1]
    prob[high_idx] <- prob[high_idx] * decline_factor
  }
  return(prob)
}

#' Combined Physiological Response Function
#' 
#' @param data Data frame with environmental variables
#' @param temp_col Name of temperature column
#' @param precip_col Name of precipitation column
#' @param temp_opt Optimal temperature
#' @param temp_tol Temperature tolerance
#' @param precip_min Minimum viable precipitation
#' @param precip_max Maximum optimal precipitation
#' @return Logit-transformed probabilities for use as GP prior mean
physiological_prior <- function(data, 
                                temp_col = "temperature", 
                                precip_col = "precipitation",
                                temp_opt = 20, 
                                temp_tol = 5,
                                precip_min = 300, 
                                precip_max = 2000) {
  
  # Extract environmental variables
  temp <- data[[temp_col]]
  precip <- data[[precip_col]]
  
  # Calculate individual responses
  temp_prob <- temp_response(temp, temp_opt, temp_tol)
  precip_prob <- precip_response(precip, precip_min, precip_max)
  
  # Combine responses (multiplicative assumption)
  combined_prob <- temp_prob * precip_prob
  
  # Constrain probabilities to avoid numerical issues
  combined_prob <- pmax(pmin(combined_prob, 0.9999), 0.0001)
  
  # Convert to logit scale
  logit_prob <- log(combined_prob / (1 - combined_prob))
  
  return(logit_prob)
}

#just temperature
temperature_prior <- function(data, 
                                temp_col = "temperature", 
                                temp_opt = 20, 
                                temp_tol = 5) {
  
  # Extract environmental variables
  temp <- data[,temp_col]
  
  # Calculate individual responses
  temp_prob <- temp_response(temp, temp_opt, temp_tol)
  
  # Combine responses (multiplicative assumption)
  combined_prob <- temp_prob
  
  # Constrain probabilities to avoid numerical issues
  combined_prob <- pmax(pmin(combined_prob, 0.9999), 0.0001)
  
  # Convert to logit scale
  logit_prob <- log(combined_prob / (1 - combined_prob))
  
  return(logit_prob)
}

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
    test_data <- data[folds == i, ]
    
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
      eval_metrics <- evaluate_model(predictions, test_data$presence)
      
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

# Example usage
# Define different prior functions to compare
mean_functions <- list(
  "Flat" = function(data) rep(0, nrow(data)),
  "Temperature_Only" = function(data) {
    temp_prob <- temp_response(data$temperature, optimal_temp = 22, tolerance = 6)
    return(log(temp_prob / (1 - temp_prob)))
  },
  "Full_Physiological" = my_prior
)

# Run cross-validation
cv_results <- cross_validate_gp(
  presence ~ temperature + precipitation,
  data = data,
  mean_functions = mean_functions,
  k = 5,
  lengthscales = c(temperature = 5, precipitation = 500)
)

# Print summary
print(cv_results$summary)

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

# Example usage
opt_result <- optimize_lengthscales(
  presence ~ temperature + precipitation,
  data = train_data,
  mean_function = my_prior,
  init_lengthscales = c(temperature = 5, precipitation = 500)
)

print(opt_result$lengthscales)

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

#use july for max
temp= brick("data/microclim/50_shade/TA1cm_soil_50_7.nc")
tmax_50= mean(temp)

#use july for max
temp= brick("data/microclim/100_shade/TA1cm_soil_100_7.nc")
tmax_100= mean(temp)

#---------------
#load presence absence
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")

spec.k<-4

#load presence absence
pa= read.csv(paste("out/presabs/PresAbs_",dat$spec[spec.k],".csv",sep=""))

#----
# Split into training and testing sets
train_idx <- sample(1:n, 0.7 * n)
train_data <- pa[train_idx, ]
test_data <- pa[-train_idx, ]

# Define physiological prior
my_prior <- function(data, temp_col) {
  temperature_prior(
    data, 
    temp_col="tmax50",
    temp_opt = 22, 
    temp_tol = 6
  )
}

# Optimize lengthscales
opt_result <- optimize_lengthscales(
  pres ~ tmax50,
  data = train_data,
  mean_function = my_prior
)

# Fit models with optimized lengthscales
# Fit GP model without physiological prior
gp_flat <- gp_sdm(
  pres ~ tmax50,
  data = train_data,
  lengthscales = opt_result$lengthscales
)

# Fit GP model with physiological prior
gp_prior <- gp_sdm(
  pres ~ tmax50,
  data = train_data,
  mean_function = my_prior,
  lengthscales = opt_result$lengthscales
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

print(results)

#---------
# Visualize model predictions

#extract values
pts= rasterToPoints(tmax50)
#pts= rasterToPoints(tmax50)
colnames(pts)=c("lon","lat","tmax50")
pts= as.data.frame(pts)

pts$pres<- NA

#predictions
pred_flat <- predict(gp_flat, newdata=pts)
pred_prior <- predict(gp_prior, newdata=pts)

#combine predictions
pts= cbind(pts, pred_flat[,1], pred_prior[,1])
colnames(pts)[(ncol(pts)-1):ncol(pts)]=c("occ_pp","occ_np")

#plot
#to long format
pts.l <- melt(pts, id=c("lon","lat","trmax"))
#presence points
pres= subset(pa, pa$pres=="1")
pres$variable="occ_pp"
#replicate pres for plotting
pres2<- pres
pres2$variable="occ_np"
pres.all= rbind(pres, pres2)

#trmax
tr.plot=ggplot(pts, aes(lat, lon, fill= trmax)) + 
  geom_tile()+scale_fill_viridis(na.value = 'grey')

#predictions
occ.plot= ggplot(pts.l, aes(lat, lon)) + 
  geom_tile(aes(fill= value))+scale_fill_viridis(na.value = 'grey') +facet_wrap(~variable, nrow=1)+
  ggtitle(dat$species[spec.k])
#add localities
occ.plot= occ.plot +geom_point(pres.all, mapping=aes(lat, lon, color="red"))

#combine plots
print(occ.plot)
#plot_grid(tr.plot, occ.plot, resp_np, resp_pp)




#&&&&&&&&&&&&&&
# Create prediction grid
temp_seq <- seq(5, 35, length.out = 50)
precip_seq <- seq(0, 3000, length.out = 50)
pred_grid <- expand.grid(temperature = temp_seq, precipitation = precip_seq)

# Make predictions on grid
pred_grid$presence <- NA
pred_grid$flat <- predict(gp_flat, pred_grid)
pred_grid$prior <- predict(gp_prior, pred_grid)
pred_grid$true <- with(pred_grid, 
                       temp_response(temperature, 22, 6) * 
                         precip_response(precipitation, 400, 1800))

# Plot results
p1 <- ggplot(pred_grid, aes(x = temperature, y = precipitation, fill = flat)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "GP Model without Prior", fill = "Probability") +
  theme_minimal()

p2 <- ggplot(pred_grid, aes(x = temperature, y = precipitation, fill = prior)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "GP Model with Physiological Prior", fill = "Probability") +
  theme_minimal()

p3 <- ggplot(pred_grid, aes(x = temperature, y = precipitation, fill = true)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "True Physiological Response", fill = "Probability") +
  theme_minimal()

# Add training data points
p1 <- p1 + geom_point(data = train_data, aes(color = factor(presence)), 
                      size = 2, alpha = 0.7) +
  scale_color_manual(values = c("white", "black"), name = "Presence")

p2 <- p2 + geom_point(data = train_data, aes(color = factor(presence)), 
                      size = 2, alpha = 0.7) +
  scale_color_manual(values = c("white", "black"), name = "Presence")

# Print plots
print(p1)
print(p2)
print(p3)



