#load libraries
library(ggplot2)
library(brms)
library(mgcv)

#RESOURCES 
#paper: 
#https://www.nature.com/articles/s41598-018-38416-3
#https://discourse.mc-stan.org/t/help-understanding-and-setting-informative-priors-in-brms/9574
#https://www.rensvandeschoot.com/tutorials/brms-priors/
#https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
#https://www.rensvandeschoot.com/tutorials/brms-started/
#https://discourse.mc-stan.org/t/space-time-models-in-stan-brms/4735

#--------------------------------
# load physiological priors from Sunday database
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/SDMpriors/out/presabs/")
dat= read.csv("SpeciesList_PresAbs.csv")

#load species data
spec.k=4

#load presence absence
pa= read.csv(paste("PresAbs_",dat$spec[spec.k],".csv",sep=""))

#----------------------------
#plot and fit environmental response fucnction

ggplot(pa, aes(x=trmax, y=pres)) +geom_point() +theme_bw()+geom_smooth()

gam1<-gam(pres ~ trmax, data=pa, family="binomial")

# fit a simple GP model #gp()
fit1 <- brm(pres ~ gp(trmax), data=pa, 
            family=bernoulli(link = "logit"), 
            #prior=set_prior("normal(0,10)", class = "b"),
            chains = 2)

prior<- get_prior(pres ~ gp(trmax), data=pa, 
    family=bernoulli(link = "logit"))

#student_t(df = 3, location = 0, scale = 1)
#inv_gamma(shape, rate)
#https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html

#-----
#Calcualte parameters for inverse beta
library(invgamma)
s <- seq(0, 50, 1)
plot(s, dinvgamma(s, 7, 10), type = 'l')

#cummulative density function
#assign 1% above and below CTmin and max

error_invgamma=function(comb, CTmin, CTmax){
  (abs(pinvgamma(CTmin, comb[1], comb[2])-0.05) + abs(pinvgamma(CTmax, comb[1], comb[2])-0.95))/2
}

error_normal=function(comb, CTmin, CTmax){
  (abs(pnorm(CTmin, comb[1], comb[2])-0.05) + abs(pnorm(CTmax, comb[1], comb[2])-0.95))/2
}

error_lnormal=function(comb, CTmin, CTmax){
  (abs(plnorm(CTmin, comb[1], comb[2])-0.05) + abs(plnorm(CTmax, comb[1], comb[2])-0.95))/2
}

param_irvgamma=function(CTmin, CTmax){
combs=expand.grid(shape=seq(1,50,1),rate=seq(1,50,1))
error= apply(combs, MARGIN=1,error_invgamma, CTmin=CTmin, CTmax=CTmax)
return(combs[which.min(error),])
}

param_normal=function(CTmin, CTmax){
  combs=expand.grid(mean=seq(1,50,1),sd=seq(1,50,1))
  error= apply(combs, MARGIN=1,error_normal, CTmin=CTmin, CTmax=CTmax)
  return(combs[which.min(error),])
}

param_lnormal=function(CTmin, CTmax){
  combs=expand.grid(mean=seq(1,50,1),sd=seq(1,50,1))
  error= apply(combs, MARGIN=1,error_lnormal, CTmin=CTmin, CTmax=CTmax)
  return(combs[which.min(error),])
}

params= param_irvgamma(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
s <- seq(0, 50, 1)
plot(s, dinvgamma(s, params$shape, params$rate), type = 'l')

params= param_normal(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
s <- seq(0, 50, 1)
plot(s, dnorm(s, params$mean, params$sd), type = 'l')

params= param_lnormal(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
s <- seq(0, 5, 0.1)
plot(s, dlnorm(s, params$mean, params$sd), type = 'l')

#------
#update prior

#params= param_irvgamma(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
params= param_normal(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
#params= param_lnormal(dat[spec.k,"tmin"],dat[spec.k,"tmax"])

#intercept
prior$prior[1]="student_t(3, 0, 10)"
#Half student_t?
#need to normalize data?

#sd of gp
prior$prior[4]="student_t(3, 0, 1)"

#length scale
#prior.in= paste("inv_gamma(",params$shape,",",params$rate,")",sep="")
prior.in= paste("normal(",params$mean,",",params$sd,")",sep="")
#prior.in= paste("lognormal(",params$mean,",",params$sd,")",sep="")
prior$prior[3] <- prior.in
#Kotta uses log gaussian

## verify that the priors indeed found their way into Stan's model code
make_stancode(pres ~ gp(trmax), data=pa, 
              family=bernoulli(link = "logit"),
              prior = prior)

# fit a simple GP model #gp()
fit1 <- brm(pres ~ gp(trmax, k=8, c=1.25), data=pa, 
            family=bernoulli(link = "logit"),
            prior=prior,
            chains = 2)

summary(fit1)
me1 <- conditional_effects(fit1, nsamples = 200, spaghetti = TRUE)
plot(me1, ask = FALSE, points = TRUE)

#predict
predict(fit1)
