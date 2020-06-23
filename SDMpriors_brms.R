#load libraries
library(ggplot2)
library(brms)
library(mgcv)
library(invgamma)
library(viridis)
library(raster)
library(reshape2)
library(cowplot)
library(ROCR)
#library(PRoC)

#RESOURCES 
#paper: 
#https://www.nature.com/articles/s41598-018-38416-3

#brms info:
#https://discourse.mc-stan.org/t/help-understanding-and-setting-informative-priors-in-brms/9574
#https://www.rensvandeschoot.com/tutorials/brms-priors/
#https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
#https://www.rensvandeschoot.com/tutorials/brms-started/
#https://discourse.mc-stan.org/t/space-time-models-in-stan-brms/4735

#--------------------------------
# load physiological priors from Sunday database
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/SDMpriors/out/presabs/")
dat= read.csv("SpeciesList_PresAbs.csv")

#DATA CAN BE DOWNLOADED TO SHARED DRIVE FOR HERE: https://figshare.com/collections/microclim_Global_estimates_of_hourly_microclimate_based_on_long_term_monthly_climate_averages/878253
#DATA WAS TOO BIG FOR ME TO TRANSFER ON MY LAPTOP
#USED 1cm Air temperature

#load all clim data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/data/microclim/0_shade/")
#use july for max
temp= brick("TA1cm_soil_0_7.nc")
tmax_0= mean(temp) #or max

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/data/microclim/50_shade/")
#use july for max
temp= brick("TA1cm_soil_50_7.nc")
tmax_50= mean(temp)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/data/microclim/100_shade/")
#use july for max
temp= brick("TA1cm_soil_100_7.nc")
tmax_100= mean(temp)

#-----
#ESTIMATE TEMPERATURE PRIORS
#set up for inverse gamma, normal, log normal

#cummulative density function
#assign 5% above and below CTmin and max

#Calculate parameters for inverse beta
#https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html
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

#==================================================
#BUILD MODELS ACROSS SPECIES

#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/SDMpriors/out/")
#pdf("priors_brms.pdf", height = 10, width = 12)

#par(mfrow=c(5,4), cex=1.2, mar=c(3, 3, 1.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

#loop species
#for(spec.k in 1:nrow(dat)){
spec.k=4

#load presence absence
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/SDMpriors/out/presabs/")
pa= read.csv(paste("PresAbs_",dat$spec[spec.k],".csv",sep=""))

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

#------
#Plot 
#plot(trmax, main=dat$spec[spec.k])
#points(pa$lon, pa$lat, pch=20, cex=0.5, col="darkgreen")

#extract values
pts= rasterToPoints(trmax)
#pts= rasterToPoints(tmax50)
colnames(pts)=c("lon","lat","trmax")
pts= as.data.frame(pts)

#------
#GP model
#update priors

#get prior
prior<- get_prior(pres ~ gp(trmax), data=pa, 
                  family=bernoulli(link = "logit"))

#intercept
prior$prior[1]="student_t(3, 0, 10)"
#Half student_t?
#need to normalize data?

#sd of gp
prior$prior[4]="student_t(3, 0, 1)"

#length scale
#params= param_irvgamma(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
params= param_normal(dat[spec.k,"tmin"],dat[spec.k,"tmax"])
#params= param_lnormal(dat[spec.k,"tmin"],dat[spec.k,"tmax"])

#prior.in= paste("inv_gamma(",params$shape,",",params$rate,")",sep="")
prior.in= paste("normal(",params$mean,",",params$sd,")",sep="")
#prior.in= paste("lognormal(",params$mean,",",params$sd,")",sep="")
prior$prior[3] <- prior.in
#Kotta uses log gaussian

## verify that the priors indeed found their way into Stan's model code
#make_stancode(pres ~ gp(trmax), data=pa, 
#              family=bernoulli(link = "logit"),
#              prior = prior)

# fit a simple GP model #gp()
fit_pp <- brm(pres ~ gp(trmax, k=5, c=1.25), data=pa, 
            family=bernoulli(link = "logit"),
            prior=prior,
            chains = 2)
#summary(fit_bp)

#---
#fit no physiological prior
#get prior
prior<- get_prior(pres ~ gp(trmax), data=pa, 
                  family=bernoulli(link = "logit"))

#intercept
prior$prior[1]="student_t(3, 0, 10)"
#Half student_t?
#need to normalize data?

#sd of gp
prior$prior[4]="student_t(3, 0, 1)"

fit_np <- brm(pres ~ gp(trmax, k=5, c=1.25), data=pa, 
              family=bernoulli(link = "logit"),
              #prior=prior,
              chains = 2)

#---------------
#PLOTS

#plot GP with phys prior
me1 <- conditional_effects(fit_pp, nsamples = 200, spaghetti = TRUE)
resp_pp= plot(me1, ask = FALSE, points = TRUE, ylim=c(0,1))

#plot GP with no phys prior
me1 <- conditional_effects(fit_np, nsamples = 200, spaghetti = TRUE)
resp_np= plot(me1, ask = FALSE, points = TRUE, ylim=c(0,1))

#predictions
pred_pp= predict(fit_pp, newdata=pts)
pred_np= predict(fit_np, newdata=pts)
#combine
pts= cbind(pts, pred_pp[,1], pred_np[,1])
colnames(pts)[(ncol(pts)-1):ncol(pts)]=c("occ_pp","occ_np")

#plot
#to long format
pts.l <- melt(pts, id=c("lon","lat","trmax"))
#presence points
pres= subset(pa, pa$pres=="1")
pres$variable="occ_pp"

#trmax
tr.plot=ggplot(pts, aes(lat, lon, fill= trmax)) + 
  geom_tile()+scale_fill_viridis(na.value = 'grey')

#predictions
occ.plot= ggplot(pts.l, aes(lat, lon)) + 
  geom_tile(aes(fill= value))+scale_fill_viridis(na.value = 'grey') +facet_wrap(~variable, nrow=1)
#add localities
occ.plot= occ.plot +geom_point(pres, mapping=aes(lat, lon, color="red"))

#combine plots
occ.plot
#plot_grid(tr.plot, occ.plot, resp_np, resp_pp)

#} #end loop species
#dev.off()

#---------------------------
#ROC assessments
#statistics
pred_pp= predict(fit_pp)
pred_np= predict(fit_np)
#make prediction object
pred_pp= prediction(pred_pp[,1], pa$pres)
pred_np= prediction(pred_np[,1], pa$pres)
#performance metrics
perf_pp <- performance(pred_pp, measure = "tpr", x.measure = "fpr")
auc_pp <- performance(pred_pp, measure = "auc")
auc_pp <- auc_pp@y.values[[1]]

perf_np <- performance(pred_np, measure = "tpr", x.measure = "fpr")
auc_np <- performance(pred_np, measure = "auc")
auc_np <- auc_np@y.values[[1]]

roc.data_pp <- data.frame(fpr=unlist(perf_pp@x.values),
                          tpr=unlist(perf_pp@y.values),
                          model="PP")
roc.data_np <- data.frame(fpr=unlist(perf_np@x.values),
                          tpr=unlist(perf_np@y.values),
                          model="NP")

#PLOT ROC
plot(roc.data_pp$fpr,roc.data_pp$tpr,type="l",col="red",ylab="TPR",xlab="FPR",main="ROC for GP vs MaxEnt",lwd=3.5)
lines(roc.data_np$fpr,roc.data_np$tpr,type="l",col="green",lwd=3.5)
legend(0.6,0.4, # places a legend at the appropriate place 
       c("PP","NP"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),col=c("red","green"))

