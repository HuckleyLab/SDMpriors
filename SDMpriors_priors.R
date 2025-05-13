#Tolerance functions

#example parameters
CTmin1=5
CTmax1=35

#PRIORS from virtual species library
library(virtualspecies)

#-----------
#betaFun(x, p1, p2, alpha, gamma)
#x a numeric value or vector. The input environmental variable.
#p1 a numeric value or vector. Lower tolerance bound for the species
#p2 a a numeric value or vector. Upper tolerance bound for the species
#alpha a numeric value or vector. Parameter controlling the shape of the curve 
#gamma a numeric value or vector. Parameter controlling the shape of the curve
#When alpha = gamma, the curve is symmetric. 
#Low values of alpha and gamma result in smooth (< 1) to plateau (< 0.01) curves.
#Higher values result in peak (> 10) curves.
#When alpha < gamma, the curve is skewed to the right. When gamma < alpha, the curve is skewed to the left.

my.betaFun= function(x, CTmin= CTmin1, CTmax= CTmax1, alpha=0.3, gamma=0.3)  betaFun(x, CTmin, CTmax, alpha, gamma)

plot(1:40, my.betaFun(1:40, CTmin1, CTmax1, 0.2, 0.2), type="l") #broad
plot(1:40, my.betaFun(1:40, CTmin1, CTmax1, 0.3, 0.3), type="l")
plot(1:40, my.betaFun(1:40, CTmin1, CTmax1, 0.5,  0.2), type="l") #skewed

#-----------
#custnorm(x, mean, diff, prob)
# x a numeric value or vector. The input environmental variable.
# mean a numeric value or vector. The optimum (mean) of the normal curve
# diff a numeric value or vector. The absolute difference between the mean and extremes.
# prob a numeric value or vector. The percentage of the area under the curve between the chosen extreme values

my.custnorm= function(x, CTmin= CTmin1, CTmax= CTmax1, prob=0.99){  
  diff= (CTmax-CTmin)/2-CTmin
  sd= -diff/qnorm(p = 1 - prob)
  custnorm(x, mean=(CTmax-CTmin)/2, diff=diff, prob=0.95)*sd/0.4
}

plot(1:40, my.custnorm(1:40, CTmin= CTmin1, CTmax= CTmax1, prob=0.99), type="l")

#-------------------
#Polynomial prior

poly.prior <- function(x, CTmin= CTmin1, CTmax= CTmax1, pmax=0.5){
  center= (CTmin+CTmax)/2
  a= pmax/((center-CTmin)*(center-CTmax))
  P1=a*(x-CTmin)*(x-CTmax)
  P1[P1<0]<- 0
  return(P1)
}

plot(1:40, poly.prior(1:40, CTmin= CTmin1, CTmax= CTmax1), type="l")

#-----------
#empirical TPC functions

#Performance Curve Function from Deutsch et al. 2008
TPC= function(T,Topt=33,CTmin=10.45, CTmax=42.62, bound="mean"){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #set negetative to zero
  F[F<0]<-0
  #set to 1 above or below min or max
  if(bound=="max") F[T<=Topt]=1
  if(bound=="min") F[T>=Topt]=1
  
  return(F)
}

TPC.prior= function(x, CTmin= CTmin1, CTmax= CTmax1){ 
  Topt= CTmin+ (CTmax-CTmin)*0.7
  P1= TPC(x, Topt, CTmin, CTmax, bound="mean")
  return(P1)
}
#Assume optimum at 70% of breadth

plot(1:40, TPC.prior(1:40, CTmin= CTmin1, CTmax= CTmax1), type="l")

#-----------

#alternative implementation of normal
norm.prior= function(x, CTmin= CTmin1, CTmax= CTmax1){ 
  sd1= (CTmax-CTmin)/6
  P1= dnorm(x, mean = (CTmax-CTmin)/2, sd = sd1)
  #scale to height 1
  P1= P1*sd1/0.4
  return(P1)
}

plot(1:40, norm.prior(1:40, CTmin= CTmin1, CTmax= CTmax1), type="l")

#----------
#threshold prior
thresh.prior <- function(x, CTmin= CTmin1, CTmax= CTmax1) ifelse(x< CTmin | x> CTmax, 0.2, 0.6)

plot(1:40, thresh.prior(1:40, CTmin= CTmin1, CTmax= CTmax1), type="l")

#----------
#sigmoidal prior

gensigmoid <- function(x, low, high, rate, v, origin) {
  # [Generalized Sigmoid function.](https://en.wikipedia.org/wiki/Generalised_logistic_function)
  return(low + ((high-low)/(1+(rate*exp((x-origin))))^(1/v)))
}

sigmoid.prior<- function(x, CTmin= CTmin1, CTmax= CTmax1, low_p=0.1, high_p=0.5, rate_p=5.5, v_p=2.5) {
  p= rep(NA, length(x))
  center= (CTmin + CTmax)/2
  p[which(x<=center)]= gensigmoid(x[which(x<=center)], low=high_p, high=low_p, rate=rate_p, v=v_p, origin=CTmin)
  p[which(x>center)]= gensigmoid(x[which(x>center)], low=low_p, high=high_p, rate=rate_p, v=v_p, origin=CTmax)
  return(p)
}  
  
plot(1:40, sigmoid.prior(1:40, CTmin= CTmin1, CTmax= CTmax1), type="l")






