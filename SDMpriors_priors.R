library(virtualspecies)

#tolerance functions
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

plot(1:60, betaFun(1:60, CTmin1, CTmax1, 0.2, 0.2)) #broad
plot(1:60, betaFun(1:60, CTmin1, CTmax1, 0.3, 0.3))
plot(1:60, betaFun(1:60, CTmin1, CTmax1, 0.5,  0.2)) #skewed

m3 <- graf(y1, as.data.frame(train[,2]), prior = my.betaFun, l=100) #opt.l = TRUE ## adjust lengthscale l = 100,
plot(m3, prior=TRUE)

#-----------
#custnorm(x, mean, diff, prob)
# x a numeric value or vector. The input environmental variable.
# mean a numeric value or vector. The optimum (mean) of the normal curve
# diff a numeric value or vector. The absolute difference between the mean and extremes.
# prob a numeric value or vector. The percentage of the area under the curve between the chosen extreme values

my.custnorm= function(x, CTmin= CTmin1, CTmax= CTmax1, prob=0.99){  
  x=as.vector(x[,1])
  diff= (CTmax-CTmin)/2-CTmin
  sd= -diff/qnorm(p = 1 - prob)
  custnorm(x, mean=(CTmax-CTmin)/2, diff=diff, prob=0.95)*sd/0.4
}

plot(1:60, my.custnorm(1:60, CTmin= CTmin1, CTmax= CTmax1, prob=0.99) )

m3 <- graf(y1, as.data.frame(train[,2]), prior = my.betaFun, l=100) #opt.l = TRUE ## adjust lengthscale l = 100,
plot(m3, prior=TRUE)
