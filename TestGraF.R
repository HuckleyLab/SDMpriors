#Modify Graf
#remove.packages("GRaF")
#devtools::install_github("HuckleyLab/GRaF",ref="test")

my.custnorm= function(x, CTmin= CTmin1, CTmax= CTmax1, prob=0.99){  
  x=as.vector(x[,1])
  diff= (CTmax-CTmin)/2-CTmin
  sd= -diff/qnorm(p = 1 - prob)
  P=custnorm(x, mean=(CTmax-CTmin)/2, diff=diff, prob=0.95)*sd/0.4
  return(P)
}

y= pa.points$Real
x= as.data.frame(pa.env)
prior= my.custnorm

error = NULL; weights = NULL
l = NULL; opt.l = FALSE; theta.prior.pars = c(log(10), 1); hessian = FALSE; opt.control = list(); verbose = FALSE; method = 'Laplace'

## Functions must be able to deal with data frames

#graf.fit.laplace
y = y; x = as.matrix(x); mn = mn; l = l; wt = weights; e = error; verbose = verbose
tol  = 10 ^ -6; itmax = 50

m3 <- graf(pa.points$Real, as.data.frame(pa.env), prior = my.custnorm, l=100) #opt.l = TRUE ## adjust lengthscale l = 100,
plot(m3, prior=TRUE)