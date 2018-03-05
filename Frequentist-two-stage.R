#This code implements the frequentist two-stage procedure by Belitser et al. (2012) Optimal
#two-stage procedures for estimating location and size of the
#maximum of a multivariate regression function, in the Annals of Statistics 44, number 6, pages 2850-2876

#The original code was supplied by one of the authors of the paper, and I have modified it for the purpose
#of comparing their procedure with our proposed Bayesian two-stage method
#################################################################################################################
#par(mai=c(0,0.1,0,0.1))

#choose optimal smoothning (bandwidth) of loess by leave-one-out cross validation
optimalband <- function(W, bandwidth){
  mse <- rep(0, nrow(W))
  for(i in 1:nrow(W)){
    testz <- W$z[i]
    testy <- W$y[i]
    testx <- W$x[i]

    trainz <- W$z[-i]  #leave one out
    trainy <- W$y[-i]
    trainx <- W$x[-i]

  m <- loess(trainz ~trainx + trainy, degree = 1, span=bandwidth)
  testobs <- data.frame(trainx = testx, trainy = testy)
  fit <- predict(m, newdata = testobs)
 
    mse[i] <- (fit - testz) ^ 2  
  }

return(mean(mse))
}

##################################################################################################################
# error standard deviation:

sigma <- 0.1


# true regression function, put in dataframe fframe:

f <- function(x, y){
   x <- 2*x - 1
   y <- 2*y - 1
   (1+exp(-5*(x)^2-2*(y)^4))*(cos(4*x)+cos(5*y))
}

#true mu
xfmax <- 0.5
yfmax <- 0.5

M0 <- f(xfmax, yfmax) 

# initial, uniformly spaced data, put in dataframe obsframe:
# total n1*n1 observations

nrep <- 1000
mseF2mu <- rep(0,nrep)
mseF2M <- rep(0, nrep)

n1 <- 30
x <- seq(0,1, length=n1)
y <- seq(0,1, length=n1)
obsframe <- expand.grid(x=x, y=y)

band <- seq(0.01,0.5,by=0.01)
mseband <- rep(0, length(band))


#set.seed(100)

#  obsframe$z <- f(obsframe$x, obsframe$y) + sigma*rnorm(n1*n1)

#for(i in 1:length(band)){
#    bandwidth <- band[i]
 
#    mseband[i] <- optimalband(W = obsframe, bandwidth = bandwidth)
#    print(i)
#}

deltacan1 <- seq(0, 0.2, by = 0.01)[-1]
deltacan2 <- seq(0, 0.2, by = 0.01)[-1]
#deltacan <- expand.grid(deltacan1, deltacan2)
deltacan <- cbind(deltacan1, deltacan2)
ndelta <- nrow(deltacan)
msedeltaF2 <-rep(0, ndelta)
 n3 <- 96 

#for(j in 1:ndelta){  #to choose optimal localization

  delta1 <- 0.06
  delta2 <- 0.06

set.seed(100)
for(i in 1:nrep){
  obsframe$z <- f(obsframe$x, obsframe$y) + sigma*rnorm(n1*n1)

# local linear regression to find stage one estimator:

  m <- loess(obsframe$z ~obsframe$x + obsframe$y, degree =1, span=0.02)
  fit <- predict(m)
  #persp(x,y, matrix(fit, n1, n1), theta = 50, phi = 20, expand = 0.5, zlab="", xlab="x", ylab="y", main="") -> res
  maxindex <- which(fit==max(fit))
  xmutilde <- obsframe$x[maxindex]
  ymutilde <- obsframe$y[maxindex]
  #points(trans3d(xmutilde,ymutilde, f(xfmax, yfmax), pmat=res), col=2, type = "b", pch = 19)


# new data

# new designpoints:

  x0 <- xmutilde
  y0 <- ymutilde
  x1 <- x0+delta1
  y1 <- y0
  x2 <- x0
  y2 <- y0+delta2
  x3 <- x0-delta1
  y3 <- y0
  x4 <- x0
  y4 <- y0-delta2
  x5 <- x0+delta1
  y5 <- y0-delta2
  x6 <- x0+delta1
  y6 <- y0+delta2
  x7 <- x0-delta1
  y7 <- y0+delta2
  x8 <- x0-delta1
  y8 <- y0-delta2

# regressors:

  p <- rep(c(x0, x1, x2, x3, x4, x5, x6, x7, x8), n3)
  q <- rep(c(y0, y1, y2, y3, y4, y5, y6, y7, y8), n3)
  r <- rep(c((x0)^2, (x1)^2, (x2)^2, (x3)^2, (x4)^2, (x5)^2, (x6)^2, (x7)^2, (x8)^2), n3)
  s <- rep(c((y0)^2, (y1)^2, (y2)^2, (y3)^2, (y4)^2, (y5)^2, (y6)^2, (y7)^2, (y8)^2), n3)
  t <- rep(c(x0*y0, x1*y1, x2*y2, x3*y3, x4*y4, x5*y5, x6*y6, x7*y7, x8*y8), n3)

# new observations:

  newz <- f(p,q) + sigma*rnorm(9*n3)
  #points(trans3d(c(x0, x1, x2, x3, x4, x5, x6, x7, x8), c(y0, y1, y2, y3, y4, y5, y6, y7, y8), newz, pmat=res), col=8)
  a <- lm(newz ~ p+q+r+s+t)$coefficients

# find muhat, the argmax of the quadratic surface

  A <- matrix(data=c(2*a[4], a[6], a[6], 2*a[5]), nrow=2, ncol=2, byrow=TRUE)
  b <- -c(a[2], a[3])
  u <- solve(A,b)
  xmuhat <- as.numeric(u[1])
  ymuhat <- as.numeric(u[2])

# plot of quadratic surface:

#g <- function(k,l) a[1] + a[2]*k + a[3]*l + a[4]*k^2+a[5]*l^2+a[6]*k*l
#gframe <- expand.grid(v=v, w=w)
#gframe$g <- g(gframe$v, gframe$w)
#persp(v,w, matrix(gframe$g, 51, 51), theta = 30, phi = 30, expand = 0.5, zlab="", xlab="x", ylab="y", main="") -> res


# plot the final estimator:

#points(trans3d(xmuhat,ymuhat, g(xmuhat, ymuhat), pmat=res), col="green", pch=19)

  mseF2mu[i] <- sqrt((xmuhat - xfmax)^2 + (ymuhat - yfmax)^2)  #produce third box-plot in Figure 2 of our paper
  #mseF2M[i] <- abs(a[1] + a[2]*xmuhat + a[3]*ymuhat + a[4]*xmuhat^2 + a[5]*ymuhat^2 + a[6]*xmuhat*ymuhat - M0) 
  print(i)
}

#msedeltaF2[j] <- mean(mseF2)
#print(j)
#}
