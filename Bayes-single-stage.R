#Author: William Weimin Yoo
#Code to implement one stage Bayesian estimation of the mode, and to draw the first box-plot
#in Figure 3 of the paper
#####################################################################################
library(splines)
library(MASS)
library(crs)

#################################################################################
#calculate sigma2tilde for empirical Bayes
Bsigma2hat <- function(y, cmiddle, B, eta){
  ydiff <- y - B %*% eta
  sigma2 = (crossprod(ydiff) - crossprod(forwardsolve(cmiddle, crossprod(B, ydiff)))) / length(y)
  return(sigma2)
}

#returns posterior mean
pfmeanf <- function(y, cmiddle, B, b, Omegainv, eta){
    ans = forwardsolve(cmiddle, crossprod(B, y) + Omegainv %*% eta)
    pmean = crossprod(t(b), backsolve(t(cmiddle), ans))
    return(pmean)
}

#calculates posterior variance
pfvarf <- function(cmiddle, b){ 
    pSigma = crossprod(forwardsolve(cmiddle, t(b)))
    return(pSigma)
}

#----------------------------------------------------------------------------------
#least squares estimator
fhatmeanf <- function(y, B, cBB, b){
  ans <- forwardsolve(cBB, crossprod(B, y))
  fmean <- crossprod(t(b), backsolve(t(cBB), ans)) 
  return(fmean)
}

#--------------------------------------------------------------------------------
#determine the best number of splines using leave-one-out cross validation (CV)
#this is not used but was included for reference
optimal <- function(W, m1, m2, Omegainv, eta){
  mse <- matrix(0, nrow(W), 2)
  for(i in 1:nrow(W)){
    testz <- W$z[i]
    testy <- W$y[i]
    testx <- W$x[i]

    trainz <- W$z[-i]  #leave one out
    trainy <- W$y[-i]
    trainx <- W$x[-i]

    Bx <- bs(x = trainx, knots = m1, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
    By <- bs(x = trainy, knots = m2, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
    B <- tensor.prod.model.matrix(list(Bx, By))
    BB = crossprod(B)
    V = BB + Omegainv
    cM <- t(chol(V))
    cBB <- t(chol(BB))
    bnewx <- predict(Bx, testx)
    bnewy <- predict(By, testy)
    bnew <- tensor.prod.model.matrix(list(bnewx, bnewy))

    pfmean = pfmeanf(y = trainz, cmiddle = cM, B = B, Omegainv = Omegainv, eta = eta, b = bnew)
    fhat = fhatmeanf(y = trainz, B = B, cBB = cBB, b = bnew)

    mse[i, 1] <- (pfmean - testz) ^ 2  #1st col Bayes

    mse[i, 2] <- (fhat - testz) ^ 2  #2nd col Frequentist
  }

return(colMeans(mse))
}

#################################################################################
#true function f0
f <- function(x, y){
   x <- 2*x - 1
   y <- 2*y - 1
   (1+exp(-5*(x)^2-2*(y)^4))*(cos(4*x)+cos(5*y))
}

#true mode
mu0 <- c(0.5, 0.5)

#true maximum
M0 <- f(mu0[1], mu0[2])

#true error varance
sigma0 <- 0.1

#Discrete uniform design
n1 <- 42
x <- seq(0, 1, length = n1)
y <- seq(0, 1, length = n1)
obsframe <- expand.grid(x = x, y = y)
nfirst <- n1 * n1  #number of first stage samples

nrep <- 1000  #number of Monte Carlo samples
mseB1mu <- rep(0, nrep)
mseB1M <- rep(0, nrep)

mutrue <- matrix(rep(mu0, nrep), nrep, 2, byrow=TRUE)

tau2 <- 1
tau2inv <- 1 / tau2
q <- 4 

#set.seed(100)
#obsframe$z <- f(obsframe$x, obsframe$y) + sigma0 * rnorm(nfirst)

####################################################################
#determine optimal number of B-splines 
#we choose the number of B-spline that maximizes its posterior 
####################################################################
#Nmax = 10
#Ncom = cbind(rep(1:Nmax, each = Nmax), rep(1:Nmax, times = Nmax))
#npair = nrow(Ncom)
#mseBF = matrix(0, npair, 2)
#piJ = rep(0, npair)

#  for(i in 1:npair){
#    J1 <- Ncom[i, 1] + q 
#    J2 <- Ncom[i, 2] + q

#    eta <- rep(0, J1 * J2)
#    Omegainv <- tau2inv * diag(rep(1, J1 * J2))

#    m1 = seq(from = 0, to = 1, length.out = Ncom[i, 1] + 2)[-c(1, Ncom[i, 1] + 2)]
#    m2 = seq(from = 0, to = 1, length.out = Ncom[i, 2] + 2)[-c(1, Ncom[i, 2] + 2)]

#    Bx <- bs(x = obsframe$x, knots = m1, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
#    By <- bs(x = obsframe$y, knots = m2, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))

#    B <- tensor.prod.model.matrix(list(Bx, By))

#    M <- crossprod(B) + Omegainv
#    cM <- t(chol(M))

#    sigma2hat <- Bsigma2hat(y = obsframe$z, cmiddle = cM, B = B, eta = eta)

#    piJ[i] <- -(nfirst/2)*log(sigma2hat) - log(prod(diag(cM)))  #marginal likelihood for number of basis J

#    #mseBF[i, ] <- optimal(W = obsframe, m1 = m1, m2 = m2, Omegainv = Omegainv, eta = eta)
#    print(i)
#  }

#---------------------------------------------------------------------------------------
N1 <- 5  #this is the optimal choice for 
N2 <- 5
J1 <- N1 + q
J2 <- N2 + q

eta <- rep(0, J1 * J2)
Omegainv <- tau2inv * diag(rep(1, J1 * J2))

m1 <- seq(from = 0, to = 1, length.out = N1 + 2)[-c(1, N1 + 2)]
m2 <- seq(from = 0, to = 1, length.out = N2 + 2)[-c(1, N2 + 2)]

Bx <- bs(x = obsframe$x, knots = m1, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
By <- bs(x = obsframe$y, knots = m2, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))

B <- tensor.prod.model.matrix(list(Bx, By))
M <- crossprod(B) + Omegainv
cM <- t(chol(M))

bx <- predict(Bx, obsframe$x)
by <- predict(By, obsframe$y)
b <- tensor.prod.model.matrix(list(bx, by))
pfvar <- pfvarf(cmiddle = cM, b = b)

nsam = 1000  #number of posterior samples

set.seed(100)

for(i in 1:nrep){
  obsframe$z <- f(obsframe$x, obsframe$y) + sigma0 * rnorm(nfirst)  #generate data

  pfmean <- pfmeanf(y = obsframe$z, cmiddle = cM, B = B, b = b, Omegainv = Omegainv, eta = eta)  #calculate posterior mean of f
  empsigma2hat <- as.numeric(Bsigma2hat(y = obsframe$z, cmiddle = cM, B = B, eta = eta))  #empirical Bayes for error variance

  fBayes <- mvrnorm(n = nsam, mu = pfmean, Sigma = empsigma2hat * pfvar)  #sample from posterior of f

  pmaxindex <- apply(fBayes, 1, function(x)which(x == max(x)))  #find maximum

  muB1 <- obsframe[pmaxindex, 1:2]  

  MaxB1 <- obsframe[pmaxindex, 3]

  muB1mean <- colMeans(muB1)  #one stage posterior mean for mode

  mseB1mu[i] <- sqrt(sum((muB1mean - mu0)^2))  #these are the values for the first box-plot in Figure 2
  mseB1M[i] <- abs(mean(MaxB1) - M0)
  print(i)
}

