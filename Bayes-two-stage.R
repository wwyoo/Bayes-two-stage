#Author: William Weimin Yoo
#Code to implement our proposed Bayesian two-stage, and to draw Figures 1, 2 and the second box-plot
#in Figure 3 of the paper
############################################################################
library(splines)
library(MASS)
library(crs)

par(mai=c(0,0.5,0,0.1))
#################################################################################
#calculate sigma2tilde for empirical Bayes
#write the inverse matrix sandwiched in the middle using another form, so as to speed up computation
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
#do least squares
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

    trainz <- W$z[-i]
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

#-----------------------------------------------------------------------------------
#calculates maximum of a quadratic polynomial (equation (7.1) in paper)
muhatf <- function(x){
  A <- matrix(c(2 * x[4], x[6], x[6], 2 * x[5]), nrow = 2, ncol = 2, byrow = TRUE)
  b <- -c(x[2], x[3])
  ans <- solve(A, b)
  return(as.numeric(ans))
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
M0 <- f(0.5,0.5)

#plot f0 (Figure (a))
v <- seq(0,1, length=51)
w <- seq(0,1, length=51)
fframe <- expand.grid(v=v, w=w)
fframe$f <- f(fframe$v, fframe$w)
persp(v,w, matrix(fframe$f, 51, 51), theta = 50, phi = 20, expand = 0.5, zlab="", xlab="x", ylab="y", main="") -> res

points(trans3d(mu0[1], mu0[2], f(mu0[1], mu0[2]), pmat=res), col=1, type = "b", pch = 19, lwd = 10)  #plot true mu

#true error varance
sigma0 <- 0.1

#Data generation
n1 <- 30
x <- seq(0, 1, length = n1)
y <- seq(0, 1, length = n1)
obsframe <- expand.grid(x = x, y = y)
nfirst <- n1 * n1  #number of first stage samples

nrep <- 1000  #number of Monte Carlo replicates
mseB2mu <- rep(0, nrep)
mseB2M <- rep(0, nrep)

#prior hyperparameters
tau2 <- 1
tau2inv <- 1 / tau2
q <- 4 

set.seed(100)
obsframe$z <- f(obsframe$x, obsframe$y) + sigma0 * rnorm(nfirst)

####################################################################
#determine optimal number of B-splines 
#we choose the number of B-spline that maximizes its posterior 
####################################################################
#Nmax = 20
#Nx <- 1:Nmax
#Ny <- 1:Nmax
#Ncom = expand.grid(x=Nx, y=Ny)
#npair = nrow(Ncom)
#mseBF = matrix(0, npair, 2)
#piJ = rep(0, npair)

#for(i in 1:npair){
# J1 <- Ncom[i, 1] + q 
# J2 <- Ncom[i, 2] + q

# eta <- rep(0, J1 * J2)
# Omegainv <- tau2inv * diag(rep(1, J1 * J2))

# m1 = seq(from = 0, to = 1, length.out = Ncom[i, 1] + 2)[-c(1, Ncom[i, 1] + 2)]
# m2 = seq(from = 0, to = 1, length.out = Ncom[i, 2] + 2)[-c(1, Ncom[i, 2] + 2)]

# Bx <- bs(x = obsframe$x, knots = m1, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
# By <- bs(x = obsframe$y, knots = m2, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))

# B <- tensor.prod.model.matrix(list(Bx, By))

# M <- crossprod(B) + Omegainv
# cM <- t(chol(M))

# sigma2hat <- Bsigma2hat(y = obsframe$z, cmiddle = cM, B = B, eta = eta)

# piJ[i] <- -(nfirst/2)*log(sigma2hat) - log(prod(diag(cM)))

 #mseBF[i, ] <- optimal(W = obsframe, m1 = m1, m2 = m2, Omegainv = Omegainv, eta = eta)  #choose N1, N2 using CV
# print(i)
#}

#plot the posterior surface of J1,J2
#grid <- cbind(Ncom, piJ)
#maxJ <- which(piJ == max(piJ))
#persp(Nx+q, Ny+q, matrix(grid[,3], Nmax, Nmax), theta = 45, phi = 20, expand = 0.5, zlab = "", xlab = expression(J1), ylab = expression(J2), main = "", ticktype = "detailed") -> res
#points(trans3d(grid[maxJ,1]+q, grid[maxJ,2]+q, grid[maxJ,3], pmat=res), col=1, type = "b", pch = 19, lwd = 5)
#text(trans3d(grid[maxJ,1]+q, grid[maxJ,2]+q, grid[maxJ,3]+100, pmat=res), labels=expression("(7,9)" ))

#---------------------------------------------------------------------------------------
N1 <- 3  #this is optimal by empirical Bayes or MAP
N2 <- 5
J1 <- N1 + q
J2 <- N2 + q

#construct basis Gram matrix and prior parameters
eta <- rep(0, J1 * J2)
Omegainv <- tau2inv * diag(rep(1, J1 * J2))

m1 <- seq(from = 0, to = 1, length.out = N1 + 2)[-c(1, N1 + 2)]
m2 <- seq(from = 0, to = 1, length.out = N2 + 2)[-c(1, N2 + 2)]

Bx <- bs(x = obsframe$x, knots = m1, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
By <- bs(x = obsframe$y, knots = m2, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))

B <- tensor.prod.model.matrix(list(Bx, By))
M <- crossprod(B) + Omegainv
cM <- t(chol(M))

u1 <- seq(0, 1, length=10)
u2 <- seq(0, 1, length=10)
newframe <- expand.grid(u1=u1, u2=u2)

bx <- predict(Bx, newframe$u1)
by <- predict(By, newframe$u2)
b <- tensor.prod.model.matrix(list(bx, by))

#second stage parameters
n2 <- 864  #number of second stage observations
#deltacan1 <- seq(0, 0.2, by = 0.01)[-1]
#deltacan2 <- seq(0, 0.2, by = 0.01)[-1]
#deltacan <- expand.grid(deltacan1, deltacan2)
#ndelta <- nrow(deltacan)
nsam <- 1000  #number of samples from the second stage posterior
mutrue <- matrix(rep(mu0, nsam), 2, nsam, byrow = FALSE)
xi <- rep(0, 6) #prior mean 
#deltap <- rep(0, ndelta)
#n <- nfirst + n2  #total samples

#---------------------------------------------------------------------------------------------------------------
#To determine how much to localize
#sample from the first stage posterior and take the maximum span

pfmean <- pfmeanf(y = obsframe$z, cmiddle = cM, B = B, b = b, Omegainv = Omegainv, eta = eta) 
sigma2first <- as.numeric(Bsigma2hat(y = obsframe$z, cmiddle = cM, B = B, eta = eta))
shrink <- pfvarf(cmiddle = cM, b = b)

fsam <- mvrnorm(n = 1000, mu = pfmean, Sigma = sigma2first * shrink)

maxind <- apply(fsam, 1, function(x){which(x==max(x))})

mu1 <- newframe[maxind, 1:2]

delta1 <- abs(max(mu1[,1]) - min(mu1[,1]) )
delta2 <- abs(max(mu1[,2]) - min(mu1[,2]) )

#-----------------------------------------------------------------------------------------------------------------
#for(j in 1:ndelta){  #uncheck this to determine best delta

# delta1 <- deltacan[j, 1]
# delta2 <- deltacan[j, 2]

#delta1 <- 0.11
#delta2 <- 0.13

#increase resolution to 1024 locations
u1 <- seq(0, 1, length=32)
u2 <- seq(0, 1, length=32)
newframe <- expand.grid(u1=u1, u2=u2)

bx <- predict(Bx, newframe$u1)
by <- predict(By, newframe$u2)
b <- tensor.prod.model.matrix(list(bx, by))

set.seed(100)
for(i in 1:nrep){   #Monte Carlo replicates to compute boxplot
 obsframe$z <- f(obsframe$x, obsframe$y) + sigma0 * rnorm(nfirst)  #generate second stage data

 pfmean <- pfmeanf(y = obsframe$z, cmiddle = cM, B = B, b = b, Omegainv = Omegainv, eta = eta) 
 sigma2first <- as.numeric(Bsigma2hat(y = obsframe$z, cmiddle = cM, B = B, eta = eta))

 pmaxindex <- which(pfmean == max(pfmean))
 pmutilde <- as.numeric(newframe[pmaxindex, 1:2])   #center of credible set 

 #plot first stage observations (Figure (a))
 points(trans3d(obsframe$x,obsframe$y, obsframe$z, pmat=res), col=1, lwd=1)

 #plot posterior mean (Figure (b))
 persp(u1, u2, matrix(pfmean, 32), theta = 50, phi = 20, expand = 0.5, zlab = "", xlab = "x", ylab = "y", main = "") -> res

###########################################################################################
#second stage
###########################################################################################
 Oinverse <- diag(c(1, delta1, delta2, delta1^2, delta2^2, delta1*delta2))^2  
 
 x2nd <- runif(n2, pmutilde[1] -  delta1, pmutilde[1] + delta1)
 y2nd <- runif(n2, pmutilde[2] -  delta2, pmutilde[2] + delta2)

 z2nd <- f(x2nd, y2nd) + sigma0 * rnorm(n2)

 #plot second stage observations (Figure (b))
 points(trans3d(x2nd, y2nd, z2nd, pmat = res), col = 1)

 x2z <- x2nd - pmutilde[1]  #must shift
 y2z <- y2nd - pmutilde[2]

 Xz <- cbind(1, x2z, y2z, x2z ^ 2, y2z ^ 2, x2z * y2z)  #original design

 Mz <- crossprod(Xz) + Oinverse
 cMz <- t(chol(Mz))

 psigma2theta <- as.numeric(Bsigma2hat(y = z2nd, cmiddle = cMz, B = Xz, eta = xi))  #empirical Bayes for error variance
 #sigma2total <- (nfirst*sigma2first + n2*psigma2theta) / n
 
 #compute posterior for thetas (polynomial coefficients)
 pthetameanz <- solve(Mz, crossprod(Xz, z2nd) + Oinverse %*% xi)
 pthetavarz <- chol2inv(t(cMz))
 
 ptheta <- mvrnorm(n = nsam, mu = pthetameanz, Sigma = psigma2theta * pthetavarz)  #draw samples from posterior of theta

 pmuhat <- pmutilde + apply(ptheta, 1, muhatf)   #translate back to original domain, samples from the second stage posterior of mu

 mseB2mu[i] <- sqrt(sum((rowMeans(pmuhat) - mu0)^2))  #compute root mse for mu
 print(i)
}

#deltap[j] <- mean(mseB2)
#print(j)
#}

#Figure 2 of paper, requires mseB1 and mseF2 from programs to compute Bayes 1 stage and Frequentist 2 stage
par(mai=c(0.5,0.5,0.1,0.1))
boxplot(mseB1,mseB2mu,mseF2mu,ylim=c(0,0.01),names=c("Bayes 1 stage", "Bayes 2 stage", "Frequentist 2 stage"))  

#plot density of mu
#pmu2 <- kde2d(pmuhat[1, ], pmuhat[2, ], n = 50, lims = c(0.38, 0.47, 0.4, 0.47))
#persp(pmu2$x, pmu2$y, matrix(pmu2$z, 50, 50), theta = 30, phi = 30, expand = 0.5, zlab = "", ylab = "", xlab = "", ticktype = "detailed")
