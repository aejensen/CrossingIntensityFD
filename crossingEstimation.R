library(mvtnorm)
library(parallel)
library(np)
library(ks)
library(fdapace)
library(boot)

#Avoid annoying messages from np
options(np.messages=FALSE)

rm(list=ls())

#Squared exponential covariance function + noise
cSE <- function(s, t, par) {
	if(length(par) != 2) stop("wrong number of parameters")
	
  tau <- abs(s-t)
	l <- par[1]
	epsilon <- par[2]

	exp(-(tau^2 / (2*l^2))) + as.numeric(tau == 0)*epsilon
}

#Rational quadratic covariance function
cRQ <- function(s, t, par) {
	if(length(par) != 3) stop("wrong number of parameters")
	
  tau <- abs(s-t)
	alpha <- par[1]
	lambda <- par[2]
	epsilon <- par[3]
  
	(1 + tau^2 / (2*alpha*lambda^2))^(-alpha) + as.numeric(tau == 0)*epsilon
}

#Fourier covariance function
cF <- function(s, t, par) {
  phi1 <- function(t) sqrt(2)*cos(2*pi*t)
  phi2 <- function(t) sqrt(2)*sin(4*pi*t)
  phi3 <- function(t) sqrt(2)*cos(4*pi*t)
	
  3*phi1(s)*phi1(t) + 2*phi2(s)*phi2(t) + 1*phi3(s)*phi3(t)
}

#Simulate functional data
simulateProcesses <- function(n, t, type) {
	if(type == 1) {
		covMat <- outer(t, t, cF)
	} else if(type == 2) {
    covMat <- outer(t, t, cSE, par=c(0.1, 0))
	} else {
		stop("wrong type")
	}
  dat <- rmvnorm(n, sigma=covMat)
}

#Get derivative of function
getDerivative <- function(t, y, spar=NULL, all.knots=TRUE) {
  s <- smooth.spline(t, y, spar=spar, all.knots=all.knots)
  predict(s, deriv=1)$y
}

#Empirical estimate
empEstimate <- function(y, thres=0) {
	sum(y[1:(length(y) - 1)] > thres & y[2:length(y)] < thres)
}

#Do FPCA analyse of functional data set
getFPCA <- function(t, y) {
  tList <- lapply(1:n, function(k) t)
  dList <- lapply(1:n, function(k) y[k,])
  pca <- fdapace::FPCA(dList, tList)
  
  t <- pca$workGrid
  mu <- pca$mu
  phi <- pca$phi
  scores <- pca$xiEst
  
  muD <- getDerivative(t, mu, 0.4, FALSE)
  phiD <- apply(phi, 2, function(f) getDerivative(t, f, 0.4, FALSE))
  
  yEst <- matrix(NA, nrow(y), length(mu))
  ydEst <- matrix(NA, nrow(y), length(mu))
  for(i in 1:nrow(y)) {
  	yEst[i,] <- mu
  	ydEst[i,] <- muD
  	for(j in 1:ncol(pca$phi)) {
  		yEst[i,] <- yEst[i,] + scores[i,j] * phi[,j]
  		ydEst[i,] <- ydEst[i,] + scores[i,j] * phiD[,j]
  	}
  }
    
  list(t = t, mu = mu, phi = phi, scores = scores,
       muD = muD, phiD = phiD, lambda=pca$lambda,
       yest = yEst, ydest = ydEst)
}

#Remove at some point
doNaive <- function(k) {
  mu <- c(mean(d[,k]), mean(dDeriv[,k]))
  cMat <- cov(cbind(d[,k], dDeriv[,k]))
  integrate(function(z) abs(z) * dmvnorm(cbind(0, z), mu, cMat), -Inf, Inf)$value
}

#Remove at some point
doNaive2 <- function(pcaRes, k) {
  mu <- c(mean(pcaRes$yest[,k]), mean(pcaRes$ydest[,k]))
  cMat <- cov(cbind(pcaRes$yest[,k], pcaRes$ydest[,k]))
  integrate(function(z) abs(z) * dmvnorm(cbind(0, z), mu, cMat), -Inf, Inf)$value
}

#Do PCA based estimation assuming normal distributed scores
doPCAnormal <- function(pcaRes, u, k, type="total") {
  muVec <- rbind(pcaRes$mu[k], pcaRes$muD[k])
  
  phiMat <- rbind(pcaRes$phi[k,], pcaRes$phiD[k,])
  #sigma <- diag(pca$lambda)
  sigma <- cov(pcaRes$scores)
  covMat <- phiMat %*% sigma %*% t(phiMat)

  if(type == "total") {
    #out <- integrate(function(z) abs(z) * dmvnorm(cbind(0, z), muVec, covMat), -Inf, Inf)$value
  	out <- pracma::quadinf(function(z) abs(z) * dmvnorm(cbind(u, z), muVec, covMat), -Inf, Inf)$Q
  } else if(type == "upper") {
  	#out <- integrate(function(z) max(0, z) * dmvnorm(cbind(0, z), muVec, covMat), -Inf, Inf)$value
  	out <- pracma::quadinf(function(z) max(0, z) * dmvnorm(cbind(u, z), muVec, covMat), -Inf, Inf)$Q
  } else if(type == "lower") {
  	#out = integrate(function(z) max(0, -z) * dmvnorm(cbind(0, z), muVec, covMat), -Inf, Inf)$value
  	out <- pracma::quadinf(function(z) max(0, -z) * dmvnorm(cbind(u, z), muVec, covMat), -Inf, Inf)$Q
  } else {
  	stop("wrong type")
  }
  out
}

doPCAkernel <- function(pcaRes, k) {
  
}

#Do non-parametric sequential estimation
doNonParametricKernel <- function(pcaRes, u, k, type="total") {
	#'arg' should be one of “cv.ml”, “cv.ls”, “normal-reference”
  kEst <- np::npudens(~x + xd, data=data.frame(x=pcaRes$yest[,k], xd=pcaRes$ydest[,k]))
  
  if(type == "total") {
  	out <- pracma::quadinf(function(z) abs(z) * predict(kEst, newdata=data.frame(x=u, xd=z)), -Inf, Inf)$Q
  } else if(type == "upper") {
  	out <- pracma::quadinf(function(z) max(0, z) * predict(kEst, newdata=data.frame(x=u, xd=z)), -Inf, Inf)$Q
  } else if(type == "lower") {
  	out <- pracma::quadinf(function(z) max(0, -z) * predict(kEst, newdata=data.frame(x=u, xd=z)), -Inf, Inf)$Q
  } else {
  	stop("wrong type")
  }
  
  out
}
 
#Main function for estimating the crossing intensity
crossingIntensity <- function(u, pcaRes, direction, method, parallel=FALSE, mc.cores) {
	t <- pcaRes$t
	
	if(method == "normal") {
    if(parallel) {
    	cross <- unlist(parallel::mclapply(1:length(t), function(k) doPCAnormal(pcaRes, u = u, k, direction), mc.cores=mc.cores))
		} else {
			cross <- sapply(1:length(t), function(k) doPCAnormal(pcaRes, u = u, k, direction))
		}
	} else if(method == "pca") {
    #
		stop("not implemented yet")
    #
	} else if(method == "nonparametric") {
    if(parallel) {
    	cross <- unlist(parallel::mclapply(1:length(t), function(k) doNonParametricKernel(pca, u = u, k, direction), mc.cores=mc.cores))
    } else {
      cross <- sapply(1:length(t), function(k) doNonParametricKernel(pca, u = u, k, direction))
    }
	} else {
		stop("wrong method")
	}
	
	list(x = t, y = cross)
}



tSeq <- seq(0, 1, length.out=100)
n <- 50
d <- simulateProcesses(n, tSeq, type=1)
#dDeriv <- t(apply(d, 1, function(y) getDerivative(tSeq, y)))

matplot(tSeq, t(d), type="l")

pca <- getFPCA(tSeq, d)
dim(pca$phi)

cross1 <- crossingIntensity(0, pca, "total", "normal", parallel=TRUE, mc.cores=6)
cross2 <- crossingIntensity(0, pca, "total", "nonparametric", parallel=TRUE, mc.cores=6)

plot(cross2, type="l")
lines(cross1, col=2)

mean(cross2$y)
mean(cross1$y)
