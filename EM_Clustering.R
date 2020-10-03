#************************************
#**************EM by Sara Mohammadi**
#************************************

require(MASS)
require(mvtnorm)

# E step
E.step <- function(X, phi, N) {
    h <-with(phi, do.call(cbind,lapply(1:N, function(i)dmvnorm(X, mu[[i]], sig[[i]]))))
    h/rowSums(h)
}

# obtain the values of mu, sigma, and alpha that maximize the likelihood
M.step <- function(X, h, N) {
    covs <- lapply(1:N, function(i) cov.wt(X, h[,i]))
    mu <- lapply(covs, "[[", "center")
    sigma <- lapply(covs, "[[", "cov")
    alpha <- colMeans(h)
    list(mu = mu, sigma = sigma, alpha = alpha)
}

#log.likelihood
log.likelihood <- function(X, phi, N) {
    probab <- with(phi, do.call(cbind,lapply(1:N, function(i)alpha[i] * dmvnorm(X, mu[[i]], sigma[[i]]))))
    sum(log(rowSums(probab)))
}

#Plot current data and visualize clusters using contours
plot.em <- function(X, phi, N){
    xpts <- seq((rr <- range(X[,1]))[1], rr[2], length.out=50)
    ypts <- seq((rr <- range(X[,2]))[1], rr[2], length.out=50)
    plot(X)
    for(i in 1:N) {
        mixture.contour <-
            with(phi, outer(xpts,ypts,function(x,y)
                alpha[i]*dmvnorm(cbind(x,y),mean=mu[[i]],sigma=sig[[i]])))
        contour(xpts, ypts, mixture.contour, nlevels=5, col=i+1, add=TRUE, lwd=3)
    }
}

run.em <- function(X, N, plot = TRUE, sleep = 0.5, max.iter = 30, save=FALSE) {
    covs <- replicate(N, list(cov.wt(X[sample(nrow(X), 30),])))
    mu <- lapply(covs, "[[", "center")
    sig <- lapply(covs, "[[", "cov")
    alpha <- rep(1/N, N)
    phi <- list(mu = mu, sig = sig, alpha = alpha)
    for(i in 1:max.iter) {
        oldphi <- phi
        h <- E.step(X, phi, N)
        phi <- M.step(X, h, N)
        if(plot) {
            plot.em(X, phi, N)
            if(save) {
                dev.copy(jpeg, sprintf("output/em-%02d-clusters-%02d.jpg", N, i))
                dev.off()
            }
            Sys.sleep(sleep)
        }
        if((log.likelihood(X, phi, N) - log.likelihood(X, oldphi, N)) < 0.01)
            break
    }
    return(list(phi = phi, aic = 2*3*N - log.likelihood(X, phi, N)))
}


# Ex.

myData <- do.call(rbind, replicate(sample(5,1), list(rmvnorm(n=1000, rnorm(2, sd=30), diag(rgamma(2, 50))))))
plot(myData)

# Run EM for 1 to 6 clusters
(aics <- sapply(1:6, function(i) run.em(myData, i, sleep=0, max.iter=30)$aic))
which.min(aics)
nrow(myData)
(aics <- sapply(1:10, function(i) min(replicate(5, run.em(myData, i, sleep=0)$aic))))
which.min(aics)
nrow(myData)
