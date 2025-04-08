# Pima Indian Women Analysis
#
# Get data
library(MASS) # Pima data in this package
library(invgamma)
head(Pima.tr)

bp <- Pima.tr$bp[Pima.tr$type=="No"]
hist(bp)

n <- length(bp)
ybar <- mean(bp)

# Histogram is fairly symmetric, proceed using
# a normal likelihood/data model.
#
# y_i ~iid N(mu, sigma2)
# mu ~ N(m, v)
# sigma2 ~ IG(a,b)
# 
# First need to select values for prior parameters

# create a gibbs sampler function wrapper
gibbs.sampler <- function(data,m,v,a,b,nsamp,nburn,nthin){

  n <- length(data)
  ybar <- mean(data)

  # empty vectors to store MCMC samples
  mu <- numeric()
  sigma2 <- numeric()

  mu[1] <- mean(data)
  sigma2[1] <- var(data)

  # number of MCMC samples
  J <- nsamp

  # Start our Gibbs Samples
  for(j in 2:J){
	
    # start by sampling mu from its complete conditional	
    mstar <- ((n/sigma2[j-1])*ybar + (1/v)*m)/(n/sigma2[j-1] + 1/v)
    vstar <- 1/(n/sigma2[j-1] + 1/v)
    mu[j] <- rnorm(1, mstar, sqrt(vstar))	

    # next sample sigma2 from its complete conditional
    astar <- a + n/2
    bstar <- b + 0.5*sum((data - mu[j])^2)	
    sigma2[j]	<- rinvgamma(1, astar, bstar)
  }
  cbind(mu[-c(1:nburn)][seq(1,(J-nburn),by=nthin)],
        sigma2[-c(1:nburn)][seq(1,(J-nburn),by=nthin)])
}

post_samp_mu_sig2_pr1 <- gibbs.sampler(data=bp, m=80, v=8, a=5, b=400,
                                      nsamp=105000, nburn=5000, nthin=10)

post_samp_mu_sig2_pr2 <- gibbs.sampler(data=bp, m=0, v=1000, a=0.01, b=0.01,
                                      nsamp=105000, nburn=5000, nthin=10)

mu_pr1 <- post_samp_mu_sig2_pr1[,1]
sigma2_pr1 <- post_samp_mu_sig2_pr1[,2]

mu_pr2 <- post_samp_mu_sig2_pr2[,1] 
sigma2_pr2 <- post_samp_mu_sig2_pr2[,2]

# check convergence for prior 1
# seems like there is no trend up or down so safe to conclude
# that the samples have converged.  
plot(mu_pr1, type='l')
plot(sigma2_pr2, type='l')

# Now check mixing.  Do this with the ACF plot.
# Appears that the mixing is good because the lag1
# autocorrelation is between the blue lines.
acf(mu_pr1)
acf(sigma2_pr1)

# Now since convergence is determined to hold and mixing
# is good, we can proceed with inference.
# Bivariate plot of the posterior distribution
library(MASS)
plot(mu_pr1, sigma2_pr1)
contour(kde2d(mu_pr1, sigma2_pr1), add=TRUE, col='red',lwd=2)

image(kde2d(mu_pr1, sigma2_pr1))
contour(kde2d(mu_pr1, sigma2_pr1), add=TRUE, col='blue',lwd=2)

# Test H0: mu <= 80 vs H1: mu > 80
mean(mu_pr1 > 80)

# Can we determine of blood pressure is too high
# or too low?  Use a credible interva
quantile(mu_pr1, c(0.025, 0.975))


# Plot the marginal posterior distribution of mu for each prior on same plot
plot(density(mu_pr2), col='red', 
     xlab=expression(mu), ylab="posterior density",
     main="", ylim=c(0,0.45), xlim=c(65,75))
lines(density(mu_pr1), col='blue')
legend(x="topleft", legend=c("Prior 1","Prior 2"), col=c("red","blue"), lty=1)
