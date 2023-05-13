## Date: Apr 27, 2023
# Author: Aminath Shausan 

###########################################
# Model is 
#  y  = beta.0 + beta.1*x + noise, 
# noise ~ N(0, sigma2)
# in this model, there are three parameters to estimate:
#     beta0, beta1 and sigma2)

#Ref: (Chapter 14) Bayesian Data Analysis, Third Edition, 3rd Edition
#Gelman, Andrew, author. Carlin, John, author. Stern, Hal, author. Dunson, David, author. Vehtari, Aki, author. Rubin, Donald, author.
#3rd edition; 2013

#The model can be whiten as: y|beta,sigma,X ~ N(X*beta, sigma2*I)
#Priors: beta ~ N(mu, tau2) ,  sigma2 ~ IG(a,b)

#posterior: (beta, sigma2|y,X) ~ N(X*beta, sigma2)N(mu, tau2)IG(a,b)

##### Gibbs sampler ###########
# posterior can be written as: 
#     p(beta, sigma2|y,X) = p(beta|sigma2,y,X)*p(sigma2|beta,y,X)

#step1: if sigma2 is given, draw posterior beta from: 
#      beta|y,sigma2,X ~ N (beta_hat, s), were
#      beta_hat = inv(tr(X)X)*tr(X)y
#      s = sigma2*inv(tr(X)X)

#step 2: if beta is given, draw posterior sigma2 from: 
#     sigma2|y, beta, X ~IG(an, bn), where 
#     an = N/2 + a,   b = tr(y-X*beta)(y-X*beta)/2
###################################################
rm(list = ls())

library(invgamma)
library(mvtnorm)
library(scatterplot3d)
library(ggplot2)


set.seed(963258)  # <-- to replicate results everytime you run this

#load data and explore structure and 1st 6 rows
load("Sim.Data.Rdata")
str(sim.dat)
head(sim.dat$dat)

#plot data
x11()
ggplot(sim.dat$dat, aes(x=x, y=y)) +
  geom_point() +
  theme_bw() +
  ggtitle("Synthetic data set")


## nor. of observations in the data (=500)
N <- length(sim.dat$dat$x)

X <- cbind(rep(1, N), sim.dat$dat$x)  # design matrix
y <- sim.dat$dat$y # response vector

##########
#precompute some expressions 
##########
trX_X <- t(X)%*%X  #XtX
beta_hat = solve(trX_X, t(X)%*%y) #this computes inv(XtX)Xty
inv_trX_X = solve(trX_X)   #

#initial values of betas and sigma2
betas = c(3,0.5)    
sigma2 = 2   

#nor. of mcmc iterations
n_iterations = 5000

#define matrices to hold  posterior samples
betas_post = matrix(data=NA, nrow=n_iterations, ncol=2)  
sigma_post = matrix(data = NA, nrow = n_iterations, ncol=1)  


for (i in 1:n_iterations){
#i =1
  betas = rmvnorm(n=1, beta_hat, sigma2 * inv_trX_X)  
  
  res = (y - X%*%t(betas))
  
  sigma2 = rinvgamma(1, N/2, t(res) %*% res * .5 ) 
  
  # save the results.
  betas_post[i,] = betas
  sigma_post[i,] = sigma2
}

#combine posterior paths for betas and sigma2
mcmc <-data.frame(cbind(betas_post, sigma_post))
colnames(mcmc)<- c("beta0", "beta1", "sigma2")
head(mcmc)
tail(mcmc)

##################################################
# MCMC diagnositics - check for convergence
# The standard checks for convergence are
# trace, density and autocorrelation plots

source("multiplot.R")

# discard burn-in values and apply thinning


burn.in<- floor(n_iterations/10)
burn.in #(=500)

thin<- 3  #(every 3rd iteration)
mcmc1<- mcmc[seq(from = burn.in, to=n_iterations, by = thin),]
dim(mcmc1)

diagnostics<- data.frame(cbind(iter = seq(1:dim(mcmc1)[1]), mcmc1))
full.results = list(diagnostics = diagnostics, sim.dat=sim.dat)
save(full.results, file = "MCMC_results2.Rdata")

## *************************************************
## Summary of posterior chains (mean and 95% CI) 
##                and comparison to the solution
## **************************************************

data.frame(parameter = c("beta0", "beta1", "sigma2"),
           posterior.mean = round(c(mean(diagnostics$beta0), mean(diagnostics$beta1), mean(diagnostics$sigma2)), 2),
           low.ci = round(c(quantile(diagnostics$beta0, probs = 0.025), quantile(diagnostics$beta1, probs = 0.025), 
                            quantile(diagnostics$sigma2, probs = 0.025)),2), 
           high.ci = round(c(quantile(diagnostics$beta0, probs = 0.975), quantile(diagnostics$beta1, probs = 0.975), 
                             quantile(diagnostics$sigma2, probs = 0.975)),2),
           Solution = c(sim.dat$B.0, sim.dat$B.1, sim.dat$noise))

## **********
## trace
## **********

trace_beta0<- ggplot(diagnostics, aes(x=iter, y=beta0)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta0")

trace_beta1<- ggplot(diagnostics, aes(x=iter, y=beta1)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta1")

trace_sigma2<- ggplot(diagnostics, aes(x=iter, y=sigma2)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace sigma2")

x11()

multiplot(trace_beta0, trace_beta1, trace_sigma2, cols = 3)


## ********
## Density
## ********

d_beta0<- ggplot(diagnostics, aes(x=beta0)) + geom_density() + 
  geom_vline(xintercept=sim.dat$B.0, col='red') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density beta0")
 

d_beta1<- ggplot(diagnostics, aes(x=beta1)) + geom_density() + 
  geom_vline(xintercept=sim.dat$B.1, col='red') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density beta1")

d_sigma2<- ggplot(diagnostics, aes(x=sigma2)) + geom_density() + 
  geom_vline(xintercept=sim.dat$noise, col='red') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density sigma2")

x11()
multiplot(d_beta0, d_beta1, d_sigma2, cols = 3)

## *****************
## Autocorrelation
## *****************

beta0_autocorr<- with(acf(diagnostics$beta0, plot=FALSE), data.frame(lag, acf))
beta0_ac<- ggplot(data = beta0_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr beta0")


beta1_autocorr<- with(acf(diagnostics$beta1, plot=FALSE), data.frame(lag, acf))
beta1_ac<- ggplot(data = beta1_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr beta1")


sigma2_autocorr<- with(acf(diagnostics$sigma2, plot=FALSE), data.frame(lag, acf))
sigma2_ac<- ggplot(data = sigma2_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr noise (sigma2)")


x11()
multiplot(beta0_ac, beta1_ac, sigma2_ac, cols = 3)

