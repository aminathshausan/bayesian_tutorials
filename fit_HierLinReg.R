## Date: May 2023
# Author: Aminath Shausan 

############ Model ########################
#Ref: 
# y_ij  = beta.0_j + beta.1_j*x + noise_j, 
 
#assume 2 (m) groups, each with 100 (n) observations
###############################
.libPaths("~/Documents/R/r-libraries-bayesian-tutes")

rm(list = ls())
set.seed(963258)


library(invgamma)
library(mvtnorm)
library(scatterplot3d)
library(ggplot2)
library(LaplacesDemon)


######## c
#### load data 
load('dataHierLinRgress.RData')

data <- data_lst$data
y1 <- data$y1 #grp1 obs.  #numeric vector (double). ### t(y1) is a 1by100 matrix
y2 <- data$y2 #grp2 obs.  #numeric vector
x1 <- data$x1
x2 <- data$x2
#plot data
x11()
#pdf(file="plots/data_hierLinReg.pdf", width = 8, height = 7)
plt_data <- ggplot(data, aes( )) +
  #geom_point(aes(x = x1, y = y1, colour = "blue")) + 
  #geom_point(aes(x = x2, y = y2, colour = "red")) +
  geom_point(aes(x = x1, y = y1, colour = "blue")) + 
  geom_point(aes(x = x2, y = y2, colour = "red")) +
  scale_color_manual(name="X", labels=c("x1","x2"),
                     values=c("blue","red")) +
  theme_bw() +
  labs(title="Synthetic data",
       x ="x", y = "y")
plt_data
#dev.off

#save plot
 

### create y and X matrices

y <- as.matrix(append( y1,  y2)) # matrix (200 by 1)
X1 <-  cbind(rep(1, length(x1)), x1) #matrix array (100 by 2)
X2 <- cbind(rep(1, length(x2)), x2) #matrix array(100 by 2)
C <- rbind(X1, X2)  #matrix 200 by 2

#X <- matrix(data=0, nrow=length(y), ncol=4)
#X[1:100, 1] <- rep(1, length(x1))
#X[101:200, 2] <- rep(1, length(x2))
#X[1:100, 3] <- x1
#X[101:200, 4] <- x2
## some useful values 

trC_C <- t(C)%*%C  #XtX. (2by2 matrix) (new = 4by4 matrix)
inv_trC_C = solve(trC_C)    #(new = 4by4)
mu0 = solve(trC_C, t(C)%*%y) #OLS beta (2 by 1 matrix)

#trX1_X1 <- t(X1)%*%X1  #XtX. (2by2 matrix) (new = 4by4 matrix)
#inv_trX1_X1 = solve(trX1_X1)    #(new = 4by4)
#mu0_1 = solve(trX1_X1, t(X1)%*%y1) #OLS beta (2 by 1 matrix)

#trX2_X2 <- t(X2)%*%X2  #XtX. (2by2 matrix) (new = 4by4 matrix)
#inv_trX2_X2 = solve(trX2_X2)    #(new = 4by4)
#mu0_2 = solve(trX2_X2, t(X2)%*%y2) #OLS beta (2 by 1 matrix)


# define parameters for zi2
n <- 100 #nor of observatios in each grp
m <- 2 # nor of grps

##initial guess
a0 <- 1
b0 <- 1
nu0 <- 1
sigma2_1 <- 20 
sigma2_2 <- 20 


#zi2 <- rgamma(1, shape =shape_zi2 ,   scale = scale_zi2) #1.44
zi2 = 1.44 #initial value of xi_sqrd

## parameters for sigma2 per grp (initial guess)
#zi_sqrd_init <- 0.87
betas_1<- matrix(c(10 , 5)) 
betas_2<- matrix(c(10 , 5))

#parameters for Sigma
n0 <- 2+2 # p=2 here
#nu0 <- 1
nu_Sigma <- n0+m
mu_betas1 <- as.matrix(c(1,1))
mu_betas2 <- as.matrix(c(1,1))
S0 <- 200*2 * inv_trC_C
#S0_1 <- 200*2 * inv_trX1_X1
#S0_2 <- 200*2 * inv_trX2_X2
#parameters for mu_betas
invLambda0 <- solve(S0)
#invLambda0_1 <- solve(S0_1)
#invLambda0_2 <- solve(S0_2)

#parameters. for betas_j

######################
#initial configuration vector:
#(betas_1 = c(10 , 5),  betas_2 =c(10 , 5),   
#mu_betas_1 = 4.530221 1.810435, mu_betas_2 = 4.530221 1.810435,  Sigma, sigma2_1, sigma2_2, zi2) =
#mu_betas_1 = (1,1), mu_betas_2 = (1,1),
#Sigma is 2 by 2 matrix with row 1 = 54.22210 23.78163
#                             raw 2 = 23.78163 10.99527
#sigma2_1 = 20, sigma2_2=20, zi2 = 1.44)

##### define matrices to hold posterior samples 
n_iterations = 5000 #nor. of mcmc iterations

betas1_post <- matrix(data=NA, nrow=n_iterations, ncol=2)
betas2_post <- matrix(data=NA, nrow=n_iterations, ncol=2) 
mu_betas_post <- matrix(data=NA, nrow=n_iterations, ncol=2)
 
sigma1_post = matrix(data = NA, nrow = n_iterations, ncol=1)
sigma2_post = matrix(data = NA, nrow = n_iterations, ncol=1)

zi2_post <- matrix(data = NA, nrow = n_iterations, ncol=1)

#### sample posteriors 

for (i in 2:n_iterations){
i=1
  
  #zi2
  sum_invsigmas <- sum(1/sigma2_1 + 1/sigma2_2) 
  shape_zi2 <- a0 + m * nu0/2 #2
  scale_zi2 <-  b0 + nu0/2 * sum_invsigmas #1.05
  zi2 <- rgamma(1, shape =shape_zi2 ,   scale = scale_zi2)
  
  #sigma2_ per grp
  #X1 is a 100 by 2 matrix.  betas_1 is numeric but t(betas_1) is matrix
  shape_sigma2 <- (nu0 + n)/2
  sum_res_sqr_1 <- sum((y1 - X1%*% betas_1)^2)
  scale_sigmma2_1 <- (nu0 * zi2 + sum_res_sqr_1)/2 #88497
  sigma2_1 <- rinvgamma(1, shape_sigma2, scale_sigmma2_1)  #1691
  
  sum_res_sqr_2 <- sum((y2 - X2%*% betas_2)^2)
  scale_sigmma2_2 <- (nu0 * zi2 + sum_res_sqr_2)/2 #88497
  sigma2_2 <- rinvgamma(1, shape_sigma2, scale_sigmma2_2)  #1077
  
  # Sigma
  res_betas_1 <- (betas_1 - mu_betas1)%*%t((betas_1 - mu_betas1))
  res_betas_2 <- (betas_2 - mu_betas2)%*%t((betas_2 - mu_betas2))
  sum_res_betas <-  res_betas_1 + res_betas_2
  S0_plus_sum_res_betas <-  S0 + sum_res_betas
  S_Sigma <- solve(S0_plus_sum_res_betas )
    
  Sigma <- rinvwishart(nu_Sigma, S_Sigma)
  
  #mu_betas
  inv_Sigma <- solve(Sigma)
  A <- invLambda0 + m * inv_Sigma 
  
  Sigma_mu_betas <- solve(A)
  sum_betas <- betas_1+betas_2
  B1 <- invLambda0 %*% mu0
  B2 <- inv_Sigma %*% sum_betas
  B <- B1 + B2
  mu_mu_betas <- Sigma_mu_betas %*% B
  
  mu_betas <- t(rmvnorm(1, mean = mu_mu_betas, sigma =  Sigma_mu_betas))
  #mu_betas <- t(mu_betas)
  
  # betas_j
  trX1_X1 <- t(X1)%*%X1 
  p1 <- inv_Sigma +  1/sigma2_1 * trX1_X1
  Sigma_betas1 <- solve(p1)
  trX1_y1 <- t(X1) %*% y1
  q11 <- inv_Sigma %*% mu_betas
  q12 <- q11 + 1/sigma2_1 *trX1_y1
  mu_betas1 <- Sigma_betas1 %*% q12
  betas_1 <- t(rmvnorm(1, mean = mu_betas1, sigma = Sigma_betas1))
  
  trX2_X2 <- t(X2)%*%X2 
  p2 <- inv_Sigma +  1/sigma2_2 * trX2_X2
  Sigma_betas2 <- solve(p2)
  trX2_y2 <- t(X2) %*% y2
  q21 <- inv_Sigma %*% mu_betas
  q22 <- q21 + 1/sigma2_2 *trX2_y2
  mu_betas2 <- Sigma_betas2 %*% q22
  betas_2 <- t(rmvnorm(1, mean = mu_betas2, sigma = Sigma_betas2))
  
  # save the results.
  betas1_post[i,] = betas_1
  betas2_post[i,] = betas_2
  
  mu_betas_post[i,] =   mu_betas
  
  sigma1_post[i,] = sigma2_1
  sigma2_post[i,] = sigma2_2
  
  zi2_post[i,] = zi2
}

#combine posterior paths for betas and sigma2
mcmc <-data.frame(cbind(betas1_post, betas2_post, sigma1_post, sigma2_post, 
                        mu_betas_post,zi2_post))
colnames(mcmc)<- c("beta0_1", "beta1_1", "beta0_2", "beta1_2","sigma2_1", "sigma2_2",
                   "beta0", "beta1", "zi2")
head(mcmc)
tail(mcmc)

##################################################
# MCMC diagnositics - check for convergence
# trace, density and autocorrelation plots

source("multiplot.R")

# discard burn-in values and apply thinning
burn.in<- floor(n_iterations/10)
burn.in #(=500)

thin<- 3  #(every 3rd iteration)
mcmc1<- mcmc[seq(from = burn.in, to=n_iterations, by = thin),]
dim(mcmc1)

diagnostics<- data.frame(cbind(iter = seq(1:dim(mcmc1)[1]), mcmc1))
full.results = list(diagnostics = diagnostics, sim_dat=data)
save(full.results, file = "MCMC_results_hiera_lin_reg.Rdata")

## *************************************************
## *************************************************
## Summary of posterior chains (mean and 95% CI) 
##                and comparison to the solution
## **************************************************

data.frame(parameter = c("beta0_1", "beta1_1", "beta0_2", "beta1_2","sigma2_1", "sigma2_2",
                         "beta0", "beta1", "zi2"),
           posterior.mean = round(c(mean(diagnostics$beta0_1), mean(diagnostics$beta1_1),
                                    mean(diagnostics$beta0_2), mean(diagnostics$beta1_2),
                                    mean(diagnostics$sigma2_1), mean(diagnostics$sigma2_2),
                                    mean(diagnostics$beta0), mean(diagnostics$beta1),
                                    mean(diagnostics$zi2)
                                    ),2),
           low.ci = round(c(quantile(diagnostics$beta0_1, probs = 0.025), quantile(diagnostics$beta1_1, probs = 0.025), 
                            quantile(diagnostics$beta0_2, probs = 0.025), quantile(diagnostics$beta1_2, probs = 0.025),
                            quantile(diagnostics$sigma2_1, probs = 0.025),quantile(diagnostics$sigma2_2, probs = 0.025),
                            quantile(diagnostics$beta0, probs = 0.025), quantile(diagnostics$beta1, probs = 0.025),
                            quantile(diagnostics$zi2, probs = 0.025)
                            ),2), 
           high.ci = round(c(quantile(diagnostics$beta0_1, probs = 0.975), quantile(diagnostics$beta1_1, probs = 0.975), 
                            quantile(diagnostics$beta0_2, probs = 0.975), quantile(diagnostics$beta1_2, probs = 0.975),
                            quantile(diagnostics$sigma2_1, probs = 0.975),quantile(diagnostics$sigma2_2, probs = 0.975),
                            quantile(diagnostics$beta0, probs = 0.975), quantile(diagnostics$beta1, probs = 0.975),
                            quantile(diagnostics$zi2, probs = 0.975)
                            ),2), 
           Solution = c(data_lst$betas_grp1[1], data_lst$betas_grp1[2], 
                        data_lst$betas_grp2[1], data_lst$betas_grp2[2],
                        data_lst$sigma2_noise_grp1, data_lst$sigma2_noise_grp2,
                        6, 1.7, 2.5
                        ))

## **********
## **********
## trace
## **********
#grp1 betas
trace_beta0_1 <- ggplot(diagnostics, aes(x=iter, y=beta0_1)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta0 grp1")

trace_beta1_1 <- ggplot(diagnostics, aes(x=iter, y=beta1_1)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta1 grp1")

#grp2 betas
trace_beta0_2 <- ggplot(diagnostics, aes(x=iter, y=beta0_2)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta0 grp2")

trace_beta1_2 <- ggplot(diagnostics, aes(x=iter, y=beta1_2)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta1 grp2")

#grp1 and grp2 sigma squared (noise)
trace_sigma2_1<- ggplot(diagnostics, aes(x=iter, y=sigma2_1)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace sigma2 grp1")

trace_sigma2_2<- ggplot(diagnostics, aes(x=iter, y=sigma2_2)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace sigma2 grp2")

#population betas
trace_beta0 <- ggplot(diagnostics, aes(x=iter, y=beta0)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta0 population")

trace_beta1 <- ggplot(diagnostics, aes(x=iter, y=beta1)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace beta1 population")

# hyper parameter zi2
trace_zi2 <- ggplot(diagnostics, aes(x=iter, y=zi2)) + geom_line() +
  theme_bw() + theme(legend.position = "none") + ggtitle("Trace hyper parameter zi2")


x11()

#pdf(file="plots/data_hierLinReg.pdf", width = 8, height = 7)
#multiplot(trace_beta0_1, trace_beta1_1, trace_beta0_2, trace_beta1_2,
  #        trace_sigma2_1, trace_sigma2_2, trace_beta0, trace_beta1, trace_zi2,
 #         cols = 4)
multiplot(trace_beta0_1, trace_beta1_1, trace_beta0_2, trace_beta1_2,
          trace_sigma2_1, trace_sigma2_2, 
          cols = 4)
#dev.off()

## ********
## Density
## ********

d_beta0_1<- ggplot(diagnostics, aes(x=beta0_1)) + geom_density() + 
  #geom_vline(xintercept=data_lst$betas_grp1[1], col='red') +
  geom_vline(xintercept=mean(diagnostics$beta0_1), col='blue') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density beta0 grp1")


d_beta1_1<- ggplot(diagnostics, aes(x=beta1_1)) + geom_density() + 
  #geom_vline(xintercept=data_lst$betas_grp1[2], col='red') +
  geom_vline(xintercept=mean(diagnostics$beta1_1), col='blue') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density beta1 grp1")

d_beta0_2<- ggplot(diagnostics, aes(x=beta0_2)) + geom_density() + 
  geom_vline(xintercept=mean(diagnostics$beta0_2), col='blue') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density beta0 grp2")

d_beta1_2<- ggplot(diagnostics, aes(x=beta1_2)) + geom_density() + 
  geom_vline(xintercept=mean(diagnostics$beta1_2), col='blue') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density beta1 grp2")

d_sigma2_1<- ggplot(diagnostics, aes(x=sigma2_1)) + geom_density() + 
  geom_vline(xintercept=mean(diagnostics$sigma2_1), col='blue') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density sigma2 (noise) grp 1")

d_sigma2_2<- ggplot(diagnostics, aes(x=sigma2_2)) + geom_density() + 
  geom_vline(xintercept=mean(diagnostics$sigma2_2), col='blue') +
  theme_bw() + theme(legend.position = "none") + ggtitle("Density sigma2 (noise) grp 2")

x11()
multiplot(d_beta0_1, d_beta1_1, d_beta0_2, d_beta1_2, 
          d_sigma2_1, d_sigma2_2,
          cols = 4)

## *****************
## Autocorrelation
## *****************

beta0_1_autocorr<- with(acf(diagnostics$beta0_1, plot=FALSE), data.frame(lag, acf))
beta0_1_ac <- ggplot(data = beta0_1_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr beta0 grp1")

beta1_1_autocorr<- with(acf(diagnostics$beta1_1, plot=FALSE), data.frame(lag, acf))
beta1_1_ac <- ggplot(data = beta1_1_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr beta1 grp2")

beta0_2_autocorr<- with(acf(diagnostics$beta0_2, plot=FALSE), data.frame(lag, acf))
beta0_2_ac <- ggplot(data = beta0_2_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr beta0 grp2")

beta1_2_autocorr<- with(acf(diagnostics$beta1_2, plot=FALSE), data.frame(lag, acf))
beta1_1_ac <- ggplot(data = beta1_2_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr beta1 grp2")

sigma2_1_autocorr<- with(acf(diagnostics$sigma2_1, plot=FALSE), data.frame(lag, acf))
sigma2_1_ac<- ggplot(data = sigma2_1_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr sigma2 grp 1")

sigma2_2_autocorr<- with(acf(diagnostics$sigma2_2, plot=FALSE), data.frame(lag, acf))
sigma2_2_ac<- ggplot(data = sigma2_1_autocorr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Auto-corr sigma2 grp 2")

x11()
multiplot(beta0_1_ac, beta1_1_ac, beta0_2_ac, beta1_1_ac, 
          sigma2_1_ac, sigma2_2_ac, cols = 4)

