## Date: May 2023
# Author: Aminath Shausan 

############ Model ########################
# y  = beta.0 + beta.1*x + noise, 
# noise ~ N(0, sigma2)
#assume 2 groups 
###############################
.libPaths("~/Documents/R/r-libraries-bayesian-tutes")

rm(list = ls())
set.seed(963258)

install.packages("LaplacesDemon", lib="~/Documents/R/r-libraries-bayesian-tutes")


library(invgamma)
library(mvtnorm)
library(scatterplot3d)
library(ggplot2)
library(LaplacesDemon)


#### simulate X data for each group ####
J <-2 #number of grps
N <-  100 # number of obs. in each grp
x1 <- rnorm(N, mean = 10, sd = 5)
mean(x1)
X1 <- cbind(rep(1, N), x1) 
x2 <- rnorm(N, mean = 5, sd = 5)
mean(x2)
X2 <- cbind(rep(1, N), x2) 

X <- rbind(X1,X2) #design matrix
#### simulate y data from priors ####



#Eg1.define hyper parameters and draw from priors
sigma_sqrd0 <- 2
a0 <- 1
b0 <- 1/sigma_sqrd0
n0 <- 2+2 # p=2 here
nu0 <- 1
g <- N*J
trX_X <- t(X)%*%X  #XtX
inv_trX_X = solve(trX_X)

Lambda0 <- g*sigma_sqrd0 * inv_trX_X 
mu0 <- c(5,1)
S0 <- Lambda0
inv_S0 <- solve(S0)

### prior distributions
si_sqrd <- rgamma(1, shape =a0,   scale = b0)
cov_betas <- rinvwishart(n0, inv_S0)
mu_betas <- rmvnorm(1, mean=mu0,   sigma = Lambda0)
sigma_sqrd_gp1 <- rinvgamma(1, shape = nu0/2, scale =  si_sqrd * nu0/2) #grp 1 [0.41]
sigma_sqrd_gp2 <- rinvgamma(1, shape = nu0/2, scale =  si_sqrd * nu0/2) #grp2 [0.13]
betas_grp1 <- rmvnorm(1, mean=mu_betas,   sigma = cov_betas) #grp1 [8.79, 6.51]
betas_grp2 <- rmvnorm(1, mean=mu_betas,   sigma = cov_betas) #grp2 [8.92, 13.1]

### simulate model from priors 
mu_y_grp1 <- X1 %*% t(betas_grp1)
mu_y_grp2 <- X2 %*% t(betas_grp2)
sigma_y_grp1  <- sigma_sqrd_gp1 *diag(length(mu_y_grp1))
sigma_y_grp2  <- sigma_sqrd_gp2 *diag(length(mu_y_grp2))

y_grp1 <- rmvnorm(1, mean=mu_y_grp1,   sigma = sigma_y_grp1 )
y_grp2 <- rmvnorm(1, mean=mu_y_grp2,   sigma = sigma_y_grp2)

######################### store data
data <- data.frame(t(y_grp1) , t(y_grp2), x1, x2)
colnames(data)<- c("y1", "y2", "x1", "x2")
data_lst <- list(data = data, betas_grp1 = betas_grp1, betas_grp2 = betas_grp2,
                 sigma_sqrd_gp1 = sigma_sqrd_gp1, sigma_sqrd_gp2)

##### plot data. (this doesn't look good)
plt_data <- ggplot(data, aes( )) +
  #geom_point(aes(x = x1, y = y1, colour = "blue")) + 
  geom_point(aes(x = x2, y = y2, colour = "red")) + 
  theme_bw() +
  ggtitle("Synthetic data set")
plt_data

#########################
save.image("data1_HierLinReg.RData")
##################################

#simulate data from a given distribution

x1 <- 