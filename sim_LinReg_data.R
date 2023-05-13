## Date: May 09, 2023
# Author: Aminath Shausan 

############ Model ########################
# y  = beta.0 + beta.1*x + noise, 
# noise ~ N(0, sigma2)

#The model can be whiten as: y|beta,sigma,X ~ N(X*beta, sigma2*I)
#Priors: beta ~ N(mu, tau2) ,  sigma2 ~ IG(a,b) 
#a shape, b = scale

#assume : mean beta0 = 5; mean beta1 = 1, mean sigma2 = 10 
###########################
.libPaths("~/Documents/R/r-libraries-bayesian-tutes")

rm(list = ls())

#install.packages("invgamma", lib="~/Documents/R/r-libraries-bayesian-tutes")
#install.packages("mvtnorm", lib="~/Documents/R/r-libraries-bayesian-tutes")
#install.packages("scatterplot3d", lib="~/Documents/R/r-libraries-bayesian-tutes")
#install.packages("ggplot2", lib="~/Documents/R/r-libraries-bayesian-tutes")

library(invgamma)
library(mvtnorm)
library(scatterplot3d)
library(ggplot2)

#Note = this simulation shows linearly decreasing trend
set.seed(963258)  # <-- to replicate results everytime you run this

N <- 100 #Number of observations

x1 <- runif(n=N, min=-10, max=10)
X <- cbind(rep(1, N), x1)

mu_betas <- c(5, 1)
sigma_betas <- 2*diag(length(mu_betas))

betas <- rmvnorm(n=1, mean=mu_betas, sigma=sigma_betas)
shape_sigma2 <- 0.95 
rate_sigma2 <- 1
scale_sigma2 <- 1
##### check inv gamma function to give mean about 10
x11() 
#sims.1<- rinvgamma(10000, shape = shape_sigma2, rate = rate_sigma2)

#sims.1<- rinvgamma(10000, shape = shape_sigma2, rate = 1, scale = 1)
#p.try1<- ggplot(data.frame(x = sims.1), aes(x=x)) + geom_density() +
#  ggtitle(paste("Package (invgamma) shape", shape_sigma2, " rate ", rate_sigma2,  
#                'mean', mean(sims.1), sep = ""))+ xlim(c(0, 3))

#p.try1
################

sigma2 <- rinvgamma(1, shape = shape_sigma2, rate = rate_sigma2, scale = scale_sigma2)
###########

#simulate data from model
mu_y <- X %*% t(betas)
sigma_y <- sigma2 *diag(length(mu_y))

y <-  rmvnorm(1, mean=mu_y,   sigma = sigma_y)

#store data points and prior values 
data <- data.frame(t(y) , x1)
colnames(data)<- c("y", "x")

data_lst <- list(data = data, beta0 = betas[1] , beta1 = betas[2], sigma2 = sigma2)

##########
#plot data 

plt_data <- ggplot(data, aes(x=x, y=y)) +
            geom_point() +
            theme_bw() +
            ggtitle("Synthetic data set")
plt_data

#####. save data
save(data_lst, file = "dataLinRgress.rData") 
