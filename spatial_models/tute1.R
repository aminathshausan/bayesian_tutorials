## Date: May 2023
# Author: Aminath Shausan 

####################################
#this program runs the 1st spatial model in section 4.7 of
#https://atlas.cancer.org.au/developing-a-cancer-atlas/Chapter_4.html
#using CARBayes
#This is a spatial Poisson model
###############################
.libPaths("~/Documents/R/r-libraries-bayesian-tutes")
install.packages("CARBayes", lib="~/Documents/R/r-libraries-bayesian-tutes")
rm(list = ls())
set.seed(963258)

library(ggplot2)
library(dplyr)

library(reshape2)
library(DiagrammeR)

library(plotly)
library(sf) # Simple Feature package has all the tools for working with 

library(maps) #to create spatial map
library(maptools) #to create spatial polygon from map
library(spdep) #to create adjacency matrix

library(CARBayes) #to fit bayesian model
####### prepare data
#data is from https://github.com/Andrew9Lawson/Bayesian-DM-code-examples/blob/master/SCcongen90_CARBayesRcode.txt

###### We first import the data and create a map (SpatialPolygons object).
SCcg90<-list(obs=c(0,7,1,5,1,1,5,16,0,17,4,0,0,1,1,7,1,3,0,0,8,2,13,7,0,8,0,3,2,4,1,11,0,1,2,3,3,8,6,14,3,11,6,0,1,5),
             expe=c(1.0807,6.3775,0.622,6.6854,0.9142,1.0744,5.6518,8.1682,0.5749,18.0989,2.174,1.6619,1.9321,1.6148,
                    1.6713,3.0819,1.7562,4.9952,0.9362,1.2001,6.1293,2.5604,15.8589,2.9437,1.0399,7.276,0.9739,2.064,2.7206,2.8275,0.9425,
                    8.828,0.3644,1.775,1.5111,1.5111,2.5321,4.5836,3.9647,15.0264,0.732,10.8292,5.9848,1.4357,1.9949,6.9807),
             pov=c(13.6,13.8,32.3,11,24.2,19.9,12.3,13.3,17,15.4,14.1,15.7,18,24.3,21.8,20.2,24.9,12,17.4,18.1,18.7,17.5,
                   10.6,13.8,22.8,13.7,21,12.6,14,14,26.9,9.4,17.8,23.1,23.1,14.4,10.8,22.1,10.1,13.6,15.9,11.2,18.3,14,26.4,10.6),
             inc=c(36.786,38.534,20.727,37.205,24.3,27.607,45.822,40.161,32.247,38.458,33.232,31.715,29.505,25.896,28.919,
                   30.776,25.552,42.886,34.297,29.96,34.009,35.008,41.658,34.109,27.65,34.654,27.117,39.04,33.698,32.32,25.144,
                   45.14,29.805,25.008,25.993,32.231,36.912,28.624,37.054,39.587,31.324,37.092,31.948,30.801,23.748,44.619))
SCcg90<-data.frame(SCcg90)

#####
# Create a map of counties in South Carolina (list of 4 items)
SC.county <- map("county", "South Carolina", plot = FALSE, fill = TRUE)
#df_SC.county <- data.frame(SC.county) #this doesn't work
SC.map <- map2SpatialPolygons(SC.county, IDs = SC.county$names)

x11()
plot(SC.map)
################################################

#Next, we need to create the spatial weights matrix which forms part of the 
#specification of the spatial random effect.

# Create binary, first-order, adjacency spatial weights matrix
#poly2nb() creates a neighbours list from the polygons based on 
#areas with contiguous boundaries. Each element of this list represents one area and contains the indices of its neighbours. For example, the second element contains the neighbours of area 2.
#nb2mat() uses neighbours list to create a spatial weights matrix.
W <- nb2mat(poly2nb(SC.map), style = "B")


################################
#Now we fit a model using CARBayes. Here we fit the BYM model.
# MCMC parameters
M.burnin <- 10000       # Number of burn-in iterations (discarded)
M <- 1000               # Number of iterations retained

# Fit model using CARBayes
#(spatial glm model, where random effects will have BYM conditional 
#autoregressive prior)
#here likelihood is: y_k ~ Poisson(mu_k)
#ln(mu_k) ~ beta0 + Offset_k + spatialRandomEffect_k
#spatialRandomEffect_k ~ BYM 
#so the parameters estimated are: 
#beta0, tau, sigma (latter 2 for spatial components)
set.seed(444)          # For reproducability
MCMC <- S.CARbym(
  formula = obs ~ 1 + offset(log(expe)), 
  data = SCcg90,
  family = "poisson",
  W = W,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  verbose = FALSE
)

# Broad summary of results
print(MCMC$summary) #converges if Geweke.diag is b/w (-1.96,1.96)
print(MCMC$summary.results)
#           Mean      2.5%  97.5%   n.sample % accept n.effective Geweke.diag
#(Intercept) 0.0319 -0.1228 0.1847     1000     44.4       339.1         1.3
#tau2        0.0159  0.0030 0.0652     1000    100.0        14.5        -1.8
#sigma2      0.0073  0.0020 0.0196     1000    100.0        29.5         1.4
print(MCMC)
MCMC$samples # A list containing the MCMC samples generated from the model, 
#where each element in the list is a matrix.
MCMC$modelfit 
#DIC           p.d          WAIC           p.w          LMPL 
#169.768575      3.056461    169.443486      2.515574    -84.778004 
#loglikelihood 
#-81.827827 

#the estimated (posterior mean) regression coefficients 
#in this model, its the intercept only.
coef(MCMC) 
#(Intercept) 
#0.03193377 

#the fitted values based on the posterior mean
#here this is a vector of 46 elements (same length as observed values)
fitted(MCMC)

#the estimated loglikelihood based on the posterior mean.
logLik(MCMC) #-81.82783

#the design matrix of covariates.
#here its 46 by 1 matrix (showing intercept only)
model.matrix(MCMC)

#The formula (as a text string) for the response, covariate and 
#offset part of the model.
formula(MCMC)
#obs ~ 1 + offset(log(expe))

#A text string describing the model that has been fitted.
MCMC$model
#"Likelihood model - Poisson (log link function)"
#"Random effects model - BYM CAR" 

#the design matrix
MCMC$X #same result as using model.matrix(MCMC)

####save MCMC results 
save(MCMC, file = "MCMC_tute1.Rdata")
