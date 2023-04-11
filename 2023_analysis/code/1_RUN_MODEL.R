# Chilkat sockeye state-space model
# authors: Sara Miller (ADF&G)
# contact: sara.miller@alaska.gov; 
# Last edited: March 2023
# must download program JAGS for this script to work

# warning: some of these packages mask commands, so need to specify the package when calling the fn
# runjags::extract vs. tidyr::extract
# coda::traceplot vs. R2jags::traceplot
# load libraries----
library(coda)
library(tidyverse)
library(R2jags)
library(rjags)
library(runjags)
library(R2OpenBUGS)
library(mcmcplots)
library(gsl)
library(stringr)
library(gdata)
library(ggplot2)
library(mcmcplots)
library(scales)
library(cowplot)
library("devtools")
devtools::install_github("commfish/fngr")
library(fngr)
library(extrafont)

# STEP 1: CHOOSE SETTINGS----
# source the model file
# this reads in a function called "mod" then writes the model to a text file to be called by JAGS if using rjags version
# if used R2Jags, can just use the "mod" object directly

jags.settings <- "full"  # "test" or "explore" or full"  # choose iteration version (see code below for details)

source("2023_analysis/code/model_source.R") 
print(mod)
model_file_loc=paste("2023_analysis/code/","chilkat_sockeye.txt", sep="") # where to write the model file
write.model(mod, model_file_loc)

# load custom functions
source('2023_analysis/code/functions.R')
source("2023_analysis/code/MCMC_CustomFunctions.R")

# choices of model runs
if(jags.settings == "test"){
  n.adapt.use <- 100 ; n.iter.use <- 1000;  n.burnin.use <- 100;   thin.use = 10
  by.use <- 10 # this is just for the progress bar
}

if(jags.settings == "explore"){
  n.adapt.use <- 10000 ; n.iter.use <- 10000; n.burnin.use <- 30000 ;   thin.use = 10
  by.use <- 100 # this is just for the progress bar
}

if(jags.settings == "full"){
  n.adapt.use <- 10000  ; n.iter.use <- 1000000    #1,000,000 per chain; 3 chains; thin by 1000
  n.burnin.use <- 10000  # consider increasing this?
  thin.use = 1000; by.use <- 1000 # this is just for the progress bar 
}

# STEP 2: READ IN DATA, MODEL, AND INITIAL VALUES----
# generates the object "dat"
source("2023_analysis/code/model_data.R")

# generates initial values
source("2023_analysis/code/model_inits.R")

# STEP 2: RUN THE MODEL AND PROCESS THE OUTPUT----
# 2 options: rjags or R2jags
# if "test" runs then do sensitivity tests with "explore", and final run with "full"
# run the R2jags package and produced outputs based on rjags
out.label <-  "R2Jags" # label to be used for the output folder (and for scenario comparisons)
package.use <- "R2jags"

# create output folder for model results
out.path <- paste0("2023_analysis/output/", out.label)
if(!exists(out.path)){dir.create(out.path)}

# This step does the actual MCMC sampling. All subsequent steps
# should just extract from "post" without rerunning the model.
# These parameters must be the same as the conv.pars for the R2jags program to run
parameters=c('S.eq.c','S.eq','S.msy.c','S.msy', 'U.msy.c','U.msy','alpha','beta','alpha.c',
             'lnalpha','lnalpha.c','phi','sigma.R','log.resid.0', 'mean.log.RO',
             'S','R','N','log.resid','mu.hbelow','pi','h.below','N.ya',
             'p','q', 'S.max','D.sum','q.weir', 'q.mr','D.scale','sigma.RO',
             'log.qS.weir', 'log.qS.mr','qS.weir', 'qS.mr',
             'S.eq.c2', 'U.msy.c2', 'S.msy.c2', 'U.max.c2')

# start the timer
# R2jags
start.jags <- proc.time()
if(package.use == "R2jags"){ # new version
  r2jags.out <- R2jags::jags(data = dat , inits = inits, 
                   parameters.to.save = parameters, model.file = mod,
                   n.chains = 3, 
                   n.iter = n.iter.use + n.burnin.use ,  
                   n.burnin = n.burnin.use, 
                   n.thin = thin.use, DIC = T)
  end.jags <- proc.time()   # store time for MCMC
  mcmc.samples <- r2jags.out$BUGSoutput$sims.matrix
  mcmc.summary <- r2jags.out$BUGSoutput$summary 
  
  # these are the same as the ones produced below
  write.csv(mcmc.samples[,c("beta","lnalpha","lnalpha.c")], file= paste0(out.path,"/coda.csv") ,row.names=FALSE)    # writes csv file
  write.csv(mcmc.summary, file= paste0(out.path,"/statsquants.csv"))    
  
  # this one is the same as coda.csv, except with all parameters 
  # - > not tracked in github
  write.csv(mcmc.samples, file= paste0(out.path,"/coda_all_parameters.csv") ,row.names=FALSE)    # writes csv file
  
# this only works for any single-value parameters
# parameters with annual values would need to be tested individually (E.g. S[1], S[2], etc)
#     -> build that later  
conv.pars <- c('S.eq.c','S.eq','S.msy.c','S.msy', 'U.msy.c','U.msy','alpha','beta',
                 'lnalpha','lnalpha.c','phi','sigma.R',
                 'S.eq.c2', 'U.msy.c2', 'S.msy.c2', 'U.max.c2','alpha.c')
  
conv.details <- checkConvergence(mcmc.out = r2jags.out, vars.check = conv.pars)

write.csv(conv.details,file=paste0(out.path,"/convergence_details.csv"), row.names=FALSE)

# for now, call it converged if gelman rubin and geweke are below critical values 
# for all the conv.pars (the acf handling is finicky)
# Note: converged = NOT flagged
conv.check <- !conv.details$Flag[conv.details$Check == "all.gelman.rubin"] & 
                   !conv.details$Flag[conv.details$Check == "all.geweke"]

if(conv.check){print("The R2jags model run converged for all the key variables!")}
if(!conv.check){print("The R2jags model run did not converge for all the key variables!")}
}

# run the rjags package and produced outputs based on rjags
out.label <-  "rjags" # label to be used for the output folder (and for scenario comparisons)
package.use <- "rjags"

# create output folder for model results
out.path <- paste0("2023_analysis/output/", out.label)
if(!exists(out.path)){dir.create(out.path)}

 if(package.use == "rjags"){
  parameters <- c('S.eq.c','S.eq','S.msy.c','S.msy', 'U.msy.c','U.msy', 'alpha','beta','alpha.c',
                  'lnalpha','lnalpha.c','phi','sigma.R','log.resid.0', 'mean.log.RO',
                  'S','R','N','log.resid','mu.hbelow','pi','h.below','N.ya',
                  'p','q', 'S.max','D.sum','q.weir', 'q.mr','D.scale','sigma.RO',
                  'log.qS.weir', 'log.qS.mr','qS.weir', 'qS.mr',
                  'S.eq.c2', 'U.msy.c2', 'S.msy.c2', 'U.max.c2')
  
  jmod <- rjags::jags.model(file='2023_analysis/code/chilkat_sockeye.txt', data=dat, n.chains=3, inits=inits, n.adapt=n.adapt.use) 
  stats::update(jmod, n.iter=n.iter.use, by=by.use, progress.bar='text', DIC=T, n.burnin=n.burnin.use) # this modifies the original object, function returns NULL
  post <- rjags::coda.samples(jmod, parameters, n.iter=n.iter.use, thin=thin.use, n.burnin=n.burnin.use)

end.jags <- proc.time()   # store time for MCMC
post.arr <- as.array(post) # convert to an accessible obj

# run the script that generates all the outputs for rjags
end.output  <- proc.time() 
print("Jags took")
print(end.jags - start.jags)
print("Output Processing took")
print(end.output - end.jags)
 }
 
    

