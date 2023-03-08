#R code to run JAGS for CHILKAT SOCKEYE SALMON STATE SPACE MODEL
#notes----
#Created by Sara Miller (sara.miller@alaska.gov)
#Originally written for Karluk Lake Chinook by Matthew Catalano (mjc0028@aubirn.edu)
#changed with the help of Hamachan and Steve Fleischman (October 2016)
#Source of equations: Fleischman et al. 2013
#RICKER STOCK-RECRUIT RELATIONSHIP WITH AR1 ERRORS;R[y] IS THE TOTAL RETURN FROM BROOD YEAR y
#Y=40; A=3; a.min=4; a.max=6;THERE ARE A TOTAL OF Y+A-1 = 41 + 3 - 1 = 43 BROOD YRS REPRESENTED IN DATA (INCL FORECAST)

#load----
rm(list=ls(all=TRUE))
library(coda)
library(emdbook)
library(MASS)
library(gtools)
library(Hmisc)
library(rbugs)
library(R2OpenBUGS)
library(rjags)
library(lattice)
library(rmarkdown)
library(boot)                                                                          
library(ggplot2)                                                                          
library(dplyr)
library(tidyr)
library(shinystan)
library(reshape2)
library(grid)
library(runjags)
library(matrixStats)
library(gdata)
#data----
rawdat<-as.data.frame(read.csv("data/Chilkat_Sockeye.csv",header=T) )
nyrs<-as.numeric(length(rawdat$year))
fyr<-min(rawdat$year)
lyr<-max(rawdat$year)
nages<-3
a.min<-4
a.max<-6
A<-3

#data clean----
year <- as.numeric(as.character(rawdat$year))
DS <- as.numeric(as.character(rawdat$DS))
DS.cv <- as.numeric(as.character(rawdat$DS.cv)) 
mr <- as.numeric(as.character(rawdat$mr))
mr.cv <- as.numeric(as.character(rawdat$mr.cv)) 
weir <- as.numeric(as.character(rawdat$weir)) 
weir.cv <- as.numeric(as.character(rawdat$weir.cv)) 
Hbelow <- as.numeric(as.character(rawdat$Hbelow))
Hbelow.cv <- as.numeric(as.character(rawdat$Hbelow.cv))
Fry <- as.numeric(as.character(rawdat$Fry))
x<-as.matrix(rawdat[,substr(colnames(rawdat), 1,1)=="x"])#age comp count data matrix
colnames(x)<-NULL
n.a<-rowSums(x)#effective sample size

#analysis----
#State Space Model Function
mod=function(){
  for (y in (A+a.min):(Y+A-1)) {
    log.R[y] ~ dnorm(log.R.mean2[y],tau.R)
    R[y] <- exp(log.R[y])
    log.R.mean1[y] <- log(S[y-a.max]) + lnalpha - beta * S[y-a.max]-gam*Fry[y-a.max] #eq.#1
    #log.R.mean1[y] <- log(S[y-a.max]) + lnalpha - beta * S[y-a.max]
    log.resid[y] <- log(R[y])-log.R.mean1[y]
  }
  log.R.mean2[A+a.min] <- log.R.mean1[A+a.min] + phi * log.resid.0 #Eq.2
  for (y in (A+a.min+1):(Y+A-1)) {
    log.R.mean2[y] <- log.R.mean1[y] + phi * log.resid[y-1]
  }
#PRIORS
    lnalpha ~ dnorm(0,1.0E-6)%_%T(0,3)#uninformative;
    #normal distribution with mean 0 and large variance; constrained to be >0 since more biologically conservative (pg. 406)
    beta ~ dnorm(0,1.0E-6)%_%T(0,)   #uninformative; normal distrib; constrained to be >0           
    phi ~ dnorm(0,1.0E-6)%_%T(-0.98,0.98)#uninformative; btw -1 and 1
    mean.log.RO ~ dnorm(0,1.0E-6)#uninformative (initial returns); mean of initial returns
    tau.RO ~ dgamma(0.001,0.001)#uninformative (initial returns); variance of initial returns
    log.resid.0 ~ dnorm(0,tau.red)#model residual with AR(1) process error
    tau.R ~ dgamma(0.001,0.001) #uninformative
    gam ~ dnorm(0,1.0E-6)%_%T(-0.98,0.98)#uninformative (for fry plants)
    
    sigma.R <- 1 / sqrt(tau.R)
    alpha <- exp(lnalpha)
    sigma.RO <- 1 / sqrt(tau.RO)#uninformative; informative (not shown) prior had large effect on parameter itself 
    #and initial R values, but not on key model quantities
    tau.red <- tau.R * (1-phi*phi)
    lnalpha.c <- lnalpha + (sigma.R * sigma.R / 2 / (1-phi*phi) )

#THE FIRST SEVERAL COHORTS ORIGINATE FROM UNMONITORED SPAWNING EVENTS
#DRAW THESE RETURNS FROM A COMMON LOGNORMAL DISTRIBUTION 
  R.O<-exp(mean.log.RO)
  for (y in 1:a.max) { 
    log.R[y] ~ dnorm(mean.log.RO,tau.RO) 
    R[y] <- exp(log.R[y]) 
  }
#REFERENCE POINTS (WITH CORRECTION FOR LOGNORMAL SKEWNESS)
  S.max <- 1 / beta #Eq.20
  alpha.c <- min(exp(lnalpha.c),1.0E4)
  #positive.lna.c <- step(lnalpha.c)
  #lnalpha.c.nonneg <- lnalpha.c * positive.lna.c
  S.eq.c <- lnalpha.c * S.max 
  #peterman.approx.c <- (0.5 - 0.65*pow(lnalpha.c.nonneg,1.27) / (8.7 +pow(lnalpha.c.nonneg,1.27)))
  U.msy.c <- lnalpha.c * (0.5-0.07*lnalpha.c)
  S.msy.c <- S.eq.c *(0.5-0.07*lnalpha.c) 
  #U.max.c <- 1 - 1 / exp(lnalpha.c.nonneg) 

  positive.lna.c <- step(lnalpha.c)
  lnalpha.c.nonneg <- lnalpha.c * positive.lna.c
  S.eq.c2 <- lnalpha.c.nonneg * S.max
  peterman.approx.c <- (0.5 - 0.65*pow(lnalpha.c.nonneg,1.27) / (8.7 +pow(lnalpha.c.nonneg,1.27)))
  U.msy.c2 <- lnalpha.c.nonneg * peterman.approx.c #peterman approximation; Eq. 17
  S.msy.c2 <- U.msy.c2 / beta  
  U.max.c2 <- 1 - 1 / exp(lnalpha.c.nonneg)
  
  
#GENERATE Y+A-1 = 42 MATURITY SCHEDULES, ONE PER BROOD YEAR USING THE DIRICHLET DISTRIB. (Eq.4-6)
  # "pi" (central tendency of "p"), and "D.scale" (dispersion of "p")
  D.scale ~ dunif(0,1)#uninformative
  D.sum <- 1 / (D.scale * D.scale)
  pi.2p ~ dbeta(1,1)#uninformative
  pi.1 ~ dbeta(1,1)#uninformative; Eq.6
  pi[1] <- pi.1
  pi[2] <- pi.2p * (1 - pi[1])
  pi[3] <- 1 - pi[1] - pi[2]

  
for (a in 1:A) {
  gamma[a] <- D.sum * pi[a]
  for (y in 1:(Y+A-1)) {                                                    
      g[y,a] ~ dgamma(gamma[a],0.001)
      p[y,a] <- g[y,a]/sum(g[y,])
    }
  }

#CALCULATE THE NUMBERS AT AGE MATRIX (Number returning to spawn at age in year y); Eq.3
#Product of the total return from brood year y-a and the prop. mature from cohort y-a returning at age a; 
  for(a in 1:A){
    for(y in a:(Y+(a-1))){
      N.ya[y-(a-1),(A+1-a)]<-p[y,(A+1-a)]*R[y]
    }
  }
  
#MULTINOMIAL SCALE SAMPLING ON TOTAL ANNUAL RETURN N; 
  for (y in 1:Y) {
    N[y] <- sum(N.ya[y,1:A])
    for (a in 1:A) {
      q[y,a] <- N.ya[y,a] / N[y]
    }
  } 
  for (t in 1:Y){  
    x[t,] ~ dmulti(q[t,],n.a[t])
  }  
   
#HARVESTS BELOW THE WEIR (No harvest above weir);
  for (y in 1:Y) {
    mu.Hbelow[y] ~ dbeta(1,1)
    H.below[y] <- mu.Hbelow[y] * N[y]
    log.Hb[y] <- log(H.below[y])
    tau.log.Hb[y] <- 1 / log(Hbelow.cv[y]*Hbelow.cv[y] + 1)
    Hbelow[y] ~ dlnorm(log.Hb[y],tau.log.Hb[y])
    S[y] <- max(N[y] - H.below[y], 1)#number of fish reaching weir = total run abundance minus harvest
    log.S[y] <- log(S[y])
  }
  
  log.q.weir ~ dnorm(0,1.0E-4)
  log.q.mr ~ dnorm(0,1.0E-4)
    
  for (y in 1:Y) { 
    tau.log.weir[y] <- 1 / log(weir.cv[y]*weir.cv[y] + 1) 
    tau.log.mr[y]<- 1 / log(mr.cv[y]*mr.cv[y] + 1) 
    log.qS.weir[y] <- log.q.weir+log.S[y]
    log.qS.mr[y] <- log.q.mr+log.S[y]
   
    weir[y]~ dlnorm(log.qS.weir[y],tau.log.weir[y]) 
    mr[y]~ dlnorm(log.qS.mr[y],tau.log.mr[y]) 
   
    tau.log.ds[y] <- 1 / log(DS.cv[y]*DS.cv[y] + 1) #direct estimates of escapement from DIDSON
    DS[y] ~ dlnorm(log.S[y],tau.log.ds[y]) #Eq.9; Spawning escapement=DIDSON count
    qS.mr[y]<-exp(log.qS.mr[y]) 
    qS.weir[y]<-exp(log.qS.weir[y]) 
  }
  q.mr<-exp(log.q.mr) 
  q.weir<-exp(log.q.weir) 
}

#write the model to a text file to be called by WinBUGS
model_file_loc=paste("code/Chilkat_Sockeye4.txt", sep="")
write.model(mod, paste("code/Chilkat_Sockeye4.txt", sep=""))

#MODEL DATA (see csv file)
dat=list(Y = nyrs, A=nages, a.min=a.min, a.max=a.max,
         x=x, DS=DS, DS.cv=DS.cv, mr=mr, mr.cv=mr.cv, weir=weir,
         weir.cv=weir.cv, Hbelow=Hbelow,Hbelow.cv=Hbelow.cv,
         n.a=n.a,Fry=Fry)


#Inital Values (these were generated from openBUGS)
log.q.weir.init = 0.34
log.q.mr.init = -8.4
D.scale.init = 0.19
beta.init = 1.055E-5
lnalpha.init = 1.6
log.resid.0.init = 0.15
mean.log.RO.init = 11.9
phi.init = 0.5
pi.1.init = 0.06
pi.2p.init = 0.4
tau.R.init = 2.5
tau.RO.init = 76
#gam.init = -0.11
log.R.inits = c(
  11.5,11.74,11.74,11.81,12.01,
  11.93,11.89,12.0,12.13,12.38,
  12.15,11.33,11.68,12.03,12.7,
  11.86,12.12,12.35,12.07,12.14,
  11.89,11.31,12.46,12.96,12.39,
  10.58,11.86,11.86,11.28,11.84,
  10.75,10.84,10.57,11.87,12.5,
  11.14,11.44,12.26,11.95,12.5,
  11.73,12.65,12)
mu.Hbelow.inits = c(
  0.38,0.42,0.5577,0.4837,0.2877,
  0.3331,0.5634,0.3847,0.3488,0.6264,
  0.8416,0.5496,0.7568,0.5665,0.7012,
  0.5651,0.52,0.2808,0.466,0.5079,
  0.7911,0.31,0.3559,0.5935,0.6256,
  0.4982,0.52,0.5602,0.3936,0.3939,
  0.2851,0.3511,0.1943,0.3446,0.4913,
  0.2395,0.29,0.368,0.5966,0.2032,0.3)
g.inits = structure(.Data = c(
  1	,	1	,	1	,	1	,	1	, 1	,	1	,	1	,	1	,	1	, 1	,	1	,	1	,	1	,	1	, 1	,	1	,	1	,	1	,	1	,
  1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,
  1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,
  1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	, 1	,	1	,	1	,	1	,	1	, 1	,	1	,	1	,	1	,	1	,
  1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,
  1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	,	1	,
  1	,	1	,	1	,	1	,	1	,1	,	1	,	1	,	1	),
  .Dim = c(43,3))

#initial parameter values to be passed to JAGS

inits1=list(
  D.scale=D.scale.init,
  beta=beta.init,
  lnalpha=lnalpha.init,
  log.resid.0=log.resid.0.init,
  mean.log.RO=mean.log.RO.init,
  phi=phi.init,
  pi.1=pi.1.init,
  pi.2p=pi.2p.init,
  tau.R=tau.R.init,
  tau.RO=tau.RO.init,
  log.q.weir=log.q.weir.init,
  log.q.mr=log.q.mr.init)# pass the initials to JAGS

inits=list(inits1,inits1,inits1)
#Define the parameters (nodes) of interest 
parameters=c('S.eq.c','S.msy.c','U.msy.c','alpha','beta',
 'lnalpha','lnalpha.c','phi','sigma.R','log.resid.0', 'mean.log.RO',
	'S','R','N','log.resid','mu.Hbelow','pi','H.below','N.ya',
	'p','q', 'S.max','D.sum','q.weir', 'q.mr','D.scale','sigma.RO',
 'log.qS.weir', 'log.qS.mr','gam','qS.weir', 'qS.mr',
  'S.eq.c2', 'U.msy.c2', 'S.msy.c2', 'U.max.c2 ')

#Run JAGS (test); if runs then run full iterations
ptm = proc.time()
jmod <- jags.model(file='code/Chilkat_Sockeye4.txt', data=dat, n.chains=3, inits=inits, n.adapt=1000) #10000
x<-update(jmod, n.iter=1000, by=10, progress.bar='text', DIC=T, n.burnin=100) #10,000,000 per chain; 3 chains; thin by 1000
endtime = proc.time()-ptm
endtime[3]
endtime[3]/60
post <- coda.samples(jmod, parameters, n.iter=1000, thin=10, n.burnin=100)
endtime = proc.time()-ptm
endtime[3]/60/60  
post.samp <- post

#Run JAGS 
ptm = proc.time()
jmod <- jags.model(file='code/Chilkat_Sockeye4.txt', data=dat, n.chains=3, inits=inits, n.adapt=10000) #10000
x<-update(jmod, n.iter=1000000, by=1000, progress.bar='text', DIC=T, n.burnin=10000) #10,000,000 per chain; 3 chains; thin by 1000
endtime = proc.time()-ptm
endtime[3]
endtime[3]/60
post <- coda.samples(jmod, parameters, n.iter=1000000, thin=1000, n.burnin=10000)
endtime = proc.time()-ptm
endtime[3]/60/60  
post.samp <- post
save(post.samp, file = "results/Model 4/coda.samples.Model4.RData")
outJG<-as.matrix(post.samp)
write.csv(outJG, file = "results/Model 4/coda.samples.Model4.csv")#save mcmc list
#load("results/Model 4/coda.samples.Model4.RData")

#results----
#Numerical summary of each parameter (mean, median, quantiles)
summary<-summary(post)                     
stats<-summary$statistics;  colnames(stats)
quants<-summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants) 
statsquants <- statsquants[,c(1,2,4,5,7,9)] #select columns of interest
write.csv(statsquants, file= paste("results/Model 4/Model4_statsquants.csv") )    

#Gelman
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold=1.10#values less than 1.2 are generally considered converged
poor <- gel[gel[,1] > poor.threshold, ]
write.csv(poor, file= paste("results/Model 4/Model4_gelman.csv") )    

#DIC
dic.pD  <-dic.samples(jmod,n.iter=1000000, thin=1000,"pD",  n.burnin=10000)
dic.popt<-dic.samples(jmod,n.iter=1000000, thin=1000,"popt",n.burnin=10000)
dev1 <- sum(dic.pD[[1]])
pD   <- sum(dic.pD[[2]])
dic.pD <- dev1 + pD
dic.pD.summary <- data.frame(dev1, pD, dic.pD)
write.csv(dic.pD.summary, file=paste("results/Model 4/Model4_DIC.csv") ) 

#Create coda samples for horsetail plots and probability plots
post2 <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=1000000, thin=1000,n.burnin=10000) #n.iter=10000 & thin=10 is 1000 samples
x <- as.array(post2)
x <- data.frame(x)
coda <- x[,1:3]
coda<- rename.vars(coda, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))
write.csv(coda, file= paste("results/Model 4/Model4_coda.csv") ,row.names=FALSE)    # writes csv file

#SMSY 15th and 65th percentile
parameters=c('S.msy.c')
post3 <- coda.samples(jmod, parameters, n.iter=1000000, thin=1000,n.burnin=10000)
x<-as.array(post3)
S.MSY<-quantile(x, probs=c(0,0.15,0.50,0.65,1))
S.MSY <- data.frame(S.MSY )
x<-as.array(post3)
x <- data.frame(x)
x <- rename.vars(x, from=c("X1","X2","X3"), to=c("Chain.1","Chain.2", "Chain.3"))
write.csv(S.MSY, file= paste("results/Model 4/Model4_Percentile_Method.csv") )
write.csv(x, file= paste("results/Model 4/Model4_SMSY_data.csv") )
png("figures/Model 4/Model4_S.MSY.png", res=600, height=4.5, width=8, units="in")
plot(post3)
dev.off()

parameters=c('alpha')
post4 <- coda.samples(jmod, parameters, n.iter=1000000, thin=1000,n.burnin=10000)
x<-as.array(post4)
x <- data.frame(x)
x <- rename.vars(x, from=c("X1","X2","X3"), to=c("Chain.1","Chain.2", "Chain.3"))
write.csv(x, file= paste("results/Model4_alpha.csv") )
png("figures/Model 4/Model4_alpha.png", res=600, height=4.5, width=8, units="in")
plot(post4)
dev.off()

parameters=c('beta')
post5 <- coda.samples(jmod, parameters, n.iter=1000000, thin=1000,n.burnin=10000)
x<-as.array(post5)
x <- data.frame(x)
x <- rename.vars(x, from=c("X1","X2","X3"), to=c("Chain.1","Chain.2", "Chain.3"))
write.csv(x, file= paste("results/Model4_beta.csv") )
png("figures/Model 4/Model4_beta.png", res=600, height=4.5, width=8, units="in")
plot(post5)
dev.off()

#Density and time series plots
post.samp <- post
nvars <- dim(post.samp[[1]])[2]
nsamps <- dim(post.samp[[1]])[1]
int <- 25
pdf("figures/Model 4/Model4_profiles.pdf",height=6, width=8)
for(j in seq(1,nvars,int)){
  par(mfrow=c(5,4),mai=c(0.3,0.3,0.2,0.2))
  for(i in 0:(int-1)){
    mindat=min(c(post.samp[[1]][,i+j],post.samp[[1]][,i+j]))
    maxdat=max(c(post.samp[[1]][,i+j],post.samp[[1]][,i+j]))
    plot(density(post.samp[[1]][,i+j]),col='blue',main=colnames(post.samp[[1]])[i+j],xlim=c(mindat,maxdat))
    lines(density(post.samp[[2]][,i+j]),col='red')
    lines(density(post.samp[[3]][,i+j]),col='green')
    
    plot(as.numeric(post.samp[[1]][,i+j]),col='blue',main=colnames(post.samp[[1]])[i+j],ylim=c(mindat,maxdat),type='l')
    lines(as.numeric(post.samp[[2]][,i+j]),col='red')
    lines(density(post.samp[[3]][,i+j]),col='green')
  }}
dev.off()
dev.off()

#Autocorrelation plot of each parameter
windows(record=T)
pdf(paste("figures/Model 4/Model4_autocorr.pdf"),onefile=T,useDingbats=F)
autocorr.plot(post, lag.max=5)
dev.off()
dev.off()
autocorr.summary<-autocorr.diag(post)
autocorr.summary<-data.frame(autocorr.summary)
write.csv(autocorr.summary, file=paste("results/Model 4/Model4_autocorr.csv") ) 



