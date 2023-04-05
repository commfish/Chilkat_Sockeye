mod=function(){
  for (y in (A+a.min):(Y+A-1)) {
    log.R[y] ~ dnorm(log.R.mean2[y],tau.R)
    R[y] <- exp(log.R[y])
    log.R.mean1[y] <- log(S[y-a.max]) + lnalpha - beta * S[y-a.max] 
    log.resid[y] <- log(R[y])-log.R.mean1[y]
  }
  log.R.mean2[A+a.min] <- log.R.mean1[A+a.min] + phi * log.resid.0 
  for (y in (A+a.min+1):(Y+A-1)) {
    log.R.mean2[y] <- log.R.mean1[y] + phi * log.resid[y-1]
  }
  
  # PRIORS
  lnalpha ~ dnorm(0,1.0E-6)%_%T(0,3)#uninformative;*
  #normal distribution with mean 0 and large variance; constrained to be >0 since more biologically conservative (pg. 406)
  beta ~ dnorm(0,1.0E-6)%_%T(0,)   #uninformative; normal distrib; constrained to be >0 *          
  phi ~ dnorm(0,1.0E-6)%_%T(-0.98,0.98)#uninformative; btw -1 and 1
  mean.log.RO ~ dnorm(0,1.0E-6)#uninformative (initial returns); mean of initial returns
  tau.RO ~ dgamma(0.001,0.001)#uninformative (initial returns); variance of initial returns
  log.resid.0 ~ dnorm(0,tau.red)#model residual with AR(1) process error
  tau.R ~ dgamma(0.001,0.001) #uninformative
  sigma.R <- 1 / sqrt(tau.R)
  alpha <- exp(lnalpha)
  sigma.RO <- 1 / sqrt(tau.RO)#uninformative; informative (not shown) prior had large effect on parameter itself 
  #and initial R values, but not on key model quantities
  tau.red <- tau.R * (1-phi*phi)
  lnalpha.c <- lnalpha + (sigma.R * sigma.R / 2 / (1-phi*phi) ) 
  
  # THE FIRST SEVERAL COHORTS ORIGINATE FROM UNMONITORED SPAWNING EVENTS
  # DRAW THESE RETURNS FROM A COMMON LOGNORMAL DISTRIBUTION 
  R.O<-exp(mean.log.RO)
  for (y in 1:a.max) { 
    log.R[y] ~ dnorm(mean.log.RO,tau.RO) 
    R[y] <- exp(log.R[y]) 
  }
  
  #REFERENCE POINTS (WITH CORRECTION FOR LOGNORMAL SKEWNESS)
  S.max <- 1 / beta 
  alpha.c <- min(exp(lnalpha.c),1.0E4)
  S.eq.c <- lnalpha.c * S.max #Eq.21
  U.msy.c <- lnalpha.c * (0.5-0.07*lnalpha.c)
  S.msy.c <- S.eq.c *(0.5-0.07*lnalpha.c)  
  # S.eq <- lnalpha * S.max #Eq.21 (no correction applied)
  # U.msy <- lnalpha * (0.5-0.07*lnalpha) (no correction applied)
  # S.msy <- S.eq *(0.5-0.07*lnalpha)  (no correction applied)
  
  positive.lna.c <- step(lnalpha.c)
  lnalpha.c.nonneg <- lnalpha.c * positive.lna.c
  S.eq.c2 <- lnalpha.c.nonneg * S.max 
  peterman.approx.c <- (0.5 - 0.65*pow(lnalpha.c.nonneg,1.27) / (8.7 +pow(lnalpha.c.nonneg,1.27)))
  U.msy.c2 <- lnalpha.c.nonneg * peterman.approx.c 
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
  
  #MULTINOMIAL SCALE SAMPLING ON TOTAL ANNUAL RETURN N; Eq.13
  for (y in 1:Y) {
    N[y] <- sum(N.ya[y,1:A])
    for (a in 1:A) {
      q[y,a] <- N.ya[y,a] / N[y]
    }
  } 
  for (t in 1:Y){  
    x[t,] ~ dmulti(q[t,],n.a[t])
  }  
  
  #HARVESTS BELOW THE WEIR (No harvest above weir);Eq.8-11
  for (y in 1:Y) {
    mu.hbelow[y] ~ dbeta(1,1)
    h.below[y] <- mu.hbelow[y] * N[y]
    log.hb[y] <- log(h.below[y])
    tau.log.hb[y] <- 1 / log(hbelow.cv[y]*hbelow.cv[y] + 1)
    hbelow[y] ~ dlnorm(log.hb[y],tau.log.hb[y])
    S[y] <- max(N[y] - h.below[y], 1)#Eq.8;number of fish reaching weir = total run abundance minus harvest
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
    
    tau.log.ds[y] <- 1 / log(DS.cv[y]*DS.cv[y] + 1) 
    DS[y] ~ dlnorm(log.S[y],tau.log.ds[y]) 
    qS.mr[y]<-exp(log.qS.mr[y]) 
    qS.weir[y]<-exp(log.qS.weir[y]) 
  }
  q.mr<-exp(log.q.mr) 
  q.weir<-exp(log.q.weir) 
}


