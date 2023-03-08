# this code gets sourced from this file, 
# and creates the "inits" object, which is the used in the main script

# initial values 
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
mu.hbelow.inits = c(
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