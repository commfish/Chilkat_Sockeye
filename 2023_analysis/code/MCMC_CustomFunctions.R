# THESE FUNCTIONS ARE FROM A DRAFT VERSION OF  A CODE PACKAGE FOR WSP METRICS
# FOR STATUS UPDATES, CONTACT
# Dr. Carrie Holt (Project Lead, Carrie.Holt"AT"dfo-mpo.gc.ca)
# Gottfried Pestal (Developer, gpestal"AT"solv.ca)

# modified here to use the  output from a R2jags::jags call 

# FUNCTION TO CALCULATE A CONVERGENCE SUMMARY FOR MCMC OUTPUT
checkConvergence <- function(mcmc.out, vars.check = c("slope","intercept","sigma")){
# function to calculate acf, acf critical value, and geweke score, then output a summary table
# mcmc.out = output from a R2jags::jags call 
# vars.check = variables for which to check convergence 

# Note: acf critical values should always be equal for all variables (same sample size),
# but calculating separately just in case

#geweke.vals <- c(abs(mcmc.out$gintercept),abs(mcmc.out$gslope),abs(mcmc.out$gsig))
#names(geweke.vals) <- list("gintercept","gslope","gsigma")

checks.vec <- c("Overall","all.acf","all.geweke","all.gelman.rubin",
				as.vector(outer(c("max.abs.acf","abs.geweke","gelman.rubin"),vars.check, paste, sep="."))
				)


conv.summary <- data.frame(Check=checks.vec,
							Values = rep(NA,length(checks.vec)),
							Trigger = rep(NA,length(checks.vec)),
							Flag = rep(NA,length(checks.vec)),
							Perc.Over = rep(NA,length(checks.vec)),
							stringsAsFactors=FALSE)

# loop through variables
for(var.use in vars.check){
print(paste("starting",var.use))

# these use the MCMC samples
acf.out <- calcACF(mcmc.par=mcmc.out$BUGSoutput$sims.matrix[,var.use],  crit.acf.buffer = 1.5,acf.plot=FALSE)
geweke.out <- calcGeweke(mcmc.par=mcmc.out$BUGSoutput$sims.matrix[,var.use])

# uses the raw output (needs individual chains)
gelman.out <- gelman.diag(as.mcmc(mcmc.out)[,var.use])$psrf[2] # extracts the Upper 95 CI

conv.summary[conv.summary$Check==paste0("max.abs.acf.",var.use) ,c("Values","Trigger")] <- round(c(acf.out[[2]],acf.out[[3]]),4)
conv.summary[conv.summary$Check==paste0("abs.geweke.",var.use) ,c("Values","Trigger")] <- c(round(geweke.out,4),2.3)
conv.summary[conv.summary$Check==paste0("gelman.rubin.",var.use) ,c("Values","Trigger")] <- c(round(gelman.out,4),1.1)

}

conv.summary[,"Flag"] <- conv.summary[,"Values"]>=conv.summary[,"Trigger"]
#conv.summary[conv.summary$Check=="all.acf","Flag"] <- any(conv.summary[conv.summary$Check %in% c("max.abs.acf.int","max.abs.acf.slope","max.abs.acf.sigma"),"Flag"],na.rm=TRUE)
#conv.summary[conv.summary$Check=="all.geweke","Flag"] <- any(conv.summary[conv.summary$Check %in% c("abs.geweke.int","abs.geweke.slope","abs.geweke.sigma"),"Flag"],na.rm=TRUE)

conv.summary[conv.summary$Check=="all.acf","Flag"] <- any(conv.summary[grepl("max.abs.acf.",conv.summary$Check) ,"Flag"],na.rm=TRUE)
conv.summary[conv.summary$Check=="all.geweke","Flag"] <- any(conv.summary[grepl("abs.geweke.",conv.summary$Check) ,"Flag"],na.rm=TRUE)
conv.summary[conv.summary$Check=="all.gelman.rubin","Flag"] <- any(conv.summary[grepl("gelman.rubin.",conv.summary$Check) ,"Flag"],na.rm=TRUE)

conv.summary[conv.summary$Check=="Overall","Flag"] <- any(conv.summary[,"Flag"],na.rm=TRUE)

conv.summary[,"Perc.Over"] <- round(100*(as.numeric(conv.summary[,"Values"])-as.numeric(conv.summary[,"Trigger"]))/as.numeric(conv.summary[,"Trigger"]))


return(conv.summary)

}


# SUB-FUNCTION TO CALCULATE ACF AND CRITICAL VALUE
calcACF <- function(mcmc.par,  crit.acf.buffer = 1,acf.plot=FALSE, plot.title = "ACF Plot"){
# mcmc.par = mcmc sample for 1 parameter
# crit.acf.buffer =   describe
# acf.plot = if TRUE, produce the default plot from acf()
# plot.title = title for plot from acf()

# critical value calc based on http://www.squaregoldfish.co.uk/2010/01/20/r-the-acf-function-and-statistical-significance/
# checked resulting value against default R plot to verify
# R default acf plot uses 0.95, but that was usually cutting off at one or the other lag in the MCMC sample
# Changing it to .90 actually lowers the critical value -> Need to work throug underlying rationale for this calc at some future point
# for now, just adding an x% tolerance around it 

acf.val <- acf(mcmc.par, main=plot.title,plot=acf.plot)
acf.crit.val <- crit.acf.buffer * qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(mcmc.par)))

return(list(acf = acf.val, max.abs.val = max(abs(acf.val$acf[-1])), crit = acf.crit.val))


}

# SUB-FUNCTION TO CALCULATE GEWEKE SCORE
calcGeweke<- function(mcmc.par){
# mcmc.par = mcmc sample for 1 parameter

gew.val <- abs(try(geweke.diag(mcmc.par), silent=TRUE)[[1]])

names(gew.val) <- "AbsGewekeScore"

return(gew.val)

}








##################################################
# Function to calculate a generational average
# apply this to raw time series! (not logged and/or smoothed values!)
# NOT USED FOR NOW
wsp.avg <- function(vec,na.rm=TRUE){

	if(na.rm){vec<-na.omit(vec)}
	if(min(vec)==0){warning("Vector constains 0. log() results in -Inf"); stop()}
	out.val <- exp(mean(log(vec)))
	return(out.val)
}



##################################################
# PLOT OF TREND FIT


plot.trend.fit <- function(vec.in , mcmc.fit,n.plot = 1000){
# vec.in is the series of abd values used to fit the trend
# mcmc.fit is an object created by calling calcPercChangeMCMC()
#  with out.type = "long" or "full"
# n.plot is the number of subsamples mcmc pars to plot to illustrate the scatter

mcmc.have <- "samples" %in% names(mcmc.fit)

yrs <- 0:(length(vec.in)-1)

det.fit <- lm(vec.in ~ yrs) # simple regression

print(vec.in)
print(det.fit)

plot(yrs,vec.in, bty="n",type="o",xlab="Year",ylab= "Abd",col="darkblue",pch=19) # Data

# plot all the mcmc fitted lines (full posterior)
if(mcmc.have){
	for(i in sample(dim(mcmc.fit$samples)[1],n.plot,replace=FALSE)){
	abline(mcmc.fit$samples[i,"intercept"], mcmc.fit$samples[i,"slope"], col="lightgrey")
		}
	}

points(yrs,vec.in,type="o",col="darkblue",pch=19) # replot all the data points
abline(det.fit,col="darkblue",lwd=2) # simple regression - plot fitted line

# independent medians
if(mcmc.have){
abline(median(mcmc.fit$samples[,"intercept"]), median(mcmc.fit$samples[,"slope"]), 
					col="red",lwd=2,lty=2)
}


# did some testing: generally almost identical (in examples so far), but
# need to discuss
# median slope with matching intercept
#med.slope <- median(mcmc.samples[,"slope"])
#med.slope.idx <- (mcmc.samples[,"slope"] > (med.slope - 0.001))  &  (mcmc.samples[,"slope"] < (med.slope +0.001) )
#sum(med.slope.idx )
#abline(median(mcmc.samples[med.slope.idx,"intercept"]), median(mcmc.samples[med.slope.idx,"slope"]), col="red",lwd=2,lty=2)


legend("top",legend=c("Data","Det Fit","MCMC Samples", "Median MCMC"),
			pch= c(19,NA,NA,NA),lty=c(NA,1,1,2),lwd=c(NA,2,1,2),col=c("darkblue","darkblue"
						,"lightgrey","red"),bty="n",cex=0.8)




} # end plot.trend.fit







