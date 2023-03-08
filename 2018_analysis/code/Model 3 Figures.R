#R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
#point estimate plots, yield plots)
#created by Ben Williams (Ben.Williams@alaska.gov);Nov 3, 2016
#Changes made by Sara Miller (Sara.Miller@alaska.gov); April 2017
#Step#1: Create three csv files from OpenBUGS or JAGS output
        #a)lnalpha, beta, and lnalpha.c called 'coda' with 1000+ values of each in columns with variable names
        #b)Parameters.csv
        #c)p_q_Nya.csv
#Step#2: Run code.
#i and z act as ways to change range of escapement based on stock size

#change to inlcude all inits????
rm(list=ls(all=T))#Remove previous variables.
LowerB<-70000 #lower bound of recommended escapement goal range
UpperB<-150000 #upper bound of recommended escapement goal range
SMSY<-98369.54#Lambert W
#load----
library(plyr)
library(reshape2)
library(extrafont)
library(tidyverse)
library(scales)
library(gsl)
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
#data----
#loadfonts(device="win") #only need to do this once; takes awhile to run!
coda <- read.csv("data/Coda.csv") #Load Data File
Parameters <- read.csv("data/Parameters.csv") #Load Data File
p_q_Nya <- read.csv("data/p_q_Nya.csv") #Load Data File

#data clean----
#Create profile parameters

coda %>% mutate(S.eq.c = lnalpha.c/beta, 
                S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
					      R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
					      MSY.c = R.msy.c-S.msy.c, 
					      Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

#analysis----
#create function for probability profiles and figures
profile <-function(i,z,xa.start, xa.end,lnalpha.c, beta){ 
xa = seq(xa.start, xa.end, by=i) 
x =(xa+i)*z
# empty dataframes
dat <- data.frame(S0=rep(1, length(coda[,1])))
dat1 <- data.frame(S0=rep(0, length(coda[,1])))
dat2 <- data.frame(S0=rep(0, length(coda[,1])))
dat3 <- data.frame(S0=rep(1, length(coda[,1])))
dat4 <- data.frame(S0=rep(0, length(coda[,1])))
dat5 <- data.frame(S0=rep(0, length(coda[,1])))
dat6 <- data.frame(S0=rep(1, length(coda[,1])))
dat7 <- data.frame(S0=rep(0, length(coda[,1])))
dat8 <- data.frame(S0=rep(0, length(coda[,1])))
dat9 <- data.frame(S0=rep(0, length(coda[,1])))
dat10 <- data.frame(S0=rep(0, length(coda[,1])))
for (i in 1:length(xa)){
  dat[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 0, ifelse(dat[,i]==0, 0,1))
  dat1[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 1,0)
  dat2[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.7*coda$Rmax), 1,0)
  dat3[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.8*coda$MSY.c), 0, ifelse(dat3[,i]==0, 0,1))
  dat4[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.8*coda$MSY.c), 1,0)
  dat5[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.8*coda$Rmax), 1,0)
  dat6[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.9*coda$MSY.c), 0, ifelse(dat6[,i]==0, 0,1))
  dat7[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.9*coda$MSY.c), 1,0)
  dat8[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.9*coda$Rmax), 1,0)
  dat9[,i+1] = x[i]*exp(coda$lnalpha.c-coda$beta*x[i])-x[i]
  dat10[,i+1] = x[i]*exp(coda$lnalpha.c-coda$beta*x[i])
  }
# Overfishing estimate ----
f.over <- function(x){
  x %>% 
    filter(complete.cases(.)) %>% 
    summarise_all(funs(mean)) %>% 
    gather() %>% 
    select(value)
}

of_0.7 <- f.over(dat)
of_0.8 <- f.over(dat3)
of_0.9 <- f.over(dat6)

# Optimal yield estimate ----
oy_0.7 <- f.over(dat1)
oy_0.8 <- f.over(dat4)
oy_0.9 <- f.over(dat7)

# Optimal recruitment ----
or_0.7 <- f.over(dat2)
or_0.8 <- f.over(dat5)
or_0.9 <- f.over(dat8)

#Bind dataframes together
Y <- cbind(of_0.7,oy_0.7,or_0.7,of_0.8,oy_0.8,or_0.8,of_0.9,oy_0.9,or_0.9, c(0, x))
names(Y) <- c('of_0.7','oy_0.7','or_0.7','of_0.8','oy_0.8','or_0.8','of_0.9','oy_0.9',
              'or_0.9','Escapement')

#Quantiles and Medians ----
summarise_all(dat9, funs(median, 
                         q95=quantile(., 0.95, na.rm=T), 
                         q90=quantile(., 0.90, na.rm=T),
                         q10=quantile(., 0.10, na.rm=T),
                         q5=quantile(., 0.05, na.rm=T))) -> mq
names(mq) <- c(rep(('Median'),length(x)+1), 
               rep(('q95'),length(x)+1), 
               rep(('q90'),length(x)+1), 
               rep(('q10'),length(x)+1), 
               rep(('q5'),length(x)+1))

qm <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
qm <- spread(qm, measure, value)
qm <- qm[c("q95", "q90", "Median","q10", "q5", "Escapement")]
Y <- Y[c("oy_0.9", "oy_0.8", "or_0.9","or_0.8", "of_0.9", "of_0.8", "oy_0.7","or_0.7","of_0.7","Escapement")]
write.csv(qm,("data/processed/QM.csv"), row.names=FALSE)
write.csv(Y,("data/processed/Y.csv"), row.names=FALSE)

#confidence intervals ----
summarise_all(dat10, funs(median, 
                         q95=quantile(., 0.95, na.rm=T), 
                         q90=quantile(., 0.90, na.rm=T),
                         q10=quantile(., 0.10, na.rm=T),
                         q5=quantile(., 0.05, na.rm=T))) -> mq
names(mq) <- c(rep(('Median'),length(x)+1), 
               rep(('q95'),length(x)+1), 
               rep(('q90'),length(x)+1), 
               rep(('q10'),length(x)+1), 
               rep(('q5'),length(x)+1))

CI <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
CI <- spread(CI, measure, value)
CI <- CI[c("q95", "q90", "Median","q10", "q5", "Escapement")]
write.csv(CI,("data/processed/CI.csv"), row.names=FALSE)

#create probability profile plots (0.7, 0.8, 0.9, 0.8 & 0.9)
Y %>% 
  dplyr::select(Escapement, oy_0.7, of_0.7,or_0.7) %>% 
  melt(., id.vars = 'Escapement') %>% 
  ggplot( aes(Escapement/1000, value, lty=variable))+geom_line()+
  xlab('Escapement (1,000)')+ylab('Probability')+
  theme(legend.justification=c(1,0), legend.position=c(1,.5), 
        legend.key = element_blank(),legend.title=element_blank())
ggsave("figures/Model 3/0.7.AR.png", dpi=200, width=8, height=5, units='in')

Y %>% 
  dplyr::select(Escapement, oy_0.8, of_0.8, or_0.8) %>% 
  melt(., id.vars = 'Escapement') %>% 
  ggplot(aes(Escapement/1000, value, lty=variable))+geom_line()+
  xlab('Escapement (1,000)')+ylab('Probability')+
  theme(legend.justification=c(1,0), legend.position=c(1,.5), 
        legend.key = element_blank(),legend.title=element_blank())
ggsave("figures/Model 3/0.8.AR.png", dpi=200, width=8, height=5, units='in')

Y %>% 
  dplyr::select(Escapement, oy_0.9, of_0.9, or_0.9) %>% 
  melt(., id.vars = 'Escapement')  %>% 
  ggplot(aes(Escapement/1000, value, lty=variable))+geom_line()+
  xlab('Escapement (1,000)')+ylab('Probability')+
  theme(legend.justification=c(1,0), legend.position=c(1,.5), 
        legend.key = element_blank(),legend.title=element_blank())
ggsave("figures/Model 3/0.9.AR.png", dpi=200, width=8, height=5, units='in')


Y <- read.csv("data/processed/Y.csv")
Y["OY0.9"] <-Y$oy_0.9 
Y["OY0.8"] <-Y$oy_0.8 
Y<-subset(Y, select=c(Escapement, OY0.9,OY0.8))
mY1 <- melt(Y, id.vars='Escapement')
mY1["sra"] <-"Yield Profile"
mY1["max_pct"] <- ifelse(grepl("OY0.8",mY1$variable), 
                         0.8,0.9)

Y <- read.csv("data/processed/Y.csv")
Y["OF0.9"] <-Y$of_0.9 
Y["OF0.8"] <-Y$of_0.8 
Y<-subset(Y, select=c(Escapement, OF0.9,OF0.8))
mY2 <- melt(Y, id.vars='Escapement')
mY2["sra"] <-"Overfishing Profile"
mY2["max_pct"] <- ifelse(grepl("OF0.8",mY2$variable), 
                         0.8,0.9)

Y <- read.csv("data/processed/Y.csv")
Y["OR0.9"] <-Y$or_0.9 
Y["OR0.8"] <-Y$or_0.8 
Y<-subset(Y, select=c(Escapement, OR0.9,OR0.8))
mY3 <- melt(Y, id.vars='Escapement')
mY3["sra"] <-"Recruitment Profile"
mY3["max_pct"] <- ifelse(grepl("OR0.8",mY3$variable), 
                         0.8,0.9)
mY4<-rbind(mY1,mY2, mY3)
mY4<-subset(mY4, select=c(Escapement, value, sra, max_pct))
mY4$Escapement<-as.numeric(mY4$Escapement)
mY4$value<-as.numeric(mY4$value)
mY4$max_pct<-as.factor(mY4$max_pct)
colnames(mY4)[2] <- "Probability"
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Fig1<-ggplot(mY4, aes(x = Escapement, y = Probability, linetype = max_pct)) 
Fig1<-Fig1+geom_rect(aes(xmin = LowerB, xmax = UpperB, ymin = 0, ymax = 1),
                     inherit.aes = FALSE, fill = "grey80", alpha = 0.3)+geom_line()+xlab('Escapement (S)')+
  scale_x_continuous(labels = comma, breaks = seq(0, 350000, 50000), limits = c(0, 350000))+
  scale_linetype_discrete(name = "Percent of Max.")+
  facet_grid(sra ~ .) +geom_vline(xintercept=SMSY, lwd=1.25)+theme(legend.key = element_blank())
 Fig1<-Fig1+ theme_bw()+theme(legend.title = element_text(size=10),
                 legend.justification=c(0,1), 
                 legend.position=c(0.7,0.5),  
                 legend.background = element_blank(),
                 legend.key = element_blank(),text=element_text(family="Times New Roman"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 options(scipen=99999)
ggsave("figures/Model 3/0.8_0.9.png", dpi=200, dev='png', width=7, height=6, units='in')


ggplot(qm, aes(Escapement, Median))+geom_line(size=1)+
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.15)+
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.15)+ xlab('Escapement (S)')+
  ylab('Expected Yield')+scale_y_continuous(labels = comma)+
  scale_x_continuous(labels = comma,breaks = seq(0, 300000, 50000), limits = c(0,300000))+
  geom_vline(xintercept = LowerB,linetype = "longdash" )+geom_vline(xintercept = UpperB ,linetype = "longdash")
ggsave("figures/Model 3/expected_sustained_yield.png", dpi=200, width=8, height=5, units='in')
}
#Run function
profile(i=10,z=500,xa.start=0, xa.end=700,lnalpha.c, beta)#can change i,z, xa.start, xa.end

################################################################################################
#Escapement, Returns, Run abundance, Residuals, and Harvest Rate by year
#Figure x.- Point estimates (posterior medians; solid lines) and 95% credibility intervals 
#(bracketed by dashed lines) of (a) spawning escapement, (b) return by brood year, (c) run abundance, 
#(d) productivity residuals, and (e) harvest rate, from state-space spawner-recruit model of Chilkat Lake  
#sockeye salmon, 1976-2015.  Posterior medians of optimal escapement SMSY and UMSY are 
#plotted as horizontal reference lines in (a) and (e).
################################################################################################
library(cowplot)
library(grid)
Parameters <- read.csv("data/Parameters.csv") #Load Data File
#Escapement-DIDSON
options(scipen=999) 
Parameters$S_val97.5pc<-as.numeric(Parameters$S_val97.5pc)
Parameters$S_median<-as.numeric(Parameters$S_median)
maxY<-max(Parameters$S_val97.5pc, na.rm=TRUE)*1.5
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Parameters$Year<-as.numeric(Parameters$Year)
Fig_a<-ggplot(Parameters, aes(x=Parameters$Year, y=(Parameters$S_median)))
Fig_a<-Fig_a+geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year")
Fig_a<-Fig_a+geom_line(aes(y=S_val2.5pc), colour="grey20", linetype="solid", size=0.1)+
  geom_ribbon(aes(ymin=(Parameters$S_val2.5pc), ymax=(Parameters$S_val97.5pc), alpha=0.001), linetype="solid", size=0.5)+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_a<-Fig_a+scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+scale_x_continuous(breaks=seq(min(Parameters$Year-5, na.rm=TRUE),
                                                                                                                                    max(Parameters$Year+1, na.rm=TRUE),5))                                              
Fig_a<- Fig_a + theme(legend.position = "none")
Fig_a<-Fig_a+geom_point(aes(y=Parameters$DIDSON), pch=8, size=2.5)

#Escapement-Weir
options(scipen=999) 
Parameters$qS_weir_val97.5pc<-as.numeric(Parameters$qS_weir_val97.5pc)
Parameters$qS_weir_median<-as.numeric(Parameters$qS_weir_median)
maxY<-max(Parameters$qS_weir_val97.5pc, na.rm=TRUE)*1.5
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Parameters$Year<-as.numeric(Parameters$Year)
Fig_b<-ggplot(Parameters, aes(x=Parameters$Year, y=(Parameters$qS_weir_median)))
Fig_b<-Fig_b+geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year")
Fig_b<-Fig_b+geom_line(aes(y=qS_weir_val2.5pc), colour="grey20", linetype="solid", size=0.1)+
  geom_ribbon(aes(ymin=(Parameters$qS_weir_val2.5pc), ymax=(Parameters$S_val97.5pc), alpha=0.001), linetype="solid", size=0.5)+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_b<-Fig_b+scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+scale_x_continuous(breaks=seq(min(Parameters$Year-5, na.rm=TRUE),
                                                                                                                                    max(Parameters$Year+1, na.rm=TRUE),5))                       
Fig_b<- Fig_b + theme(legend.position = "none")
Fig_b<-Fig_b+geom_point(aes(y=Parameters$weir), pch=8, size=2.5)


#Escapement-mark-recapture
options(scipen=999) 
Parameters$qS_mr_val97.5pc<-as.numeric(Parameters$qS_mr_val97.5pc)
Parameters$qS_mr_median<-as.numeric(Parameters$qS_mr_median)
maxY<-max(Parameters$qS_weir_val97.5pc, na.rm=TRUE)*1.5
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Parameters$Year<-as.numeric(Parameters$Year)
Fig_c<-ggplot(Parameters, aes(x=Parameters$Year, y=(Parameters$qS_mr_median)))
Fig_c<-Fig_c+geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year")
Fig_c<-Fig_c+geom_line(aes(y=qS_mr_val2.5pc), colour="grey20", linetype="solid", size=0.1)+
  geom_ribbon(aes(ymin=(Parameters$qS_mr_val2.5pc), ymax=(Parameters$qS_mr_val97.5pc), alpha=0.001), linetype="solid", size=0.5)+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_c<-Fig_c+scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+scale_x_continuous(breaks=seq(min(Parameters$Year-5, na.rm=TRUE),
                                                                                                                                    max(Parameters$Year+1, na.rm=TRUE),5))                       
Fig_c<- Fig_c + theme(legend.position = "none")
Fig_c<-Fig_c+geom_point(aes(y=Parameters$mr), pch=8, size=2.5)

png(file='figures/Model 3/Escapement.png', res=200, width=8, height=9, units ="in") 
plot_grid(Fig_a, Fig_b,Fig_c, labels = c("A", "B", "C"), ncol = 1, align="v",hjust=-7,
          vjust=2, label_size=14)
dev.off()

#Escapement
options(scipen=999) 
Parameters$S_val97.5pc<-as.numeric(Parameters$S_val97.5pc)
Parameters$S_median<-as.numeric(Parameters$S_median)
maxY<-max(Parameters$S_val97.5pc, na.rm=TRUE)*1.5
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Parameters <- read.csv("data/Parameters.csv") #Load Data File
Parameters$Year<-as.numeric(Parameters$Year)
Fig_a<-ggplot(Parameters, aes(x=Parameters$Year, y=(Parameters$S_median)))
Fig_a<-Fig_a+geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year")
Fig_a<-Fig_a+geom_line(aes(y=S_val2.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_ribbon(aes(ymin=(Parameters$S_val2.5pc), ymax=(Parameters$S_val97.5pc), alpha=0.15), colour="grey70", linetype="dotted", size=0.5)+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_a<-Fig_a+scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+scale_x_continuous(breaks=seq(min(Parameters$Year, na.rm=TRUE),
                                                                                                                                    max(Parameters$Year+1, na.rm=TRUE),5))                       
Fig_a<- Fig_a + theme(legend.position = "none")
Fig_a<-Fig_a+geom_line(aes(y=Parameters$SMSY), colour="grey70", size=1, linetype=2)

#Returns
maxY<-max(Parameters$R_val97.5pc, na.rm=TRUE)*1.5
Fig_b<-ggplot(Parameters, aes(x=Parameters$Year, y=(Parameters$R_median)))+geom_line(size=0.75)+ geom_point (size=2)+ylab("Recruitment (R)")+xlab("Year")
Fig_b<-Fig_b+geom_line(aes(y=R_val2.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_line(aes(y=Parameters$R_val97.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_ribbon(aes(ymin=Parameters$R_val2.5pc, ymax=R_val97.5pc, alpha=0.15))+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_b<-Fig_b+scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+scale_x_continuous(breaks=seq(min(Parameters$Year, na.rm=TRUE),
max(Parameters$Year+1, na.rm=TRUE),5)) 
Fig_b<- Fig_b + theme(legend.position = "none")
                                                                                                                                    
maxY<-max(Parameters$N_val97.5pc, na.rm=TRUE)*1.5
Fig_c<-ggplot(Parameters, aes(x=Parameters$Year, y=Parameters$N_median))+geom_line(size=0.75)+ geom_point (size=2)+ylab("Total Run Abundance (N)")+xlab("Year")
Fig_c<-Fig_c+geom_line(aes(y=Parameters$N_val2.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_ribbon(aes(ymin=Parameters$N_val2.5pc, ymax=N_val97.5pc, alpha=0.15))+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_c<-Fig_c+scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+scale_x_continuous(breaks=seq(min(Parameters$Year, na.rm=TRUE),
max(Parameters$Year+1, na.rm=TRUE),5))                  
Fig_c<- Fig_c + theme(legend.position = "none")  

#Ricker Productivity Residuals
scaleFUN <- function(x) sprintf("%.2f", x)
maxY<-max(Parameters$log.resid_val97.5pc, na.rm=TRUE)+0.5
minY<-min(Parameters$log.resid_val2.5pc, na.rm=TRUE)-0.5
Fig_d<-ggplot(Parameters, aes(x=Parameters$Year, y=Parameters$log.resid_median))+geom_line(size=0.75)+geom_point (size=2)+ ylab("Productivity Residuals")+xlab("Year")
Fig_d<-Fig_d+geom_line(aes(y=Parameters$log.resid_val2.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_line(aes(y=Parameters$log.resid_val97.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_ribbon(aes(ymin=Parameters$log.resid_val2.5pc, ymax=Parameters$log.resid_val97.5pc, alpha=0.15))+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_d<-Fig_d+geom_line(aes(y=0), colour="black", size=0.5)
Fig_d<-Fig_d+scale_y_continuous(breaks = seq(-2, 2, 0.5), limits = c(-2, 2))+scale_x_continuous(breaks=seq(min(Parameters$Year, na.rm=TRUE),
      max(Parameters$Year+1, na.rm=TRUE),5))                       
Fig_d<- Fig_d + theme(legend.position = "none")

#Harvest Rate
Fig_e<-ggplot(Parameters, aes(x=Parameters$Year, y=Parameters$mu.HB_median))+geom_line(size=0.75)+geom_point (size=2)+ ylab("Harvest Rate")+xlab("Year")
Fig_e<-Fig_e+geom_line(aes(y=Parameters$mu.HB_val2.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_line(aes(y=Parameters$mu.HB_val97.5pc), colour="grey70", linetype="dotted", size=0.5)+
  geom_ribbon(aes(ymin=Parameters$mu.HB_val2.5pc, ymax=Parameters$mu.HB_val97.5pc, alpha=0.15))+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_e<-Fig_e+coord_cartesian(ylim=c(0,1))+scale_x_continuous(breaks=seq(min(Parameters$Year, na.rm=TRUE),
        max(Parameters$Year+1, na.rm=TRUE),5))                       
Fig_e<-Fig_e+geom_line(aes(y=Parameters$UMSY), colour="grey70", size=1, linetype=2)
Fig_e<- Fig_e + theme(legend.position = "none")
png(file='figures/Model 3/Point Estimates.png', res=200, width=8, height=11, units ="in") 
plot_grid(Fig_a, Fig_b,Fig_c,Fig_d, Fig_e, labels = c("A", "B", "C","D", "E"), ncol = 1, align="v",hjust=-7,
          vjust=2, label_size=14)
dev.off()


################################################################################################
#Mean Age at Maturity, Age Composition, Total Run by Age Proportions
#Figure x.  Area graph of  mean age-at-maturity proportions (p) of Chilkat Lake sockeye salmon 
#by brood years 1970-2011. Spaces between the solid lines are posterior medians of proportions 
#of fish returning by age for the cohort specified on the horizontal axis. Horizontal lines are 
#posterior medians of age-at-maturity central tendency proportions pa.
################################################################################################
library(cowplot)
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

p_q_Nya <- read.csv("data/p_q_Nya.csv") #Load Data File
Parameters<-subset(p_q_Nya, select = c(Year,Age, p, q, age_comp, N.ya))
newdata<-Parameters[order(Parameters$Year, Parameters$Age),]
newdata$Year<-as.numeric(newdata$Year)
newdata$age_comp<-as.numeric(newdata$age_comp)
newdata['variable']<-paste("Age",newdata$Age)
#newdata['variable'] <- ifelse(newdata$Age!='5',"Ages","Age")
#newdata['variable']<-paste(newdata$variable,newdata$Age)
newdata$Year<-as.numeric(newdata$Year)
Fig1<-ggplot(newdata,aes(x=newdata$Year, y=newdata$p, fill=factor(newdata$variable))) +
  geom_area(position=position_stack(reverse=TRUE))+scale_fill_grey(start=0.2, end=0.8)
Fig1<-Fig1 + theme(legend.title=element_blank())
Fig1<-Fig1 +  ylab("Age-at-Maturity Proportions")+xlab("Year")
Fig1<-Fig1+scale_x_continuous(breaks = seq(1970, 2020, 5), limits = c(1970, 2020))    

Fig2<-ggplot(newdata,aes(x=newdata$Year, y=newdata$q, fill=factor(newdata$variable))) +
  geom_area(position=position_stack(reverse=TRUE))+scale_fill_grey(start=0.2, end=0.8)
Fig2<-Fig2 +  ylab("Age Composition Proportions")+xlab("Year")
#Fig2<-Fig2+geom_point(data=subset(newdata,newdata$variable!='Age 6'),mapping=aes(x=Year, y=age_comp,fill=factor(variable)),
#                      position="stack", color="black", size=1.6)
Fig2<-Fig2 + theme(legend.title=element_blank())
Fig2<-Fig2+scale_x_continuous(breaks = seq(1970, 2020, 5), limits = c(1970, 2020))

Fig3<-ggplot(newdata,aes(x=newdata$Year, y=newdata$N.ya, fill=factor(newdata$variable))) +
  geom_area(position=position_stack(reverse=TRUE))+scale_fill_grey(start=0.2, end=0.8)
Fig3<-Fig3 +  ylab("Total Run by Age")+xlab("Year")
Fig3<-Fig3 + guides(fill = guide_legend(reverse=FALSE)) 
Fig3<-Fig3 + theme(legend.title=element_blank())
Fig3<-Fig3+scale_y_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000))
Fig3<-Fig3+scale_x_continuous(breaks = seq(1970, 2020, 5), limits = c(1970, 2020))

png(file='figures/Model 3/Proportions.png', res=200, width=8, height=11, units ="in") 
plot_grid(Fig1, Fig2, Fig3, labels = c("A", "B", "C"), ncol = 1, align="v",hjust=-7,
          vjust=2, label_size=14)
dev.off()
################################################################################################
#Horesetail Plots
#Figure x.- Graphical summary of knowledge of spawner-recruitment relationship for Chilkat Lake 
#sockeye salmon as derived from age-structured state space model fitted to abundance, harvest, 
#and age data 1976-2015. Symbols are posterior medians of R and S; error bars bracket 95% credibility 
#intervals. Heavy dashed line is Ricker relationship constructed from ln(a) and b posterior medians,  
#Ricker relationships are also plotted for 50 paired values of ln(a) and b sampled from the posterior 
#probability distribution, representing plausible Ricker relationships that could have generated the 
#observed data.  Diagonal line is replacement line (R=S).
################################################################################################
coda <- read.csv("data/Coda.csv") 
Parameters <- read.csv("data/Parameters.csv")
QM <- read.csv("data/processed/QM.csv")
CI<- read.csv("data/processed/CI.csv")
coda <- subset(coda, select = c(lnalpha.c, beta))
coda<-as.data.frame(coda[1:50,])#select first 50 rows of dataframe
QM<-subset(QM, select = c(Escapement))
num<-nrow(QM)
coda<-data.frame(lnalpha.c=rep(coda$lnalpha.c,each=num), beta=rep(coda$beta, each=num))
dataset<-cbind(coda,QM) #lnalpha.c, beta, and S
dataset['Recruitment']<-dataset$Escapement*exp(dataset$lnalpha.c-dataset$beta*dataset$Escapement)
dataset['Variable']<-data.frame(dataset=rep(1:50,each=num))

myvars <- c("lnalpha.c", "beta")
Parameters<-Parameters[myvars]
Parameters<-data.frame(lnalpha.c=rep(Parameters$lnalpha.c,each=num), beta=rep(Parameters$beta, each=num))
row.has.na <- apply(Parameters, 1, function(x){any(is.na(x))})
sum(row.has.na)
final.filtered <- Parameters[!row.has.na,]
final.filtered <-cbind(final.filtered,QM)
final.filtered['Recruitment']<-final.filtered$Escapement*exp(final.filtered$lnalpha.c-final.filtered$beta*final.filtered$Escapement)
final.filtered['Variable']<-data.frame(final.filtered=rep(51,each=num))
dataset<-rbind(dataset, final.filtered)
dataset['Year']<-'NA'
dataset['R_val2.5pc']<-'0'
dataset['R_val97.5pc']<-'0'
dataset['S_val2.5pc']<-'0'
dataset['S_val97.5pc']<-'0'
Parameters <- read.csv('./data/Parameters.csv') #Load Data File
Parameters<-subset(Parameters, Parameters$Year>1975) 
Parameters<-subset(Parameters, Parameters$Year<2013) 
Parameters<-subset(Parameters, select = c(S_median, R_median,Year,R_val2.5pc,R_val97.5pc,S_val2.5pc,S_val97.5pc))
Parameters['lnalpha.c']<-'NA'
Parameters['beta']<-'NA'
Parameters['Escapement']<-Parameters$S_median
Parameters['Recruitment']<-Parameters$R_median
Parameters['Variable']<-52
Parameters<-subset(Parameters, select=c(lnalpha.c, beta, Escapement, Recruitment, Variable,Year,R_val2.5pc,R_val97.5pc,S_val2.5pc,S_val97.5pc))
dataset<-rbind(dataset, Parameters)
dataset$R_val2.5pc<-as.numeric(dataset$R_val2.5pc)
dataset$R_val97.5pc<-as.numeric(dataset$R_val97.5pc)
dataset$S_val2.5pc<-as.numeric(dataset$S_val2.5pc)
dataset$S_val97.5pc<-as.numeric(dataset$S_val97.5pc)
dataset <- dataset[order(dataset$Variable, dataset$Escapement),] 
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Fig<-ggplot(data=dataset, aes(x=Escapement, y=Recruitment, group=Variable))+
  geom_line(data=subset(dataset,dataset$Variable<52),linetype="solid", size=0.5, color="grey80")+
  scale_y_continuous(labels = comma,breaks = seq(0, 650000, 50000), limits = c(0, 650000))+
  scale_x_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000))+ylab("Recruits (R)")+xlab("Spawners (S)")+
  geom_line(data=dataset, aes(x=Escapement, y=Escapement, group=1),linetype="solid", size=1)+#replacement line
  geom_line(data=subset(dataset,dataset$Variable==51),colour = "black", lty=2, size=2)+
  geom_text(data=subset(dataset,dataset$Variable==52), aes(x=Escapement, y=Recruitment, label=Year,family="Times"))+
  geom_errorbar(aes(ymax = R_val97.5pc, ymin=R_val2.5pc), width=0.20,linetype = 2, colour="grey50")+
  geom_errorbarh(aes(xmax = S_val97.5pc, xmin=S_val2.5pc), height=0.20,linetype = 2,colour="grey50")
png(file='figures/Model 3/Horsetail_Plot.png', res=200, width=6, height=4, units ="in")  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(Fig,vp=vplayout(1,1:1)) 
dev.off()

dataset <- subset(dataset,Variable == 52) 
dataset['Escapement1']<-dataset$Escapement
dataset<-subset(dataset, select = -c(lnalpha.c, beta, Escapement))
dataset['Escapement']<-'NA'
CI['Year']<-'NA'
CI['R_val2.5pc']<-'NA'
CI['R_val97.5pc']<-'NA'
CI['S_val2.5pc']<-'NA'
CI['S_val97.5pc']<-'NA'
CI['Variable']<-51
CI['Recruitment']<-'NA'
CI['Escapement1']<-'NA'
dataset['Median']<-'NA'
dataset['q95']<-'NA'
dataset['q90']<-'NA'
dataset['q10']<-'NA'
dataset['q5']<-'NA'
dataset1<-rbind(dataset, CI)
dataset1$Escapement<-as.numeric(dataset1$Escapement)
dataset1$Median<-as.numeric(dataset1$Median)
dataset1$q5<-as.numeric(dataset1$q5)
dataset1$q95<-as.numeric(dataset1$q95)
dataset1$q10<-as.numeric(dataset1$q10)
dataset1$q90<-as.numeric(dataset1$q90)
dataset1$Recruitment<-as.numeric(dataset1$Recruitment)
dataset1$Escapement1<-as.numeric(dataset1$Escapement1)
dataset1$R_val2.5pc<-as.numeric(dataset1$R_val2.5pc)
dataset1$R_val97.5pc<-as.numeric(dataset1$R_val97.5pc)
dataset1$S_val2.5pc<-as.numeric(dataset1$S_val2.5pc)
dataset1$S_val97.5pc<-as.numeric(dataset1$S_val97.5pc)

Fig1<-ggplot(data=dataset1, aes(x=Escapement, y=Median, group=Variable))+geom_line(size=2, lty=2, group=51)+
  geom_ribbon(aes(ymin = q5, ymax = q95, group=51), alpha=.15)+
  geom_ribbon(aes(ymin = q10, ymax = q90, group=51), alpha=.15)+
  xlab('Spawners (S)')+
  ylab('Recruits (R)')+scale_y_continuous(labels = comma)+
  scale_x_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0,350000))+
  geom_line(aes(x=Escapement, y=Escapement, group=51),linetype="solid", size=1)+
  geom_text(size=3, data=dataset1, aes(x=Escapement1, y=Recruitment, group=52, label=Year,family="Times", 
                                       hjust = -0.1, vjust= -0.4))
Fig1<-Fig1+
  geom_point(data=dataset1, aes(x=Escapement1, y=Recruitment, group=52),pch=16, size=1)+
  geom_errorbar(data=dataset1, aes(x=Escapement1, ymax=R_val97.5pc, ymin=R_val2.5pc, group=52), width=0.2,linetype = 1, colour="grey50")
Fig1<-Fig1+
  geom_point(data=dataset1, aes(x=Escapement1, y=Recruitment, group=52),pch=16, size=1)+
  geom_errorbarh(data=dataset1, aes(x=Escapement1, y=Recruitment, xmax=S_val97.5pc, xmin=S_val2.5pc),na.rm=T,  linetype = 1, colour="grey50")

png(file='figures/Model 3/Horsetail_Plot_Reconfig.png', res=200, width=6, height=4, units ="in")  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(Fig1,vp=vplayout(1,1:1)) 
dev.off()

################################################################################################
#ESCAPEMENT ESTIMATES WITH REFERENCE LINES
################################################################################################
Parameters <- read.csv("data/Parameters.csv") #Load Data File
options(scipen=999) 
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Fig_a<-ggplot(Parameters, aes(x=Year, y=S_median))+geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year")
Fig_a<-Fig_a+theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))
Fig_a<-Fig_a+scale_y_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000))+scale_x_continuous(breaks=seq(min(Parameters$Year, na.rm=TRUE),
       max(Parameters$Year+4, na.rm=TRUE),5))                       
Fig_a<- Fig_a + theme(legend.position = "topright")
Fig_a<-Fig_a+geom_line(aes(y=SMAX), colour="grey40", size=1, linetype=2)+
geom_errorbar(aes(ymin=S_val2.5pc, ymax=S_val97.5pc),size=0.5, linetype = "solid",colour="black", width=0.02)
Fig_a<-Fig_a+geom_line(aes(y=SEQ), colour="grey50", size=1, linetype=1)+
annotate("rect", xmin = min(Parameters$Year), xmax = max(Parameters$Year)+1, ymin = LowerB, ymax = UpperB,
         alpha = .3)
png(file='figures/Model 3/Historical Escapement Plot.png', res=200, width=6, height=4, units ="in")  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(Fig_a,vp=vplayout(1,1:1)) 
dev.off()
################################################################################################
#AGE COMPOSITION RESIDUALS-TEST FIGURE
################################################################################################
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

p_q_Nya <- read.csv("data/p_q_Nya.csv") #Load Data File
p_q_Nya ['residual']<-p_q_Nya$age_comp-p_q_Nya$q
Fig<-ggplot(data=p_q_Nya, aes(x=Year, y=q, size=abs(p_q_Nya$residual))) +
  geom_point(aes(shape=factor(sign(residual)), colour = factor(residual)))+ 
  scale_shape_discrete(name="Mean Is", breaks=c(-1, 1), labels=c("Negative", "Positive"))+
theme(legend.position = "none")

png(file='figures/Model 3/test.png', res=200, width=6, height=4, units ="in")  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(Fig,vp=vplayout(1,1:1)) 
dev.off()

################################################################################################
#Escapement Goal Ranges Across State Comparison Figure (after Fleishmand and Reimer 2017;Appendix D1)
################################################################################################
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

Data <- read.csv("data/EG_Sockeye.csv") #Load Data File
Data$UB<-as.numeric(Data$UB)
Data$LB<-as.numeric(Data$LB)
Fig<-ggplot() + 
  geom_segment(data=Data, mapping=aes(x=LB, y=System, xend=UB, yend=System))+
  xlab(bquote(paste("Multiples of ", italic("S"[MSY]))))+
  geom_vline(xintercept=1, lwd=0.75, linetype=1, color="grey50")+
  scale_x_continuous(breaks = seq(0, 2, 0.25), limits = c(0,2))
      
png(file='figures/Model 3/Goal_Comparison.png', res=200, width=6, height=4, units ="in")  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(Fig,vp=vplayout(1,1:1)) 
dev.off()