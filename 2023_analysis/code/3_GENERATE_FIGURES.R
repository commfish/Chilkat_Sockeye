# THIS SCRIPT IS RUN AFTER 1_RUN_MODELS.R and 2b_GENERATE_OUTPUTS.r
# i and z act as ways to change range of escapement based on stock size

# input values below based on stats output
lowerB <- 70000  #lower bound of recommended escapement goal range
upperB <- 150000 #upper bound of recommended escapement goal range
SMSY <- 43857  #Lambert W from lambert file
UMSY <- 0.75  #median from staquants file
SMAX <- 59145  #median from staquants file
SEQ <- 124106 #median from staquants file
lnalpha.c <-  2.1115 #median from staquants file
beta <-1.69076E-05  #median from staquants file

# load----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(gsl)
library(scales)
library(dplyr)
library(extrafont)
library(grid)
library(FField)
library("devtools")
devtools::install_github("commfish/fngr")
library(fngr)
source('2023_analysis/code/functions.r')
#extrafont::font_import()


windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_report(base_size = 14))
if(!dir.exists(file.path("2023_analysis", "output", "rjags", "processed"))){dir.create(file.path("2023_analysis", "output", "rjags", "processed"))}

# data----
# loadfonts(device="win") #only need to do this once; takes awhile to run!
coda <- read.csv("2023_analysis/output/rjags/coda.csv") 
coda  %>%
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
         MSY.c = R.msy.c-S.msy.c, 
         Umsy = (1-lambert_W0(exp(1-lnalpha.c))),
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

# analysis----
# create function for probability profiles and figures
profile(i=10, z=500, xa.start=0, xa.end=700,lnalpha.c, beta) #can change i,z, xa.start, xa.end

p_q_Nya <- read.csv("2023_analysis/output/rjags/p_q_Nya.csv")
parameters <- read.csv("2023_analysis/output/rjags/parameters.csv")
read.csv("2023_analysis/data/chilkat_sockeye.csv") %>%
merge(., parameters, by=c("year"), all=TRUE)-> parameters
xaxis = tickr(parameters, year, 4)
read.csv("2023_analysis/output/rjags/CI.csv")-> CI

# escapement-DIDSON
options(scipen=999) 
parameters$S_val97.5pc <- as.numeric(parameters$S97.5.)
parameters$S_val2.5pc <- as.numeric(parameters$S2.5.)
parameters$S_median <- as.numeric(parameters$S50.)
maxY<-max(parameters$S97.5., na.rm=TRUE)*1.5
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
parameters$year<-as.numeric(parameters$year)
ggplot(parameters, aes(x=year, y=(S_median))) +
geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year") +
geom_line(aes(y=S_val2.5pc), colour="grey20", linetype="solid", size=0.1) +
  geom_ribbon(aes(ymin=(S_val2.5pc), ymax=(S_val97.5pc)), alpha=0.3, linetype="solid", size=0.5)+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black")) +
scale_y_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+
theme(legend.position = "none") +
geom_point(aes(y=DS), pch=21, size=3)-> plot1

# escapement-weir
parameters$qS_weir_val97.5pc<-as.numeric(parameters$qS.weir97.5.)
parameters$qS_weir_val2.5pc<-as.numeric(parameters$qS.weir2.5.)
parameters$qS_weir_median<-as.numeric(parameters$qS.weir50.)
maxY<-max(parameters$qS_weir_val97.5pc, na.rm=TRUE)*1.5
ggplot(parameters, aes(x=year, y=(qS_weir_median))) +
geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year") +
geom_line(aes(y=qS_weir_val2.5pc), colour="grey20", linetype="solid", size=0.1) +
  geom_ribbon(aes(ymin=(qS_weir_val2.5pc), ymax=(qS_weir_val97.5pc)), alpha=0.3, linetype="solid", size=0.5) +
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
scale_y_continuous(labels = comma,breaks = seq(0, 300000, 50000), limits = c(0, 300000))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+
theme(legend.position = "none") +
geom_point(aes(y=weir), pch=21, size=3)-> plot2

# escapement-mark-recapture
parameters$qS_mr_val97.5pc<-as.numeric(parameters$qS.mr97.5.)
parameters$qS_mr_val2.5pc<-as.numeric(parameters$qS.mr2.5.)
parameters$qS_mr_median<-as.numeric(parameters$qS.mr50.)
maxY<-max(parameters$qS_weir_val97.5pc, na.rm=TRUE)*1.5
ggplot(parameters, aes(x=year, y=(qS_mr_median))) +
geom_line(size=0.75)+ geom_point(size=2)+ylab("Escapement (S)")+xlab("Year") +
geom_line(aes(y=qS_mr_val2.5pc), colour="grey20", linetype="solid", size=0.1) +
  geom_ribbon(aes(ymin=(qS_mr_val2.5pc), ymax=(qS_mr_val97.5pc)), alpha=0.3, linetype="solid", size=0.5)+
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+                  
theme(legend.position = "none") +
geom_point(aes(y=mr), pch=21, size=3) -> plot3

png(file='2023_analysis/figures/escapement.png', res=200, width=8, height=9, units ="in") 
plot_grid(plot1, plot2, plot3, labels = c("A", "B", "C"), ncol = 1, align="v",hjust=-7,
          vjust=2, label_size=14)
dev.off()

# escapement
parameters$S_val97.5pc<-as.numeric(parameters$S97.5.)
parameters$S_val2.5pc<-as.numeric(parameters$S2.5.)
parameters$S_median<-as.numeric(parameters$S50.)
maxY<-max(parameters$S_val97.5pc, na.rm=TRUE)*1.5
ggplot(parameters, aes(x=year, y=(S_median))) +
geom_line(size=0.75)+ geom_point (size=2)+ylab("Escapement (S)")+xlab("Year") +
geom_line(aes(y=S_val2.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_ribbon(aes(ymin=(S_val2.5pc), ymax=(S_val97.5pc)), alpha=0.3, colour="grey70", linetype="dotted", size=0.5) +
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
scale_y_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+                      
theme(legend.position = "none") +
geom_line(aes(y=upperB), colour="grey70", size=1, linetype=2) +
geom_line(aes(y=lowerB), colour="grey70", size=1, linetype=2) + 
geom_line(aes(y=SMSY), colour="grey70", size=1, linetype=3) -> plot4

png(file='2023_analysis/figures/spawners.png', res=200, width=6, height=5, units ="in") 
plot_grid(plot4, ncol = 1, align="v",hjust=-7,
          vjust=2, label_size=14)
dev.off()

# returns
parameters$R_val97.5pc<-as.numeric(parameters$R97.5.)
parameters$R_val2.5pc<-as.numeric(parameters$R2.5.)
parameters$R_median<-as.numeric(parameters$R50.)
maxY<-max(parameters$R_val97.5pc, na.rm=TRUE)*1.5
ggplot(parameters, aes(x=year, y=(R_median)))+geom_line(size=0.75)+ geom_point (size=2)+ylab("Recruitment (R)")+xlab("Year") +
geom_line(aes(y=R_val2.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_line(aes(y=R_val97.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_ribbon(aes(ymin=R_val2.5pc, ymax=R_val97.5pc), alpha=0.3) +
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+
theme(legend.position = "none") -> fig_a

# abundance
parameters$N_val97.5pc<-as.numeric(parameters$N97.5.)
parameters$N_val2.5pc<-as.numeric(parameters$N2.5.)
parameters$N_median<-as.numeric(parameters$N50.)
maxY<-max(parameters$N_val97.5pc, na.rm=TRUE)*1.5
ggplot(parameters, aes(x=year, y=N_median))+geom_line(size=0.75)+ geom_point (size=2)+ylab("Total Run Abundance (N)")+xlab("Year") +
geom_line(aes(y=parameters$N_val2.5pc), colour="grey70", linetype="dotted", size=0.3) +
  geom_ribbon(aes(ymin=N_val2.5pc, ymax=N_val97.5pc), alpha=0.15) +
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+               
theme(legend.position = "none") -> fig_b

# Ricker productivity residuals
scaleFUN <- function(x) sprintf("%.2f", x)
parameters$log.resid_val97.5pc<-as.numeric(parameters$log.resid97.5.)
parameters$log.resid_val2.5pc<-as.numeric(parameters$log.resid2.5.)
parameters$log.resid_median<-as.numeric(parameters$log.resid50.)
maxY<-max(parameters$log.resid_val97.5pc, na.rm=TRUE)+0.5
minY<-min(parameters$log.resid_val2.5pc, na.rm=TRUE)-0.5
ggplot(parameters, aes(x=year, y=log.resid_median))+geom_line(size=0.75)+geom_point (size=2)+ ylab("Productivity Residuals")+xlab("Year") +
geom_line(aes(y=log.resid_val2.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_line(aes(y=log.resid_val97.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_line(aes(y=log.resid_val2.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_ribbon(aes(ymin=log.resid_val2.5pc, ymax=log.resid_val97.5pc), alpha=0.3) +
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
  geom_line(aes(y=0), colour="grey70", size=1, linetype=2) +
scale_y_continuous(breaks = seq(-2.5, 2, 0.5), limits = c(-2.5, 2)) +
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+                    
theme(legend.position = "none") -> fig_c

# harvest rate
parameters$mu.hbelow_val97.5pc<-as.numeric(parameters$mu.hbelow97.5.)
parameters$mu.hbelow_val2.5pc<-as.numeric(parameters$mu.hbelow2.5.)
parameters$mu.hbelow<-as.numeric(parameters$mu.hbelow50.)
ggplot(parameters, aes(x=year, y=parameters$mu.hbelow))+geom_line(size=0.75) + geom_point (size=2)+ ylab("Harvest Rate")+xlab("Year") +
geom_line(aes(y=mu.hbelow_val2.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_line(aes(y=mu.hbelow_val97.5pc), colour="grey70", linetype="dotted", size=0.5) +
  geom_ribbon(aes(ymin=mu.hbelow_val2.5pc, ymax=parameters$mu.hbelow_val97.5pc), alpha=0.3) +
  theme(panel.background = element_rect(colour="white")) +
  theme(panel.border = element_rect(colour = "black")) +
coord_cartesian(ylim=c(0,1))+
scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020))+                     
geom_line(aes(y=UMSY), colour="grey70", size=1, linetype=2) +
theme(legend.position = "none") -> fig_d
png(file='2023_analysis/figures/point_estimates.png', res=200, width=8, height=11, units ="in") 
plot_grid(fig_a, fig_b, fig_c, fig_d, labels = c("A", "B", "C","D"), ncol = 1, align="v",hjust=-7,
          vjust=2, label_size=14)
dev.off()

# proportions
# mean age at maturity (p), age composition (q), terminal run by age proportions (Nya)
p_q_Nya %>%
    mutate(age_comp = as.numeric(age_comp),
           year = as.numeric(year),
           Age = factor(age, ordered = TRUE, 
                        levels = c( "Ages 2-4", "Age 5", "Ages 6-8"),
                        labels = c("Ages 2-4", "Age 5", "Ages 6-8"))) %>%
    mutate(age_comp = ifelse(Age == 'Ages 2-4', NA, age_comp)) -> data
  
ggplot(data,aes(x=year, y=p, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) + scale_fill_grey(start=0.1, end=0.8) +
    theme(legend.title=element_blank(),legend.position=c(0.88,0.80)) +
    ylab("Age-at-Maturity Proportions") + xlab("") + 
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot1  

  
ggplot(data,aes(x=year, y=q, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) +
    scale_fill_grey(start=0.1, end=0.8) + 
    ylab("Age Composition Proportions") + xlab("") +
    theme(legend.title=element_blank(), legend.position="none") + geom_point(aes(x=year, y=age_comp), position='stack') +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot2
  
ggplot(data,aes(x=year, y=Nya, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) + scale_fill_grey(start=0.1, end=0.8) +
    ylab("Terminal Run by Age")+xlab("Year") + 
    guides(fill = guide_legend(reverse=TRUE)) + 
    theme(legend.title=element_blank(), legend.position="none") +
    scale_y_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot3

png(file='2023_analysis/figures/proportions.png', res=200, width=8, height=11, units ="in") 
plot_grid(plot1, plot2, plot3, labels = c("A", "B", "C"), ncol = 1, align="v",hjust=-7,
            vjust=2, label_size=14)
dev.off()

# horsetail (spawner recruit) plots
ggplot(data=CI, aes(x=Escapement, y=Median)) +
  geom_line(size=0.75, lty=2) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.08) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.08) +
  xlab("Spawners (S)") +
  ylab("Recruits (R)") +
  geom_vline(xintercept = SMSY, color ="gray70", lty=2)+
  scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000)) +
  geom_line(aes(x=Escapement, y=Escapement),linetype="solid", size=0.75, color ="gray80") +
  geom_point(data=parameters, aes(x=S_median, y=R_median),pch=1, size=2) +
  geom_text(size=3, data=parameters, aes(x=S_median, y=R_median, label=year,family="Times",
                                             hjust = -0.1, vjust= -0.4)) -> plot1
cowplot::plot_grid(plot1,  align = "v", nrow = 1, ncol=1)
ggsave("2023_analysis/figures/SR_curve.png", dpi = 500, height = 5, width = 8, units = "in")

# density plots
ggplot(coda, aes(x=S.msy.c, fill=Smsy, color = S.msy.c)) +
  geom_density(fill ="#999999", alpha=0.5) + 
  scale_color_manual(values=c("#999999")) +
  scale_fill_manual(values=c("#999999")) +
  geom_vline(xintercept = 43857,linetype = "longdash" ) +
  labs(x="Smsy",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                                    'Arial')+
                                           theme(panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank())) + theme(legend.position="none")+ 
  scale_x_continuous(labels = comma,breaks = seq(0, 200000, 25000), limits = c(0, 200000)) -> plot1

ggplot(coda, aes(x=Umsy, fill=Umsy)) +
  geom_density(fill ="#999999", alpha=0.5)+ 
  scale_color_manual(values=c("#999999"))+
  scale_fill_manual(values=c("#999999"))+geom_vline(xintercept = 0.75,linetype = "longdash" ) +
  labs(x="Umsy",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                                    'Arial')+
                                           theme(panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank())) + theme(legend.position="none")+ 
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) -> plot2

png(file='2023_analysis/figures/density.png', res=200, width=8, height=11, units ="in") 
plot_grid(plot1, plot2, labels = c("A", "B"), ncol = 1, align="v",hjust=-9,
          vjust=2, label_size=14)
dev.off()

