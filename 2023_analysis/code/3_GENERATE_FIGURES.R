# THIS SCRIPT IS RUN AFTER 1_RUN_MODELS.R and 2b_GENERATE_OUTPUTS.r
# i and z act as ways to change range of escapement based on stock size

# input values below based on stats output
LowerB <- 49500  #lower bound of recommended escapement goal range
UpperB <- 55700 #upper bound of recommended escapement goal range
SMSY <- 43857  #Lambert W from lambert file
UMSY <- 0.75  #median from staquants file
SMAX <- 59145  #median from staquants file
SEQ <- 124106 #median from staquants file
lnalpha.c <-  2.1115 #median from staquants file
beta <-1.69076E-05  #median from staquants file
SGEN<-5873

# load----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(gsl)
library(FField)
library(scales)
library(dplyr)
library(extrafont)
library("devtools")
devtools::install_github("commfish/fngr")
library(fngr)
source('state_space_model/code/functions.r')
extrafont::font_import()

# read in the custom dictionary (english to french translation)
terms.use <- read.csv("state_space_model/data/translation_terms.csv",stringsAsFactors = FALSE)

# test the translation function
translate("Overfishing Profile", from = "english", to = "french", 
          terms = terms.use, allow_missing = FALSE)

windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_report(base_size = 14))
if(!dir.exists(file.path("state_space_model", "output", "rjags_Explore_Basecase", "processed"))){dir.create(file.path("state_space_model", "output", "rjags_Explore_Basecase", "processed"))}

# data----
# loadfonts(device="win") #only need to do this once; takes awhile to run!
coda <- read.csv("state_space_model/output/rjags_Explore_Basecase/coda.csv") 
coda  %>%
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
         MSY.c = R.msy.c-S.msy.c, 
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

# analysis----
# create function for probability profiles and figures
profile(i=10, z=50, xa.start=0, xa.end=8000,lnalpha.c, beta) #can change i,z, xa.start, xa.end
profile_french(i=10, z=50, xa.start=0, xa.end=8000,lnalpha.c, beta) #can change i,z, xa.start, xa.end

parameters <- read.csv("state_space_model/output/rjags_Explore_BaseCase/parameters.csv") 
Taku_sockeye <- read.csv("state_space_model/data/Taku_sockeye.csv") 
ir <- read.csv("state_space_model/data/ir.csv") 
p_q_Nya<- read.csv("state_space_model/output/rjags_Explore_BaseCase/p_q_Nya.csv") 
QM <- read.csv("state_space_model/output/rjags_Explore_BaseCase/processed/QM.csv")
CI<- read.csv("state_space_model/output/rjags_Explore_BaseCase/processed/CI.csv")
coda <- read.csv("state_space_model/output/rjags_Explore_BaseCase/coda.csv") 
parameters <- read.csv("state_space_model/output/rjags_Explore_BaseCase/parameters.csv") 
pubestimates <- read.csv("state_space_model/data/pubestimates.csv") 

# escapement, returns, run abundance, and residuals by year
parameters %>%
  mutate (inriver.run97.5. = as.numeric(inriver.run97.5.),
          inriver.run50. = as.numeric(inriver.run50.),
          year = as.numeric(year),
          S97.5. =as.numeric(S97.5.),
          S50. =as.numeric(S50.)) -> parameters
maxY<-max(parameters$inriver.run97.5., na.rm=TRUE)*1.5
data <- merge(parameters, Taku_sockeye, by=c("year"), all=TRUE)
data <- merge(data, ir, by=c("year"), all=TRUE)
xaxis = tickr(data, year, 4)

# POINT ESTIMATES----
for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/point_est_",lang.use,".png")
  ylab.use <- paste(translate("Inriver Run", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  # inriver run   
  ggplot(data, aes(x=year, y=(inriver.run50.))) +
    geom_line(size=0.75)+ geom_point(size=2)+ylab(ylab.use)+xlab("") +
    geom_ribbon(aes(ymin=(inriver.run2.5.), ymax=(inriver.run97.5.)), alpha=0.20) +
    scale_y_continuous(labels = comma,breaks = seq(0, 250000, 50000), limits = c(0, 250000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    annotate("text",x = 1974, y=250000, label="a)", family="Arial" ,size=6) +
    theme(legend.position = "none") + geom_point(aes(x=year, y=ir), pch=8, size=3,colour="grey40") -> plot1 
  
  # escapement
  ylab.use <- paste(translate("Spawners", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE),"(S)")
  maxY<-max(parameters$S97.5., na.rm=TRUE)*1.5
  ggplot(parameters, aes(x=year, y=(S50.))) +
    geom_line(size=0.75)+ geom_point (size=2)+ylab(ylab.use) + xlab("") +
    geom_ribbon(aes(ymin=(parameters$S2.5.), ymax=(parameters$S97.5.)), alpha=0.20) +
    scale_y_continuous(labels = comma,breaks = seq(0, 150000, 50000), limits = c(0, 150000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") +
    geom_line(aes(y=SMSY), colour="grey40", size=1, linetype=2)+
    annotate("text",x = 1974, y=150000, label="b)", family="Arial" ,size=6) +
    geom_line(aes(y=SGEN), colour="grey40", size=1, linetype=3)-> plot2
  
  # terminal run abundance   
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  
  ylab.use <- paste(translate("Terminal Run Abundance", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE),"(N)")
  maxY<-max(parameters$N97.5., na.rm=TRUE)*1.5
  ggplot(parameters, aes(x=year, y=N50.))+geom_line(size=0.75) + 
    geom_point (size=2) + xlab(xlab.use) +
    ylab(ylab.use) +annotate("text",x = 1974, y=400000, label="c)", family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=N2.5., ymax=N97.5.), alpha=0.20) +
    scale_y_continuous(labels = comma,breaks = seq(0, 400000, 100000), limits = c(0, 400000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") -> plot3
  cowplot::plot_grid(plot1,plot2, plot3,  align = "v", nrow = 3, ncol=1) 
  ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")}

for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/point_est_",lang.use,".png")
  ylab.use <- paste(translate("Inriver Run", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  # inriver run   
  ggplot(data, aes(x=year, y=(inriver.run50.))) +
    geom_line(size=0.75)+ geom_point(size=2)+ylab(ylab.use)+xlab("") +
    geom_ribbon(aes(ymin=(inriver.run2.5.), ymax=(inriver.run97.5.)), alpha=0.20) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE),breaks = seq(0, 250000, 50000), limits = c(0, 250000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    annotate("text",x = 1974, y=250000, label="a)", family="Arial" ,size=6) +
    theme(legend.position = "none") + geom_point(aes(x=year, y=ir), pch=8, size=3,colour="grey40") -> plot1 
  
  # escapement
  ylab.use <- paste(translate("Spawners", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE),"(S)")
  maxY<-max(parameters$S97.5., na.rm=TRUE)*1.5
  ggplot(parameters, aes(x=year, y=(S50.))) +
    geom_line(size=0.75)+ geom_point (size=2)+ylab(ylab.use) + xlab("") +
    geom_ribbon(aes(ymin=(parameters$S2.5.), ymax=(parameters$S97.5.)), alpha=0.20) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 150000, 50000), limits = c(0, 150000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") +
    geom_line(aes(y=SMSY), colour="grey40", size=1, linetype=2)+
    annotate("text",x = 1974, y=150000, label="b)", family="Arial" ,size=6) +
    geom_line(aes(y=SGEN), colour="grey40", size=1, linetype=3)-> plot2
  
  # terminal run abundance   
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  
  ylab.use <- paste(translate("Terminal Run Abundance", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE),"(N)")
  maxY<-max(parameters$N97.5., na.rm=TRUE)*1.5
  ggplot(parameters, aes(x=year, y=N50.))+geom_line(size=0.75) + 
    geom_point (size=2) + xlab(xlab.use) +
    ylab(ylab.use) +annotate("text",x = 1974, y=400000, label="c)", family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=N2.5., ymax=N97.5.), alpha=0.20) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 400000, 100000), limits = c(0, 400000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") -> plot3
  cowplot::plot_grid(plot1,plot2, plot3,  align = "v", nrow = 3, ncol=1) 
  ggsave(out.file, dpi = 500, height = 10, width = 9, units = "in")}


# POINT ESTIMATES RECRUITS----
for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/point_est_rec_",lang.use,".png")
  ylab.use <- paste(translate("Recruits", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE, "(R)"))
  # returns
  maxY<-max(parameters$R97.5., na.rm=TRUE)*1.5
  ggplot(parameters, aes(x=year, y=(R50.))) + geom_line(size=0.75) + 
    geom_point(size=2)+ylab(ylab.use) + xlab("") + annotate("text",x = 1974, y=400000, label="a)", family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=R2.5., ymax=R97.5.), alpha=0.20) +
    scale_y_continuous(labels = comma,breaks = seq(0, 400000, 100000), limits = c(0, 400000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") -> plot1
  
  # Ricker productivity residuals
  maxY<-max(parameters$log.resid97.5., na.rm=TRUE)+0.5
  minY<-min(parameters$log.resid2.5., na.rm=TRUE)-0.5
  ylab.use <- paste(translate("Productivity Residuals", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE, "(R)"))
  xlab.use <- paste(translate("Brood Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE, "(R)"))
  ggplot(parameters, aes(x=year, y=log.resid50.))+geom_line(size=0.75)+geom_point (size=2) + 
    ylab(ylab.use)+xlab(xlab.use) + annotate("text",x = 1974, y=1.00, label="b)", family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=parameters$log.resid2.5., ymax=parameters$log.resid97.5.), alpha=0.20) +
    geom_line(aes(y=0), colour="black", size=0.5) +
    scale_y_continuous(breaks = seq(-1, 1, 0.25), limits = c(-1, 1)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") -> plot2
  cowplot::plot_grid(plot1,plot2,  align = "v", nrow = 2, ncol=1) 
  ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")}

for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/point_est_rec_",lang.use,".png")
  ylab.use <- paste(translate("Recruits", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE, "(R)"))
  # returns
  maxY<-max(parameters$R97.5., na.rm=TRUE)*1.5
  ggplot(parameters, aes(x=year, y=(R50.))) + geom_line(size=0.75) + 
    geom_point(size=2)+ylab(ylab.use) + xlab("") + annotate("text",x = 1974, y=400000, label="a)", family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=R2.5., ymax=R97.5.), alpha=0.20) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 400000, 100000), limits = c(0, 400000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") -> plot1
  
  # Ricker productivity residuals
  maxY<-max(parameters$log.resid97.5., na.rm=TRUE)+0.5
  minY<-min(parameters$log.resid2.5., na.rm=TRUE)-0.5
  ylab.use <- paste(translate("Productivity Residuals", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE, "(R)"))
  xlab.use <- paste(translate("Brood Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE, "(R)"))
  ggplot(parameters, aes(x=year, y=log.resid50.))+geom_line(size=0.75)+geom_point (size=2) + 
    ylab(ylab.use)+xlab(xlab.use) + annotate("text",x = 1974, y=1.00, label="b)", family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=parameters$log.resid2.5., ymax=parameters$log.resid97.5.), alpha=0.20) +
    geom_line(aes(y=0), colour="black", size=0.5) +
    scale_y_continuous(breaks = seq(-1, 1, 0.25), limits = c(-1, 1)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    theme(legend.position = "none") -> plot2
  cowplot::plot_grid(plot1,plot2,  align = "v", nrow = 2, ncol=1) 
  ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")}


# HARVEST RATES WILD (english)----
for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/harvest_rates_wild_",lang.use,".png")
  ylab.use <- paste(translate("Naturally Spawned Harvest Rate", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  label.use <-paste(translate("a) Below Border", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameters, aes(x=year, y=mu.hbelow50.))+geom_line(size=0.75)+geom_point (size=2)+ ylab(ylab.use) +
    xlab("") +
    geom_ribbon(aes(ymin=mu.hbelow_wild2.5., ymax=mu.hbelow_wild97.5.), alpha=0.15) +
    coord_cartesian(ylim=c(0,1)) + annotate("text",x = 1978, y=1.00, label=label.use, family="Arial", size=6) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    geom_line(aes(y=UMSY), colour="grey70", size=1, linetype=2) + theme(legend.position = "none") -> plot1
  
  ylab.use <- paste(translate("Naturally Spawned Harvest Rate", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  label.use <-paste(translate("b) Above Border", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameters, aes(x=year, y=mu.habove_wild50.))+geom_line(size=0.75)+geom_point (size=2)+ ylab(ylab.use) +
    xlab(xlab.use) +  annotate("text",x = 1978, y=1.00, label=label.use, family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=mu.habove_wild2.5., ymax=mu.habove_wild97.5.), alpha=0.15) +
    coord_cartesian(ylim=c(0,1)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    geom_line(aes(y=UMSY), colour="grey70", size=1, linetype=2) + theme(legend.position = "none") -> plot2
  
  cowplot::plot_grid(plot1,plot2,  align = "v", nrow = 2, ncol=1) 
  ggsave(out.file, dpi = 500, height = 10, width = 9, units = "in")}

# HARVEST RATES WILD (french)----
for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/harvest_wild_",lang.use,".png")
  ylab.use <- paste(translate("Naturally Spawned Harvest Rate", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  label.use <-paste(translate("a) Below Border", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameters, aes(x=year, y=mu.hbelow50.))+geom_line(size=0.75)+geom_point (size=2)+ ylab(ylab.use) +
    xlab("") +
    geom_ribbon(aes(ymin=mu.hbelow_wild2.5., ymax=mu.hbelow_wild97.5.), alpha=0.15) +
    coord_cartesian(ylim=c(0,1)) + annotate("text",x = 1981.4, y=1.00, label=label.use, family="Arial", size=6) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    geom_line(aes(y=UMSY), colour="grey70", size=1, linetype=2) + theme(legend.position = "none") -> plot1
  
  ylab.use <- paste(translate("Naturally Spawned Harvest Rate", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  label.use <-paste(translate("b) Above Border", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameters, aes(x=year, y=mu.habove_wild50.))+geom_line(size=0.75)+geom_point (size=2)+ ylab(ylab.use) +
    xlab(xlab.use) +  annotate("text",x = 1981, y=1.00, label=label.use, family="Arial" ,size=6) +
    geom_ribbon(aes(ymin=mu.habove_wild2.5., ymax=mu.habove_wild97.5.), alpha=0.15) +
    coord_cartesian(ylim=c(0,1)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
    geom_line(aes(y=UMSY), colour="grey70", size=1, linetype=2) + theme(legend.position = "none") -> plot2
  
  cowplot::plot_grid(plot1,plot2,  align = "v", nrow = 2, ncol=1) 
  ggsave(out.file, dpi = 500, height = 10, width = 9, units = "in")}

# PROPORTIONS (english)----
for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/props_",lang.use,".png")
  ylab.use <- paste(translate("Age-at-Maturity Proportions", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
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
    theme(legend.title=element_blank(),legend.position=c(0.90,0.86)) +
    ylab(ylab.use) + xlab("") + annotate("text",x = 1973, y=1.00, label= "a)", family="Arial" ,size=6) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot1  
  
  ylab.use <- paste(translate("Age Composition Proportions", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  
  ggplot(data,aes(x=year, y=q, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) +
    scale_fill_grey(start=0.1, end=0.8) + annotate("text",x = 1973, y=1.00, label="b)", family="Arial" ,size=6) +
    ylab(ylab.use) + xlab("") +
    theme(legend.title=element_blank(), legend.position="none") + geom_point(aes(x=year, y=age_comp), position='stack') +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot2
  
  ylab.use <- paste(translate("Terminal Run by Age", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(data,aes(x=year, y=Nya, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) + scale_fill_grey(start=0.1, end=0.8) +
    ylab(ylab.use)+xlab(xlab.use) + annotate("text",x = 1973, y=350000, label="c)", family="Arial" ,size=6) +
    guides(fill = guide_legend(reverse=TRUE)) + 
    theme(legend.title=element_blank(), legend.position="none") +
    scale_y_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot3
  cowplot::plot_grid(plot1,plot2,plot3, align = "v", nrow = 3, ncol=1) 
  ggsave(out.file, dpi = 500, height = 11, width = 8, units = "in")}

# PROPORTIONS (french)----
for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/props_",lang.use,".png")
  ylab.use <- paste(translate("Age-at-Maturity Proportions", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  # mean age at maturity (p), age composition (q), terminal run by age proportions (Nya)
  p_q_Nya %>%
    mutate(age_comp = as.numeric(age_comp),
           year = as.numeric(year),
           Age = factor(age, ordered = TRUE, 
                        levels = c( "Ages 2-4", "Age 5", "Ages 6-8"),
                        labels = c("Âges 2-4", "Âge 5", "Âges 6-8"))) %>%
    mutate(age_comp = ifelse(Age == 'Âges 2-4', NA, age_comp))-> data
  
  ggplot(data,aes(x=year, y=p, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) + scale_fill_grey(start=0.1, end=0.8) +
    theme(legend.title=element_blank(), legend.position=c(0.90,0.86)) +
    ylab(ylab.use) + xlab("") + annotate("text",x = 1973, y=1.00, label= "a)", family="Arial" ,size=6) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot1  
  
  ylab.use <- paste(translate("Age Composition Proportions", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(data,aes(x=year, y=q, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) +
    scale_fill_grey(start=0.1, end=0.8) + annotate("text",x = 1973, y=1.00, label="b)", family="Arial" ,size=6) +
    ylab(ylab.use) + xlab("") +
    theme(legend.title=element_blank(), legend.position="none")  + geom_point(aes(x=year, y=age_comp), position='stack') +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot2
  
  ylab.use <- paste(translate("Terminal Run by Age", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(data,aes(x=year, y=Nya, fill=as.factor(Age))) +
    geom_area(position=position_stack(reverse=FALSE)) + scale_fill_grey(start=0.1, end=0.8) +
    ylab(ylab.use)+xlab(xlab.use) + annotate("text",x = 1973, y=350000, label="c)", family="Arial" ,size=6) +
    guides(fill = guide_legend(reverse=FALSE)) + 
    theme(legend.title=element_blank(), legend.position="none") + 
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 350000, 50000), limits = c(0, 350000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1973, 2020)) -> plot3
  cowplot::plot_grid(plot1,plot2,plot3, align = "v", nrow = 3, ncol=1) 
  ggsave(out.file, dpi = 500, height = 11, width = 8, units = "in")}

# HORSETAIL PLOTS----
parameters <- read.csv("state_space_model/output/rjags_Explore_BaseCase/parameters.csv") 
coda <- read.csv("state_space_model/output/rjags_Explore_Basecase/coda.csv") 
coda  %>%
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
         MSY.c = R.msy.c-S.msy.c, 
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda
QM <- read.csv("state_space_model/output/rjags_Explore_BaseCase/processed/QM.csv")
num <- nrow(QM)
QM %>%
  dplyr::select(c(Escapement)) -> x
coda %>%
  dplyr::select(c(lnalpha.c, beta)) %>%
  filter(row_number()==1:50) %>%
  slice(rep(1:n(), each = num)) %>%
  cbind(., x) %>%#lnalpha.c, beta, and S
  mutate(Recruitment = Escapement*exp(lnalpha.c-beta*Escapement),
         variable = rep(1:50,each=num)) -> dataset

QM %>%
  dplyr::select(c(Escapement)) %>%
  mutate (lnalpha.c = lnalpha.c,
          beta = beta,
          Recruitment = Escapement*exp(lnalpha.c-beta*Escapement),
          variable = 51) %>%
  rbind(., dataset) %>%
  mutate(year = 'NA',
         R2.5. = 0,
         R97.5. = 0,
         S2.5. = 0,
         S97.5. = 0) -> dataset

parameters %>%
  filter (year %in% c(1980:2014)) %>%
  dplyr::select(c(year, S50., R50., R2.5., R97.5., S2.5., S97.5.)) %>%
  mutate(lnalpha.c = 'NA',
         beta = 'NA',
         Escapement = S50.,
         Recruitment = R50.,
         variable = 52) %>%
  dplyr::select(year, lnalpha.c, beta, Escapement, Recruitment, variable, R2.5., R97.5., S2.5., S97.5.) %>%
  rbind(., dataset) -> dataset

dataset %>%
  filter (variable %in% c(52)) %>%
  mutate(Escapement1 = Escapement) %>%
  dplyr::select(-c(lnalpha.c, beta, Escapement)) %>%
  mutate(Escapement = 'NA',
         Median = 'NA',
         q95 = 'NA',
         q90 ='NA',
         q10 ='NA',
         q5 = 'NA') -> dataset
CI %>%
  mutate(year = 'NA',
         R2.5.= 'NA',
         R97.5. = 'NA',
         S2.5. ='NA',
         S97.5. ='NA',
         variable = 51,
         Recruitment = 'NA',
         Escapement1 ='NA') %>%
  rbind(., dataset) %>%
  mutate_if(is.character, as.numeric) -> dataset1

x.fact <- 100/max(dataset1$Escapement1) 
y.fact <- 100/max(dataset1$Recruitment)
coords <- FFieldPtRep(coords = cbind(dataset1$Escapement1 * x.fact, dataset1$Recruitment * y.fact), rep.fact = 40)
x.t <- coords$x/x.fact 
y.t <- coords$y/y.fact

for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/horsetail_",lang.use,".png")
  ylab.use <- paste(translate("Recruits", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Spawners", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(data=dataset1, aes(x=Escapement, y=Median, group=variable)) + 
    geom_line(size=1, lty=2, group=51) +
    geom_ribbon(aes(ymin = q5, ymax = q95, group=51), alpha=.08) +
    geom_ribbon(aes(ymin = q10, ymax = q90, group=51), alpha=.08) +
    xlab(xlab.use) +
    ylab(ylab.use) +
    scale_y_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
    scale_x_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
    geom_line(aes(x=Escapement, y=Escapement, group=51),linetype="solid", size=1) +
    geom_point(data=dataset1, aes(x=x.t, y=y.t, group=52),pch=16, size=1) +
    geom_errorbar(data=dataset1, aes(x=Escapement1, ymax=R97.5., ymin=R2.5., group=52), width=0.2,linetype = 1, colour="grey70") +
    geom_point(data=dataset1, aes(x=Escapement1, y=Recruitment, group=52),pch=16, size=1) +
    geom_errorbarh(data=dataset1, aes(x=Escapement1, y=Recruitment, xmax=S97.5., xmin=S2.5.),na.rm=T,  linetype = 1, colour="grey70") +
    geom_text(size=3, data=dataset1, aes(x=Escapement1, y=Recruitment, group=52, label=year,family="Times", 
                                         hjust = -0.1, vjust= -0.4)) 
  ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")}

for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/horsetail_",lang.use,".png")
  ylab.use <- paste(translate("Recruits (R)", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Spawners (S)", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(data=dataset1, aes(x=Escapement, y=Median, group=variable)) + 
    geom_line(size=1, lty=2, group=51) +
    geom_ribbon(aes(ymin = q5, ymax = q95, group=51), alpha=.08) +
    geom_ribbon(aes(ymin = q10, ymax = q90, group=51), alpha=.08) +
    xlab(xlab.use) +
    ylab(ylab.use) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
    scale_x_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
    geom_line(aes(x=Escapement, y=Escapement, group=51),linetype="solid", size=1) +
    geom_point(data=dataset1, aes(x=x.t, y=y.t, group=52),pch=16, size=1) +
    geom_errorbar(data=dataset1, aes(x=Escapement1, ymax=R97.5., ymin=R2.5., group=52), width=0.2,linetype = 1, colour="grey70") +
    geom_point(data=dataset1, aes(x=Escapement1, y=Recruitment, group=52),pch=16, size=1) +
    geom_errorbarh(data=dataset1, aes(x=Escapement1, y=Recruitment, xmax=S97.5., xmin=S2.5.),na.rm=T,  linetype = 1, colour="grey70") +
    geom_text(size=3, data=dataset1, aes(x=Escapement1, y=Recruitment, group=52, label=year,family="Times", 
                                         hjust = -0.1, vjust= -0.4)) 
  ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")}

# ESCAPEMENT ESTIMATES
parameters <- read.csv("state_space_model/output/rjags_Explore_BaseCase/parameters.csv") 
parameters %>%
  filter (year %in% c(1980:2018)) -> parameter_set

for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/hist_esc_",lang.use,".png")
  ylab.use <- paste(translate("Escapement", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameter_set, aes(x=year, y=S50.)) + 
    geom_line(size=0.75) + geom_point (size=2) + ylab(ylab.use) + xlab(xlab.use) +
    scale_y_continuous(labels = comma,breaks = seq(0, 150000, 25000), limits = c(0, 150000)) + 
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1974, 2020)) +                    
    theme(legend.position = "topright") +
    geom_hline(yintercept = SMAX, colour="grey40", size=1, linetype=4) +
    geom_hline(yintercept = SMSY, colour="grey40", size=1, linetype=2) +
    geom_errorbar(aes(ymin=S2.5., ymax=S97.5.),size=0.5, linetype = "solid", colour="grey40", width=0.02) +
    geom_hline(yintercept=SEQ, colour="grey50", size=1, linetype=1)+
    geom_hline(yintercept=SGEN, colour="grey40", size=1, linetype=3)
  ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")}

for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/hist_esc_",lang.use,".png")
  ylab.use <- paste(translate("Escapement", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameter_set, aes(x=year, y=S50.)) + 
    geom_line(size=0.75) + geom_point (size=2) + ylab(ylab.use) + xlab(xlab.use) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 150000, 25000), limits = c(0, 150000)) + 
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1974, 2020)) +                    
    theme(legend.position = "topright") +
    geom_hline(yintercept = SMAX, colour="grey40", size=1, linetype=4) +
    geom_hline(yintercept = SMSY, colour="grey40", size=1, linetype=2) +
    geom_errorbar(aes(ymin=S2.5., ymax=S97.5.),size=0.5, linetype = "solid", colour="grey40", width=0.02) +
    geom_hline(yintercept=SEQ, colour="grey50", size=1, linetype=1)+
    geom_hline(yintercept=SGEN, colour="grey40", size=1, linetype=3)
  ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")}

# ln(R/S)
parameters <- read.csv("state_space_model/output/rjags_Explore_BaseCase/parameters.csv") 
parameters %>%
  filter (year %in% c(1980:2014)) %>%
  mutate(lnRS = log(R50./S50.))-> parameters

for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/lnRS_",lang.use,".png")
  xlab.use <- paste(translate("Brood Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameters, aes(x=year, y=lnRS))+
    geom_point(size=4, color = "black", shape = 18) +
    xlab(xlab.use) + 
    ylab("ln(R/S)") +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1980, 2020))+
    geom_label((aes(label=year)))+ geom_smooth(method = lm, formula = y ~ splines::bs(x, 10), se = TRUE, colour = "black") 
  ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")}

for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/lnRS_",lang.use,".png")
  xlab.use <- paste(translate("Brood Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  ggplot(parameters, aes(x=year, y=lnRS))+
    geom_point(size=4, color = "black", shape = 18) +
    xlab(xlab.use) + 
    ylab("ln(R/S)") +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1980, 2020))+
    geom_label((aes(label=year)))+ geom_smooth(method = lm, formula = y ~ splines::bs(x, 10), se = TRUE, colour = "black") 
  ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")}

# PREVIOUSLY PUBLISHED ESTIMATES
for(lang.use in c("english")){
  changeLangOpts(L="e")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/prev_publ_",lang.use,".png")
  ylab.use <- paste(translate("Capture-recapture abundance estimates", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  pubestimates %>%
    filter(type == "PP") %>%
    mutate(CI_lower = mean - (sd *2),
           CI_upper = mean + (sd *2))-> data
  pubestimates %>%
    mutate(type = factor(type, ordered = TRUE, 
                         levels = c("PP", "PubEstAdj", "PubEst"),
                         labels = c("PP", "Previously Published (adjusted)", "Previously Published"))) %>%
    ggplot(., aes(x=year, y=mean, fill = type, colour=type, lty=type, shape =type)) + 
    geom_line(size=1) + geom_point (size = 3) +
    xlab(xlab.use) + 
    ylab(ylab.use) +  
    scale_shape_manual(values = c(16,NA,NA)) +
    scale_color_manual(values=c("grey70", "black","black")) +
    scale_linetype_manual(values = c("solid", "dotted", "solid")) +
    scale_y_continuous(labels = comma,breaks = seq(0, 225000, 25000), limits = c(0, 225000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1980, 2018)) +
    geom_ribbon(data=data, aes(ymin=CI_lower, ymax=CI_upper), alpha=0.10, colour ="grey95") +
    theme(legend.title=element_blank(), legend.position=c(0.70,0.10), legend.key = element_blank(), legend.background = element_rect(color = "white")) +
    guides(color=guide_legend(override.aes=list(fill=NA))) -> plot1
  cowplot::plot_grid(plot1, align = "v", nrow = 1, ncol=1) 
  ggsave(out.file, dpi = 500, height = 4, width = 6, units = "in")}

for(lang.use in c("french")){
  changeLangOpts(L="f")
  out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/processed/prev_publ_",lang.use,".png")
  ylab.use <- paste(translate("Capture-recapture abundance estimates", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  xlab.use <- paste(translate("Year", from = "english", to = lang.use, 
                              terms = terms.use, allow_missing = FALSE))
  pubestimates %>%
    filter(type == "PP") %>%
    mutate(CI_lower = mean - (sd *2),
           CI_upper = mean + (sd *2))-> data
  pubestimates %>%
    mutate(type = factor(type, ordered = TRUE, 
                         levels = c("PP", "PubEstAdj", "PubEst"),
                         labels = c("PP", "Publié avec Prévisé (ajusté)", "Publié avec Prévisé"))) %>%
    ggplot(., aes(x=year, y=mean, fill = type, colour=type, lty=type, shape =type)) + 
    geom_line(size=1) + geom_point (size = 3) +
    xlab(xlab.use) + 
    ylab(ylab.use) +  
    scale_shape_manual(values = c(16,NA,NA)) +
    scale_color_manual(values=c("grey70", "black","black")) +
    scale_linetype_manual(values = c("solid", "dotted", "solid")) +
    scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, 225000, 25000), limits = c(0, 225000)) +
    scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1980, 2018)) +
    geom_ribbon(data=data, aes(ymin=CI_lower, ymax=CI_upper), alpha=0.10, colour ="grey95") +
    theme(legend.title=element_blank(), legend.position=c(0.70,0.10), legend.key = element_blank(), legend.background = element_rect(color = "white")) +
    guides(color=guide_legend(override.aes=list(fill=NA))) -> plot1
  cowplot::plot_grid(plot1, align = "v", nrow = 1, ncol=1) 
  ggsave(out.file, dpi = 500, height = 4, width = 6, units = "in")}