# THIS SCRIPT IS SOURCED FROM INSIDE 1_RUN_MODEL.R
#(if package.use = "rjags" and location of file is rjags_base_case) 

# - post.samp and post 
      # these are the mcmc.list objects created by coda.samples
      # you can access individual variables like this: post[,"var.name"]
# - post.arr 
# - coda
# create coda file
#Gelman statistic
# Brooks and Gelman (1997) have suggested, if Rc<1.2 for all model parameters, 
# one can be fairly confident that convergence has been reached. Otherwise, longer chains or other means for improving the convergence may be needed
# Brooks, S. P., and A. Gelman. 1997. General Methods for Monitoring Convergence of Iterative Simulations. 
# Journal of Computational and Graphical Statistics 7: 434–455.

loadfonts(device="win")
windowsFonts(Times=windowsFont("Arial"))
theme_set(theme_report(base_size = 14))

parameters <- c('S.eq.c','S.msy.c','U.msy.c','alpha','beta','alpha.c',
                  'lnalpha','lnalpha.c','phi','sigma.R','log.resid.0', 'mean.log.RO',
                  'S','R','N','log.resid','mu.hbelow','pi','h.below','N.ya',
                  'p','q', 'S.max','D.sum','q.weir', 'q.mr','D.scale','sigma.RO',
                  'log.qS.weir', 'log.qS.mr','qS.weir', 'qS.mr',
                  'S.eq.c2', 'U.msy.c2', 'S.msy.c2', 'U.max.c2')

# Numerical summary of each parameter (mean, median, quantiles of posteriers)----
summary<-summary(post)  
stats<-summary$statistics;  colnames(stats)
quants<-summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants) 
statsquants %>% 
  as.data.frame() %>% 
  dplyr::select(Mean, SD, 'Time-series SE', '2.5%', '50%', '97.5%') %>%
  dplyr::rename(time_series_se = 'Time-series SE') %>%
  rownames_to_column('variable') %>%
  mutate (mc_error = time_series_se/SD,
  converge = ifelse(mc_error < 0.05, "true", "false")) %>%
  write.csv(., file= paste0(out.path,"/statsquants.csv"))    

# Gelman Diagnostics----
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold = 1.2 #values less than 1.2 are generally considered converged
gel %>%
  rownames_to_column('variable') %>%
  rename(point_estimate = "Point est.") %>%
  mutate (converge = ifelse(point_estimate < poor.threshold, "true", "false")) %>%
  write.csv(., file= paste0(out.path,"/gelman.csv") )   

# Geweke Diagnostics----
#Examine convergence of the Markov chains using the Geweke's convergence diagnostic
#a convergence diagnostic for Markov chains based on a test for equality 
#of the means of the ﬁrst and last part of a Markov chain (by default the ﬁrst 10% and the last 50%).
#As for a cut-off you can compare to the standard normal critical values z α/2 where α=0.05
#null hypothesis is rejected if Z large (used to determine the burn-in period)
# Geweke diagnostic
l <- geweke.diag(post)
df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T),stringsAsFactors=FALSE)
df <- t(df) %>%
  as.data.frame() %>%
  rownames_to_column('variable1') -> df
colnames(df) <- c("variable1", "chain1", "chain2", "chain3")
names <-read.csv(file = paste0(out.path,"/statsquants.csv"))
names  %>%
  dplyr::select(variable) %>% 
  add_row(variable = "frac1") %>% #0.1 and 0.5 default
  add_row(variable = "frac2")  -> names
x <- cbind(names, df)
x  %>% 
  dplyr::select(-variable1) %>%
  write.csv(., paste0(out.path,"/geweke.csv")) 

pdf("2023_analysis/output/rjags_base_case/geweke.pdf",height=10, width=8,onefile=T)
geweke.plot(post)
dev.off()

# 90th and 95th percentiles of output quants----
parameters=c('alpha')
x <- post.arr[,parameters,]
alpha<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
alpha <- data.frame(alpha)

parameters=c('lnalpha')
x <- post.arr[,parameters,]
lnalpha<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
lnalpha <- data.frame(lnalpha)

parameters=c('lnalpha.c')
x <- post.arr[,parameters,]
lnalpha.c<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
lnalpha.c <- data.frame(lnalpha.c)

parameters=c('beta')
x <- post.arr[,parameters,]
beta<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
beta <- data.frame(beta)

parameters=c('phi')
x <- post.arr[,parameters,]
phi<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
phi <- data.frame(phi)

parameters=c('sigma.R')
x <- post.arr[,parameters,]
sigma.R<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
sigma.R <- data.frame(sigma.R)

parameters=c('S.eq.c')
x <- post.arr[,parameters,]
S.eq.c<-quantile(x, probs=c(0,0.025, 0.05, 0.5, 0.95, 0.975, 1))
S.eq.c <- data.frame(S.eq.c)

parameters=c('S.max')
x <- post.arr[,parameters,]
S.max<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
S.max <- data.frame(S.max)

parameters=c('D.sum')
x <- post.arr[,parameters,]
D.sum<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
D.sum <- data.frame(D.sum)

parameters=c('pi[1]')
x <- post.arr[,parameters,]
pi<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
pi1 <- data.frame(pi)

parameters=c('pi[2]')
x <- post.arr[,parameters,]
pi<-quantile(x, probs=c(0,0.025, 0.05, 0.5, 0.95, 0.975, 1))
pi2 <- data.frame(pi)

parameters=c('pi[3]')
x <- post.arr[,parameters,]
pi<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
pi3 <- data.frame(pi)

parameters=c('alpha.c')
x <- post.arr[,parameters,]
alpha.c<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
alpha.c <- data.frame(alpha.c)

step1 <- cbind(alpha, lnalpha)
step2 <- cbind(step1, lnalpha.c)
step3 <- cbind(step2, beta)
step4 <- cbind(step3, phi)
step5 <- cbind(step4, sigma.R)
step6 <- cbind(step5, S.eq.c)
step7 <- cbind(step6, S.max)
step8 <- cbind(step7, D.sum)
step9 <- cbind(step8, pi1)
step10 <- cbind(step9, pi2)
step11 <- cbind(step10, pi3)
step12 <- cbind(step11, alpha.c)
step12 %>% 
  t()%>%
write.csv(., file= paste0(out.path,"/percentiles.csv")) 

# p_q_Nya.csv file----
rawdat<-as.data.frame(read.csv("state_space_model/data/Taku_sockeye.csv",header=T) ) #get age_comp data
rawdat  %>% 
  dplyr::select (c(year, x4, x5, x6)) %>%
  mutate(age4 = x4/(x4+x5+x6),
         age5 = x5/(x4+x5+x6),
         age6 = x6/(x4+x5+x6)) %>%
  dplyr::select (c(year, age4, age5, age6)) %>%
  gather(key="age", value="age_comp", -year) %>%
  mutate(age = ifelse(age == 'age4', "Ages 2-4",
                      ifelse(age== 'age5', "Age 5","Ages 6-8"))) -> age_data
summary <- summary(post)# extract p
stats <- summary$statistics;  colnames(stats)
stats <- stats[,c(1)] #select columns of interest
data.frame(stats) %>% 
  rownames_to_column('variable') -> stats
stats %>% 
  filter(str_detect(variable, "p")) %>% 
  filter(!str_detect(variable, "alpha")) %>%
  filter(!str_detect(variable, "lnalpha")) %>% 
  filter(!str_detect(variable, "lnalpha.c")) %>%
  filter(!str_detect(variable, "phi")) %>%
  filter(!str_detect(variable, "pi")) %>%
  rename_at(vars(starts_with("stats")), 
            funs(str_replace(., "stats", "p"))) %>%
  dplyr::select(-c(variable))-> df
new.df <- data.frame(
  year = c(1974:2014, 1974:2014, 1974:2014),
  age = rep(c("Ages 2-4", "Age 5", "Ages 6-8"), each = 41))
p_data <- cbind(df, new.df)

stats %>% 
  filter(str_detect(variable, "q")) %>% 
  filter(!str_detect(variable, "S.eq.c")) %>%
  filter(!str_detect(variable, " S.eq.c2")) %>%
  rename_at(vars(starts_with("stats")), 
            funs(str_replace(., "stats", "q"))) %>%
  dplyr::select(-c(variable))-> df
new.df <- data.frame(
  year = c(1980:2018, 1980:2018, 1980:2018),
  age = rep(c("Ages 2-4", "Age 5", "Ages 6-8"), each = 39))
q_data <- cbind(df, new.df)

stats %>% 
  filter(str_detect(variable, "N.ya")) %>% 
  rename_at(vars(starts_with("stats")), 
            funs(str_replace(., "stats", "Nya"))) %>%
  dplyr::select(-c(variable))-> df
Nya_data <- cbind(df, new.df)
step1 <- merge(age_data, p_data, by=c("year", "age"), all=TRUE)
step2 <- merge(step1, q_data, by=c("year", "age"), all=TRUE)
p_q_Nya <- merge(step2, Nya_data, by=c("year", "age"), all=TRUE)
p_q_Nya <- p_q_Nya[order(p_q_Nya$age, p_q_Nya$year),]
p_q_Nya %>% mutate(age = as.character(age))
write.csv(p_q_Nya, file= paste0(out.path,"/p_q_Nya.csv")) 

# parameters file----
data.frame(quants) %>% 
  rownames_to_column('variable') -> quants
quants %>%
  filter(str_detect(variable, "h.above")) %>% 
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "h.above"))) %>%
  dplyr::select(-c(variable))-> df
new.df <- data.frame(
  year = rep(c(1980:2018),1))
df1 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "h.below")) %>% 
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "h.below"))) %>%
  dplyr::select(-c(variable))-> df
df2 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "mu.hbelow_wild")) %>% 
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "mu.hbelow_wild"))) %>%
  dplyr::select(-c(variable))-> df
df3 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "mu.habove_wild")) %>% 
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "mu.habove_wild"))) %>%
  dplyr::select(-c(variable))-> df
df4 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "mu.hbelow")) %>% 
  filter(!str_detect(variable, "mu.hbelow_wild")) %>%
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "mu.hbelow"))) %>%
  dplyr::select(-c(variable))-> df
df5 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "mu.habove")) %>% 
  filter(!str_detect(variable, "mu.habove_wild")) %>%
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "mu.habove"))) %>%
  dplyr::select(-c(variable))-> df
df6 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "inriver.run")) %>% 
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "inriver.run"))) %>%
  dplyr::select(-c(variable))-> df
df7 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "S")) %>% 
  filter(!str_detect(variable, "S.eq.c")) %>%
  filter(!str_detect(variable, " S.eq.c2")) %>%
  filter(!str_detect(variable, "S.max")) %>%
  filter(!str_detect(variable, "S.msy.c")) %>%
  filter(!str_detect(variable, "S.msy.c2")) %>%
  filter(!str_detect(variable, "S.msy.c.80")) %>%
  filter(!str_detect(variable, "S.msy.c2.80")) %>%
  filter(!str_detect(variable, "lnRS")) %>%
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "S"))) %>%
  dplyr::select(-c(variable))-> df
df8 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "N")) %>% 
  filter(!str_detect(variable, "N.ya")) %>%
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "N"))) %>%
  dplyr::select(-c(variable))-> df
df9 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "log.resid")) %>% 
  filter(!str_detect(variable, "log.resid.0")) %>%
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "log.resid"))) %>%
  dplyr::select(-c(variable))-> df
new.df <- data.frame(
  year = rep(c(1980:2014),1))
df10 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "R")) %>%
  filter(!str_detect(variable, "mean.log.RO")) %>%
  filter(!str_detect(variable, "sigma.R")) %>%
  filter(!str_detect(variable, "sigma.R0")) %>%
  filter(!str_detect(variable, "lnRS")) %>%
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "R"))) %>%
  dplyr::select(-c(variable))-> df
new.df <- data.frame(
  year = rep(c(1974:2014),1))
df11 <- cbind(df, new.df)

quants %>%
  filter(str_detect(variable, "lnRS")) %>% 
  rename_at(vars(starts_with("X")), 
            funs(str_replace(., "X", "lnRS"))) %>%
  dplyr::select(-c(variable))-> df
new.df <- data.frame(
  year = rep(c(1980:2018),1))
df12 <- cbind(df, new.df)

step1 <- merge(df1, df2, by=c("year"), all=TRUE)
step2 <- merge(step1, df3, by=c("year"), all=TRUE)
step3 <- merge(step2, df4, by=c("year"), all=TRUE)
step4 <- merge(step3, df5, by=c("year"), all=TRUE)
step5 <- merge(step4, df6, by=c("year"), all=TRUE)
step6 <- merge(step5, df7, by=c("year"), all=TRUE)
step7 <- merge(step6, df8, by=c("year"), all=TRUE)
step8 <- merge(step7, df9, by=c("year"), all=TRUE)
step9 <- merge(step8, df10, by=c("year"), all=TRUE)
step10 <- merge(step9, df11, by=c("year"), all=TRUE)
step11 <- merge(step10, df12, by=c("year"), all=TRUE)
write.csv(step11, file= paste0(out.path,"/parameters.csv")) 

# lambert calc----
parameters=c("lnalpha", "beta", "lnalpha.c")
x <- as.data.frame(post.arr[,parameters,])
coda1 <- x[,1:3]
coda2 <- x[,4:6]
coda3 <- x[,7:9]
coda1 %>% 
  dplyr::rename(beta = beta.1,
         lnalpha = lnalpha.1,
         lnalpha.c = lnalpha.c.1)-> coda1
coda2 %>% 
  dplyr::rename(beta = beta.2,
         lnalpha = lnalpha.2,
         lnalpha.c = lnalpha.c.2) -> coda2
coda3 %>% 
  dplyr::rename(beta = beta.3,
         lnalpha = lnalpha.3,
         lnalpha.c = lnalpha.c.3) -> coda3
coda<-rbind(coda1,coda2,coda3)
write.csv(coda, file= paste0(out.path,"/coda.csv") ,row.names=FALSE)  
coda %>% 
  mutate(Smsy_lambert = (1-lambert_W0(exp(1-lnalpha.c)))/beta,
         Umsy_lambert = (1-lambert_W0(exp(1-lnalpha.c))),
         Smsy_lambert80 = Smsy_lambert *0.80)  %>%
  as.data.frame() %>%
  dplyr::select(Smsy_lambert,Umsy_lambert, Smsy_lambert80) -> coda
coda %>% 
  apply(., 2, sd) %>%
  as.data.frame()%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('variable') %>%
  mutate(variable = ifelse(variable == '.', "sd", "sd"))-> sd
coda %>% 
  apply(., 2, mean) %>%
  as.data.frame() %>%
  t()%>%
  as.data.frame() %>%
  rownames_to_column('variable') %>%
  mutate(variable = ifelse(variable == '.', "mean", "mean"))-> mean
rbind(sd,mean) -> x1
q1<-apply(coda,2,quantile,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
q1 %>%
  as.data.frame() %>%
  rownames_to_column('variable') -> x2
rbind(x1,x2) %>%
  t() %>%
write.csv(., file= paste0(out.path,"/quantiles_lambert.csv"))    

# lambert density plot----
parameters=c("lnalpha", "beta", "lnalpha.c")
x <- as.data.frame(post.arr[,parameters,])
coda1 <- x[,1:3]
coda2 <- x[,4:6]
coda3 <- x[,7:9]
coda1 %>% 
  dplyr::rename(beta = beta.1,
                lnalpha = lnalpha.1,
                lnalpha.c = lnalpha.c.1)-> coda1
coda2 %>% 
  dplyr::rename(beta = beta.2,
                lnalpha = lnalpha.2,
                lnalpha.c = lnalpha.c.2) -> coda2
coda3 %>% 
  dplyr::rename(beta = beta.3,
                lnalpha = lnalpha.3,
                lnalpha.c = lnalpha.c.3) -> coda3
coda<-rbind(coda1,coda2,coda3)
coda %>% 
  mutate(Smsy = (1-lambert_W0(exp(1-lnalpha.c)))/beta,
         Umsy = (1-lambert_W0(exp(1-lnalpha.c))),
         Smsy80 = Smsy *0.80)  %>%
  as.data.frame() -> coda
coda %>% 
  dplyr::select(Smsy) -> Smsy
coda %>% 
  dplyr::select(Smsy80) -> Smsy80
coda %>% 
  dplyr::select(Umsy) -> Umsy

# density plot
out.file <- paste0("state_space_model/output/rjags_Explore_BaseCase/density.png")
options(scipen=999)
ggplot(Smsy, aes(x=Smsy, fill=Smsy, color = Smsy)) +
  geom_density(fill ="#999999", alpha=0.5)+ annotate("text",x = 0, y=0.00004, label="a)", family="Arial" ,size=6) +
  scale_color_manual(values=c("#999999"))+
  scale_fill_manual(values=c("#999999"))+
  geom_vline(xintercept = 43857,linetype = "longdash" ) +
  labs(x="Smsy",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                             'Arial')+
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())) + theme(legend.position="none")+ 
 scale_x_continuous(labels = comma,breaks = seq(0, 150000, 25000), limits = c(0, 150000)) -> plot1

ggplot(Smsy80, aes(x=Smsy80, fill=Smsy80)) +
  geom_density(fill ="#999999", alpha=0.5)+ annotate("text",x = 0, y=0.00005, label="b)", family="Arial" ,size=6) +
  scale_color_manual(values=c("#999999"))+geom_vline(xintercept = 35086,linetype = "longdash" ) +
  scale_fill_manual(values=c("#999999"))+
  labs(x="Smsy80",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                                    'Arial')+
                                           theme(panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank())) +theme(legend.position="none")+ 
  scale_x_continuous(labels = comma,breaks = seq(0, 150000, 25000), limits = c(0, 150000)) -> plot2

ggplot(Umsy, aes(x=Umsy, fill=Umsy)) +
  geom_density(fill ="#999999", alpha=0.5)+ annotate("text",x = 0, y=5, label="c)", family="Arial" ,size=6) +
  scale_color_manual(values=c("#999999"))+
  scale_fill_manual(values=c("#999999"))+geom_vline(xintercept = 0.75,linetype = "longdash" ) +
  labs(x="Umsy",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                                      'Arial')+
                                             theme(panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank())) + theme(legend.position="none")+ 
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) -> plot3
cowplot::plot_grid(plot1,plot2,plot3,  align = "v", nrow = 3, ncol=1) 
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

 #trace and density plots----
parameters <- c('alpha', 'lnalpha', 'lnalpha.c','beta','phi','sigma.R','S.eq.c','S.max','D.sum', 'pi')
png("state_space_model/output/rjags_Explore_BaseCase/density_other.png",res=500, width=8, height=9, units ="in")
denplot(post, parms = c(parameters), style="plain", greek=TRUE, col="gray30", collapse =T)
dev.off()
dev.off()

parameters <- c('S')
pdf("state_space_model/output/rjags_Explore_BaseCase/density2.pdf",height=10, width=8,onefile=F)
denplot(post, parms = c(parameters))
dev.off()
pdf("state_space_model/output/rjags_Explore_BaseCase/trace2.pdf",height=10, width=8,onefile=F)
traplot(post, parms = c(parameters))
dev.off()
parameters <- c('R')
pdf("state_space_model/output/rjags_Explore_BaseCase/density3.pdf",height=10, width=8,onefile=F)
denplot(post, parms = c(parameters))
dev.off()
pdf("state_space_model/output/rjags_Explore_BaseCase/trace3.pdf",height=10, width=8,onefile=F)
traplot(post, parms = c(parameters))
dev.off()
parameters <- c('N')
pdf("state_space_model/output/rjags_Explore_BaseCase/density4.pdf",height=10, width=8,onefile=F)
denplot(post, parms = c(parameters))
dev.off()
pdf("state_space_model/output/rjags_Explore_BaseCase/trace4.pdf",height=10, width=8,onefile=F)
traplot(post, parms = c(parameters))
dev.off()
parameters <- c('p')
pdf("state_space_model/output/rjags_Explore_BaseCase/density5.pdf",height=10, width=8,onefile=F)
denplot(post, parms = c(parameters))
dev.off()
pdf("state_space_model/output/rjags_Explore_BaseCase/trace5.pdf",height=10, width=8,onefile=F)
traplot(post, parms = c(parameters))
dev.off()
parameters <- c('q')
pdf("state_space_model/output/rjags_Explore_BaseCase/density6.pdf",height=10, width=8,onefile=F)
denplot(post, parms = c(parameters))
dev.off()
pdf("state_space_model/output/rjags_Explore_BaseCase/trace6.pdf",height=10, width=8,onefile=F)
traplot(post, parms = c(parameters))
dev.off()

# autocorrelation plots----
windows(record=T)
pdf("state_space_model/output/rjags_Explore_BaseCase/autocorr.pdf",height=6, width=8,onefile=T,useDingbats=F)
autocorr.plot(post, lag.max=5)
dev.off()
dev.off()
autocorr.summary<-autocorr.diag(post)
autocorr.summary<-data.frame(autocorr.summary)
write.csv(autocorr.summary, file= paste0(out.path,"/autocorr.csv")) 

#density and time series plots----
post.samp <- post
nvars <- dim(post.samp[[1]])[2]
nsamps <- dim(post.samp[[1]])[1]
int <- 25
pdf("state_space_model/output/rjags_Explore_BaseCase/profiles.pdf",height=6, width=8)
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