#COMPLETE ANALYSIS OF 2016 WILDTRACKS DATA, FEBRUARY 2018

#set working directory

setwd("~/Spider_monkey_stats")

#install packages
install.packages("readxl")
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("nortest")
install.packages("lattice")
install.packages("nlme")
install.packages("lme4")
install.packages("lmerTest")
install.packages("multcomp")
install.packages("igraph")
install.packages("ggplot2")

#load useful packages

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(nortest)
library(lattice)
library(nlme)
library(lme4)
library(lmerTest)
library(multcomp)
library(igraph)
library(ggplot2)

#import code

focals <- read.table(file="focals_Wildtracks.csv", sep=",", header=T)#alternate import code
focals$n_encl<-factor(focals$n_encl, levels=c("1","2"))
focals$month<-factor(focals$month, levels=c("jun-jul","jul-aug","aug-sep"))

#create summarized table

sum_focals1<-focals %>%
  dplyr::filter(actor_id!="FG") %>% 
  dplyr::filter(actor_id!="HUMAN")
  
sum_focals1$actor_id<-factor(sum_focals1$actor_id, levels=c("CL","DU","FY","MA","ME","PA","PE","PO","PP","RK","TR"))
str(sum_focals1)  

sum_focals<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  dplyr::select(duration_sec, occurrences, month, focal_id, time_meal, n_encl, location, activity, focal_sex, focal_group, focal_age) %>%
  group_by(focal_id, focal_sex, focal_group, focal_age, month, time_meal, n_encl, location, activity) %>%
  complete(activity, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))

sum_focals$nest <- with(sum_focals, factor(paste(focal_group,focal_id)))

View(sum_focals)

#Look at mixed effects model
#start without random factor

sum_abnormal<-sum_focals %>%
  filter(activity=="abnormal")

M0_abnormal<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_abnormal, method="ML")

#add random factor

M1_abnormal<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_abnormal, method="ML")


#add nesting

M2_abnormal<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

#anova
anova(M0_abnormal, M1_abnormal, M2_abnormal)
anova(M0_abnormal, M2_abnormal)

#M2 is lowest AIC

#Analyze base model residuals

E2<-residuals(M2_abnormal)
str(E2)
E2
summary(E2)

plot(sum_abnormal$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_abnormal$location, E2, xlab="Location", ylab="Residuals")
plot(sum_abnormal$month, E2, xlab="Month", ylab="Residuals")
plot(sum_abnormal$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_abnormal$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_abnormal$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_abnormal) 

plot(sum_abnormal$prop_time, E2)

#log normalize proportions

sum_abnormal$l.prop_time=log(sum_abnormal$prop_time+1)

#repeat model attempts with normalized data

M0_abnormal<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_abnormal, method="ML")

M1_abnormal<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_abnormal, method="ML")

M2_abnormal<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

anova(M0_abnormal, M1_abnormal, M2_abnormal)
anova(M0_abnormal, M2_abnormal)

E2<-residuals(M2_abnormal)
str(E2)
E2
summary(E2)

plot(sum_abnormal$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_abnormal$location, E2, xlab="Location", ylab="Residuals")
plot(sum_abnormal$month, E2, xlab="Month", ylab="Residuals")
plot(sum_abnormal$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_abnormal$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_abnormal$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_abnormal) 

plot(sum_abnormal$l.prop_time, E2)

#arcsine transformation

asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation

sum_abnormal$as.prop_time=asinTransform(sum_abnormal$prop_time)

#repeat model attempts with arcsine normalized data

M0_abnormal<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_abnormal, method="ML")

M1_abnormal<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_abnormal, method="ML")

M2_abnormal<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

anova(M0_abnormal, M1_abnormal, M2_abnormal)
anova(M0_abnormal, M2_abnormal)

E2<-residuals(M2_abnormal)
str(E2)
E2
summary(E2)

plot(sum_abnormal$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_abnormal$location, E2, xlab="Location", ylab="Residuals")
plot(sum_abnormal$month, E2, xlab="Month", ylab="Residuals")
plot(sum_abnormal$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_abnormal$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_abnormal$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_abnormal) 

plot(sum_abnormal$as.prop_time, E2)

#check for autocorrelation

acf(E2, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#little autocorrelation in the first lag

#backwards selection to optimize model in Chapter 5 of Zuur

M.full<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                       random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

summary(M.full)

#can drop focal_age due to high P value

M.full.a<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex, 
            random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

#compare M.full and M.full.a

anova(M.full,M.full.a)

summary(M.full.a)



M.full.b<-lme(as.prop_time ~ n_encl + location + month + time_meal, 
              random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

anova(M.full.a,M.full.b)

summary(M.full.b)

M.full.c<-lme(as.prop_time ~ location + month + time_meal, 
              random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

anova(M.full.b,M.full.c)

summary(M.full.c)

#####################################################
#Get Full Model Statistics and Make Graph
#####################################################

#best overall model was arcsine transformed
M.full<-lme(as.prop_time ~ location + month + time_meal, 
              random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="REML")
     #need to change method to REML or drop the method term since it defaults to REML

anova(M.full)#this is the summary of your full model

#post hoc tests
model.matrix.gls <- function(M.full, ...){
  model.matrix(terms(M.full), data = getData(M.full), ...)  
}
model.frame.gls <- function(M.full, ...){
  model.frame(formula(M.full), data = getData(M.full), ...)  
}
terms.gls <- function(M.full, ...){
  terms(model.frame(M.full),...)  
}

multCompTukey1 <- glht(M.full, linfct = mcp(location = "Tukey")) 
summary(multCompTukey1)

multCompTukey2 <- glht(M.full, linfct = mcp(month = "Tukey")) 
summary(multCompTukey2)

multCompTukey3 <- glht(M.full, linfct = mcp(time_meal = "Tukey")) 
summary(multCompTukey3)


#make a graphic

x <- group_by(sum_abnormal, month) %>% 
  summarize(time.avg = mean(as.prop_time, na.rm = TRUE), 
            time.sd=sd(as.prop_time, na.rm = TRUE),  
            n = sum(!is.na(as.prop_time)), 
            time.se=time.sd/sqrt(n))

ggplot(data=x, 
       aes(x=month, y=time.avg)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=time.avg, ymax=time.avg+time.se), width=0.2, 
                position=position_dodge(0.9)) + 
  #scale_fill_manual(values=c("black","gray","white")) +
  xlab("Time") +
  ylab(expression(Abnormal~Behavior~(arcsin~sqrt~proportion))) +
  #ylim(0,0.085) +
  #labs(fill="Amendment") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=6),
        legend.key=element_blank(),
        legend.position=c(0.5,0.95),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        legend.direction="horizontal",
        legend.key.size=unit(0.3, "cm"),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8))

ggsave('Figures/abnormal_month.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")
