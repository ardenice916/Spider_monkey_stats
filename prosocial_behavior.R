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

sum_prosocial<-sum_focals %>%
  filter(activity=="prosocial")

M0_prosocial<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_prosocial, method="ML")

#add random factor

M1_prosocial<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_prosocial, method="ML")


#add nesting

M2_prosocial<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_prosocial, method="ML")

#anova
anova(M0_prosocial, M1_prosocial, M2_prosocial)
anova(M0_prosocial, M2_prosocial)
anova(M1_prosocial, M2_prosocial)

#M2 has lowest AIC

#Analyze base model residuals

E2<-residuals(M2_prosocial)
str(E2)
E2
summary(E2)

plot(sum_prosocial$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_prosocial$location, E2, xlab="Location", ylab="Residuals")
plot(sum_prosocial$month, E2, xlab="Month", ylab="Residuals")
plot(sum_prosocial$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_prosocial$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_prosocial$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_prosocial) 

plot(sum_prosocial$prop_time, E2)

#log normalize proportions
sum_prosocial$l.prop_time=log(sum_prosocial$prop_time+1)

#repeat model attempts with normalized data

M0_prosocial<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_prosocial, method="ML")

M1_prosocial<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_prosocial, method="ML")

M2_prosocial<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_prosocial, method="ML")

anova(M0_prosocial, M1_prosocial, M2_prosocial)
anova(M1_prosocial, M2_prosocial)

#M2 still lowest AIC

E2<-residuals(M2_prosocial)
str(E2)
E2
summary(E2)

plot(sum_prosocial$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_prosocial$location, E2, xlab="Location", ylab="Residuals")
plot(sum_prosocial$month, E2, xlab="Month", ylab="Residuals")
plot(sum_prosocial$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_prosocial$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_prosocial$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_prosocial) 

plot(sum_prosocial$prop_time, E2)

#arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation
sum_prosocial$as.prop_time=asinTransform(sum_prosocial$prop_time)

#repeat model attempts with arcsine normalized data

M0_prosocial<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_prosocial, method="ML")

M1_prosocial<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_prosocial, method="ML")

M2_prosocial<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_prosocial, method="ML")

anova(M0_prosocial, M1_prosocial, M2_prosocial)
anova(M0_prosocial, M2_prosocial)
anova(M0_prosocial, M1_prosocial)

#M2 still best

E2<-residuals(M2_prosocial)
str(E2)
E2
summary(E2)

plot(sum_prosocial$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_prosocial$location, E2, xlab="Location", ylab="Residuals")
plot(sum_prosocial$month, E2, xlab="Month", ylab="Residuals")
plot(sum_prosocial$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_prosocial$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_prosocial$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_prosocial) 

plot(sum_prosocial$prop_time, E2)

#check for autocorrelation
acf(E2, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#potential autocorrelation in the first lag or 4-5, 19?

#even arcsine transformation doesn't normalize data... plotting residuals vs. model on line 192 looks baaaad


