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

sum_not_visible<-sum_focals %>%
  filter(activity=="not_visible")

M0_not_visible<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_not_visible, method="ML")

#add random factor

M1_not_visible<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_not_visible, method="ML")


#add nesting

M2_not_visible<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_not_visible, method="ML")

#anova
anova(M0_not_visible, M1_not_visible, M2_not_visible)
anova(M0_not_visible, M2_not_visible)
anova(M1_not_visible, M2_not_visible)

#M0 has lowest AIC

#Analyze base model residuals

E0<-residuals(M0_not_visible)
str(E0)
E0
summary(E0)

plot(sum_not_visible$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(sum_not_visible$location, E0, xlab="Location", ylab="Residuals")
plot(sum_not_visible$month, E0, xlab="Month", ylab="Residuals")
plot(sum_not_visible$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(sum_not_visible$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(sum_not_visible$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_not_visible) 

plot(sum_not_visible$prop_time, E0)

#log normalize proportions
sum_not_visible$l.prop_time=log(sum_not_visible$prop_time+1)

#repeat model attempts with normalized data

M0_not_visible<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_not_visible, method="ML")

M1_not_visible<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_not_visible, method="ML")

M2_not_visible<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_not_visible, method="ML")

anova(M0_not_visible, M1_not_visible, M2_not_visible)
anova(M1_not_visible, M2_not_visible)

#M0 still lowest AIC

E0<-residuals(M0_not_visible)
str(E0)
E0
summary(E0)

plot(sum_not_visible$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(sum_not_visible$location, E0, xlab="Location", ylab="Residuals")
plot(sum_not_visible$month, E0, xlab="Month", ylab="Residuals")
plot(sum_not_visible$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(sum_not_visible$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(sum_not_visible$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_not_visible) 

plot(sum_not_visible$l.prop_time, E0)

#arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation
sum_not_visible$as.prop_time=asinTransform(sum_not_visible$prop_time)

#repeat model attempts with arcsine normalized data

M0_not_visible<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_not_visible, method="ML")

M1_not_visible<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_not_visible, method="ML")

M2_not_visible<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_not_visible, method="ML")

anova(M0_not_visible, M1_not_visible, M2_not_visible)
anova(M0_not_visible, M2_not_visible)
anova(M0_not_visible, M1_not_visible)

#M0 still lowest

E0<-residuals(M0_not_visible)
str(E0)
E0
summary(E0)

plot(sum_not_visible$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(sum_not_visible$location, E0, xlab="Location", ylab="Residuals")
plot(sum_not_visible$month, E0, xlab="Month", ylab="Residuals")
plot(sum_not_visible$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(sum_not_visible$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(sum_not_visible$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_not_visible) 

plot(sum_not_visible$as.prop_time, E0)

#check for autocorrelation
acf(E0, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#no autocorrelation

#arcsine transformation helped, but residual boxplots still look uneven

#too many zeros? very difficult to fit a model.

