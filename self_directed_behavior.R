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

sum_self_directed<-sum_focals %>%
  filter(activity=="self_directed")

M0_self_directed<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_self_directed, method="ML")

#add random factor

M1_self_directed<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_self_directed, method="ML")


#add nesting

M2_self_directed<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_self_directed, method="ML")

#anova
anova(M0_self_directed, M1_self_directed, M2_self_directed)
anova(M0_self_directed, M2_self_directed)
anova(M1_self_directed, M2_self_directed)

#M1 has lowest AIC

#Analyze base model residuals

E1<-residuals(M1_self_directed)
str(E1)
E1
summary(E1)

plot(sum_self_directed$n_encl, E1, xlab="# Enclosures", ylab="Residuals")
plot(sum_self_directed$location, E1, xlab="Location", ylab="Residuals")
plot(sum_self_directed$month, E1, xlab="Month", ylab="Residuals")
plot(sum_self_directed$time_meal, E1, xlab="Meal Status", ylab="Residuals")
plot(sum_self_directed$focal_sex, E1, xlab="Focal Sex", ylab="Residuals")
plot(sum_self_directed$focal_age, E1, xlab="Focal Age", ylab="Residuals")

qqnorm(E1)
qqline(E1)
ad.test(E1)

plot(M1_self_directed) 

plot(sum_self_directed$prop_time, E1)

#log normalize proportions
sum_self_directed$l.prop_time=log(sum_self_directed$prop_time+1)

#repeat model attempts with normalized data

M0_self_directed<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_self_directed, method="ML")

M1_self_directed<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_self_directed, method="ML")

M2_self_directed<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_self_directed, method="ML")

anova(M0_self_directed, M1_self_directed, M2_self_directed)
anova(M1_self_directed, M2_self_directed)

#M1 still lowest AIC

E1<-residuals(M1_self_directed)
str(E1)
E1
summary(E1)

plot(sum_self_directed$n_encl, E1, xlab="# Enclosures", ylab="Residuals")
plot(sum_self_directed$location, E1, xlab="Location", ylab="Residuals")
plot(sum_self_directed$month, E1, xlab="Month", ylab="Residuals")
plot(sum_self_directed$time_meal, E1, xlab="Meal Status", ylab="Residuals")
plot(sum_self_directed$focal_sex, E1, xlab="Focal Sex", ylab="Residuals")
plot(sum_self_directed$focal_age, E1, xlab="Focal Age", ylab="Residuals")

qqnorm(E1)
qqline(E1)
ad.test(E1)

plot(M1_self_directed) 

plot(sum_self_directed$prop_time, E1)

#arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation
sum_self_directed$as.prop_time=asinTransform(sum_self_directed$prop_time)

#repeat model attempts with arcsine normalized data

M0_self_directed<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_self_directed, method="ML")

M1_self_directed<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_self_directed, method="ML")

M2_self_directed<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_self_directed, method="ML")

anova(M0_self_directed, M1_self_directed, M2_self_directed)
anova(M0_self_directed, M2_self_directed)
anova(M0_self_directed, M1_self_directed)

#M2 now best

E2<-residuals(M2_self_directed)
str(E2)
E2
summary(E2)

plot(sum_self_directed$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_self_directed$location, E2, xlab="Location", ylab="Residuals")
plot(sum_self_directed$month, E2, xlab="Month", ylab="Residuals")
plot(sum_self_directed$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_self_directed$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_self_directed$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_self_directed) 

plot(sum_self_directed$prop_time, E2)

#check for autocorrelation
acf(E2, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#no autocorrelation except 1st lag

#M1 w/ log-transformed data looks slightly better?
#M2 Q-Q plot looks better w/ arcsine transformation?


