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

#view and summarize

View(focals)
summary(focals)
str(focals$actor_id)
str(focals$focal_id)

#convert durations to minutes

focals$dur_m <- focals$duration_sec/60
focals$dur_h <- focals$dur_m/60

sum(focals$duration_sec)
sum(focals$dur_m)
sum(focals$dur_h)
#199.67 focal hours collected

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
anova(M0_abnormal, M1_abnormal, M2_abnormal, M3_abnormal)

#M1 is lowest AIC

#Analyze base model residuals

E1<-residuals(M1_abnormal)

plot(sum_abnormal %>% dplyr::select(focal_group), E1, xlab="Focal Group", ylab="Residuals")
plot(sum_abnormal %>% dplyr::select(n_encl),
     E1, xlab="# Enclosures", ylab="Residuals")
plot(sum_abnormal %>% dplyr::select(location),
     E1, xlab="Location", ylab="Residuals")
plot(sum_abnormal %>% dplyr::select(month),
     E1, xlab="Month", ylab="Residuals")
plot(sum_abnormal %>% dplyr::select(time_meal),
     E1, xlab="Meal Status", ylab="Residuals")
plot(sum_abnormal %>% dplyr::select(focal_sex),
     E1, xlab="Focal Sex", ylab="Residuals")
plot(sum_abnormal %>% dplyr::select(focal_age),
     E1, xlab="Focal Age", ylab="Residuals")

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))

plot(M1) 

#

x<-tf$col.volume.mm[!is.na(tf$col.volume.mm)]#removes na values from column
E2<-residuals(M2,type="normalized")
plot(x, E2)
#residuals are linear
