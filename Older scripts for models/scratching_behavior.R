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

#view and summarize

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
  dplyr::select(duration_sec, occurrences, month, focal_id, time_meal, n_encl, pref_encl, location, behavior, focal_sex, focal_group, focal_age) %>%
  group_by(focal_id, focal_sex, focal_group, focal_age, month, time_meal, n_encl, pref_encl, location, behavior) %>%
  complete(behavior, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))

sum_focals$nest <- with(sum_focals, factor(paste(focal_group,focal_id)))

#Look at mixed effects model
#start without random factor

sum_scratch<-sum_focals %>%
  filter(behavior=="sdb_scratching")

M0_scratch<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                      na.action=na.omit, data=sum_scratch, method="ML")

#add random factor

M1_scratch<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                      random = ~1|focal_group, na.action=na.omit, data=sum_scratch, method="ML")


#add nesting

M2_scratch<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                      random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

#anova
anova(M0_scratch, M1_scratch, M2_scratch)
anova(M0_scratch, M2_scratch)
anova(M1_scratch, M2_scratch)

#M2 has lowest AIC

#Analyze base model residuals

E2<-residuals(M2_scratch)
str(E2)
E2
summary(E2)

plot(sum_scratch$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_scratch$location, E2, xlab="Location", ylab="Residuals")
plot(sum_scratch$month, E2, xlab="Month", ylab="Residuals")
plot(sum_scratch$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_scratch$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_scratch$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_scratch) 

plot(sum_scratch$prop_time, E2)

#log normalize proportions
sum_scratch$l.prop_time=log(sum_scratch$prop_time+1)

#repeat model attempts with normalized data

M0_scratch<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                      na.action=na.omit, data=sum_scratch, method="ML")

M1_scratch<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                      random = ~1|focal_group, na.action=na.omit, data=sum_scratch, method="ML")

M2_scratch<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                      random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

anova(M0_scratch, M1_scratch, M2_scratch)
anova(M1_scratch, M2_scratch)

#M2 still lowest AIC

E2<-residuals(M2_scratch)
str(E2)
E2
summary(E2)

plot(sum_scratch$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_scratch$location, E2, xlab="Location", ylab="Residuals")
plot(sum_scratch$month, E2, xlab="Month", ylab="Residuals")
plot(sum_scratch$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_scratch$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_scratch$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_scratch) 

plot(sum_scratch$prop_time, E2)

#arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation
sum_scratch$as.prop_time=asinTransform(sum_scratch$prop_time)

#repeat model attempts with arcsine normalized data

M0_scratch<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                      na.action=na.omit, data=sum_scratch, method="ML")

M1_scratch<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                      random = ~1|focal_group, na.action=na.omit, data=sum_scratch, method="ML")

M2_scratch<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                      random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

anova(M0_scratch, M1_scratch, M2_scratch)
anova(M0_scratch, M2_scratch)
anova(M0_scratch, M1_scratch)
anova(M1_scratch, M2_scratch)

#M1 now best, and significantly better than M0 but not M2

E2<-residuals(M2_scratch)
str(E2)
E2
summary(E2)

plot(sum_scratch$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_scratch$location, E2, xlab="Location", ylab="Residuals")
plot(sum_scratch$month, E2, xlab="Month", ylab="Residuals")
plot(sum_scratch$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_scratch$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_scratch$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_scratch) 

plot(sum_scratch$prop_time, E2)

#check for autocorrelation
acf(E2, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#no autocorrelation except 1st lag

#M1 w/ arcsine transformed data looks slightly better?
#M2 not significantly different, not a big change in residuals

#make M2 w/ arcsine sqrt transformation the new full model

M2.full<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
             random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")
summary(M2.full)

#drop month

M2.full.a<-lme(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
             random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

anova(M2.full,M2.full.a)

#improved insignificantly

summary(M2.full.a)

#drop focal age

M2.full.b<-lme(as.prop_time ~ n_encl + location + time_meal + focal_sex, 
               random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

anova(M2.full.a,M2.full.b)

#improved insignificantly

summary(M2.full.b)

#drop focal_sex

M2.full.c<-lme(as.prop_time ~ n_encl + location + time_meal, 
               random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

anova(M2.full.b,M2.full.c)

#improved insignificantly

summary(M2.full.c)

#drop time_meal

M2.full.d<-lme(as.prop_time ~ n_encl + location, 
               random = ~1|nest, na.action=na.omit, data=sum_scratch, method="ML")

anova(M2.full.c,M2.full.d)

#got a little bit better, insignificantly

summary(M2.full.d)

#best overall model was arcsine transformed
M.full<-lme(as.prop_time ~ n_encl + location, 
            random = ~1|nest, na.action=na.omit, data=sum_scratch, method="REML")
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

multCompTukey2 <- glht(M.full, linfct = mcp(n_encl = "Tukey")) 
summary(multCompTukey2)
