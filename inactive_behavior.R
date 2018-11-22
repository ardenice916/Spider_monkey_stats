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


#Look at mixed effects model
#start without random factor

sum_inactive<-sum_focals %>%
  filter(activity=="inactive")

M0_inactive<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_inactive, method="ML")

#add random factor

M1_inactive<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_inactive, method="ML")


#add nesting

M2_inactive<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

#anova
anova(M0_inactive, M1_inactive, M2_inactive)
anova(M0_inactive, M2_inactive)
anova(M1_inactive, M2_inactive)

#M2 has lowest AIC

#Analyze base model residuals

E2<-residuals(M2_inactive)
str(E2)
E2
summary(E2)

plot(sum_inactive$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_inactive$location, E2, xlab="Location", ylab="Residuals")
plot(sum_inactive$month, E2, xlab="Month", ylab="Residuals")
plot(sum_inactive$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_inactive$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_inactive$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_inactive) 

plot(sum_inactive$prop_time, E2)

#log normalize proportions
sum_inactive$l.prop_time=log(sum_inactive$prop_time+1)

#repeat model attempts with normalized data

M0_inactive<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_inactive, method="ML")

M1_inactive<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_inactive, method="ML")

M2_inactive<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

anova(M0_inactive, M1_inactive, M2_inactive)
anova(M1_inactive, M2_inactive)

#M2 still lowest AIC

E2<-residuals(M2_inactive)
str(E2)
E2
summary(E2)

plot(sum_inactive$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_inactive$location, E2, xlab="Location", ylab="Residuals")
plot(sum_inactive$month, E2, xlab="Month", ylab="Residuals")
plot(sum_inactive$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_inactive$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_inactive$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_inactive) 

plot(sum_inactive$l.prop_time, E0)

#arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation
sum_inactive$as.prop_time=asinTransform(sum_inactive$prop_time)

#repeat model attempts with arcsine normalized data

M0_inactive<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_inactive, method="ML")

M1_inactive<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_inactive, method="ML")

M2_inactive<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

anova(M0_inactive, M1_inactive, M2_inactive)
anova(M0_inactive, M2_inactive)
anova(M0_inactive, M1_inactive)

#M2 still lowest

E2<-residuals(M2_inactive)
str(E2)
E2
summary(E2)

plot(sum_inactive$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(sum_inactive$location, E2, xlab="Location", ylab="Residuals")
plot(sum_inactive$month, E2, xlab="Month", ylab="Residuals")
plot(sum_inactive$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_inactive$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_inactive$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_inactive) 

plot(sum_inactive$as.prop_time, E2)

#check for autocorrelation
acf(E2, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#potential autocorrelation in the first lag 

#unclear if transformations helped... M2 looks best regardless

#make M2 the new full model

M2.full<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
             random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

summary(M2.full)

#drop n_encl

M2.full.a<-lme(as.prop_time ~ location + month + time_meal + focal_sex + focal_age, 
             random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

anova(M2.full,M2.full.a)

#improved insignificantly

summary(M2.full.a)

#drop month

M2.full.b<-lme(as.prop_time ~ location + time_meal + focal_sex + focal_age, 
               random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

anova(M2.full.a,M2.full.b)

#improved insignificantly

summary(M2.full.b)

#drop focal sex

M2.full.c<-lme(as.prop_time ~ location + time_meal + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

anova(M2.full.b,M2.full.c)

#improved insignificantly

summary(M2.full.c)

#drop focal age

M2.full.d<-lme(as.prop_time ~ location + time_meal, 
                 random = ~1|nest, na.action=na.omit, data=sum_inactive, method="ML")

anova(M2.full.c,M2.full.d)

#got significantly worse! M2.full.c (w/ focal age) is best fit

#best overall model was arcsine transformed
M.full<-lme(as.prop_time ~ location + time_meal + focal_age, 
            random = ~1|nest, na.action=na.omit, data=sum_inactive, method="REML")
#need to change method to REML or drop the method term since it defaults to REML

anova(M.full)#this is the summary of your full model

#####################################################
#Get Full Model Statistics and Make Graph
#####################################################

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

multCompTukey2 <- glht(M.full, linfct = mcp(time_meal = "Tukey")) 
summary(multCompTukey2)

multCompTukey3 <- glht(M.full, linfct = mcp(focal_age = "Tukey")) 
summary(multCompTukey3)
