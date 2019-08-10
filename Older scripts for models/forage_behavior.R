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

sum_forage<-sum_focals %>%
  filter(activity=="forage")

M0_forage<-gls(prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_forage, method="ML")

#add random factor

M1_forage<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_forage, method="ML")


#add nesting

M2_forage<-lme(prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")

#anova
anova(M0_forage, M1_forage, M2_forage)
anova(M0_forage, M2_forage)

#M0 has lowest AIC

#Analyze base model residuals

E0<-residuals(M0_forage)
str(E0)
E0
summary(E0)

plot(sum_forage$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(sum_forage$location, E0, xlab="Location", ylab="Residuals")
plot(sum_forage$month, E0, xlab="Month", ylab="Residuals")
plot(sum_forage$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(sum_forage$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(sum_forage$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_forage) 

plot(sum_forage$prop_time, E0)

#log normalize proportions
sum_forage$l.prop_time=log(sum_forage$prop_time+1)

#repeat model attempts with normalized data

M0_forage<-gls(l.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_forage, method="ML")

M1_forage<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_forage, method="ML")

M2_forage<-lme(l.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")

anova(M0_forage, M1_forage, M2_forage)
anova(M0_forage, M2_forage)

#M0 still lowest AIC, better fit than before

E0<-residuals(M0_forage)
str(E0)
E0
summary(E0)

plot(sum_forage$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(sum_forage$location, E0, xlab="Location", ylab="Residuals")
plot(sum_forage$month, E0, xlab="Month", ylab="Residuals")
plot(sum_forage$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(sum_forage$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(sum_forage$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_forage) 

plot(sum_forage$l.prop_time, E0)

#arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }

#execute arcsine transformation
sum_forage$as.prop_time=asinTransform(sum_forage$prop_time)

#repeat model attempts with arcsine normalized data

M0_forage<-gls(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_forage, method="ML")

M1_forage<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|focal_group, na.action=na.omit, data=sum_forage, method="ML")

M2_forage<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
                 random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")

anova(M0_forage, M1_forage, M2_forage)
anova(M0_forage, M2_forage)
anova(M0_forage, M1_forage)

#M0 still lowest, not sig diff

E0<-residuals(M0_forage)
str(E0)
E0
summary(E0)

plot(sum_forage$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(sum_forage$location, E0, xlab="Location", ylab="Residuals")
plot(sum_forage$month, E0, xlab="Month", ylab="Residuals")
plot(sum_forage$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(sum_forage$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(sum_forage$focal_age, E0, xlab="Focal Age", ylab="Residuals")

#bartlett test compares variance b/w groups of residuals
#you can use this to see if variance differs by category
#in these residual tests

bartlett.test(E0~sum_forage$month)

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_forage) 

plot(sum_forage$as.prop_time, E0)

#check for autocorrelation
acf(E0, na.action=na.pass,
    main="Auto-correlation plot for residuals")
#little autocorrelation in the first lag

#arcsine looks best... M0 best fitting technically, but M2 would be acceptable too

#use M2 as new full model

M.full<-lme(as.prop_time ~ n_encl + location + month + time_meal + focal_sex + focal_age, 
            random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")
summary(M.full)

#drop month

M.full.a<-lme(as.prop_time ~ n_encl + location + time_meal + focal_sex + focal_age, 
              random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")
anova(M.full,M.full.a)

#much better... 

summary(M.full.a)

#drop n_encl

M.full.b<-lme(as.prop_time ~ location + time_meal + focal_age + focal_sex, 
              random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")
anova(M.full.a,M.full.b)

#somewhat better, not significant

summary(M.full.b)

#drop focal_sex

M.full.c<-lme(as.prop_time ~ location + time_meal + focal_age,
              random = ~1|nest, na.action=na.omit, data=sum_forage, method="ML")
anova(M.full.b,M.full.c)

#somewhat better, albeit not significantly

summary(M.full.c)

#best overall model was arcsine transformed
M.full<-lme(as.prop_time ~ location + time_meal + focal_age,
            random = ~1|nest, na.action=na.omit, data=sum_forage, method="REML")
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
