#COMPLETE ANALYSIS OF 2016 WILDTRACKS DATA, FEBRUARY 2018

#set working directory

setwd("~/Spider_monkey_stats")

#install packages

install.packages("TMB", type = 'source')
install.packages("glmmTMB")
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
install.packages("MASS")
install.packages("pscl")
install.packages("GLMMadaptive")

#load useful packages

library(TMB)
library(glmmTMB)
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
library(MASS)
library(pscl)
library(GLMMadaptive)

####Enter Focal Data Set########################################

focals <- read.table(file="focals_Wildtracks.csv", sep=",", header=T)
focals$n_encl<-factor(focals$n_encl, levels=c("1","2"))
focals$month<-factor(focals$month, levels=c("jun-jul","jul-aug","aug-sep"))
focals$time_meal<-factor(focals$time_meal, levels=c("before","after","none"))
focals$dur_m <- focals$duration_sec/60
focals$dur_h <- focals$dur_m/60

sum_focals1<-focals %>%
  dplyr::filter(actor_id!="FG") %>% 
  dplyr::filter(actor_id!="HUMAN")

sum_focals1$actor_id<-factor(sum_focals1$actor_id, levels=c("CL","DU","FY","MA","ME","PA","PE","PO","PP","RK","TR"))

#1: Enclosure Preferences######################################################

sum_cp<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  dplyr::filter(n_encl=="2") %>%
  group_by(focal_id, focal_sex, focal_group, 
           focal_age, time_meal, month, location) %>% 
  complete(location, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur)) %>%
  dplyr::filter(location=="center")

summary(sum_cp$prop_time)

histogram(sum_cp$prop_time)

plot(sum_cp$prop_time)

plot(sum_cp$prop_time ~ sum_cp$focal_id)
plot(sum_cp$prop_time ~ sum_cp$focal_group)
plot(sum_cp$prop_time ~ sum_cp$focal_age)
plot(sum_cp$prop_time ~ sum_cp$focal_sex)
plot(sum_cp$prop_time ~ sum_cp$month)

M_cp1<-lme(prop_time ~ time_meal + focal_sex + focal_age, 
               random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

Ecp<-residuals(M_cp)
str(Ecp)
Ecp
summary(Ecp)

plot(sum_cp$month, Ecp, xlab="Month", ylab="Residuals")
plot(sum_cp$time_meal, Ecp, xlab="Meal Status", ylab="Residuals")
plot(sum_cp$focal_sex, Ecp, xlab="Focal Sex", ylab="Residuals")
plot(sum_cp$focal_age, Ecp, xlab="Focal Age", ylab="Residuals")

qqnorm(Ecp)
qqline(Ecp)
ad.test(Ecp)

plot(M_cp) 

plot(sum_cp$prop_time, Ecp)

M_cp2<-lme(prop_time ~ focal_age, 
        random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

anova(M_cp2)

model.matrix.gls <- function(M_cp2, ...){
  model.matrix(terms(M_cp2), data = getData(M_cp2), ...)  
}
model.frame.gls <- function(M_cp2, ...){
  model.frame(formula(M_cp2), data = getData(M_cp2), ...)  
}
terms.gls <- function(M_cp2, ...){
  terms(model.frame(M_cp2),...)  
}

multCompTukey1 <- glht(M_cp2, linfct = mcp(focal_age = "Tukey")) 
summary(multCompTukey1)

plot_cp_id<-ggplot(sum_cp, aes(focal_id, prop_time)) + geom_boxplot() +
  ggtitle("Time Spent in the Center Enclosure (by Choice)") + 
  labs(x="Focal Individual", y="Proportion of Time") +
  theme_classic()
plot_cp_id

plot_cp_age<-ggplot(sum_cp, aes(focal_age, prop_time)) + geom_boxplot() +
  ggtitle("Time Spent in the Center Enclosure (by Choice)") +
  xlab("Focal Age") + ylab("Proportion of Time") + theme_classic() +
  scale_x_discrete(labels = c('Adults','Subadults'))
plot_cp_age

#2: Activity Budgets #################################################

#2A: Check which activities were most prevalent

sum_overall<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  group_by(activity) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))

sum_overall

#drop agonistic, caregiver, other, parental, prosocial
#drop sociosexual, solitary play, vocalization, not visible
#keep/analyze abnormal, forage, inactive, movement, 
#keep/analyze self_directed/scratching

#2B: Summarize focal data set

sum_focals<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  dplyr::select(duration_sec, occurrences, month, focal_id, time_meal, n_encl, location, activity, focal_sex, focal_group, focal_age) %>%
  group_by(focal_id, focal_sex, focal_group, focal_age, month, time_meal, n_encl, location, activity) %>%
  complete(activity, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))

sum_focals$nest <- with(sum_focals, factor(paste(focal_group,focal_id)))

#2A: Activity - ABNORMAL

sum_abnormal<-sum_focals %>%
  filter(activity=="abnormal")

asinTransform <- function(p) { asin(sqrt(p)) }

sum_abnormal$as.prop_time=asinTransform(sum_abnormal$prop_time)

M_ab1<-lme(as.prop_time ~ n_encl + location + month +
              time_meal + focal_sex + focal_age, 
            random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="ML")

anova(M_ab1)
summary(M_ab1)

M_ab2<-lme(as.prop_time ~ location + month + time_meal, 
            random = ~1|nest, na.action=na.omit, data=sum_abnormal, method="REML")

model.matrix.gls <- function(M_ab2, ...){
  model.matrix(terms(M_ab2), data = getData(M_ab2), ...)  
}
model.frame.gls <- function(M_ab2, ...){
  model.frame(formula(M_ab2), data = getData(M_ab2), ...)  
}
terms.gls <- function(M_ab2, ...){
  terms(model.frame(M_ab2),...)  
}

multCompTukey1 <- glht(M_ab2, linfct = mcp(location = "Tukey")) 
summary(multCompTukey1)

multCompTukey2 <- glht(M_ab2, linfct = mcp(month = "Tukey")) 
summary(multCompTukey2)

multCompTukey3 <- glht(M_ab2, linfct = mcp(time_meal = "Tukey")) 
summary(multCompTukey3)

plot_ab_id<-ggplot(sum_abnormal, aes(focal_id, prop_time)) + geom_boxplot() +
  ggtitle("Time Engaged in Abnormal Behavior") + 
  labs(x="Focal Individual", y="Proportion of Time") +
  theme_classic()
plot_ab_id

plot_ab_loc<-ggplot(sum_abnormal, aes(location, prop_time)) + geom_boxplot() +
  ggtitle("Time Engaged in Abnormal Behavior by Location") + 
  labs(x="Location", y="Proportion of Time") +
  theme_classic() +
  scale_x_discrete(labels = c('Center Enclosure','Satellite Enclosures'))
plot_ab_loc

plot_ab_mon<-ggplot(sum_abnormal, aes(month, prop_time)) + geom_boxplot() +
  ggtitle("Time Engaged in Abnormal Behavior by Month") + 
  labs(x="Month", y="Proportion of Time") +
  theme_classic() +
  scale_x_discrete(labels = c('June-July','July-August','August-September'))
plot_ab_mon

plot_ab_ms<-ggplot(sum_abnormal, aes(time_meal, prop_time)) + geom_boxplot() +
  ggtitle("Time Engaged in Abnormal Behavior by Meal Status") + 
  labs(x="Meal Status", y="Proportion of Time") +
  theme_classic() +
plot_ab_ms

#: MODELS FOR SCAN DATA

scans<-read.table(file="scans_Wildtracks.csv", sep=",", header=T)#alternate import code
scans$n_encl<-factor(scans$n_encl, levels=c("1","2"))
scans$location_focal<-factor(scans$location_focal, levels=c("satellite","center"))
scans$month<-factor(scans$month, levels=c("jun_jul","jul_aug","aug_sep"))
scans$time_meal<-factor(scans$time_meal, levels=c("none","before","after"))
scans$observer_id<-factor(scans$observer_id, levels=c("916","1231"))

