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

#import code

scans <- read.table(file="scans_Wildtracks.csv", sep=",", header=T)#alternate import code
scans$n_encl<-factor(scans$n_encl, levels=c("1","2"))
scans$location_focal<-factor(scans$location_focal, levels=c("satellite","center"))
scans$month<-factor(scans$month, levels=c("jun_jul","jul_aug","aug_sep"))
scans$time_meal<-factor(scans$time_meal, levels=c("none","before","after"))
scans$observer_id<-factor(scans$observer_id, levels=c("916","1231"))

scans <- scans %>% filter(n_prox_1m!="NA")

#view and summarize

View(scans)
summary(scans)

scans$nest <- with(scans, factor(paste(focal_group,focal_id)))

#check n_prox_1m variable
str(scans)
str(scans$n_prox_1m)
summary(scans$n_prox_1m)
print(scans$n_prox_1m)

x=sum(scans$scans)
x

y=sum(scans$scans[scans$n_prox_1m==0])
y

y/x

#no individuals in proximity 93% of scans

histogram(scans$n_prox_1m)

plot(scans$n_prox_1m)

ggplot(scans, aes(x = n_prox_1m)) + geom_bar() + facet_wrap(location_focal ~ n_encl)
ggplot(scans, aes(x = n_prox_1m)) + geom_bar() + facet_wrap(focal_group ~ time_meal)
ggplot(scans, aes(x = n_prox_1m)) + geom_bar() + facet_wrap(focal_group ~ month)

#try with glmmTMB

M0_p1<-glmmTMB(n_prox_1m ~ n_encl + location_focal + 
                 time_meal + focal_sex + focal_age + month + 
                 (1|focal_group/focal_id),data= scans, ziformula=~1, family = poisson)

M1_p1<-glmmTMB(n_prox_1m ~ n_encl + location_focal + 
                 time_meal + focal_sex + focal_age + month + 
                 (1|focal_group/focal_id),data= scans, ziformula=~1, family = nbinom1(link="log"))

M2_p1<-glmmTMB(n_prox_1m ~ n_encl + location_focal + 
                 time_meal + focal_sex + focal_age + month + 
                 (1|focal_group/focal_id),data= scans, ziformula=~1, family = nbinom2(link="log"))

anova(M0_p1,M1_p1,M2_p1)

#ZINB2 M2 looks better than ZIP M0 or other ZINB1 M1
#use M2_p1

#examine residuals

E2<-residuals(M2_p1)
str(E2)
E2
summary(E2)

plot(scans$n_encl, E2, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E2, xlab="Location", ylab="Residuals")
plot(scans$month, E2, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E2, xlab="Focal Age", ylab="Residuals")

qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_p1) 

plot(scans$n_prox_1m, E2)

#those residuals look atrocious... that model clearly did not work haha

#Let's examine the others?

E0<-residuals(M0_p1)
str(E0)
E0
summary(E0)

plot(scans$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E0, xlab="Location", ylab="Residuals")
plot(scans$month, E0, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_p1) 

plot(scans$n_prox_1m, E0)

#what is going on?

#########repeat with log-transformed data#########????

scans$n_prox_1m<-log(scans$n_prox_1m+1)

summary(scans$n_prox_1m)
print(scans$n_prox_1m)

histogram(scans$n_prox_1m)

plot(scans$n_prox_1m)

ggplot(scans, aes(x = n_prox_1m)) + geom_bar() + facet_wrap(location_focal ~ n_encl)

#try with glmmTMB

M0_p1<-glmmTMB(n_prox_1m ~ n_encl + location_focal + 
                 time_meal + focal_sex + focal_age + month + 
                 (1|focal_group/focal_id),data= scans, ziformula=~1, family = poisson)

M1_p1<-glmmTMB(n_prox_1m ~ n_encl + location_focal + 
                 time_meal + focal_sex + focal_age + month + 
                 (1|focal_group/focal_id),data= scans, ziformula=~1, family = nbinom1(link="log"))

M2_p1<-glmmTMB(n_prox_1m ~ n_encl + location_focal + 
                 time_meal + focal_sex + focal_age + month + 
                 (1|focal_group/focal_id),data= scans, ziformula=~1, family = nbinom2(link="log"))

anova(M0_p1,M1_p1,M2_p1)

#Only ZIP (M0) worked

#Let's examine the residuals

E0<-residuals(M0_p1)
str(E0)
E0
summary(E0)

plot(scans$n_encl, E0, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E0, xlab="Location", ylab="Residuals")
plot(scans$month, E0, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E0, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E0, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E0, xlab="Focal Age", ylab="Residuals")

qqnorm(E0)
qqline(E0)
ad.test(E0)

plot(M0_p1) 

plot(scans$n_prox_1m, E0)

#residuals look terrible...