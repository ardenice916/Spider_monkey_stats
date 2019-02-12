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
install.packages("MASS")
install.packages("pscl")
install.packages("GLMMadaptive")

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

histogram(scans$n_prox_1m)

plot(scans$n_prox_1m)

ggplot(scans, aes(x = n_prox_1m)) + geom_bar() + facet_wrap(location_focal ~ n_encl)

#begin to fit models

M0_n_prox_1m<-glmer(n_prox_1m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                  (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

M1_n_prox_1m<-glmer(n_prox_1m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group), 
                  na.action=na.omit, data=scans, family = poisson)

M2_n_prox_1m<-zeroinfl(n_prox_1m ~ n_encl + location_focal + time_meal + 
                         focal_sex + focal_age + month + 1|focal_group, 
                       data=scans, na.action=na.omit)

M3_n_prox_1m<-zeroinfl(n_prox_1m ~ n_encl + location_focal + time_meal + 
                         focal_sex + focal_age + month + 1|focal_id, 
                       data=scans, na.action=na.omit)


anova(M0_n_prox_1m,M1_n_prox_1m,M2_n_prox_1m)
#can't compare across packages?

summary(M2_n_prox_1m,M3_n_prox_1m)




#M2_n_prox_1m<-zeroinfl(n_prox_1m ~ n_encl + location_focal +
#                        time_meal + focal_sex + focal_age +
#                         month +(1|focal_group/focal_id), data=scans, na.action=na.omit)
#zeroinfl() doesn't accept random factors?

#try GLMMadaptive package:

M2_n_prox_1m<- mixed_model(n_prox_1m ~ n_encl + location_focal + time_meal + focal_sex + 
                             focal_age + month, random = ~ 1 | focal_group/focal_id, 
                           zi_fixed = ~ n_encl, data = scans, family = zi.poisson())

M3_n_prox_1m<- mixed_model(n_prox_1m ~ n_encl + location_focal + time_meal + focal_sex + 
                             focal_age + month, random = ~ 1 | focal_group/focal_id, 
                           zi_fixed = ~ n_encl, data = scans, family = zi.negative.binomial())

#anova
M1_p1<-M1_n_prox_1m
M2_p1<-M2_n_prox_1m
M3_p1<-M3_n_prox_1m
M0_p1<-M0_n_prox_1m

anova(M0_p1,M1_p1)
#M0 better than M1 as expected

anova(M2_p1,M3_p1)
#ZINB M3 looks better than ZIP M2

anova(M0_p1, M3_p1)
#anova function not comparing M0 to M3? glmer and GLMMadaptive not compatible?

#M0 and M3 have lowest AIc

E0<-residuals(M0_p1, type = "pearson")
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

#those residuals look atrocious... that model clearly did not work haha

# ZERO-INFLATED... 


E3<-residuals(M3_p1)
str(E3)
E3
summary(E3)

plot(scans$n_encl, E3, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E3, xlab="Location", ylab="Residuals")
plot(scans$month, E3, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E3, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E3, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E3, xlab="Focal Age", ylab="Residuals")

qqnorm(E3)
qqline(E3)
ad.test(E3)

plot(M3_p1) 

plot(scans$n_prox_1m, E3)

#what is going on?

#let's continue on with M3 as if everything is normal...

overdisp_fun <- function(M3_p1) {
  rdf <- df.residual(M0_n_prox_1m)
  rp <- residuals(M0_n_prox_1m,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(M3_p1)

#drop factors

M3_p1<- mixed_model(n_prox_1m ~ n_encl + location_focal + time_meal + focal_sex + 
                             focal_age + month, random = ~ 1 | focal_group/focal_id, 
                           zi_fixed = ~ n_encl, data = scans, family = zi.negative.binomial())
summary(M3_p1)

M3a_p1<-mixed_model(n_prox_1m ~ n_encl + location_focal + focal_sex + focal_age + month, 
                    random = ~ 1 | focal_group/focal_id, zi_fixed = ~ n_encl, 
                    data = scans, family = zi.negative.binomial())

summary(M3a_p1)

anova(M3_p1,M3a_p1)

