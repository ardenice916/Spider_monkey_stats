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

#import code

scans <- read.table(file="scans_Wildtracks.csv", sep=",", header=T)#alternate import code
scans$n_encl<-factor(scans$n_encl, levels=c("1","2"))
scans$month<-factor(scans$month, levels=c("jun_jul","jul_aug","aug_sep"))
scans$time_meal<-factor(scans$time_meal, levels=c("none","before","after"))
scans$observer_id<-factor(scans$observer_id, levels=c("916","1231"))
scans$n_prox_5m<-as.numeric(scans$n_prox_5m)

scans <- scans %>% filter(n_prox_5m!="NA")

#view and summarize

View(scans)
summary(scans)

scans$nest <- with(scans, factor(paste(focal_group,focal_id)))

#check n_prox_5m variable
str(scans)
str(scans$n_prox_5m)
summary(scans$n_prox_5m)

histogram(scans$n_prox_5m)


M0_n_prox_5m<-glmer(n_prox_5m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                  (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

M1_n_prox_5m<-glmer(n_prox_5m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group), 
                  na.action=na.omit, data=scans, family = poisson)

M0_n_prox_5m
summary(M0_n_prox_5m)

#anova
anova(M0_n_prox_5m, M1_n_prox_5m)

anova(M0_n_prox_5m, M1_n_prox_5m, test = "Chi")

#M0 has lowest AIc

E0<-residuals(M0_n_prox_5m, type = "pearson")
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

plot(M0_n_prox_5m) 

plot(scans$n_prox_5m, E0)

# not too bad... can we try a transformation?


overdisp_fun <- function(M0_n_prox_5m) {
  rdf <- df.residual(M0_n_prox_5m)
  rp <- residuals(M0_n_prox_5m,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(M0_n_prox_5m)

#significantly lower than 1... underdispersed?