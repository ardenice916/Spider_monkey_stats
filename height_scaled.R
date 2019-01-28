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

#view and summarize

View(scans)
summary(scans)

scans$nest <- with(scans, factor(paste(focal_group,focal_id)))

#check height_scaled variable
str(scans)
str(scans$height_scaled)
summary(scans$height_scaled)


M0_h_scaled<-glmer(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                  (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

M1_h_scaled<-glmer(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group), 
                  na.action=na.omit, data=scans, family = poisson)



M0_h_scaled
summary(M0_h_scaled)

#anova
anova(M0_h_scaled, M1_h_scaled)

anova(M0_h_scaled, M1_h_scaled, test = "Chi")

#M0 has lowest AIc

E0<-residuals(M0_h_scaled, type = "pearson")
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

plot(M0_h_scaled) 

plot(scans$height_scaled, E0)

# test for overdispersion? quasipoisson does not seem to be an option in lme4 package anymore
#coding for random factors didn't work unless I used lme4
#can I use a negative binomial distribution to correct for potential overdispersion instead?


overdisp_fun <- function(M0_h_scaled) {
  rdf <- df.residual(M0_h_scaled)
  rp <- residuals(M0_h_scaled,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(M0_h_scaled)

#underdispersion? ratio is 0.19..., significantly different from 1

# "reasons for underdispersion are the model is fitting a couple of outliers rather too well 
# or there are too many explanatory variables or interactions in the model.
# If this is not the case, then the consensus is NOT to correct for underdispersion"

#so.... do I correct?

