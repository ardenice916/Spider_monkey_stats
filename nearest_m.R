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

scans<-read.table(file="scans_Wildtracks.csv", sep=",", header=T)#alternate import code
scans$n_encl<-factor(scans$n_encl, levels=c("1","2"))
scans$month<-factor(scans$month, levels=c("jun_jul","jul_aug","aug_sep"))
scans$time_meal<-factor(scans$time_meal, levels=c("none","before","after"))
scans$observer_id<-factor(scans$observer_id, levels=c("916","1231"))
scans$nearest_m<-as.numeric(scans$nearest_m)

scans<-scans %>%
  filter(nearest_m != "NA")

#view and summarize

View(scans)
summary(scans)

scans$nest <- with(scans, factor(paste(focal_group,focal_id)))

#check nearest_m variable
str(scans)
str(scans$nearest_m)
summary(scans$nearest_m)

histogram(scans$nearest_m)

M0_nearest_m<-glmer(nearest_m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                  (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

M1_nearest_m<-glmer(nearest_m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group), 
                  na.action=na.omit, data=scans, family = poisson)

M2_nearest_m<-lmer(nearest_m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                      (1|focal_group/focal_id), na.action=na.omit, data=scans)

scans$l.nearest_m=log(scans$nearest_m+1)

histogram(scans$l.nearest_m)

M3_nearest_m<-lmer(l.nearest_m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                     (1|focal_group), na.action=na.omit, data=scans)

#anova
anova(M0_nearest_m, M1_nearest_m, M2_nearest_m, M3_nearest_m)

anova(M0_nearest_m, M1_nearest_m, test = "Chi")

#M0 and M1 have similar AIc, M2 is worse

#M0 poisson residuals

E0<-residuals(M0_nearest_m, type = "pearson")
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

plot(M0_nearest_m) 

plot(scans$nearest_m, E0)

#residuals M2 gaussian

E2<-residuals(M2_nearest_m, type = "pearson")
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

plot(M2_nearest_m) 

plot(scans$nearest_m, E2)

#residuals M3 gaussian

E3<-residuals(M3_nearest_m, type = "pearson")
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

plot(M3_nearest_m) 

plot(scans$nearest_m, E3)

# test for overdispersion? quasipoisson does not seem to be an option in lme4 package anymore
#coding for random factors didn't work unless I used lme4

overdisp_fun <- function(M0_nearest_m) {
  rdf <- df.residual(M0_nearest_m)
  rp <- residuals(M0_nearest_m,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(M0_nearest_m)

#no over- or under-dispersion, but the residuals look like they fit terribly... maybe we need to transform?

#use M0_nearest_m

M0_nm<-glmer(nearest_m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                      (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

anova(M0_nearest_m)

M1_nm<-glmer(nearest_m ~ n_encl + location_focal + focal_sex + focal_age + month +
               (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

M.full<-M0_nm

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

multCompTukey2 <- glht(M.full, linfct = mcp(month = "Tukey")) 
summary(multCompTukey2)

multCompTukey3 <- glht(M.full, linfct = mcp(time_meal = "Tukey")) 
summary(multCompTukey3)