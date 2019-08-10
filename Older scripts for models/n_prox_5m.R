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

plot(scans$n_prox_5m)

ggplot(scans, aes(x = n_prox_5m)) + geom_bar() + facet_wrap(location_focal ~ n_encl)

#begin to fit models

M0_n_prox_5m<-glmer(n_prox_5m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                  (1|focal_group/focal_id), na.action=na.omit, data=scans, family = poisson)

M1_n_prox_5m<-glmer(n_prox_5m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group), 
                  na.action=na.omit, data=scans, family = poisson)

#try a log transformation?

scans$l.n_prox_5m<-log(scans$n_prox_5m+1)

histogram(scans$l.n_prox_5m)

M3_n_prox_5m<-lmer(l.n_prox_5m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                     (1|focal_group/focal_id), na.action=na.omit, data=scans)

M0_n_prox_5m
summary(M0_n_prox_5m)
anova(M0_n_prox_5m)

#anova
anova(M0_n_prox_5m, M1_n_prox_5m, M3_n_prox_5m)

anova(M0_n_prox_5m, M1_n_prox_5m, test = "Chi")

#M3 has lowest AIc, but might not work. M0 might be preferable

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

#examine residuals of M3

E3<-residuals(M3_n_prox_5m, type = "pearson")
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

plot(M3_n_prox_5m) 

plot(scans$n_prox_5m, E3)

#test for overdispersion

overdisp_fun <- function(M0_n_prox_5m) {
  rdf <- df.residual(M0_n_prox_5m)
  rp <- residuals(M0_n_prox_5m,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(M0_n_prox_5m)

overdisp_fun <- function(M3_n_prox_5m) {
  rdf <- df.residual(M3_n_prox_5m)
  rp <- residuals(M3_n_prox_5m,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(M3_n_prox_5m)

#not significantly lower than 1...?

#use M3, rename M0_p5
#drop factors

M0_p5<-lmer(l.n_prox_5m ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                     (1|focal_group/focal_id), na.action=na.omit, data=scans)

summary(M0_p5)

#drop meal status

M1_p5<-lmer(l.n_prox_5m ~ n_encl + location_focal + focal_sex + focal_age + month +
              (1|focal_group/focal_id), na.action=na.omit, data=scans)

anova(M0_p5, M1_p5)

#no significant change

summary(M1_p5)

#drop focal age

M2_p5<-lmer(l.n_prox_5m ~ n_encl + location_focal + focal_sex + month +
              (1|focal_group/focal_id), na.action=na.omit, data=scans)

anova(M1_p5, M2_p5)

#no significant change

summary(M2_p5)

#all remaining factors appear somewhat significant

M.full<-M2_p5

anova(M.full)

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

multCompTukey1 <- glht(M.full, linfct = mcp(n_encl = "Tukey")) 
summary(multCompTukey1)

multCompTukey2 <- glht(M.full, linfct = mcp(location_focal = "Tukey")) 
summary(multCompTukey2)

multCompTukey3 <- glht(M.full, linfct = mcp(focal_sex = "Tukey")) 
summary(multCompTukey3)

multCompTukey4 <- glht(M.full, linfct = mcp(month= "Tukey")) 
summary(multCompTukey4)

#examine relationships visually

plot(scans$n_prox_5m ~ scans$month)
plot(scans$n_prox_5m ~ scans$n_encl)
plot(scans$n_prox_5m ~ scans$location_focal)
plot(scans$n_prox_5m ~ scans$focal_sex)

ggplot(scans, aes(x = n_prox_5m)) + geom_bar() + facet_wrap(focal_group ~ focal_id)

scans_jun_jul <- scans %>%
  filter(month=="jun_jul")
summary(scans_jun_jul$n_prox_5m)

scans_jul_aug <- scans %>%
  filter(month=="jul_aug")
summary(scans_jul_aug$n_prox_5m)

scans_aug_sep <- scans %>%
  filter(month=="aug_sep")
summary(scans_aug_sep$n_prox_5m)

scans_male <- scans %>%
  filter(focal_sex=="male")
summary(scans_male$n_prox_5m)

scans_female <- scans %>%
  filter(focal_sex=="female")
summary(scans_female$n_prox_5m)

