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

scans <- scans %>% filter(height_scaled!="NA")

#view and summarize

View(scans)
summary(scans)

scans$nest <- with(scans, factor(paste(focal_group,focal_id)))

#check height_scaled variable
str(scans)
str(scans$height_scaled)
summary(scans$height_scaled)
print(scans$height_scaled)

histogram(scans$height_scaled)

plot(scans$height_scaled)

ggplot(scans, aes(x = height_scaled)) + geom_bar() + facet_wrap(focal_group ~ n_encl)
ggplot(scans, aes(x = height_scaled)) + geom_bar() + facet_wrap(focal_group ~ time_meal)
ggplot(scans, aes(x = height_scaled)) + geom_bar() + facet_wrap(focal_group ~ month)

#across conditions and categories, monkeys preferred upper half of enclosures?

M0_h_scaled<-lmer(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                  (1|focal_group/focal_id), na.action=na.omit, data=scans)

M1_h_scaled<-glmer(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group/focal_id), 
                  na.action=na.omit, data=scans, family = poisson)

M2_h_scaled<-lmer(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                    (1|focal_group), na.action=na.omit, data=scans)

M3_h_scaled<-glmer(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month +
                     (1|focal_group), 
                   na.action=na.omit, data=scans, family = poisson)

M4_h_scaled<-glm(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month, na.action=na.omit, data=scans)

M5_h_scaled<-glm(height_scaled ~ n_encl + location_focal + time_meal + focal_sex + focal_age + month,
                   na.action=na.omit, data=scans, family = poisson)

M6_h_scaled<-glm(height_scaled ~ n_encl + time_meal + month + focal_sex + focal_age, na.action=na.omit, data=scans)

M7_h_scaled<-lmer(height_scaled ~ n_encl + location_focal + 
                    time_meal + month + focal_sex + focal_age + (1|focal_id), 
                  na.action=na.omit, data=scans)


#anova
anova(M0_h_scaled, M1_h_scaled, M2_h_scaled, M3_h_scaled, M4_h_scaled, M5_h_scaled, M6_h_scaled)

#M0 has lowest AIc but only M4, M5, M6 were not overfitted

#let's examine some residuals...

E0<-residuals(M0_h_scaled)
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

E4<-residuals(M4_h_scaled)
str(E4)
E4
summary(E4)

plot(scans$n_encl, E4, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E4, xlab="Location", ylab="Residuals")
plot(scans$month, E4, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E4, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E4, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E4, xlab="Focal Age", ylab="Residuals")

qqnorm(E4)
qqline(E4)
ad.test(E4)

plot(M4_h_scaled) 

plot(scans$height_scaled, E4)

E5<-residuals(M5_h_scaled)
str(E5)
E5
summary(E5)

plot(scans$n_encl, E5, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E5, xlab="Location", ylab="Residuals")
plot(scans$month, E5, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E5, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E5, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E5, xlab="Focal Age", ylab="Residuals")

qqnorm(E5)
qqline(E5)
ad.test(E5)

plot(M5_h_scaled) 

plot(scans$height_scaled, E5)

#residuals and QQ lines look bad... data looks heteroscedastic...

E6<-residuals(M6_h_scaled)
str(E6)
E6
summary(E6)

plot(scans$n_encl, E6, xlab="# Enclosures", ylab="Residuals")
plot(scans$month, E6, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E6, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E6, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E6, xlab="Focal Age", ylab="Residuals")

qqnorm(E6)
qqline(E6)
ad.test(E6)

plot(M6_h_scaled) 

plot(scans$height_scaled, E6)

E7<-residuals(M7_h_scaled)
str(E7)
E7
summary(E7)

plot(scans$n_encl, E7, xlab="# Enclosures", ylab="Residuals")
plot(scans$location_focal, E7, xlab="# Enclosures", ylab="Residuals")
plot(scans$month, E7, xlab="Month", ylab="Residuals")
plot(scans$time_meal, E7, xlab="Meal Status", ylab="Residuals")
plot(scans$focal_sex, E7, xlab="Focal Sex", ylab="Residuals")
plot(scans$focal_age, E7, xlab="Focal Age", ylab="Residuals")
plot(scans$focal_id, E7, xlab="Focal ID", ylab="Residuals")

qqnorm(E7)
qqline(E7)
ad.test(E7)

plot(M7_h_scaled) 

plot(scans$height_scaled, E7)

#let's try a transformation?########################3

scans$l.height_scaled<-log(scans$height_scaled+1)
summary(scans$l.height_scaled)

histogram(scans$l/height_scaled)

plot(scans$l/height_scaled)

l.M0_h_scaled<-lmer(l/height_scaled ~ n_encl + location_focal + 
                    time_meal + month + focal_sex + focal_age + (1|focal_id), 
                  na.action=na.omit, data=scans)

E0<-residuals(l.M0_h_scaled)
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

plot(l.M0_h_scaled) 

plot(scans$l.height_scaled, E0)

#looks better before transformation!

#let's try backwards selection...

M.full<-lmer(height_scaled ~ n_encl + location_focal + 
                    time_meal + month + focal_sex + focal_age + (1|focal_id), 
                  na.action=na.omit, data=scans)
summary(M.full)
anova(M.full)

#drop focal age
M1<-lmer(height_scaled ~ n_encl + location_focal + 
               time_meal + month + focal_sex + (1|focal_id), 
             na.action=na.omit, data=scans)
anova(M.full,M1)

#slight, insignificant difference

anova(M1)

#drop month

M2<-lmer(height_scaled ~ n_encl + location_focal + 
           time_meal + focal_sex + (1|focal_id), 
         na.action=na.omit, data=scans)

anova(M1,M2)

#slight, insignificant difference

anova(M2)

#drop meals

M3<-lmer(height_scaled ~ n_encl + location_focal + focal_sex + (1|focal_id), 
         na.action=na.omit, data=scans)

anova(M3,M2)

#slight, insignificant difference

anova(M3)

#drop sex

M4<-lmer(height_scaled ~ n_encl + location_focal + (1|focal_id), 
         na.action=na.omit, data=scans)

anova(M3,M4)

#M3 not singificantly different, but has slightly lower AIC (p = .06)

anova(M4)

#drop n_encl

M5<-lmer(height_scaled ~ location_focal + (1|focal_id), 
         na.action=na.omit, data=scans)

anova(M5,M4)

anova(M5)

#keep M3 for Tukey HSd post hoc tests

M.full<-M5

summary(M.full)

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

multCompTukey2 <- glht(M.full, linfct = mcp(location_focal = "Tukey")) 
summary(multCompTukey2)

ggplot(scans, aes(x = n_encl, y = height_scaled)) + geom_boxplot()

plot(scans$height_scaled ~ scans$n_encl)


