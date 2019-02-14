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

#view and summarize

View(focals)
summary(focals)
str(focals$actor_id)
str(focals$focal_id)

#convert durations to minutes

focals$dur_m <- focals$duration_sec/60
focals$dur_h <- focals$dur_m/60

sum(focals$duration_sec)
sum(focals$dur_m)
sum(focals$dur_h)
#199.67 focal hours collected

#create summarized table

sum_focals1<-focals %>%
  dplyr::filter(actor_id!="FG") %>% 
  dplyr::filter(actor_id!="HUMAN")
  
sum_focals1$actor_id<-factor(sum_focals1$actor_id, levels=c("CL","DU","FY","MA","ME","PA","PE","PO","PP","RK","TR"))
str(sum_focals1)  

#spatial analysis - preferences for Center enclosure

sum_cp<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  dplyr::filter(n_encl=="2") %>%
  group_by(focal_id, focal_sex, focal_group, 
           focal_age, time_meal, month, location) %>% 
  complete(location, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur)) %>%
  dplyr::filter(location=="center")

View(sum_cp)
summary(sum_cp$prop_time)

histogram(sum_cp$prop_time)

plot(sum_cp$prop_time)

plot(sum_cp$prop_time ~ sum_cp$focal_id)
plot(sum_cp$prop_time ~ sum_cp$focal_age)
plot(sum_cp$prop_time ~ sum_cp$focal_sex)
plot(sum_cp$prop_time ~ sum_cp$month)

#4 monkeys preferred satellites across conditions; 7 preferred Center

#let's make some models

M0_center<-gls(prop_time ~ time_meal + focal_sex + focal_age, 
                 na.action=na.omit, data=sum_cp, method="ML")

#add random factor

M1_center<-lme(prop_time ~ time_meal + focal_sex + focal_age, 
               random = ~1|focal_group, na.action=na.omit, data=sum_cp, method="ML")

#add nesting

M2_center<-lme(prop_time ~ time_meal + focal_sex + focal_age, 
               random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

#anova
anova(M0_center, M1_center, M2_center)

#M2 looks good... let's examine residuals

E2<-residuals(M2_center)
str(E2)
E2
summary(E2)

plot(sum_cp$month, E2, xlab="Month", ylab="Residuals")
plot(sum_cp$time_meal, E2, xlab="Meal Status", ylab="Residuals")
plot(sum_cp$focal_sex, E2, xlab="Focal Sex", ylab="Residuals")
plot(sum_cp$focal_age, E2, xlab="Focal Age", ylab="Residuals")


qqnorm(E2)
qqline(E2)
ad.test(E2)

plot(M2_center) 

plot(sum_cp$prop_time, E2)

#M2 looks good...

#let's drop factors?

M.full<-lme(prop_time ~ time_meal + focal_sex + focal_age + month, 
               random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

anova(M.full)

#drop meal status

M1<-lme(prop_time ~ month + focal_sex + focal_age, 
            random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

anova(M.full,M1)

#improved!

anova(M1)

#drop focal sex

M2<-lme(prop_time ~ month + focal_age, 
        random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

anova(M2,M1)

#model improved slightly, insignificantly

anova(M2)

#drop month

M3<-lme(prop_time ~ focal_age, 
        random = ~1|focal_group/focal_id, na.action=na.omit, data=sum_cp, method="ML")

anova(M3,M2)

#M3 looks better, insignificantly

#leave focal age in there?
  
#stick with M3...

M.full<-M3

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

multCompTukey1 <- glht(M.full, linfct = mcp(focal_age = "Tukey")) 
summary(multCompTukey1)

anova(M.full)


