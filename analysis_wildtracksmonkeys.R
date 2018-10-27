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

sum_focals<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  dplyr::select(duration_sec, occurrences, month, focal_id, time_meal, n_encl, location, activity, focal_sex, focal_group, focal_age) %>%
  group_by(focal_id, focal_sex, focal_group, focal_age, month, time_meal, n_encl, location, activity) %>%
  complete(activity, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))

sum_focals$nest <- with(sum_focals, factor(paste(focal_group,focal_id)))

View(sum_focals)

#look at overall activity budgets

sum_overall<-sum_focals1 %>%
  dplyr::filter(actor_id==focal_id) %>%
  group_by(activity) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))
View(sum_overall)
sum_overall

#drop agonistic, caregiver, other, parental, sociosexual, solitary play, vocalization
#keep/analyze abnormal, forage, inactive, movement, not_visible, prosocial, self_directed