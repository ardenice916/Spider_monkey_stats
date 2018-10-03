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
install.packages("ggplots")

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

#import

focals <- read.table(file="focals_Wildtracks.csv", sep=",", header=T)#alternate import code

focals <- read_csv("focals_Wildtracks.csv", 
                              col_types = cols(condition = col_skip(), 
                                               date = col_date(format = "%m/%d/%Y"), 
                                               month = col_factor(levels = c("jun-jul", "jul-aug", "aug-sep")), 
                                               activity = col_factor(levels = c("abnormal","agonistic","caregiver",
                                                                                "forage","inactive","movement","not_visible",
                                                                                "other","parental","prosocial","self_directed",
                                                                                "sociosexual","solitary_play","vocalization")),
                                               duration_sec = col_number(), enclosure_access = col_skip(), 
                                               location = col_factor(levels = c("satellite", 
                                                                                "center")), 
                                               n_encl = col_factor(levels = c("1","2")), 
                                               observer_id = col_character(), 
                                               occurrences = col_integer(), record = col_skip(), 
                                               time = col_time(format = "%H:%M:%S"), 
                                               focal_sex = col_factor(levels = c("male","female")),
                                               focal_group = col_factor(levels = c("sat1","sat2","sat3")),
                                               focal_age = col_factor(levels = c("subadult","adult")),
                                               time_meal = col_factor(levels = c("none", 
                                                                                 "before", "after")), 
                                               trim_ws = FALSE))

#view and summarize

View(focals)
summary(focals)
str(focals)

#convert durations to minutes

focals$dur_m <- focals$duration_sec/60
focals$dur_h <- focals$dur_m/60

sum(focals$duration_sec)
sum(focals$dur_m)
sum(focals$dur_h)
#199.67 focal hours collected

#create summarized table

sum_focals<-focals %>%
  filter(focal_id==actor_id) %>%
  dplyr::select(duration_sec, occurrences, month, focal_id, time_meal, n_encl, location, activity, focal_sex, focal_group, focal_age) %>%
  group_by(focal_id, focal_sex, focal_group, focal_age, month, time_meal, n_encl, location, activity) %>%
  complete(focal_id, activity, fill = list(duration_sec = 0)) %>%
  summarize(sum_dur=sum(duration_sec)) %>%
  mutate(prop_time=sum_dur/sum(sum_dur))

View(sum_focals)

#Look at mixed effects model
#start without random factor

sum_abnormal<-sum_focals %>%
  filter(activity=="abnormal")

M0_abnormal<-gls(prop_time ~ n_encl*location + focal_id + time_meal + month, data=sum_abnormal, method="ML")

