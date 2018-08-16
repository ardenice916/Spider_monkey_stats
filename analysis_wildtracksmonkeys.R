#COMPLETE ANALYSIS OF 2016 WILDTRACKS DATA, FEBRUARY 2018

#set working directory

setwd("~/Spider_monkey_stats")

#load useful packages

library(dplyr)
library(igraph)
library(ggplot2)
library(readxl)
library(readr)
library(nortest)

#import

focals <- read_csv("focals_Wildtracks.csv", 
                              col_types = cols(condition = col_skip(), 
                                               date = col_date(format = "%m/%d/%Y"), 
                                               duration_sec = col_number(), enclosure_access = col_skip(), 
                                               location = col_factor(levels = c("satellite", 
                                                                                "center")), n_encl = col_factor(levels = c("1", 
                                                                                                                           "2")), observer_id = col_character(), 
                                               occurrences = col_integer(), record = col_skip(), 
                                               time = col_time(format = "%H:%M:%S"), 
                                               time_meal = col_factor(levels = c("none", 
                                                                                 "before", "after")), X18 = col_skip(), X19 = col_skip(), 
                                               X20 = col_skip(), X21 = col_skip()), trim_ws = FALSE)

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
  select(focal_id, time_meal, n_encl, location, activity) %>%
  group_by(focal_id, time_meal, n_encl, location, activity)

View(sum_focals)


