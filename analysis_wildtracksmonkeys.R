#COMPLETE ANALYSIS OF 2016 WILDTRACKS DATA, FEBRUARY 2018##############################

#Github

##############set working directory####################################################

setwd("~/Wildtracks - Spider Monkeys/Publication Work/Re-Analysis")

##############load packages############################################################

library(dplyr)
library(igraph)
library(ggplot2)
library(readxl)
library(readr)

##############import###################################################################

focals <- read_excel("~/Wildtracks - Spider Monkeys/Publication Work/Re-Analysis/focals_R.xlsx", 
                       col_types = c("skip", "text", "date", 
                                     "date", "skip", "skip", "text", 
                                     "text", "text", "numeric", "numeric", 
                                     "text", "text", "text", "skip"))
#change $time to proper format

focals$time<-format(focals$time, "%H:%M:%S")

#view and summarize

View(focals)
summary(focals)

##############convert durations to minutes#############################################

focals$dur_m <- focals$duration_sec/60
focals$dur_h <- focals$dur_m/60

sum(focals$dur_m)
sum(focals$dur_h)

##############create subset table by id################################################
summary_id1<-focals %>% 
  group_by(focal_id, behavior) %>%
  filter(focal_id==actor_id) %>% 
  summarize(total_dur_m=sum(dur_m), total_dur_h=sum(dur_h),
            total_dur_cbc=sum(dur_m, focals$condition=="cbc"), 
            total_dur_sbc=sum(dur_m, focals$condition=="sbc"),
            total_dur_sat=sum(dur_m, focals$condition=="satellite"),
            total_dur_cen=sum(dur_m, focals$condition=="center"))

View(summary_id1)

##############convert to rates/rel dur##########################################