#####################################################
#R Script to Analyze Water Data Set
#from Bailey's and Izak's Theses
#
#19-Sep-17
#####################################################

#####################################################
#Load Packages
#####################################################
install.packages("nlme")
install.packages("lme4")
install.packages("lmerTest")
install.packages("dplyr")
install.packages("nortest")
install.packages("ggplot2")
install.packages("multcomp")

library(nlme)
library(lme4)
library(lmerTest)
library(dplyr)
library(nortest)
library(ggplot2)
library(multcomp)

#####################################################
#Load and Assess Throughfall Data
#####################################################

tf = read.table(file="tf.summary.csv", header=T, sep=",")
str(tf)
unique(tf$date)
unique(tf$site.name)
unique(tf$sample.event)
unique(tf$site)

     #ammonium and nitrate concentrations have NA so they come across
     #as factors.  need to change if we analyze concentrations

#set sample event and rep from integer to factor
tf$f.sample.event=as.factor(tf$sample.event)
tf$f.rep=as.factor(tf$rep)
str(tf)

#create a factor equivalent to nesting
tf$nest <- with(tf, factor(paste(site,f.rep)))

#remove rainfall samples
tf<-subset(tf, type=="TF")

#####################################################
#Analyze NH4 flux base models (units, kg N ha-1)
#####################################################

#Look at mixed effects model
#start without random factor
M0<-gls(col.volume.mm ~ budworm*f.sample.event, 
        na.action=na.omit, data=tf, method="ML")
     #adding site as a blocking factor produces a singularity and no solution

#add random factor - refer to chapter 5 of zuur
M1<-lme(col.volume.mm ~ budworm*f.sample.event,  
        random = ~1|site, na.action=na.omit, data=tf, method="ML")

#try nesting
M2<-lme(col.volume.mm ~ budworm*f.sample.event,
        random = ~1|nest, na.action=na.omit, data=tf, method="ML")

anova(M0,M1,M2)
     #we'll use M2

#Analyze base model residuals

E2<-residuals(M2)

plot(filter(tf, !is.na(col.volume.mm)) %>% dplyr::select(site),
     E2, xlab="Site", ylab="Residuals")
plot(filter(tf, !is.na(col.volume.mm)) %>% dplyr::select(budworm),
     E2, xlab="Budworm", ylab="Residuals")

qqnorm(residuals(M2))
qqline(residuals(M2))
ad.test(residuals(M2))
     
plot(M2) 

x<-tf$col.volume.mm[!is.na(tf$col.volume.mm)]#removes na values from column
E2<-residuals(M2,type="normalized")
plot(x, E2)
     #residuals are linear

#####################################################
#Analyze models with alternate variance structures
#####################################################

#base model from above
M2<-lme(col.volume.mm ~ budworm*f.sample.event,  
        random = ~1|nest, na.action=na.omit, data=tf)#run with method="REML" default for comparison

#alternate variance structures
vf1=varIdent(form = ~1|site)
vf2=varPower(form = ~fitted(.))
vf3=varExp(form = ~fitted(.))
vf4=varConstPower(form = ~fitted(.))
vf5=varPower(form = ~fitted (.)|site)
vf6=varExp(form = ~fitted(.)|site)
vf7=varConstPower(form = ~fitted(.)|site)
vf8=varIdent(form = ~1|f.sample.event)

#alternate models
M2.1<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf1)

M2.2<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf2)

M2.3<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf3)

M2.4<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf4)

M2.5<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf5)
     
M2.6<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf6)
     
M2.7<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf7)
     #singularity
     
M2.8<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf8)

M2.9<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=varComb(vf1,vf2))
     #singularity

M2.10<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=varComb(vf1,vf3))

M2.11<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=varComb(vf1,vf4))     
     #singularity

anova(M2,M2.1,M2.2,M2.3,M2.4,M2.5,M2.6,M2.8,M2.10)
     #M2.2 is the best, alternate variance structure power of fits

#Analyze residuals of best model with alternate variance structure

E2.2<-residuals(M2.2)

plot(filter(tf, !is.na(col.volume.mm)) %>% dplyr::select(site),
     E2.2, xlab="Site", ylab="Residuals")
plot(filter(tf, !is.na(col.volume.mm)) %>% dplyr::select(budworm),
     E2.2, xlab="Budworm", ylab="Residuals")

qqnorm(residuals(M2.2))
qqline(residuals(M2.2))
ad.test(residuals(M2.2))
     
plot(M2.2) 
     
x<-tf$col.volume.mm[!is.na(tf$col.volume.mm)]#removes na values from column
E2.2<-residuals(M2.2,type="normalized")
plot(x, E2.2)
     
#check for autocorrelation
acf(E2.2, na.action=na.pass,
    main="Auto-correlation plot for residuals")
     #little autocorrelation in the first lag 

#####################################################
#Check Autocorrelation Models
#####################################################

#new base model from above
M2.2<-lme(col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf2)

#autocorrelation models
M3<-lme(col.volume.mm ~ budworm*f.sample.event,
        random = ~1|nest, na.action=na.omit, data=tf, weights=vf2,
        correlation=corCompSymm(form=~f.sample.event))

cs1<-corARMA(c(0.2), p=1, q=0)

M5<-lme(col.volume.mm ~ budworm*f.sample.event,
        random = ~1|nest, na.action=na.omit, data=tf, weights=vf2,
        correlation=cs1)

anova(M2.2,M3,M5)
     #autocorrelation does not improve model

#Analyze Residuals of Autocorrelation Model

#E3<-residuals(M3)

#plot(filter(tf, !is.na(col.volume.mm)) %>% dplyr::select(site),
#     E3, xlab="Site", ylab="Residuals")
#plot(filter(tf, !is.na(col.volume.mm)) %>% dplyr::select(budworm),
#     E3, xlab="Budworm", ylab="Residuals")

#qqnorm(residuals(M3))
#qqline(residuals(M3))
#ad.test(residuals(M3))

#plot(M3) 
     
#x<-tf$col.volume.mm[!is.na(tf$col.volume.mm)]#removes na values from column
#E3<-residuals(M3,type="normalized")
#plot(x, E3)
     #this is about the best they've looked

#####################################################
#Analyze NH4 flux base models with NORMALIZED data (units, kg N ha-1)
#####################################################

tf$l.col.volume.mm=log(tf$col.volume.mm+1)

#Look at mixed effects model
#start without random factor
M0<-gls(l.col.volume.mm ~ budworm*f.sample.event, 
        na.action=na.omit, data=tf, method="ML")
     #adding site as a blocking factor produces a singularity and no solution

#add random factor - refer to chapter 5 of zuur
M1<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
        random = ~1|site, na.action=na.omit, data=tf, method="ML")

#try nesting
M2<-lme(l.col.volume.mm ~ budworm*f.sample.event,
        random = ~1|nest, na.action=na.omit, data=tf, method="ML")

anova(M0,M1,M2)
     #M2 is best

#look at residuals of model with normalized data
E2<-residuals(M2)

plot(filter(tf, !is.na(l.col.volume.mm)) %>% dplyr::select(site),
     E2, xlab="Site", ylab="Residuals")
plot(filter(tf, !is.na(l.col.volume.mm)) %>% dplyr::select(budworm),
     E2, xlab="Budworm", ylab="Residuals")

qqnorm(residuals(M2))
qqline(residuals(M2))
ad.test(residuals(M2))
     
plot(M2) 
     #residuals v model fit look passable

x<-tf$l.col.volume.mm[!is.na(tf$l.col.volume.mm)]#removes na values from column
E2<-residuals(M2,type="normalized")

plot(x, E2)
     
#####################################################
#Analyze models using NORMALIZED DATA with alternate variance structures 
#####################################################

#base model from above
M2<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
        random = ~1|nest, na.action=na.omit, data=tf)#run with method="REML" default for comparison

#alternate variance structures
vf1=varIdent(form = ~1|site)
vf2=varPower(form = ~fitted(.))
vf3=varExp(form = ~fitted(.))
vf4=varConstPower(form = ~fitted(.))
vf5=varPower(form = ~fitted (.)|site)
vf6=varExp(form = ~fitted(.)|site)
vf7=varConstPower(form = ~fitted(.)|site)
vf8=varIdent(form = ~1|f.sample.event)

#alternate models
M2.1<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf1)

M2.2<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf2)
     
M2.3<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf3)

M2.4<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf4)

M2.5<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf5)
     #singularity

M2.6<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf6)

M2.7<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf7)

M2.8<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=vf8)

M2.9<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
          random = ~1|nest, na.action=na.omit, data=tf, weights=varComb(vf1,vf2))

M2.10<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
           random = ~1|nest, na.action=na.omit, data=tf, weights=varComb(vf1,vf3))

M2.11<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
           random = ~1|nest, na.action=na.omit, data=tf, weights=varComb(vf1,vf4))     

anova(M2,M2.1,M2.2,M2.3,M2.4,M2.6,M2.7,M2.8,M2.9,M2.10,M2.11)    
     #M2.7 is best

#look at model residuals
E2.7<-residuals(M2.7)

plot(filter(tf, !is.na(l.col.volume.mm)) %>% dplyr::select(site),
     E2.7, xlab="Site", ylab="Residuals")
plot(filter(tf, !is.na(l.col.volume.mm)) %>% dplyr::select(budworm),
     E2.7, xlab="Budworm", ylab="Residuals")

qqnorm(residuals(M2.7))
qqline(residuals(M2.7))
ad.test(residuals(M2.7))
     
plot(M2.7) 
     
x<-tf$col.volume.mm[!is.na(tf$col.volume.mm)]#removes na values from column
E2.7<-residuals(M2.7,type="normalized")
plot(x, E2.7)
     #residuals v observation are not good

#check for autocorrelation
acf(E2.7, na.action=na.pass,
    main="Auto-correlation plot for residuals")
     #there is some weak autocorrelation in the second lag, Zuur Ch 6

#####################################################
#Check Autocorrelation Models using NORMALIZED DATA
#####################################################

M3<-lme(l.col.volume.mm ~ budworm*f.sample.event,
        random = ~1|nest, na.action=na.omit, data=tf,
        correlation=corCompSymm(form=~f.sample.event), weights=vf7)

cs1<-corARMA(c(0.2), p=1, q=0)

M5<-lme(l.col.volume.mm ~ budworm*f.sample.event,
        random = ~1|nest, na.action=na.omit, data=tf, weights=vf7,
        correlation=cs1)

anova(M2.7,M3,M5)
     #autocorrelation models are not better

#look at model residuals
#E5<-residuals(M5)

#plot(filter(tf, !is.na(l.col.volume.mm)) %>% dplyr::select(site),
#     E5, xlab="Site", ylab="Residuals")
#plot(filter(tf, !is.na(l.col.volume.mm)) %>% dplyr::select(budworm),
#     E5, xlab="Budworm", ylab="Residuals")

#qqnorm(residuals(M5))
#qqline(residuals(M5))
#ad.test(residuals(M5))

#plot(M5) 
     #residuals v fit looks OK

#x<-tf$l.col.volume.mm[!is.na(tf$l.col.volume.mm)]#removes na values from column
#E5<-residuals(M5,type="normalized")
#plot(x, E5)
     
#####################################################
#Get Full Model Statistics and Make Graph
#####################################################

#best overall model was log normalized, nested, and no alternate variance
M.full<-lme(l.col.volume.mm ~ budworm*f.sample.event,  
            random = ~1|nest, na.action=na.omit, data=tf)

anova(M.full)

#post hoc test
model.matrix.gls <- function(M.full, ...){
  model.matrix(terms(M.full), data = getData(M.full), ...)  
}
model.frame.gls <- function(M.full, ...){
  model.frame(formula(M.full), data = getData(M.full), ...)  
}
terms.gls <- function(M.full, ...){
  terms(model.frame(M.full),...)  
}

multCompTukey <- glht(M.full, linfct = mcp(f.sample.event = "Tukey")) 
summary(multCompTukey)

#create table x to interpret post hoc test
x <- group_by(tf, f.sample.event) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(col.volume.mm.mean = mean(col.volume.mm, na.rm = TRUE), # na.rm = TRUE to remove missing values
            col.volume.mm.sd=sd(col.volume.mm, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(col.volume.mm)), # of observations, excluding NAs. 
            col.volume.mm.se=col.volume.mm.sd/sqrt(n))

sort(x$col.volume.mm.mean, index.return=T)

#create new table x to graph by categorical time for graphing purposes
x <- group_by(tf, budworm, f.sample.event) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(col.volume.mm.mean = mean(col.volume.mm, na.rm = TRUE), # na.rm = TRUE to remove missing values
            col.volume.mm.sd=sd(col.volume.mm, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(col.volume.mm)), # of observations, excluding NAs. 
            col.volume.mm.se=col.volume.mm.sd/sqrt(n))
#this code defines graphing.interval as date format, but we won't use it for now
#x$graph.interval <-as.Date(as.character(x$graph.interval), format="%m/%d/%Y")

#make a new vector with the categorical times
cat.time<-c("11Sep15", "11Oct15", "29Oct15", "8Nov15", "8May16", "22May16", "21Jun16", "13Jul16", "21Jul16", "19Sep16")
#force the new vector to be characters
x$cat.time<-as.character(cat.time)
#force the new vector to be ordered in the order you gave it instead of alphabetical
x$cat.time<-factor(x$cat.time, levels=unique(x$cat.time))

pd=position_dodge(0.1)

ggplot(x, aes(x=cat.time, y=col.volume.mm.mean)) + 
  geom_errorbar(aes(ymin=col.volume.mm.mean-col.volume.mm.se, ymax=col.volume.mm.mean+col.volume.mm.se), color="black", width=0.1, position=pd) + 
  geom_line(position=pd, color="black", aes(group=budworm)) +
  geom_point(size=3, pch=21, aes(fill=budworm)) +
  xlab("Sample Date") +
  ylab(expression(Throughfall~water~transmission~(mm))) +
  scale_fill_manual(name="Budworm Activity", values=c("white", "black")) +
  theme_bw() +
  theme(legend.justification=c(0.03,0.6),
        legend.position=c(0.65,0.88),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size= 12),
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=1, y=11, label="bc", size=4) +
  annotate("text", x=2, y=7, label="ab", size=4) +
  annotate("text", x=2.7, y=14, label="c", size=4) +
  annotate("text", x=3.8, y=59, label="d", size=4) +
  annotate("text", x=5.1, y=8, label="a", size=4) +
  annotate("text", x=6, y=19, label="c", size=4) +
  annotate("text", x=7, y=11, label="bc", size=4) +
  annotate("text", x=7.9, y=6, label="a", size=4) +
  annotate("text", x=8.7, y=28, label="c", size=4) +
  annotate("text", x=10.2, y=6, label="a", size=4) +
  annotate("text", x=0.8, y=60, label="budworm,", size=4, hjust=0) +
  annotate("text", x=0.8, y=57, label="p=0.74", size=4, hjust=0) +
  annotate("text", x=0.8, y=54, label="sample event,", size=4, hjust=0) +
  annotate("text", x=0.8, y=51, label="p<0.0001", size=4, hjust=0) +
  annotate("text", x=0.8, y=48, label="budworm*", size=4, hjust=0) +
  annotate("text", x=1, y=45, label="sample date,", size=4, hjust=0) +
  annotate("text", x=0.8, y=42, label="p<0.0001", size=4, hjust=0)


#this will save the file
ggsave('figures/TFWaterFlux.tiff',
       units="in",
       width=5.5,
       height=4.5,
       dpi=1200,
       compression="lzw")

