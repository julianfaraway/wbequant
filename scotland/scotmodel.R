#' ---
#' title: "Modelling for Scotland"
#' author: "Julian Faraway"
#' output:
#'   github_document:
#'      toc: true
#' ---
#+ echo=FALSE
knitr::opts_chunk$set(comment=NA, 
                      echo = TRUE,
                      fig.path="figs/", 
                      dev = 'svglite',  
                      fig.ext = ".svg",
                      warning=FALSE, 
                      message=FALSE)
ggplot2::theme_set(ggplot2::theme_bw())
par(mgp=c(1.5,0.5,0), mar=c(3.1,3.1,3.1,0), pch=20)

#' # Setup
#' 
#' Set some graphics options.
#' The here package ensures files are found in the expected places.

library(svglite)
library(here)
library(tidyverse)
library(lubridate)
library(mgcv)
library(broom)

#' Load previously processed data

load(here("data/scot.rda"))

#' Restrict to one year

startdate = ymd("2020-09-01")
enddate = ymd("2021-08-31")

allscot %>% 
  filter(between(Date, startdate,enddate)) -> allscot

#' Merge in population data and create cases per 100,000 population variable.
#' Need to add +1 to deal with zeroes on log scale

left_join(allscot, scotpop, by = c("hab" = "HBName")) %>% 
  mutate(rollper = (rollcases+1)*1e5/hbpop) %>% 
  select(Date,hab,N1,cases,hbpop,rollper,rollcases) -> allscot
allhabs = unique(allscot$hab)
#' Glasgow as the reference level
allscot$hab = fct_relevel(allscot$hab,"Greater Glasgow and Clyde")
#'
#' # Data
#' 
#' Some information about the data
#' 
#' Population and WWTP coverage fraction
#' 
scotpop %>% select(HBName,hbpop,wwpop,covfrac) %>% 
  knitr::kable(digits=2,
               col.names = c("Health Area", "Population",
                             "WWTP Population", "Fraction"))
#'
#' WWTP sites information
#' 
scotrna %>% 
  filter(between(Date, startdate,enddate)) %>% 
  group_by(hab) %>% 
  summarise(popsize = sum(as.numeric(unique(pop))),
            nobs = n(),
            nsites = length(unique(site))) %>% 
  knitr::kable(col.names = c("Health Area", "WW Population",
                             "WW observations", "Number of sites"))
#'
#' # Comparison to Infection Survey
#' 
isd = read_csv(here(file.path("data","isd.csv")))
#' pick out and group PH data
scotcc %>% group_by(Date) %>% summarise(cases = sum(DailyPositive)) -> scc
#' compute weekly rolling sums
scc %>% mutate(wcases = zoo::rollsumr(scc$cases, k=7, fill="extend")) -> scc
#' join
left_join(isd, scc, by=c("end" = "Date")) -> phinf
#' 
#+ frac
ggplot(phinf, aes(x=end, y=wcases/ninf)) + geom_point() +
  geom_smooth(color='black') + 
  ylab("Ratio") + xlab("Date") + xlim(as_date(c("2020-10-31","2021-09-25")))
#' - Ratio peaks around 0.35
#' - Underreporting and shedding profiles mean ratio will not be 1
#' - Ratio is increasing in time which has several possible explanations
#' 
#' # Glasgow only
#' 
#' Only consider Glasgow (highest population, most measurements)
#' 
allscot %>% 
  filter(hab == "Greater Glasgow and Clyde") -> glasgow
#+ glasn1inc
glasgow %>% 
  ggplot(aes(x=N1+20, y = rollper)) + 
  geom_point(size=0.5,alpha=0.25) +
  scale_x_continuous(name = "RNA gc/l", trans = "log10") +
  scale_y_continuous(name = "COVID-19 case rate", trans = "log10") +
  geom_smooth(method="lm", color='black')
#' Looks fine on a log scale but now consider untransformed
#+ glasn1incnolog
glasgow %>% 
  ggplot(aes(x=N1+20, y = rollper)) + 
  geom_point(size=0.5,alpha=0.25) +
  scale_x_continuous(name = "RNA gc/l") +
  scale_y_continuous(name = "COVID-19 case rate") +
  geom_smooth(color='black') 
#' Can see non-linearity on the untransformed scale
#' 
#' Fit the a linear model with the variables on a log scale
glasmod = lm(log(rollper) ~ log(N1), glasgow)
summary(glasmod)
#'
#' # Level of Aggregation
#'
#' Create a whole of Scotland dataset 
left_join(allscot, scotpop, by = c("hab" = "HBName")) %>% 
  mutate(wbyn1 = N1*wwpop) %>% 
  group_by(Date) %>% 
  summarise(N1 = sum(wbyn1)/sum(scotpop$wwpop),
            cases = sum(rollper)*1e5/sum(scotpop$hbpop)) %>% 
  ungroup() -> aggscot
#+ aggscot
ggplot(aggscot,aes(N1, cases)) +
  geom_point(size=0.5,alpha=0.25) +
  scale_x_continuous(name = "RNA gc/l", trans = "log10") +
  scale_y_continuous(name = "COVID-19 case rate", trans = "log10") +
  geom_smooth(method="lm", color='black')    
#' A stronger relationship is seen compared to the individual areas
#' 
#' Now fit a model
lmodagg = lm(log(cases) ~ log(N1), aggscot)
summary(lmodagg)
exp(predict(lmodagg,data.frame(N1=mean(aggscot$N1)),interval = "prediction"))
#' compared to a Glasgow only model
glasmod = lm(log(rollper) ~ log(N1), glasgow)
summary(glasmod)
exp(predict(glasmod,data.frame(N1=mean(aggscot$N1)),interval = "prediction"))
#' 
#' 
#' # Relationships across Scotland 
#' 
#' Four selected areas as seen in the paper
#' 
selset = c("Borders","Lothian","Lanarkshire","Orkney")
#+ relhab4
allscot %>% 
  filter(hab %in% selset) %>% 
  ggplot(aes(N1,rollper)) + geom_point(alpha=0.25,size=0.25) +
  ylab("COVID-19 case rate") + xlab("RNA gc/l") +
  facet_wrap(~ hab, scales = "free") +
  geom_smooth(se=FALSE,color='black') +
  scale_x_continuous(labels = scales::scientific)
#' Now all of them
#+ allrelhab, height=11
allscot %>% 
  ggplot(aes(N1,rollper)) + geom_point(alpha=0.25,size=0.25) +
  ylab("COVID-19 case rate") + xlab("RNA gc/l") +
  facet_wrap(~ hab, scales = "free") +
  geom_smooth(se=FALSE,color='black') +
  scale_x_continuous(labels = scales::scientific) +
  theme(axis.text.x = element_text(angle = 30,hjust=1))
#' 
#' # Varying Coefficient models
#' 
#' Assemble data
scotplus = allscot %>% 
  filter(between(Date, as.Date("2020-09-01"),as.Date("2021-08-31"))) %>% 
  mutate(y = log10(rollper),
         x = log10(N1 + 20),
         timet = as.numeric(Date)-18505,
         hab = factor(hab))
#' Fit a model where each HAB has its own time varying intercept but
#' the slope is in common for all
#' $$ y = \alpha_{i}(t) + \beta(t)x $$ 
gmod = gam(y ~ s(timet,by=hab) + s(timet,by=x), data=scotplus)
summary(gmod)
#' Plot the coefficient functions, gray lines at Jan1, Apr1 and Jul1
#+ allhabgam
par(mgp=c(1.5,0.5,0), mar=c(3.1,3.1,3.1,0), pch=20, mfrow=c(2,2))
for(i in 1:4){
  pgam = plot(gmod, rug=FALSE, xlab="Time", ylab='a(t)',
       main=allhabs[i], select=i)
  abline(v=c(90,180,270)+30,col=gray(0.75))
}
for(i in 5:8){
  plot(gmod, rug=FALSE, xlab="Time", ylab='a(t)',
       main=allhabs[i], select=i)
  abline(v=c(90,180,270)+30,col=gray(0.75))
}
for(i in 9:12){
  plot(gmod, rug=FALSE, xlab="Time", ylab='a(t)',
       main=allhabs[i], select=i)
  abline(v=c(90,180,270)+30,col=gray(0.75))
}
par(mfrow=c(1,2))
for(i in 13:14){
  plot(gmod, rug=FALSE, xlab="Time", ylab='a(t)',
       main=allhabs[i], select=i)
  abline(v=c(90,180,270)+30,col=gray(0.75))
}
par(mfrow=c(1,1))
#' - Most HABs have a similar shape
#' - Number of WWTPs sampled increased from around 30 to 90 early in 2021
#+ slopegam
plot(gmod, rug=FALSE, xlab="Time", ylab='b(t)',
     select=15,scale=0)
abline(v=c(90,180,270)+30,col=gray(0.75))
#' Improvement on slope over time but not reaching the desired target
#' of one. Maybe lab procedures are improving
#' 
#' Generate nicer plots
#+ allgamplot, fig.show='hide'
pgam = plot(gmod, pages = 1)
pgdf = data.frame(x=unlist(map(pgam,"x")),
                  fit=unlist(map(pgam,"fit")),
                  se=unlist(map(pgam,"se")),
                  hab = rep(c(allhabs,"beta"),each=100))
pgdf$ub = pgdf$fit+2*pgdf$se
pgdf$lb = pgdf$fit-2*pgdf$se
#+ pgam4
pgdf %>% filter(hab %in% allhabs[1:4]) %>% 
  ggplot(aes(x=x)) + 
  geom_ribbon(aes(ymin = lb, ymax = ub),fill="gray90") +
  geom_line(aes(y=fit)) +
  facet_wrap(~ hab) +
  ylab("a(t)") +
  scale_x_continuous(name="Date",breaks=c(90,180,270),
                     labels=c("Jan1","Apr1","Jul1"))
#+ pgamb
pgdf %>% filter(hab == "beta") %>% 
  ggplot(aes(x=x)) + 
  geom_ribbon(aes(ymin = lb, ymax = ub),fill="gray90") +
  geom_line(aes(y=fit)) +
  ylab("b(t)") +
  scale_x_continuous(name="Date",breaks=c(90,180,270),
                     labels=c("Jan1","Apr1","Jul1"))
#' Look at change in sampling over time
#+ wwtpmon
scotrna %>% mutate(Mon = floor_date(Date, "month")) %>% 
  group_by(Mon) %>% count() %>% 
  ggplot(aes(Mon,n)) + geom_step() +
  xlab("Date") + ylab("Samples per Month")
#' 
#' # Lagged Models for Prediction
#' 
#' Can we predict the future? Try lagging the N1 variable:
#' 
timewindow = 7
ldf = tibble(lag = 0:timewindow, r2 = rep(NA, timewindow+1),
             se = rep(NA,timewindow+1))

for(i in 0:timewindow){
  allscot %>% group_by(hab) %>% mutate(N1lag = lag(N1,i)) %>% ungroup() -> scotlag
  lmodlag = lm(log(rollper) ~ log(N1lag+20)*hab, scotlag)
  gd = glance(lmodlag)
  ldf$r2[i+1] = gd$r.squared
  ldf$se[i+1] = gd$sigma
}

ldf %>% knitr::kable(digits = 3)
#' We can see that the fit of the model decreases as the lag is increased
#' while the standard error increases. While we might obtain better predictions
#' by using more than one lagged variable, this does not indicate any
#' strong precognitive ability for N1.
#' 
#' Consider trying to predict a week ahead:

allscot %>% group_by(hab) %>% mutate(rollperlag = lag(rollper,7)) %>% 
  ungroup() -> scotlag
lmodlag7 = lm(log(rollper) ~ log(rollperlag)*hab, scotlag)
glance(lmodlag7)[c('r.squared','sigma')]

#' This model predicts a week ahead using current cases. We can
#' see this works better than using N1. Now suppose
#' we have the current N1 count and we add this information to the model:
#' 
#' 
allscot %>% group_by(hab) %>% 
  mutate(rollperlag = lag(rollper,7),
         N1lag7 = lag(N1,7)) %>% 
  ungroup() -> scotlag
lmodlagcase7 = lm(log(rollper) ~ (log(rollperlag) + log(N1lag7+20))*hab, scotlag)
glance(lmodlagcase7)[c('r.squared','sigma')]

#' We can see this is only fractionally better than using cases alone i.e N1
