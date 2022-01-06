#' ---
#' title: "US RNA and Covid case modelling"
#' author: "Julian Faraway"
#' output:
#'   github_document
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
par(mgp=c(1.5,0.5,0), mar=c(3.1,3.1,0.1,0), pch=20)
#'
#' # Setup
#'
library(svglite)
library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(lubridate)
library(mgcv)
library(broom)

load(here("data/usa.rda"))

#' Plot of all counties
#' 
#+ allN1cases
usa %>% ggplot(aes(N1+1,cases)) + 
  geom_point(size=0.5) +
  scale_x_continuous(name="N1",trans="log10") +
  scale_y_continuous(name= "Cases", trans="log10") +
  facet_wrap(~ county) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6))
#' Correlations between N1 and incidence
usa %>% 
  group_by(county) %>% 
  summarise(incidence = cor(log(N1+1),log(cases+1),use = "complete.obs")) %>% 
  knitr::kable(digits=3)
#' Indiana county as seen in the paper
ssite = "Indiana"
#+ indianalog
usa %>% filter(county == ssite) %>% 
  ggplot(aes(N1+1,cases+1)) + geom_point() + 
  geom_smooth(method="lm",color='black') +
  xlab("RNA gc/l") + ylab("Case rate") +
  scale_y_log10() + scale_x_log10() 
#+ indianauntran
usa %>% filter(county == ssite) %>% 
  ggplot(aes(N1,cases)) + geom_point() + 
  geom_smooth(color='black') +
  xlab("RNA gc/l") + ylab("Case rate") 
#' # Models
#' 
#' Linear fit in each county on log scales
#' 
#' Need +1 for zero values
#' 
lmodpop = lm(log(cases+1) ~ log(N1 + 1) * county, usa)
anova(lmodpop)
#'
summary(lmodpop)
#'
#' - Intercepts vary greatly by county (more so than in Scotland) indicating
#' large heterogeneity in sewage systems
#' - Slopes are less than one (sometimes substantially) which shows the same
#' less than linear (on the original scale) behaviour as also seen in the Scotland
#' data
#' 
#' 
#' # Time Varying Models
#' 
#' Create a week number for the time index
#' 
#' Date start on
min(usa$date)
#'
usa %>% 
  mutate(timet = (as.numeric(date) - as.numeric(min(usa$date)))/7,
         y = log(cases+1),
         x = log(N1+1),
         county = factor(county)) -> usa
#'
#' Fit the model
#' 
#' 
gmod = gam(y ~ s(timet,by=county) + s(timet,by=x), data=usa)
summary(gmod)
#'
#' First plot the intercept coefficient function $\alpha_i (t)$.
#' Grey vertical lines represent quarters (start Jul20, Oct20, Jan21, Apr21, 
#' Jul21, Oct21)
#+ iceptusa
allcounties = levels(usa$county)
par(mgp=c(1.5,0.5,0), mar=c(1.5,1.5,1.1,0), pch=20, mfrow=c(2,3))
for(i in 1:6){
  plot(gmod, rug=FALSE, xlab="", ylab="",
       main=allcounties[i], select=i,ylim=c(-5,5))
  abline(v=c(17,30,43,56,69, 82),col=gray(0.75))
}
for(i in 7:12){
  plot(gmod, rug=FALSE, xlab="", ylab="",
       main=allcounties[i], select=i,ylim=c(-5,5))
  abline(v=c(17,30,43,56,69,82),col=gray(0.75))
}
for(i in 13:18){
  plot(gmod, rug=FALSE, xlab="", ylab="",
       main=allcounties[i], select=i,ylim=c(-5,5))
  abline(v=c(17,30,43,56,69, 82),col=gray(0.75))
}
for(i in 19:24){
  plot(gmod, rug=FALSE, xlab="", ylab="",
       main=allcounties[i], select=i,ylim=c(-5,5))
  abline(v=c(17,30,43,56,69,82),col=gray(0.75))
}
for(i in 25:30){
  plot(gmod, rug=FALSE, xlab="", ylab="",
       main=allcounties[i], select=i,ylim=c(-5,5))
  abline(v=c(17,30,43,56,69, 82),col=gray(0.75))
}
par(mfrow=c(1,1))
#' - Some counties have missing data, particularly at the beginning of 
#' the period which leads to greater uncertainty (wide bands)
#' 
#' Nicer plots for publication
#+ fig.show='hide'
pgam = plot(gmod)
pgdf = data.frame(x=unlist(map(pgam,"x")),
                  fit=unlist(map(pgam,"fit")),
                  se=unlist(map(pgam,"se")),
                  county = rep(c(allcounties,"beta"),each=100))
pgdf$ub = pgdf$fit+2*pgdf$se
pgdf$lb = pgdf$fit-2*pgdf$se
#+ usaalpha
pgdf %>% filter(county == ssite) %>% 
  ggplot(aes(x=x)) + 
  geom_ribbon(aes(ymin = lb, ymax = ub),fill="gray90") +
  geom_line(aes(y=fit)) +
  ylab("a(t)") +
  scale_x_continuous(name="Date",
                     breaks=c(17,30,43,56,69,82),
                     labels=c("Jul20", "Oct20", "Jan21", "Apr21", 
                              "Jul21", "Oct21"))
#+ usabeta
pgdf %>% filter(county == "beta") %>% 
  ggplot(aes(x=x)) + 
  geom_ribbon(aes(ymin = lb, ymax = ub),fill="gray90") +
  geom_line(aes(y=fit)) +
  ylab("b(t)") +
  scale_x_continuous(name="Date",
                     breaks=c(17,30,43,56,69,82),
                     labels=c("Jul20", "Oct20", "Jan21", "Apr21", 
                              "Jul21", "Oct21"))
#'
#' # Lagged Models for Prediction
#' 
#' Can we predict the future? Fit the concurrent model
#' 
lmodpop = lm(log(cases+1) ~ log(N1 + 1) * county, usa)
glance(lmodpop)[c('sigma','adj.r.squared')]
#' 
#' Try lagging the N1 variable by one week
#' 
lmodlag1 = lm(log(cases+1) ~ lag(log(N1 + 1),1) * county, usa)
glance(lmodlag1)[c('sigma','adj.r.squared')]
#'
#' Worse than the concurrent model
#' 
#' Now suppose we use case data to predict one week ahead:
#' 
lmodnoww = lm(log(cases+1) ~ lag(log(cases + 1),1) * county, usa)
glance(lmodnoww)[c('sigma','adj.r.squared')]
#'
#' Works better than using WW even though it is a very simplistic model
#' 
lmodboth = lm(log(cases+1) ~ lag(log(cases + 1),1) + lag(log(N1+1),1) * county, usa)
glance(lmodboth)[c('sigma','adj.r.squared')]
#'
#' Adding WW info results in no improvement to the predictive
#' ability of the model 

