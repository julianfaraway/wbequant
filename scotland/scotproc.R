#Process Scottish cases and RNA counts

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(zoo)
library(here)

dataloc = "data"

# Read in RNA data

vnames=c("hab","site","Date","popband","pop","resdesc","N1","dss")

sc = read_csv(file.path(dataloc,"scotrna.csv"),col_names=vnames,skip=1)
              

# Fixup the data
sc$Date = dmy(sc$Date)
sc$dss = parse_number(sc$dss)

# Remove small number of rows with empty site

sc %>% filter(hab != "(Empty)") -> scotrna

#
# Covid Cases
#

cc = read_csv(file.path(dataloc,"scotcases.csv"))
cc$Date = ymd(cc$Date)
# Make same names as RNA data
cc$HBName = str_replace(str_remove(cc$HBName,"NHS "),"&","and")
# Don't need all Scotland data
cc %>% filter(HBName != "Scotland") -> scotcc

HBs = unique(scotcc$HBName)

# Assemble a data frame

pidf = NULL
for(hb in HBs){
  ii = which(hb == scotcc$HBName)
  dfd = data.frame(Date=scotcc$Date[ii],
                   DailyPositive=scotcc$DailyPositive[ii],
                   HBName = scotcc$HBName[ii],
                   HB = scotcc$HB[ii])
  dfd$rollcases = rollmeanr(dfd$DailyPositive,k=7,fill="extend")
  pidf = bind_rows(pidf,dfd)
}
scotcc=pidf

#
# Population estimates
#

sp = read_csv(file.path(dataloc,"scotpop.csv"))

scotcc %>% select(HB,HBName) %>% 
  group_by(HB) %>% summarise(HBName=first(HBName)) -> shb

# Only use 2020 and aggregate pop not broken down by age and sex
# S92000003 is all of Scotland

sp %>% 
  filter(Year == 2020 & Sex == "All" & HB != "S92000003") %>% 
  select(HB,AllAges) %>% 
  rename(hbpop = AllAges) %>% 
  left_join(shb) -> scotpop

# Merge the RNA and case data estimates of population

scotrna %>% group_by(site) %>% 
  summarise(hab=first(hab), pop=as.numeric(first(pop))) %>% 
  group_by(hab) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup() %>% 
  rename(HBName = hab, wwpop = pop) -> wwpop

left_join(scotpop,wwpop) %>% 
  mutate(covfrac = wwpop/hbpop) -> scotpop

#' Start and end dates

mindate = as_date("2020-08-01")
maxdate = max(scotrna$Date)

#' Fill out df with all days in the range 
#' Use LOCF to avoid use of future values with interpolation
#' values within each site

scotrna %>% group_by(site) %>% 
  select(site,Date,N1) %>% 
  filter(Date >= mindate) %>% 
  complete(Date = seq.Date(mindate, maxdate, by="day")) %>% 
  mutate(N1 = zoo::na.locf(N1, na.rm=FALSE)) %>% 
  ungroup() -> scotfill

#' Merge in the HAB and population info

scotrna %>% group_by(site) %>% 
  summarise(hab=first(hab), pop=first(pop)) %>% 
  ungroup() -> habpop

scotfill = left_join(scotfill, habpop)
scotfill$pop = as.numeric(scotfill$pop)

#'
#' Compute weighted averages
#'

#' weighted average function
wavf = function(x,w){
  if(all(is.na(x))) return(NA)
  sum(x*w, na.rm = TRUE)/sum(as.numeric(!is.na(x))*w, na.rm = TRUE)
}

habnames = unique(scotrna$hab)

rnahab = NULL

for(nhab in habnames){
  scotfill %>% filter(hab == nhab) %>% 
    select(-c(hab,pop)) %>% 
    pivot_wider(names_from = site, values_from = N1) -> nw
  
  scotfill %>% filter(hab == nhab) %>% 
    group_by(site) %>% 
    summarise(pop=first(pop)) -> pophab
  
  wn1 = rep(NA, nrow(nw))
  for(i in 1:nrow(nw)){
    wn1[i] = wavf(nw[i,-1],pophab$pop)
  }
  
  rnahab = bind_rows(rnahab,tibble(Date = nw$Date, hab = nhab, N1 = wn1 ))
}

#' Join the rna and case data

scotcc %>% select(Date,HBName,DailyPositive,rollcases) %>% 
  rename(hab = HBName, cases = DailyPositive) -> sc

left_join(rnahab, sc) -> allscot

#' Save the results

save(scotrna,scotcc,scotpop,allscot,file = here("data/scot.rda"))

