# Challenges in realising the potential of wastewater-based epidemiology to quantitatively monitor and predict the spread of disease

This repository provides code and links to data supporting the published
article <https://doi.org/10.2166/wh.2022.020>

## Abstract

Researchers around the world have demonstrated correlations between measurements of SARS-CoV-2 RNA in wastewater and case rates derived from direct testing of individuals. This has raised hopes that WBE methods might be used to quantify the spread of disease, perhaps faster than direct testing, and with less expense and intrusion. We illustrate, using data from Scotland and the USA, the issues regarding the construction of effective predictive models for disease case rates. We discuss the effects of variation in, and the problem of aligning, public health reporting and wastewater measurements. We investigate time-varying effects in public health-reported case rates and their relationship to wastewater measurements. We show the lack of proportionality of wastewater measurements to case rates with associated spatial heterogeneity. We illustrate how the level of aggregation chosen affects the precision of predictions. We determine whether public health or wastewater measurements are the leading indicator of disease and how they may be used in conjunction to produce predictive models. The prospects of using wastewater based predictive models with or without ongoing public health data are discussed.

## Data

- Scottish RNA data can be found at
<https://informatics.sepa.org.uk/RNAmonitoring
- Case data can be downloaded from
<https://www.opendata.nhs.scot/dataset/weekly-covid-19-statistical-data-in-scotland>
- ONS Scottish infection survey data from <https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/covid19infectionsurveyscotland/2021>
- USA data from <https://github.com/biobotanalytics/covid19-wastewater-data>
and displayed at <https://biobot.io/data/>.

CSVs derived from these are in the [data directory](data)

## Code

- Process Scottish cases and RNA counts: [R script](scotland/scotproc.R)
- Modelling Scottish data: [Rscript](scotland/scotmodel.R) and [output](scotland/scotmodel.md)
- Process the USA data [R script](usa/usaproc.R)
- Model the USA data [R script](usa/usamodel.R) and 
[output](usa/usamodel.md)


