# Max's COVID-19 Prediction Script
This script attempts to predict the impact of COVID-19, based off of currently-available data.

## Requirements
* MATLAB
* Internet connection

## Instructions
1. Download the latest COVID-19 statistics from [opendata](https://opendata.ecdc.europa.eu/covid19/casedistribution/json)
In a terminal, type
```
wget -O covid.json https://opendata.ecdc.europa.eu/covid19/casedistribution/json/
```
2. Run PlotCovidStats.m in MATLAB

## TODO
* Make this run in Octave and MATLAB
* Incorporate the wget into the MATLAB script
* Find a suitable `goodness-of-fit' statistic. Currently using r-squared, which is not great for non-linear models
* Further improve the cost functional minimizer
