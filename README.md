# Max's COVID-19 Prediction Script
This script attempts to predict a lower bound on the impact of COVID-19, based off of currently-available data.

# Method
Predictions are made by fitting a sigmoid function to the reported death cases. Fitting a sigmoid function to data is a non-linear optimization problem. This problem 
is solved in two steps:
1. An exhaustive search over the three parameters describing the sigmoid function.
2. A conjugate-gradient optimization of those three parameters, using the result from the exhaustive search as an initial guess

## Disclaimer
This script assumes that people will maintain a constant level of infection-control measures. These measured include
* closing schools and workplaces
* closing restaurants and bars
* minimizing public trips and social events

That assumption is not true. All of these measures will eventually stop. Therefore, the predicted death tolls reported here are a **lower bound** on the total number of deaths that will be observed.

## Graphs
See graphs for all qualifying countries [here](./Latest)
![Cases Since Threshold](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/SinceThresh.png)
![Deaths By Date](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/ConfDeaths.png)

![United_States_of_America](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/United_States_of_America.png)
![Italy](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Italy.png)
![France](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/France.png)
![China](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/China.png)
![ChinaAdjusted](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/China_Adjusted.png)

![Australia](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Australia.png)
![Belgium](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Belgium.png)
![Canada](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Canada.png)
![India](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/India.png)
![Iran](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Iran.png)
![Japan](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Japan.png)
![Norway](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Norway.png)
![Portugal](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Portugal.png)
![South_Korea](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/South_Korea.png)
![Spain](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Spain.png)
![Switzerland](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Switzerland.png)
![Turkey](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/Turkey.png)
![United_Kingdom](https://raw.githubusercontent.com/hughsonm/covid-prediction/master/Latest/United_Kingdom.png)

## Requirements
* MATLAB
* Internet connection

## Instructions
1. Download the latest COVID-19 statistics from [opendata](https://opendata.ecdc.europa.eu/covid19/casedistribution/json)
In a terminal, type
```bash
wget -O covid.json https://opendata.ecdc.europa.eu/covid19/casedistribution/json/
```
2. Run PlotCovidStats.m in MATLAB



## TODO
* Make this run in Octave and MATLAB
* Incorporate the wget into the MATLAB script
* Find a suitable 'goodness-of-fit' statistic. Currently using r-squared, which is not great for non-linear models
* Further improve the cost functional minimizer
