
# Toolbox of modelling COVID-19 tools 


[![GitHub license](https://img.shields.io/badge/License-Creative%20Commons%20Attribution%204.0%20International-blue)](https://github.com/ec-jrc/COVID-19/blob/master/LICENSE)
[![GitHub commit](https://img.shields.io/github/last-commit/ec-jrc/COVID-19)](https://github.com/ec-jrc/COVID-19/commits/master)

## Introduction
The objective of this section is to illustrate the various elements of the scenario modelling that have been provided in the report:  'Scenarios for targeted and local COVID-19 Non Pharmaceutical Intervention Measures' , developed by JRC and provided in this folder. 

## Tools and folder structure
The following files are present in the main folder of the toolbox
* scenarioReport.docx - scenario report v.1.0  (** to be included)


The following folders are contained in this folder

- [scenario forecast](https://github.com/ec-jrc/COVID-19/tree/master/Modelling%20Scenarios/scenario%20forecast%20folders) folders that contains
  - [graphFinals.xlsx](https://github.com/ec-jrc/COVID-19/tree/master/Modelling%20Scenarios/graphFinal.xlsx) -  contains an EXCEL file with all the plots and data contained in the scenario report
  - [Epidemic Modelling v0710](https://github.com/ec-jrc/COVID-19/tree/master/Modelling%20Scenarios/scenario%20forecast%20folders/forecasts)  - calibration program in python
  - [foreProg](https://github.com/ec-jrc/COVID-19/tree/master/Modelling%20Scenarios/scenario%20forecast%20folders/foreProg)  - forecast programs in python to estimate the behaviour in 6 months. It uses the output created by the Epidemic Modelling v0710 programmes
  - [Forecasts](https://github.com/ec-jrc/COVID-19/tree/master/Modelling%20Scenarios/scenario%20forecast%20folders/forecasts) - all forecasts folders produced for the scenario report
- Regional Analyses (to be compiled by E.2)
- Mobility Analysis Tools  (to be compiled by E.6)

## Scenario Forecast Model Description
In order to understand the possible evolution of the current situation, a SIR model has been setup, calibrated for the period 27 July-27 August 2020 and with the conditions frozen.
For the simulation a simple SIR model was adopted (Susceptible, Infected, Recovery), in which the only parameters to be calibrated are the Rt and Trecovery:
>       dSdt=-Rt/Trecov*S*I/N
>       dIdt=Rt/Trecov*S*I/N - 1/Trecov*I
>       dRdt=1/Trecov*I
where N is the overall population, S is the susceptible population, I are the infected individuals and R are the recovered. Rt is the Reproduction number and Trecov the recovery tim, both obtaned by the calibration process.
In this simulation the following hypotheses have been adopted
-	Initial situation based on the calibration of a SIR model applied to each region in Europe, based on the 1 month of data between 26 July and 26 Aug 2020
-	R0 and Trecov maintained constant in the following 6 months for each region, unless the conditions for lockdown or partial lockdown are met
-	When the release conditions are reached the Rt returns to the value had during the calibration phase

The calibration for 30 days allows to identify the current parameters of the SIR model. As an example, Figure 1 shows the calibration of Italy: dotted blue line are the observation, the orange line are the prediction. The model has been applied to each region in Italy, Spain and Germany while a country model was used for France because the regional epidemiological data are not available. An example of fitting quality is shown below for Italy.

(Include image for Italy)

As it can be seen from the plots below during the calibration period (27/7-27/8) the quality of the forecast is rather good.  In the following period, after 27/8 in some cases the curve respect quite well the observed behaviour which means that no difference in the mitigation of the outbreak.

## Regional Analyses 
(to be compiled by E.2)

## Mobility Analysis Tools  
(to be compiled by E.6)
