
# Toolbox of modelling COVID-19 tools 


[![GitHub license](https://img.shields.io/badge/License-Creative%20Commons%20Attribution%204.0%20International-blue)](https://github.com/ec-jrc/COVID-19/blob/master/LICENSE)
[![GitHub commit](https://img.shields.io/github/last-commit/ec-jrc/COVID-19)](https://github.com/ec-jrc/COVID-19/commits/master)

## Introduction
The objective of this section is to illustrate the various elements of the scenario modelling that have been provided in the report:  'Scenarios for targeted and local COVID-19 Non Pharmaceutical Intervention Measures' , developed by JRC and provided in this folder. 

## Tools and folder structure
The following files are present in the main folder of the toolbox
* scenarioReport.docx - scenario report v.1.0  (** to be included)


The following folders are contained in this folder

- [scenario forecast](https://github.com/ec-jrc/COVID-19/tree/master/Modelling%20Scenarios/scenario%20forecast%20folders) folder that contains
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

![Example for Italy](https://github.com/ec-jrc/COVID-19/blob/master/Modelling%20Scenarios/scenario%20forecast%20folders/exampleCalibration.png?raw=true)

As it can be seen from the plots below during the calibration period (27/7-27/8) the quality of the forecast is rather good.  In the following period, after 27/8 in some cases the curve respect quite well the observed behaviour which means that no difference in the mitigation of the outbreak.  It should be noted that the case considered below is without any restriction action after the calibration period (e.g. no control case).
![forecast Italy](https://github.com/ec-jrc/COVID-19/blob/master/Modelling%20Scenarios/scenario%20forecast%20folders/forecasts/1000000_1000000_control_at_0_0.95_BYREGION/000/Italy_newInfe.jpg?raw=true)

## Regional Analyses 
(to be compiled by E.2)

## Mobility Functional Areas (MFAs):
MFAs are built starting from the daily Origin-Destination-Matrix (ODM) by considering the relative distribution of movements from each given origin to all possible destinations. This distribution is truncated so that only relative movements above 15% (mode = 1, or 30% mode = 2) remain. All positive entries are set to 1 and this new matrix P represents a proximity matrix of connected regions. Then, an adjacency matrix A is obtained through the formula A = 0.5 * (P * P’). Each element of the matrix A[i,j] is either 0, if regions i and j are not connected, A[I,j] = 0.5 if there are movements only from i to j or from j to i, and, A[I,j] =1 if the connection between i and j is bilateral. 

The A matrix is then used to create a network on a graph of connectivity and the daily MFAs are the result of a clustering algorithm on this graph. As the MFAs can vary from day to day, we look for a common set of clustering for the days before and during the lockdown. These two sets are used to identify the persistent (or pre-lockdown) MFAs and the lockdown MFAs taking all the areas which appears at least 50% of the times in the same cluster. More details can be found here.

### Anomaly Detection Dashboard:
The Anomaly detection system is an on-line tool and dashboard to detect abrupt changes of mobility from or to regions. It is based on a simple change point analysis strategy. Anomalies are detected as follows. Each cell of the ODM matrix, say X(t) = OD[i,j](t), observed through time, represents a time series for the connectivity from region i to region j. Looking at the moving average ma[i,j] of the previous values of the X(t) at week t-1, t-2, t-3 and t-4, any excess of X(t) of more than two standard deviations from the moving average m[i,j] is considered an anomaly. Excesses are then classified as signal of level 1, 2 and 3 if the relative (absolute) increment of X(t) with respect to ma[i,j] is below 50%, from 50 to 100%, above 100%.

### Mobility Visualisation Platform
The Mobility Visualisation Platform provides access to mobility data products over 22 MSs and Norway. In particular:
-	“Mobility Indicators”<sup>1</sup>  to quantify the impact of adopting/lifting restrictive measures through mobility variations in almost real-time;
-	“Connectivity”<sup>2</sup>  between regions predict/understand dynamics in early phases of the epidemic and monitor the effects of outward/inward mobility restriction measures;
-	“Mobility Functional Areas”<sup>3</sup>  of highly interconnected regions that can be the target of selective intelligent measures, providing a balance between epidemiological effects and socio-economic impacts.

Access to the platform and products (i.e. indicators, connectivity and Mobility Functional Area) is provided upon request to local, regional and national practitioners and policymakers in the Member States and the Commission. 

<sub>
(1) Santamaria, C., Sermi, F., Spyratos, S., Iacus, S. M., Annunziato, A., Tarchi, D., and Vespe, M. (2020). Measuring the impact of covid-19 confinement measures on human mobility using mobile positioning data. a european regional analysis. Safety Science, 132:104925.

(2) Iacus, S. M., Santamaria, C., Sermi, F., Spyratos, S., Tarchi, D., and Vespe, M. (2020a). Human mobility and covid-19 initial dynamics. Nonlinear Dynamics

(3) Iacus, S. M., Santamaria, C., Sermi, F., Spyratos, S., Tarchi, D., and Vespe, M. (2020b). Mapping mobility functional areas (MFA) using mobile positioning data to inform covid-19 policies, JRC121299.
</sub>
