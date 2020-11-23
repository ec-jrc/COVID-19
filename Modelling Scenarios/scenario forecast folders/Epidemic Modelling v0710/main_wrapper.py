#!/usr/bin/env python
# coding: utf-8
'''
~~~~~~~~~~~~
JRC Ispra
COVID-19 project, June 2020
~~~~~~~~~~~~
[WORK IN PROGRESS]

Python implementation of a family of epidemiological models and associated analysis instruments.

We currently include:
- /_data (demographic and epidemiological datasets) => see also init.py (data fetching utility)
- /_geospatial (geospatial datasets)
- /_report (output files useful for reporting)
- /data_handlers (for the manipulation of input data, containing e.g. epidemiological and demographic information)
- /doc (utilities for documentation generation)
- /models (an array of epidemiological models)
- /optimization (currently devoted to model fitting)
- /util (general-purpose utilities, e.g. for plotting)
'''

# todo automatically manage regions/countries wherein R is lacking

import datetime
import itertools
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_handlers.areanameGetter import areanameGetter
from doc.docgenerator import docgenerator
from doc.get_metadata import get_study_filestring, get_study_metadata
from doc.get_report_entry import get_report_entry
from main import main_fitting, main_forecast

from CONF import DIR_REPORTING, COUNTRIES_SUBSET, COUNTRIES_ALL, FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES,FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS,PWD
import os
import sys,getopt

##########################
##########################
# [START OF] SET ANALYSIS PARAMETERS

STUDY_TYPE = "FITTING_ALLREGIONS"
MODEL_TYPE = "SIR"
ANALYSIS_HORIZON = (-30, -30+29)
#ANALYSIS_HORIZON = (-260, -260+29)

# OTHER STUDY_TYPE SUPPORTED


# [END OF] SET ANALYSIS PARAMETERS
##########################
##########################

# defaults
COUNTRY_FILTER = ""
DEFAULT_PARAMS = {}
VAR_MAP_NONDEFAULT = {}
VARIABLES_TO_PLOT = []

AREA_TYPE = "Country"
# get list of iterables to combine a produce the analysis batch; iterable_var lists the names of variables to scan, if they are given as lists
iterable_var = ["STUDY_TYPE", "MODEL_TYPE", "ANALYSIS_HORIZON", "REFERENCE_AREA"]
detach = False
makereport = False #True
# download latest data
url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-country/jrc-covid-19-all-days-by-country.csv'

if not os.path.exists(FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES):
    cmd='curl -o "'+ FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES+ '" '+url
    print (cmd)
    os.system(cmd)

url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
if not os.path.exists(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS):
    cmd='curl -o "'+ FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS+ '" '+url
    print (cmd)
    os.system(cmd)
#for firstDay in range(-260, 0, 7):   #to have all the forecasts since march


if len(sys.argv)>1:
    try:
        opts, other = getopt.getopt(sys.argv[1:], "c:s:", ["calibration=",30,"start=",''])

    except getopt.GetoptError:
        print
        'test.py -cp <calibrationPeriod> -o <outputfile>'
        sys.exit(2)
    
    for opt,arg in opts:
        if opt=='-c':
            calibrationPeriod=int(arg)
        elif opt=='-s':
            start=int(arg)

else:
    calibrationPeriod=15  #  portare a 30
    start=0

print("Calibration period=",calibrationPeriod,' start=',start)


for firstDay in range(-calibrationPeriod+start,-calibrationPeriod+start+1):

    ANALYSIS_HORIZON = (firstDay, firstDay+calibrationPeriod-1)
    
    lastDay=datetime.datetime.utcnow()+datetime.timedelta(days=ANALYSIS_HORIZON[1]+1)
    print("Start "+format(ANALYSIS_HORIZON[0]),lastDay)
    DIR_REPORTING = PWD+ "/_report/"+lastDay.strftime('%Y-%m-%d')
    if not os.path.isdir(DIR_REPORTING):
        os.makedirs(DIR_REPORTING)

    for STUDY_TYPE in ["FITTING_SELECTEDCOUNTRIES","FITTING_ALLCOUNTRIES","FITTING_ALLREGIONS"]:
    #for STUDY_TYPE in ["FITTING_SELECTEDCOUNTRIES"]:
        fname1=DIR_REPORTING+'/FITTING_ALLCOUNTRIES/data/FITTING_ALLCOUNTRIES_SIR_(-'+format(calibrationPeriod)+' -1).csv'
        fname2=DIR_REPORTING+'/FITTING_ALLREGIONS/data/FITTING_ALLREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).csv'
    
        print(fname1,fname2)
        if (STUDY_TYPE=='FITTING_ALLCOUNTRIES' or STUDY_TYPE=='FITTING_SELECTEDCOUNTRIES') and os.path.exists(fname1): continue
        if STUDY_TYPE=='FITTING_ALLREGIONS' and os.path.exists(fname2): continue
    
        if STUDY_TYPE == "FITTING_SELECTEDCOUNTRIES":  # "COUNTRIES_SUBSET":
            fname3=DIR_REPORTING+'/FITTING_SELECTEDCOUNTRIES/data/FITTING_SELECTEDCOUNTRIES_SIR_(-'+format(calibrationPeriod)+' -1).csv'
            if os.path.exists(fname3): 
                continue
            # REFERENCE_AREA (string/list) => leave [] for all countries, or specify e.g. ["Italy","Germany"]
            REFERENCE_AREA = COUNTRIES_SUBSET

        elif STUDY_TYPE == "FITTING_SELECTEDREGIONS":

            AREA_TYPE = "Region"
            # region-wide study => either specify a list of REFERENCE_AREA, or specify a country filter
            REFERENCE_AREA = []  #["Schleswig Holstein"]
            COUNTRY_FILTER = "" #"COUNTRIES_SUBSET"  #"" #"Germany" #"Spain"  # "Spain"

        elif STUDY_TYPE == "FITTING_ALLCOUNTRIES" or STUDY_TYPE == "FITTING_ALLREGIONS":

            SUP_ANALYSIS = "FITTING_SELECTEDCOUNTRIES"

            if STUDY_TYPE == "FITTING_ALLCOUNTRIES":
                REFERENCE_AREA = COUNTRIES_ALL
                COUNTRY_FILTER="Belgium"

            elif STUDY_TYPE == "FITTING_ALLREGIONS":
                AREA_TYPE = "Region"
                REFERENCE_AREA = []
                COUNTRY_FILTER ="" # "Germany"

            # load dataset from (previously-run) same-type analysis for COUNTRIES_SUBSET
            try:
                filename_in = DIR_REPORTING + "/" + SUP_ANALYSIS + "/data/" + \
                              get_study_filestring(study_type=SUP_ANALYSIS, model_type=MODEL_TYPE,
                                                   analysis_horizon=ANALYSIS_HORIZON, reference_area=COUNTRIES_SUBSET) + ".csv"
                df = pd.read_csv(filename_in)
            except:
                raise Exception("Corresponding FITTING_SELECTEDCOUNTRIES/FITTING_ALLREGIONS-level analysis should be run before")
            # attribute non-default configuration to model, in order to bypass native variable mapping et cetera
            DEFAULT_PARAMS = df.mean().to_dict()
            DEFAULT_PARAMS.pop('r0', None)  # todo generalize if relevant
            VAR_MAP_NONDEFAULT = {"S": ["Susceptible"], "I+R": ["CumulativePositive"]}
            VARIABLES_TO_PLOT = ["S", "delta_S", "I+R"]

        elif STUDY_TYPE == "FORECASTING_COUNTRIES" or STUDY_TYPE == "FORECASTING_REGIONS":

            FORECAST_HORIZON = 180
            FORECAST_FACTORS = [{"r0": x} for x in 0.5 + np.arange(11) / 10]

            if STUDY_TYPE == "FORECASTING_COUNTRIES":
                SUP_ANALYSIS = "FITTING_ALLCOUNTRIES"
                FILENAME_IN_REFERENCEAREA = COUNTRIES_ALL
            elif STUDY_TYPE == "FORECASTING_REGIONS":
                SUP_ANALYSIS = "FITTING_ALLREGIONS"
                FILENAME_IN_REFERENCEAREA = "Germany" # todo generalize

            FILENAME_IN = DIR_REPORTING + "/" + SUP_ANALYSIS + "/data/" + \
                          get_study_filestring(study_type=SUP_ANALYSIS, model_type=MODEL_TYPE,
                                               analysis_horizon=ANALYSIS_HORIZON,
                                               reference_area=FILENAME_IN_REFERENCEAREA) + ".json"

            # we include FORECAST_FACTORS among the set of variables to loop over
            iterable_var += ["FORECAST_FACTORS"]

            # REFERENCE_AREA is not significant here as it is derived from baseline scenarios used for forecasting
            REFERENCE_AREA = None
            # detach = True means that this STUDY_TYPE will produce separate output files for each analysis
            detach = True
            makereport = False

        else:

            raise Exception("Select supported STUDY_TYPE.")


        REFERENCE_AREA, SUPER_REFERENCE_AREA = areanameGetter(REFERENCE_AREA, area_type=AREA_TYPE,
                                                              country_filter=COUNTRY_FILTER)

        # ************
        print("")
        print("* " * 12)
        print("STARTING THE ANALYSIS BATCH FOR")
        for iv in iterable_var:
            print("\t" + iv + ": " + str(eval(iv)) + ".")
        # print("\tSTUDY TYPE: " + str(STUDY_TYPE) +".")
        # print("\tMODEL TYPE: " + str(MODEL_TYPE) +".")
        # print("\tREFERENCE AREA: " + str(REFERENCE_AREA) +".")


        # init list useful to collect reports of the analysis process
        report = {}
        report["metadata"] = get_study_metadata(STUDY_TYPE, MODEL_TYPE, ANALYSIS_HORIZON, REFERENCE_AREA, SUPER_REFERENCE_AREA)

        # my_locals = locals()
        # iter_vars = [x for x in iterable_var if x in my_locals and isinstance(eval(x), list)]
        my_locals = locals()
        iter_vars = [x for x in iterable_var if x in my_locals and isinstance(eval(x), list)]


        if len(iter_vars) == 1:
            iter_vars_comb = eval(iter_vars[0])
        else:
            iter_vars_vals = [eval(x) for x in iter_vars]
            iter_vars_comb = list(itertools.product(*iter_vars_vals))

        for entry in iter_vars_comb:

            # init metadata
            entry_metadata = {}
            analysis_outcomes = {}
            analysis_output = {}

            # assign variables based on currently selected combination of iterables
            if len(iter_vars) == 1:
                entry_metadata[iter_vars[0]] = entry
                if isinstance(entry, str):
                    curr_entry = '"' + entry + '"'
                else:
                    curr_entry = str(entry.replace("'","\'"))
                exec(iter_vars[0] + " = " + curr_entry)
            else:
                for ii in range(len(entry)):
                    entry_metadata[iter_vars[ii]] = entry[ii]
                    if isinstance(entry[ii], str):
                        curr_entry = "'" + entry[ii] + "'"
                    else:
                        curr_entry = str(entry[ii])
                    exec(iter_vars[ii] + " = " + curr_entry)

            print("")
            print("*" * 5)
            print(entry_metadata)

            # #  append map descriptors to metadata
            # # todo move geodata fetching to reporting modules
            # entry_metadata["map"] = geodataGetter(location=REFERENCE_AREA, super_location=SUPER_REFERENCE_AREA).plot()

            # execute the core routine
            if "FITTING" in STUDY_TYPE:

                epidemiology, mymodel, myoptim = main_fitting(
                    area_type=AREA_TYPE, reference_area=REFERENCE_AREA, model_type=MODEL_TYPE,
                    analysis_horizon=ANALYSIS_HORIZON,
                    make_plots=False, default_params=DEFAULT_PARAMS, var_map_nondefault=VAR_MAP_NONDEFAULT,
                    variables_to_plot=VARIABLES_TO_PLOT
                )

                try:
                    # reporting information, if fitting produced a significant output
                    analysis_outcomes["optimum"] = myoptim.get_optimum_parameters()
                    analysis_outcomes["optimum"].pop('_ivp', None) #json.dumps(myoptim.get_optimum_parameters().pop('_ivp', None), indent=4)
                    # analysis_outcomes["figures"] = {}
                    analysis_outcomes["fig_fitting"] = myoptim.model.plot_intime()
                    # close all figures to avoid overloading
                    plt.close('all')
                    analysis_output["json"] = myoptim.model.get_params()
                    analysis_output["csv"] = myoptim.get_optimum_parameters()
                    analysis_output["csv"].pop('_ivp', None)
                except:
                    pass

            elif "FORECASTING" in STUDY_TYPE:

                filename_suffix = str(FORECAST_FACTORS)

                # load reference scenarios to start from;
                # restart at every loop from the original dataset
                try:
                    with open(FILENAME_IN) as f:
                        SCENARIOS = json.load(f)
                except:
                    raise Exception("Corresponding FITTING_ALLCOUNTRIES-level analysis should be run before")

                VARS_TO_STORE = ["S", "I", "R", "delta_S", "I+R"]
                df_ = {var: pd.DataFrame() for var in VARS_TO_STORE} #mymodel.VARIABLES_TO_PLOT}

                # prepare all scenarios
                for scenario in SCENARIOS.keys():

                    # elaborate time horizon => extend forward for forecasting purposes
                    SCENARIOS[scenario]["time"] = [datetime.datetime.strptime(x, '%Y/%m/%d') for x in
                                                   SCENARIOS[scenario]["time"]]
                    base = SCENARIOS[scenario]["time"][-1]
                    arr = np.array([base + datetime.timedelta(days=i) for i in range(1, FORECAST_HORIZON)]).tolist()
                    SCENARIOS[scenario]["time"] = SCENARIOS[scenario]["time"] + arr

                    # modify the baseline scenario loaded so far by applying multiplicative factors defined by FORECAST_FACTORS
                    for key in FORECAST_FACTORS.keys():
                        SCENARIOS[scenario][key] = SCENARIOS[scenario][key] * FORECAST_FACTORS[key]

                    # run scenario
                    model_name = SCENARIOS[scenario]["MODEL_NAME"]
                    reference_area = SCENARIOS[scenario]["reference_area"]
                    parameters_model = SCENARIOS[scenario]  # we rely on the fact that the model will only pick parameter-named entries
                    parameters_ivp = SCENARIOS[scenario]["_ivp"]
                    time = SCENARIOS[scenario]["time"]
                    population = SCENARIOS[scenario]["population"]

                    mymodel = main_forecast(model_name, reference_area, population, parameters_model, parameters_ivp, time)

                    mymodel.VARIABLES_TO_PLOT = VARS_TO_STORE # ["S", "I", "R", "delta_S", "I+R"]
                    # todo check if this is what we need

                    # # init of the data-collection dict of dataframes
                    # if "df_" not in locals():
                    #     df_ = {}
                    #     for var in mymodel.VARIABLES_TO_PLOT:
                    #         # here we assume that all iterations will cover the same compartments or subsets of the first one
                    #         df_[var] = pd.DataFrame()

                    # df_ = {var:pd.DataFrame() for var in mymodel.VARIABLES_TO_PLOT}


                    # appending values
                    for var in mymodel.VARIABLES_TO_PLOT:
                        df_[var][reference_area] = mymodel.predictions[var]

                # append totals
                for var in mymodel.VARIABLES_TO_PLOT:
                    df_[var]["total"] = df_[var].sum(axis=1) #.sum(axis=1)

                # reporting information
                try:
                    # reporting information, if fitting produced a significant output
                    analysis_outcomes = {}
                    analysis_output["xlsx"] = df_
                except:
                    pass

            # append data produced by current iteration to report, if significant information was produced
            if analysis_outcomes or analysis_output:
                loop_id = str(entry)  # todo generalize
                report[loop_id] = get_report_entry(entry_metadata, analysis_outcomes, analysis_output)

        # output report_entries, generating a report
        docgenerator(DIR_REPORTING, report, detach, makereport)

    print("Done "+format(ANALYSIS_HORIZON[0]))
print('****  Completed ****')