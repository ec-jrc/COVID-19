#!/usr/bin/env python
# coding: utf-8

import numpy as np

# todo add support for population normalization

from data_handlers.dhWebcritech import dhWebcritech, get_geofilter_Webcritech
from models.SIR import SIR
from models.SEIR import SEIR
from models.SIR1 import SIR1
from models.SIRD import SIRD
from models.SIRD1 import SIRD1
from models.SIRD2 import SIRD2
from models.SEIRD import SEIRD
from optimization.optiLMFIT import optiLMFIT

d = {"SIR": SIR1, "SIRD": SIRD1, "SIRD2": SIRD2, "SEIRD": SEIRD}


def main_fitting(area_type, reference_area, model_type, analysis_horizon, make_plots,
                 default_params={}, var_map_nondefault={}, variables_to_plot=[]):
    '''
    Performs model fitting in three main steps: gathering of epidemiological data; model setup; model fitting.
    :param file_epidemiology:           filename for epidemiological information source
    :param reference_area:              name of the analysis area, to be matched with file_demography and file_epidemiology
    :param model_type:                  string defining the analysis model to exploit, e.g."SIR"
    :param backward_analysis_horizon:   backward time horizon to consider for fitting, starting from the most recent epidemiological data
    :param make_plots:                  flag indicating whether or not to plot the results of the processing steps
    :param default_params:              parameters of a model we want to assign a priori and exclude from the optimization problem

    :return:

    Typical analysis pipeline:
    [0] load epidemiological data (as an object), e.g.
    >>> epidemiology = dhWebcritech(area_type, filters={"CountryName":"Italy"})
    for loading a Webcritech dataset (source: file_epidemiology), filtering information based on field CountryName, and
    attaching demographic information associated to the reference area
    [1] set up the epidemiological model (as an object), e.g.
    >>> PM={'beta':0.3,'gamma':1/7}
    >>> PIVP={"I":100,"R":0}
    >>> mymodel = SIR(parameters_model=PM,parameters_ivp=PIVP,population=epidemiology.population,observations=epidemiology.data)
    where
    - PM describes the model parameters
        (for SIR, for instance, these are described by class attribute PARAMETERS in class definition, see SIR.py)
    - PIVP are the initial conditions associated to the initial value problem, specified in "named variables"
        (i.e. tagged with populations names, e.g. "S"; see class attribute VARIABLES, e.g. in SIR.py)
    - SIR is the model we selected for our analysis, see subfolder /models for alternatives
    Observe that PIVP should not include what - in class definition - we denoted as VARIABLE_PIVOT, e.g. VARIABLE_PIVOT="S" in SIR.
    This represents the subpopulation that will be initialized as (total population - sum of the rest of subpopulations).
    * This pivot mechanism allows to auto-adjust population counts when specifying initial conditions via the optimization problem.
    [2] Model fitting, e.g.
    >>> VTF = ["I","R"]
    >>> myoptim=optiLMFIT(model=mymodel,variables_to_fit=VTF)
    where
    - VTF defines the variables with respect to which we want to perform our fitting
    - optiLMFIT is the optimization utility we prepared for fitting
    '''

    mymodel = []
    myoptim = []

    ##########################
    # Epidemiology information
    # Observations:
    # - epidemiology dataset is filtered based on analysis_horizon
    # - we support dataset smoothing
    print("Epidemiological data is processed.")

    epidemiology = dhWebcritech(area_type=area_type, filters={get_geofilter_Webcritech(area_type): reference_area})

    if 'Susceptible' not in epidemiology.data:
        print("Epidemiological/demographic information insufficient to proceed further.")
        return epidemiology, mymodel, myoptim

    epidemiology.data = epidemiology.data.iloc[analysis_horizon[0]:analysis_horizon[1]]

    ##########################
    # Model initialization
    # Observations:
    # - we use a model mapping, specifying the source class used per model type
    # - inside object mymodel,
    #       observations equal epidemiology.data converted into the named variables of the model (e.g. "S", "I", "R")
    #       predictions are the ODE-IVP solutions converted into named variables; model solution is auto-triggered when accessing mymodel.predictions
    # - at init, it is not necessary to fully specify model parameters and ivp, as these will be handled by the fitting routine

    # d = {"SIR": SIR1, "SIRD": SIRD1, "SIRD2": SIRD2, "SEIRD": SEIRD}

    print("Model is set up.")
    print("Model in use is " + model_type + ", exploiting mapping: " + str(d))

    mymodel = d[model_type](reference_area=epidemiology.reference_area, population=epidemiology.population,
                            observations=epidemiology.data, parameters_model_static=default_params,
                            var_map_nondefault=var_map_nondefault, variables_to_plot=variables_to_plot)

    ##########################
    # Model fitting
    # Observations:
    # - optimization variables are specified inside mymodel (and exclude the pivot variable)
    # - variables_to_fit defines the terms included in residual calculation, and must refer to terms represented in
    #   mymodel.observations and mymodel.predictions
    # - optimal solution will be embedded in myoptim.model object

    # # based on discussion with Alessandro:
    # # variables_to_fit = ["S", "delta_S", "D"]

    # # Alessandro's latest version:
    # if "D" in model_type:
    #     variables_to_fit = ["delta_S", "S", "D"]
    # else:
    #     variables_to_fit = ["delta_S", "S", "R"]

    # Alternative
    # Here we assume that delta_S is always in (new cases) and that we consider all the mapped variables
    #variables_to_fit = ["delta_S"] + list(mymodel.VAR_MAP.keys())
    variables_to_fit = ["delta_S","S"] #+ list(mymodel.VAR_MAP.keys())

    # we assign weights for the residual terms
    fitting_weights = {var: 1 for var in variables_to_fit}

    print("Model is fitted to epidemiological data.")
    print('variables_to_fit=', variables_to_fit)

    # # round 1 optimization
    myoptim = optiLMFIT(model=mymodel, variables_to_fit=variables_to_fit, fitting_weights=fitting_weights,
                        optimization_method="least_squares")
    # # further optimization rounds, do e.g.
    # # mymodel.set_limits(limit_entry="value", input_values=myoptim.get_optimum_parameters())
    # # myoptim = optiLMFIT(model=mymodel, variables_to_fit=variables_to_fit, optimization_method="dual_annealing")
    # # NOTE: list_of_methods = ["leastsq","least_squares","differential_evolution","brute","basinhopping","ampgo","`nelder","lbfgsb","powell","cg","cobyla","bfgs","tnc",

    if make_plots:
        myoptim.model.plot_intime()

    return epidemiology, mymodel, myoptim


def main_forecast(model_name, reference_area, population, parameters_model, parameters_ivp, time):
    '''
    Performs forecasting based on a specified model and IVP configuration
    :param model_name:
    :param reference_area:
    :param population:
    :param parameters_model:
    :param parameters_ivp:
    :param time:
    :return:
    '''
    model_class = eval(model_name)
    mymodel = model_class(reference_area=reference_area, population=population, parameters_model=parameters_model,
                          parameters_ivp=parameters_ivp, time=time)

    return mymodel


if __name__ == "__main__":

    # MODEL_TYPE = "SIR"
    # ANALYSIS_HORIZON = (-10, -1)
    #
    # MODEL_TYPE = "SIR"
    # ANALYSIS_HORIZON = [-30, -1]
    # AREA_TYPE = "Country"
    # REFERENCE_AREA = "Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden".split(
    #     ', ')
    #
    # for a in REFERENCE_AREA:
    #     print('Analysing ' + a)
    #
    #     epidemiology, mymodel, myoptim = main_fitting(
    #         area_type=AREA_TYPE, reference_area=a,
    #         model_type=MODEL_TYPE, analysis_horizon=ANALYSIS_HORIZON, make_plots=True)
    #
    # # mymodel = SIRD(population=60000000, time=np.arange(20),
    # #                parameters_model={"beta": 0.17931595659932845, "gamma": 0.019688397759320023,
    # #                                  "mu": 0.015577455120231655},
    # #                parameters_ivp={"I": 2521.700187106392, "R": 69.79346598667506, "D": 1.7660613590794605})
    # # mymodel.plot_intime()

    # simple numerical experiment
    import datetime

    base = datetime.datetime(2000, 1, 1)
    time = np.array([base + datetime.timedelta(days=i) for i in range(366)])

    mymodel = SIR1(population=60000000, time=time,
                   parameters_model={"r0": 2.5, "Trecov": 20},
                   parameters_ivp={"I": 100, "R": 0})
    mymodel.VARIABLES_TO_PLOT = ["S", "delta_S", "I", "R"]
    mymodel.plot_intime()

    print("DONE")
