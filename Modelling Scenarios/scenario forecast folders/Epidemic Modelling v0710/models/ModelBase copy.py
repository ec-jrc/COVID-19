from pylab import *
from scipy.integrate import odeint
import warnings
import pandas as pd
from util.plotting import plot_intime_comparison


class ModelBase:
    '''
    Implements a baseline class for epidemiological models, including
    - variable check and remapping
    - ODE IVP solver
    - metrics computation
    - basic visualization

    Class parameters:
    PARAMETERS          list of model parameter names
    VARIABLES           list of model variable names
    VARIABLES_TO_PLOT   list of names of model variables to be included in plots
    VAR_MAP             dict mapping entries in VARIABLES to fields in property 'observations' useful to remap observations to model variables
    VARIABLE_PIVOT      name of the variable that is evaluated, at initialization, as difference between total population and rest of the compartments
    LIMITS              reference value/min/max of each variable, useful for fitting purposes
    '''

    PARAMETERS = []
    VARIABLES = []
    VARIABLES_TO_PLOT = []
    VAR_MAP = {}
    VARIABLE_PIVOT = []
    LIMITS = []

    VAR_OPTIM = []

    def __init__(self, verbose: bool = True, parameters_model: dict = {}, parameters_ivp: dict = {}, reference_area: str="", population=0,
                 observations=None, time=None, parameters_model_static={}):  # , init={}, t=[], t_dates=[]):
        '''
        Inits attributes
        :param verbose:             defines verbosity level, used e.g. in issuing warnings based on parameter checks
        :param parameters_model:    contains model parameters, dict keys should match PARAMETERS for a complete parameterization
        :param parameters_ivp:      contains the specification of the initial value problem, including
                                    - initial conditions, keyed as indicated in VARIABLES
                                    - the array of time instants at which IVP solutions should be returned, keyed as 't'
        :param parameters_ivp:      Pandas Dataframe contains the raw observations;
                                    its columns will be automatically remapped to model variables based on VAR_MAP matching
        '''

        # parameter checks
        self.checks(params=list(parameters_model.keys()) + list(parameters_ivp.keys()))
        # parameters_model => extract to individual attributes, preserving key names
        # vars(self).update(parameters_model)
        for key in self.PARAMETERS:
            if key in parameters_model.keys():
                setattr(self, key, parameters_model[key])
            else:
                setattr(self, key, None)

        self.VAR_OPTIM = self.PARAMETERS + ["_ivp"]

        for key in self.PARAMETERS:
            if key in parameters_model_static.keys():
                setattr(self, key, parameters_model_static[key])
                self.VAR_OPTIM.remove(key)


        # parameters_ivp => extract to individual attributes, preserving key names
        # vars(self).update(parameters_ivp)
        # self._ivp = parameters_ivp
        self._ivp = {}

        for key in self.VARIABLES:
            if key != self.VARIABLE_PIVOT:
                if key in parameters_ivp.keys():
                    self._ivp[key] = parameters_ivp[key]
                else:
                    self._ivp[key] = None

        # for key in self.VARIABLES:
        #     if key != self.VARIABLE_PIVOT:
        #         if key in parameters_ivp.keys():
        #             self._ivp[key] = parameters_ivp[key]
        #         else:
        #             self._ivp[key] = None

            # if key in parameters_ivp.keys() and key != self.VARIABLE_PIVOT:
            #     self._ivp[key] = parameters_ivp[key]
            # else:
            #     # todo check
            #     # self._ivp[key] = None
            #     pass
        # reference area
        self.reference_area = reference_area
        # population
        self._population = population
        # observations
        self._observations = observations

        # model name
        self.MODEL_NAME = self.getName()
        # indexes used in constructing the DataFrame for model predictions
        # self._indexes_for_predictions = []
        # absolute time indexes, useful e.g. for time-dependent models
        # self.time = []
        self._time = time

        # verbosity
        self._verbose = verbose

        # default checks
        if self._verbose:
            self.checks()


    def checks(self, params=None):
        '''
        Performs checks on data structures. In particular,
        - is params==None, performs standard checks on attributes;
        - if params are provided in input, checks their internal consistency and consistency with other variables.
        Throws errors or warnings on inconsistencies.
        '''

        # lambdas for checking conditions and issuing warnings/errors if needed
        check = lambda test, message: 0 if test == True else warnings.warn("check " + message + " failed")
        check_hard = lambda test, message: 0 if test == True else exec(
            'raise(Exception(' + "check " + message + " failed" + '))')

        if params is None:
            # name collisions
            check(len(self.PARAMETERS + self.VARIABLES) == len(set(self.PARAMETERS + self.VARIABLES)),
                  "PARAMETERS AND VARIABLES NO NAME COLLISIONS")
            # parameters
            check(all([hasattr(self, attr) for attr in self.PARAMETERS]), "ALL MODEL PARAMETERS ENTERED")
            # variables
            check(set(self.VARIABLES_TO_PLOT) <= set(self.VARIABLES), "VARIABLES_TO_PLOT \subseteq VARIABLES")
            check(set(list(self.VAR_MAP.keys())) <= set(self.VARIABLES), "VAR_MAP.keys() \subseteq VARIABLES")
            check(isinstance(self.VARIABLE_PIVOT, str), "VARIABLE_PIVOT IS STRING")
            # the pivot variable must not be present in _init, as it is computed from total population and rest of the subpopulations
            # check(not (self.VARIABLE_PIVOT in self._init), "VARIABLE_PIVOT NOT IN _init")
            # # the non-pivot variables and total population must be in _init
            # check(not (self.VARIABLE_PIVOT in self.VARIABLES_NONPIVOT), "VARIABLES_NONPIVOT IN _init")
        else:
            # name collisions
            check_hard(len(params) == len(set(params)), "INPUT DICTIONARIES NO NAME COLLISIONS")
            avoid_these = ["PARAMETERS", "VARIABLES", "VARIABLES_TO_PLOT", "VAR_MAP", "VARIABLE_PIVOT", "verbose",
                           "observations", "properties", "_observations", "_properties", "_population", "MODEL_NAME"]
            check_hard(len(params + avoid_these) == len(set(params + avoid_these)),
                       "INPUT DICTIONARIES NO NAME COLLISIONS")
            # # consistency with expected attributes
            # check(set(params) <= set(self.PARAMETERS + self.VARIABLES, ["time"]), "COMPLETENESS OF MODEL&IVP PARAMETERIZATION")

            # todo add check on N>= sum of other assigned variables in input

    def getName(self):
        '''
        Returns the class name.
        '''

        return self.__class__.__name__

    def model(self, x, t):
        '''
        Defines the model structure.
        '''

        pass

    @property
    def _x0(self):
        '''
        Returns the array-form initial state for the IVP.
        Entries are taken from self._ivp, with the except of VARIABLE_PIVOT variable, which is derived.
        '''
        # todo if incomplete default to observations

        _x0 = np.zeros(len(self.VARIABLES))
        pivot_index = self.VARIABLES.index(self.VARIABLE_PIVOT)

        # vars_non_pivot = self.VARIABLES[:]
        # vars_non_pivot.remove(self.VARIABLE_PIVOT)
        # entry_pivot = self._population - sum([self._ivp[key] for key in vars_non_pivot])
        #
        # _x0= [self._ivp[key] for key in vars_non_pivot]
        # #     = self._ivp[key]

        _x0[pivot_index] = 0

        for key in self.VARIABLES:
            if key != self.VARIABLE_PIVOT:
                _x0[self.VARIABLES.index(key)] = self._ivp[key]
                _x0[pivot_index] -= self._ivp[key]
            else:
                _x0[pivot_index] += self._population

        if _x0[pivot_index] < 0:
            raise Exception("Unfeasible choice of initial conditions for the IVP.")

        return _x0

    @property
    def _t(self):
        '''
        Returns the array of sampling times for the IVP. If time is not specified in _ivp, it is collected from
        _observations.
        '''
        # todo check if _observations has those attributes
        # todo improve _ivp structure

        if self._time is not None: #in locals(): #self._ivp.keys():
            _t = self._time #self._ivp["time"]
            self._t_dfindex = _t
            # self._indexes_for_predictions = _t
            # self.time = _t
        elif isinstance(self._observations, pd.DataFrame) and "time" in self._observations:
            _t = self._observations["time"]
            self._t_dfindex = self._observations.index
            # self._indexes_for_predictions = self._observations.index
            # self.time = self._observations["time"]
        else:
            # _t = [0]
            raise Exception("Unfeasible choice of time indexes for the IVP.")

        return _t

    @property
    def _state(self):
        '''
        Returns the array-form system's state found as solution of the ODE-IVP.
        '''

        return odeint(self.model, self._x0, self._t)

    @property
    def predictions(self):
        '''
        Returns the named-variables-form of the system's state remapped from _state.
        '''

        # map system's states to named variables
        predictions = pd.DataFrame(data=self._state, index=self._t_dfindex, columns=self.VARIABLES)
        # append differences, fill na with zeros in order to allow use in optimization (same as observations)
        for variable in self.VARIABLES:
            predictions["delta_" + variable] = predictions[variable].diff().fillna(0)
            predictions["delta_" + variable].values[0]=predictions["delta_" + variable].values[1]

        import itertools

        # get list of column names
        list_fields = list(predictions.columns.values)
        # get order two combinations of column names
        comb_list = list(itertools.combinations(list_fields, 2))

        for comb in comb_list:
            # for each combination, append to dataframe a combination column obtained as sum of the respective column names

            comb_name = '+'.join(comb)
            comb_entries = predictions[list(comb)].sum(axis=1)

            predictions[comb_name] = comb_entries

        return predictions

    @property
    def predictions_normalized(self):
        '''
        Returns normalized predictions
        '''
        return self.predictions.div(self._population)

    @property
    def observations(self):
        '''
        This property remaps observations_raw (raw epidemiological data)
        to the named state variables of the model (e.g. "S", "I", "R"), based on mapping VAR_MAP
        '''

        try:
            # map observations to named variables, using sums for combinations of original dataset columns
            observations = pd.DataFrame(index=self._observations.index)  # TODO check here index=self.t_dates)
            for key in self.VAR_MAP.keys():
                observations[key] = self._observations[self.VAR_MAP[key]].sum(axis=1)
            # append differences, fill na with zeros in order to allow use in optimization (same as predictions)
            for variable in self.VAR_MAP.keys():
                observations["delta_" + variable] = observations[variable].diff().fillna(0)
                observations["delta_" + variable].values[0]=observations["delta_" + variable].values[1]
        except:
            observations = []

        return observations

    @property
    def observations_normalized(self):
        '''
        Returns normalized predictions
        '''
        return self.observations.div(self._population)

    @property
    def VAR_OPTIM_ATTRIBUTES(self):
        '''
        Returns reference init value and ranges for the model parameters involved in parameter fitting.
        '''
        # todo this has to be generalized and extended, as well as better placed

        VAR_OPTIM_ATTRIBUTES = {}

        # attributes for model parameters
        for key in self.PARAMETERS:
            if key in self.VAR_OPTIM:
                if key in self.LIMITS:
                    VAR_OPTIM_ATTRIBUTES[key] = {"value": self.LIMITS[key][0], "min": self.LIMITS[key][1],
                                                 "max": self.LIMITS[key][2]}
                else:
                    VAR_OPTIM_ATTRIBUTES[key] = {"value": 0.25, "min": 0, "max": 55.5}

        # attributes for initial conditions
        for key in self._ivp:
            try:
                VAR_OPTIM_ATTRIBUTES["_ivpVARIABLE" + key] = {"value": max(self.observations[key]) / 2, "min": 0,
                                                              "max": max(1, max(self.observations[key]))}
                # todo manage issue as some regions may have zeros for long, avoid min-max overlap
            except:
                VAR_OPTIM_ATTRIBUTES["_ivpVARIABLE" + key] = {"value": 0, "min": 0, "max": 10000000}
                # todo derive from rest of ranges

        return VAR_OPTIM_ATTRIBUTES

    def set_limits(self, limit_entry, input_values):
        '''
        Sets self.LIMITS entries based on limit_entry (vqlue, min, max) and input_values.
        Useful to iterate over initializations of the fitting problem.
        :param limit_entry:
        :param input_values:
        :return:
        '''
        if limit_entry == "value":
            index = 0
        elif limit_entry == "min":
            index = 1
        elif limit_entry == "max":
            index = 2

        for key in input_values.keys():
            try:
                self.LIMITS[key][index] = input_values[key]
            except:
                pass

    @property
    def model_metrics(self):
        '''
        Derives key epidemiological metrics associated to the MODEL.
        TODO do the same for the observations

        From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1732833/pdf/v058p00538.pdf,
        - PREVALENCE
            is the proportion of people in a population who have some attribute or condition
            at a given point in time or during a specified time period
        - INCIDENCE (INCIDENCE RATE)
            is the number of new events (for example, new cases of a disease)
            in a defined population, occurring within a specified period of time
        - INCIDENCE PROPORTION (CUMULATIVE INCIDENCE)
            is the proportion of people who develop a condition within a fixed time period

        From https://www.cdc.gov/csels/dsepd/ss1978/lesson3/section2.html [MORBIDITY FREQUENCY MEASURES, see also other measures therein],
        - INCIDENCE PROPORTION (RISK)
            (number of new cases of disease or injury during specified period)/(size of population at start of period)
        - INCIDENCE RATE
            (number of new cases of disease or injury during specified period)/(time each person was observed, totaled for all persons)
        - PREVALENCE OF DISEASE
            [(all new and pre-existing cases during a given time period)/(population during the same time period)]*10^n
            where the value of 10 n is usually 1 or 100 for common attributes
        - PREVALENCE OF AN ATTRIBUTE
            [(persons having a particular attribute during a given time period)/(population during the same time period)]*10^n
            where the value of 10 n is usually 1 or 100 for common attributes

        Further references:
        https://courses.lumenlearning.com/microbiology/chapter/the-language-of-epidemiologists/
        https://en.wikipedia.org/wiki/Incidence_(epidemiology)
        https://en.wikipedia.org/wiki/Cumulative_incidence

        On the use of metrics for model fitting, see also
        King, A. A., Domenech de Cellès, M., Magpantay, F. M., & Rohani, P. (2015).
        Avoidable errors in the modelling of outbreaks of emerging pathogens, with special reference to Ebola.
        Proceedings of the Royal Society B: Biological Sciences, 282(1806), 20150347.
        Key observations from there:
            - "Recognizing that quantification of uncertainty is prerequisite
            to reliable forecasting, we computed parameter estimate
            confidence intervals and investigated their accuracy."
            - "Using the raw incidence data, one recovers the true observation variability."
            - "When a deterministic model is fit to cumulative incidence
            data, the net result is a potentially quite over-optimistic estimate
            of precision, for three reasons. First, failure to account
            for the non-independence of successive measurement errors
            leads to an underestimate of parameter uncertainty (figure
            1c). Second, as seen in figure 1b, the variance of measurement
            noise will be substantially underestimated. Finally, because
            the model ignores environmental and demographic stochasticity—
            treating the unfolding outbreak as a deterministic
            process—forecast uncertainty will grow unrealistically slowly
            with the forecast horizon.We elaborate on the last point in §4."

        todo consider also computing compartmental model fractions where relevant, see e.g.
        https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model
        todo add R0, Rt et cetera
        '''

        # todo implement
        model_metrics = {}
        model_metrics["prevalence"] = 0
        model_metrics["incidence_rate"] = 0
        model_metrics["cumulative_incidence"] = 0

        return model_metrics

    def plot_intime(self,plot_predictions=True, plot_observations=True, SinglePlot=False):
        '''
        Plots predictions and observations in time
        '''

        return plot_intime_comparison(cou=self.reference_area,df_reference=self.predictions, df_comparison=self.observations,
                                      fields_to_plot=self.VARIABLES_TO_PLOT,
                                      labelsuffix_df_reference="predictions", labelsuffix_df_comparison="observations")

    def get_optimization_variables(self):
        '''
        Returns a dict describing the model variables to exploit in optimization problems
        '''

        return {key: getattr(self, key) for key in self.VAR_OPTIM}
