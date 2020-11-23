'''
SIR MODEL

PARAMETERS:
beta (infection rate); gamma (recovery rate)
VARIABLES:
S (susceptible); I (infectious); R: (removed)
REFERENCES:
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
'''

from models.ModelBase import ModelBase


class SIR_betahat(ModelBase):
    PARAMETERS = ["gamma"]
    LIMITS = {"gamma": [0.1, 0, 0.5]}  # ref, min, max for each parameter
    VARIABLES = ["S", "I", "R"]
    VARIABLES_TO_PLOT = ["S", "I", "R", "delta_S", "delta_R"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeDeceased", "CumulativeRecovered"]}
    VARIABLE_PIVOT = "S"

    def __init__(self, verbose: bool = True, parameters_model: dict = {}, parameters_ivp: dict = {}, population=0,
                 observations=None, time=None, betahat=None):
        super().__init__(verbose=verbose,parameters_model=parameters_model,parameters_ivp=parameters_ivp,population=population,observations=observations,time=time)
        self.betahat = betahat

    def model(self, x, t):
        S, I, R = x
        N = S + I + R
        if I > 0:
            beta = self.betahat * ((I+R)/I)
        else:
            beta = 0

        dxdt = [
            -beta * I * S / N,
            beta * I * S / N - self.gamma * I,
            self.gamma * I
        ]
        return dxdt
