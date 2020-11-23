'''
TIME-DEPENDENT SIR MODEL

PARAMETERS:
beta (contact rate, time-dependent); gamma (inverse of average infectious period, time-dependent)
VARIABLES:
S (susceptible); I (infective); R: (removed)
REFERENCES:
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
'''

from models.ModelBase import ModelBase
import numpy as np


class SIR_td(ModelBase):
    PARAMETERS = ["beta", "gamma"]
    VARIABLES = ["S", "I", "R"]
    VARIABLES_TO_PLOT = ["I", "R"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeDeceased", "CumulativeRecovered"]}
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, I, R = x
        N = S + I + R
        # select the closes lower entry index for parameters
        entry = min(np.argmax(1 / (t - self.time)), len(self.beta) - 1)  # todo
        dxdt = [
            -self.beta[entry] * S * I / N,
            self.beta[entry] * S * I / N - self.gamma[entry] * I,
            self.gamma[entry] * I
        ]
        return dxdt
