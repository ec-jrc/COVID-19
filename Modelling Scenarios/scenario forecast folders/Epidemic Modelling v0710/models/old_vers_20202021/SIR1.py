'''
SIR1 MODEL

PARAMETERS:
r0 (reprod number); Tinf (infection time); Trecov (time of recovery)
VARIABLES:
S (susceptible);  I (infectious); R: (recovered)
REFERENCES:
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
'''

from models.ModelBase import ModelBase


class SIR1(ModelBase):
    PARAMETERS = ["r0", "Tinf", "Trecov"]
    # LIMITS = {"r0": [2.5, 1.1, 3.5], "Tinf": [7, 5.6, 8.5], "Trecov": [14, 5, 30]}  # ref, min, max for each parameter
    LIMITS = {"r0": [2.5, 0, 10], "Tinf": [7, 1, 15], "Trecov": [7, 1, 50]}  # ref, min, max for each parameter
    # LIMITS     = {"r0":[2.5,1.1,3.5],"Tinf":[7,6.9,7.1],"Trecov":[14,5,30]}  # ref, min, max for each parameter
    VARIABLES = ["S", "I", "R"]
    VARIABLES_TO_PLOT = ["I", "R"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeDeceased", "CumulativeRecovered"]}

    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, I, R = x
        N = S + I + R
        dSdt = - self.r0 / (self.Tinf) * S / N * I
        dIdt = self.r0 / (self.Tinf) * S / N * I - 1 / self.Trecov * I
        dRdt = 1 / self.Trecov * I

        dxdt = [dSdt, dIdt, dRdt]
        return dxdt

# -*- coding: utf-8 -*-
