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


class SIR(ModelBase):
    PARAMETERS = ["beta", "gamma"]
    LIMITS = {"beta": [0.1, 0, 0.5], "gamma": [0.1, 0, 0.5]}  # ref, min, max for each parameter
    VARIABLES = ["S", "I", "R"]
    VARIABLES_TO_PLOT = ["S", "I", "R", "delta_S", "delta_R"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeDeceased", "CumulativeRecovered"]}
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, I, R = x
        N = S + I + R
        dxdt = [
            -self.beta * I * S / N,
            self.beta * I * S / N - self.gamma * I,
            self.gamma * I
        ]
        return dxdt
