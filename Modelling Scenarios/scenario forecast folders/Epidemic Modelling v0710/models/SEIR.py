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


class SEIR(ModelBase):
    PARAMETERS = ["beta", "gamma", "mu"]
    LIMITS = {"beta": [0.1, 0, 0.5], "gamma": [0.1, 0, 0.5], "mu": [0.1, 0, 0.5]}  # ref, min, max for each parameter
    VARIABLES = ["S", "E", "I", "R"]
    VARIABLES_TO_PLOT = ["E", "I", "R"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeDeceased", "CumulativeRecovered"]}
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, E, I, R = x
        N = S + E + I + R
        a = 7  # average incubation period
        dxdt = [
            self.mu * N - self.mu * S - self.beta * I * S / N,
            self.beta * I * S / N - (self.mu + a) * E,
            a * E - (self.gamma + self.mu) * I,
            self.gamma * I - self.mu * R
        ]
        return dxdt
