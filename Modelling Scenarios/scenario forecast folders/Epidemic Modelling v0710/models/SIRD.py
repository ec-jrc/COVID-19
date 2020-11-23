'''
SIRD MODEL

PARAMETERS:
beta (infection rate); gamma (recovery rate); mu (mortality rate)
VARIABLES:
S (susceptible); I (infectious); R: (recovered); D (deceased)
REFERENCES:
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
'''

from models.ModelBase import ModelBase


class SIRD(ModelBase):
    PARAMETERS = ["beta", "gamma", "mu"]
    LIMITS     = {"beta":[0.5,0.,1],"gamma":[0.5,0.,1],"mu":[0.5,0.,1]}  # ref, min, max for each parameter
    VARIABLES = ["S", "I", "R", "D"]
    VARIABLES_TO_PLOT = ["I", "R", "D"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeRecovered"],
               "D": ["CumulativeDeceased"]}
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, I, R, D = x
        N = S + I + R + D
        dxdt = [
            -self.beta * I * S / N,
            self.beta * I * S / N - self.gamma * I - self.mu * I,
            self.gamma * I,
            self.mu * I
        ]
        return dxdt
