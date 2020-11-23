'''
SEIRD MODEL

PARAMETERS:
r0 (reprod number); Tinf (infection time); mu (mortality rate) ; Tinc (time of incubation)
VARIABLES:
S (susceptible); E (exposed); I (infectious); R: (recovered); D (deceased)
REFERENCES:
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
'''

from models.ModelBase import ModelBase


class SEIRD(ModelBase):
    PARAMETERS = ["r0", "Tinf", "mu", "Tinc", "Trecov"]
    VARIABLES = ["S", "E", "I", "R", "D"]
    VARIABLES_TO_PLOT = ["E", "I", "R", "D"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeRecovered"],
               "D": ["CumulativeDeceased"]}
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, E, I, R, D = x
        N = S + E + I + R + D
        dSdt = - self.r0 / self.Tinc * S / N * I
        dEdt = self.r0 / self.Tinc * S / N * I - 1 / self.Tinc * E
        dIdt = 1 / self.Tinf * E - ((1 - self.mu) * 1 / self.Trecov + self.mu * 1 / self.Trecov) * I
        dRdt = (1 - self.mu) * 1 / self.Trecov * I
        dDdt = (self.mu) * 1 / self.Trecov * I

        dxdt = [dSdt, dEdt, dIdt, dRdt, dDdt]
        return dxdt
