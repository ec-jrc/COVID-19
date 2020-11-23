'''
SIRD1 MODEL

PARAMETERS:
r0 (reprod number); Tinf (infection time); Trecov (time of recovery); mu (fatality rate)
VARIABLES:
S (susceptible);  I (infectious); R: (recovered) ; D: (deceased)
REFERENCES:
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
'''

from models.ModelBase import ModelBase


class SIRD1(ModelBase):
    PARAMETERS = ["r0", "Tinf",  "Trecov", "Tdeath", "mu"]
    LIMITS     = {"r0":[0.7,0.,3.5],"Tinf":[7,5.6,8.5],"Trecov":[14,5,30],"Tdeath":[14,5,30], "mu":[0.02,0.0,0.15]}  # ref, min, max for each parameter
    # LIMITS     = {"r0":[0.7,0.,3.5],"Tinf":[7,6.9,7.1],"Trecov":[7,5.6,8.5],"Tdeath":[7,5.6,8.5], "mu":[0.02,0.0,0.15]}  # ref, min, max for each parameter
    VARIABLES = ["S",  "I", "R" , "D"]
    VARIABLES_TO_PLOT = ["I", "R", "D"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "R": ["CumulativeRecovered"],
               "D": ["CumulativeDeceased"]}
               
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S,  I, R, D = x
        N = S + I + R  + D
        dSdt = - self.r0 /(self.Tinf) * S / N *I
        dIdt =   self.r0/(self.Tinf)  * S / N *I - 1/self.Trecov * I *(1-self.mu) - 1/self.Tdeath * I *self.mu
        dRdt =   1/self.Trecov * I *(1-self.mu)
        dFdt =   1/self.Tdeath * I *(self.mu)
        

        dxdt =  [dSdt,dIdt,dRdt, dFdt]
        return dxdt
    
# -*- coding: utf-8 -*-



