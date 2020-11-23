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
    PARAMETERS = ["r0",  "Trecov", "Tdeath", "mu"]

    LIMITS     = {"r0":[0.7,0.5,3.5],"Trecov":[7,5.5,8],"Tdeath":[14,8,30], "mu":[0.02,0.0,0.15]}  # ref, min, max for each parameter
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
        Teff=(1-self.mu)*1/self.Trecov+self.mu*1/self.Tdeath
        dSdt = - self.r0 /(Teff) * S / N *I        
        dIdt =   self.r0/(Teff)  * S / N *I - 1/self.Trecov * I *(1-self.mu) - 1/self.Tdeath * I *self.mu
        dRdt =   1/self.Trecov * I *(1-self.mu)
        dFdt =   1/self.Tdeath * I *(self.mu)
        

        dxdt =  [dSdt,dIdt,dRdt, dFdt]
        return dxdt
    
# -*- coding: utf-8 -*-



