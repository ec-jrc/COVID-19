'''
MODIFIED SIRD MODEL, FORMULATION PROPOSED BY DIEGO CACCAVO, SEE REFERENCE

PARAMETERS:
beta0; beta1; tau_beta; delta0; delta1; tau_delta; gamma0

VARIABLES:
S (susceptible); I (infectious); R: (recovered); D (deceased)
REFERENCES:
@article{caccavo2020chinese,
  title={Chinese and Italian COVID-19 outbreaks can be correctly described by a modified SIRD model},
  author={Caccavo, Diego},
  journal={medRxiv},
  year={2020},
  publisher={Cold Spring Harbor Laboratory Press}
}
'''

from models.ModelBase import ModelBase
import numpy as np


class SIRD_modified_Caccavo(ModelBase):
    PARAMETERS = ["beta0", "beta1", "tau_beta", "gamma0", "gamma1", "tau_gamma", "delta0", "delta1", "tau_delta"]
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
        beta = self.beta0 * np.exp(-t/self.tau_beta) + self.beta1
        gamma = self.gamma0 + self.gamma1/(1 + np.exp(-t + self.tau_gamma))
        delta = self.delta0 * np.exp(-t/self.tau_delta) + self.delta1
        dxdt = [
            -beta * I * S / N,
            beta * I * S / N - gamma * I - delta * I,
            gamma * I,
            delta * I
        ]
        return dxdt
