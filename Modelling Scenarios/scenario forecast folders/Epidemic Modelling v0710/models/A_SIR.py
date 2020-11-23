'''
A-SIR MODEL (SIR with asymptomatic infectives) [\cite{gaeta2020simple}, formula (30)]

PARAMETERS:
alpha, beta, eta, xi
VARIABLES:
S (susceptible); I (symptomatic); J (asymptomatic); R (symptomatic removed); U (asymptomatic removed)
REFERENCES:
@article{gaeta2020simple,
  title={A simple SIR model with a large set of asymptomatic infectives},
  author={Gaeta, Giuseppe},
  journal={arXiv preprint arXiv:2003.08720},
  year={2020}
}
'''

from models.ModelBase import ModelBase


class A_SIR(ModelBase):
    PARAMETERS = ["alpha", "beta", "eta", "xi"]
    VARIABLES = ["S", "I", "J", "R", "U"]
    VARIABLES_TO_PLOT = ["I", "J", "R", "U"]
    VAR_MAP = {"S": ["Susceptible"],
               "I": ["CurrentlyPositive"],
               "J": [],
               "R": ["CumulativeDeceased", "CumulativeRecovered"],
               "U": []}
    VARIABLE_PIVOT = "S"

    def model(self, x, t):
        S, I, J, R, U = x
        dxdt = [
            -self.alpha * S * (I + J),
            self.alpha * self.xi * S * (I + J) - self.beta * I,
            self.alpha * (1 - self.xi) * S * (I + J) - self.eta * J,
            self.beta * I,
            self.eta * J
        ]
        return dxdt
