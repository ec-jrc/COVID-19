import numpy as np

class optiBase:
    '''
    Implements a baseline class for model optimization.
    '''

    def __init__(self, model, variables_to_fit, fitting_weights, verbose=True, optimization_method=None):  # ,fitted_variables):
        self.model = model
        self.variables_to_fit = variables_to_fit
        self.fitting_weights = fitting_weights
        self.verbose = verbose

        self.optimize(optimization_method = optimization_method)

    @property
    def residuals(self):
        '''
        Returns residuals between predictions and observations.
        '''

        residuals = []

        for key in self.variables_to_fit:
            residuals = np.append(residuals,self.fitting_weights[key]*(self.model.predictions[key].to_numpy() - self.model.observations[key].to_numpy()))

        return residuals

    def optimize(self):
        '''
        Solves the optimization problem.
        '''
        pass
