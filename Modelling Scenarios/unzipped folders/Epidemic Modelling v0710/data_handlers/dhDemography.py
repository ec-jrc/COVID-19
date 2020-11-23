import numpy as np
from data_handlers.dhBase import dhBase
from util.df_matcher import df_matcher_unique

from CONF import FILE_DEMOGRAPHY

# todo check we skip "Repatriierte" and the like

class dhDemography(dhBase):
    '''
    Handles the EUROSTAT-sourced demographic dataset demo_r_pjangrp3.
    '''

    _ENCODING = 'latin'

    def __init__(self, filters={}):

        super().__init__(filename=FILE_DEMOGRAPHY)

        for key in filters.keys():
            if key == "GEO_LABEL" and filters[key] == "Czech Republic":
                filters[key] = "Czechia"

        self.data = df_matcher_unique(self.data,filters)

    @property
    def population(self):
        if self.data.Value.size > 0:
            # (as a last resort against multiple instances,) take the maximum
            return max(self.data.Value.str.replace(',', '').replace(':', '0').to_numpy().astype(np.int))
        else:
            return None
