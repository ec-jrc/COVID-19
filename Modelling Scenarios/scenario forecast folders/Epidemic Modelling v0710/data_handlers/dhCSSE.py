from data_handlers.dhBase import dhBase
import pandas as pd
from CONF import FILE_EPIDEMIOLOGY_CSSE_RECOVERED


class dhCSSE(dhBase):
    '''
    Handles the CSSE epidemiological dataset and jointly retrieves demographic information
    '''

    _ENCODING = 'utf-8'
    _FIELDS_TO_PLOT = ["CumulativeDeceased", "CumulativePositive", "CumulativeRecovered", "CurrentlyPositive",
                       "Hospitalized", "IntensiveCare"]

    def __init__(self, file_epidemiology, filters={}):
        super().__init__(filename=file_epidemiology, filters=filters)

        aa = pd.read_csv(file_epidemiology, encoding=self._ENCODING).pivot(columns=0, values=1)

    ############
    # WIP


if __name__ == "__main__":
    file_epidemiology = FILE_EPIDEMIOLOGY_CSSE_RECOVERED
    aa = pd.read_csv(file_epidemiology, encoding='UTF-8', index_col=0, header=None).T

    csse = dhCSSE(file_epidemiology=file_epidemiology, file_demography="")
