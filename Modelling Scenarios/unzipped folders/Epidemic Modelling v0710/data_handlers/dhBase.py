import numpy as np
import pandas as pd
from util.plotting import plot_intime_comparison
# from Levenshtein import distance
from util.df_matcher import df_matcher

class dhBase:
    """
    Basic class for dataset handling and cleaning.
    """
    _ENCODING = 'utf-8'
    _FIELDS_TO_PLOT = []

    def __init__(self, filename=None, filters={}):

        # _data collects the original dataset
        self._data = pd.read_csv(filename, encoding=self._ENCODING)  # .compute()
        # data is the working copy, based e.g. on filtering or further manipulations
        self.data = self._data.copy()

        self.data = df_matcher(self.data, filters)

        # # apply filters to data
        # for key in filters.keys():
        #     if isinstance(filters[key], str):
        #         # if string, try perfect matching => Levenshtein-bounded matching => substring matching
        #         strict_matches = self.data[self.data[key] == filters[key]]
        #         if strict_matches.shape[0] > 0:
        #             self.data = strict_matches
        #         else:
        #             quasistrict_matches = self.data[[distance(x,filters[key]) <= 1 for x in self.data[key]]]
        #             if quasistrict_matches.shape[0] > 0:
        #                 self.data = quasistrict_matches
        #             else:
        #                 self.data = self.data[self.data[key].str.contains(filters[key], na=False)]
        #     else:
        #         # if not string, search for exact matches
        #         self.data= self.data[self.data[key] == filters[key]]
        #         # self.select(key=key, value=filters[key])

        self.reference_area = ""

    def select(self, key, value):
        '''
        Filters dataframe self.data by column column_name having entry value.
        '''

        relevant_rows = self.data[key] == value
        self.data = self.data[relevant_rows]

    @property
    def data_smoothed(self, filtering_method="Default"):
        '''
        Returns the filtered time-series (excluding "lat" and "lon" fields)
        '''

        if filtering_method is None:
            data_smoothed = []
        else:
            # deep copy self.data
            data_smoothed = self.data.copy()
            # find all numeric fields and fit them
            numeric_fields = data_smoothed.select_dtypes(include=np.number).columns.tolist()
            for entry in numeric_fields:
                if entry not in ["lat", "lon"]:
                    # todo here we assume that time is a column in data_smoothed; ensure that this is true not only for Webcritech
                    data_smoothed[entry] = self.smooth(data_smoothed.time.to_numpy(), data_smoothed[entry].to_numpy())

        return data_smoothed

    def smooth(self, x, y):
        '''
        Returns an smoothed array.
        '''
        # todo generalize fitting_method and parameters
        # todo handle data gaps

        # https://math.stackexchange.com/questions/677699/how-can-i-find-a-non-negative-interpolation-function
        # output_array = savgol_filter(y, 5, 3)
        # output_array = savgol_filter(np.square(y), 5, 3)
        # output_array = np.sqrt(output_array)
        output_array = []
        # todo

        return output_array

    def plot_intime(self):
        '''
        Plots original and filtered datasets
        '''

        # todo reintroduce when data smooting is done
        # return plot_intime_comparison(cou=self.reference_area,df_reference=self.data_smoothed, df_comparison=self.data,
        #                               fields_to_plot=self._FIELDS_TO_PLOT,
        #                               labelsuffix_df_reference="filtered", labelsuffix_df_comparison="raw")

        return plot_intime_comparison(cou=self.reference_area, df_reference=self.data, df_comparison=self.data,
                                      fields_to_plot=self._FIELDS_TO_PLOT,
                                      labelsuffix_df_reference="filtered", labelsuffix_df_comparison="raw")
