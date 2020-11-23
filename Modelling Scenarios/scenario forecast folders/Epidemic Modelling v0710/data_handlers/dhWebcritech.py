from data_handlers.dhBase import dhBase
from data_handlers.dhDemography import dhDemography
import pandas as pd
from CONF import FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS, FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES,FILE_DEMOGRAPHY_B


def get_filename_Webcritech(area_type):
    if area_type == "Region":
        file_epidemiology = FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS
    elif area_type == "Country":
        file_epidemiology = FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES

    return file_epidemiology


def get_geofilter_Webcritech(area_type):
    '''
    Returns the fieldname useful to filter Webcritech information on a geographic basis (e.g. by country or by region).
    '''
    if area_type == "Country":
        geofilter = "CountryName"
    elif area_type == "Region":
        geofilter = "Region"
    else:
        raise Exception("Unsupported area_type.")

    return geofilter


class dhWebcritech(dhBase):
    '''
    Handles the Webcritech epidemiological dataset and jointly retrieves demographic information
    '''

    _ENCODING = 'utf-8'
    _FIELDS_TO_PLOT = ["CumulativeDeceased", "CumulativePositive", "CumulativeRecovered", "CurrentlyPositive",
                       "Hospitalized", "IntensiveCare"]

    def __init__(self, area_type, filters={}):

        super().__init__(filename=get_filename_Webcritech(area_type), filters=filters)

        # todo check => we fill nan entries with zeros
        # self.data = self.data.fillna(0)
        #
        # # self.data = self.data.interpolate(method='linear', axis=0)

        # convert Date to datetime => extract time indexes from Date => make Date the dataframe index
        self.data['Date'] = pd.to_datetime(self.data['Date'], format='%Y-%m-%d')
        # self.data['time'] = [(x - self.data.Date.iloc[0]).days for x in self.data.Date] # todo move
        self.data.set_index('Date', inplace=True)

        try:
            # reindex in order to include all dates, including gaps
            all_days = pd.date_range(self.data.index.min(), self.data.index.max(), freq='D')
            self.data = self.data.reindex(all_days)

            # interpolate missing values
            # self.data = self.data.fillna(0)
            self.data = self.data.interpolate(method='linear', axis=0)

        except:
            # reindexing may be impossible if filters are not enough to have uniqueness; in this case we skip
            pass

        # introduce an integer time index, useful for model integration todo consider moving this to ModelBase
        self.data['time'] = [(x - self.data.index[0]).days for x in self.data.index]

        # store filters applied in data processing
        self.filters_epidemiology = filters

        # if relevant filters are present, introduce them for population evaluation; give pre-eminence to Region (lower scale)
        geofilter = None
        if "CountryName" in self.filters_epidemiology.keys():
            geofilter = "CountryName"
        if "Region" in self.filters_epidemiology.keys():
            geofilter = "Region"

        if geofilter is not None:
            self.reference_area = self.filters_epidemiology[geofilter]
            self.filters_demography = {'TIME': 2019, 'SEX_LABEL': 'Total', 'AGE_LABEL': 'Total',
                                       'GEO_LABEL': self.reference_area}

        else:
            self.reference_area = ""
            self.filters_demography = {'TIME': 2019, 'SEX_LABEL': 'Total', 'AGE_LABEL': 'Total'}

        self.population = self.get_population()
        if self.population is not None:
            # derive the susceptible population, here assuming that it equals total population - cumulative positive
            self.data['Susceptible'] = self.population - self.data['CumulativePositive']

    def get_population1(self):

      	import json
      	ff=open(FILE_DEMOGRAPHY_B,'r',encoding='UTF-8')
      	data=json.load(ff)
      	ff.close()

      	region=self.reference_area
      	country=self.data.CountryName.values[0]
      	for feat in data['features']:
      	    pro=feat['properties']
      	    cou=pro['Cnt_Name']
      	    reg=pro['Region']
      	    if reg=='': reg=cou
      	    ghsl_sum=pro['ghsl_sum']
      	    if reg==region and cou==country:
      	         print('ALT_DEMO,'+cou+','+reg+','+format(ghsl_sum))
      	         return ghsl_sum
      	return None

    def get_population(self):
        '''
        Returns population size
        '''
        
        demography = dhDemography(filters=self.filters_demography)
        if demography.population is not None:
            return demography.population
        else:
            return self.get_population1()

