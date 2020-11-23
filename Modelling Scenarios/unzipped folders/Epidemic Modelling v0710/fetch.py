#!/usr/bin/env python
# coding: utf-8

'''
Imports datasets to subfolder specified by DATA_DIR.
'''

import dload
import shutil

DATA_DIR = './_data'

# [0] Epidemiological data

# webcritech
shutil.rmtree(DATA_DIR + '/webcritech', ignore_errors=True)
EPIDEMIOLOGY_WEBCRITECH_COUNTRY = 'https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
EPIDEMIOLOGY_WEBCRITECH_REGION = 'https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
dload.save_multi([EPIDEMIOLOGY_WEBCRITECH_COUNTRY, EPIDEMIOLOGY_WEBCRITECH_REGION], DATA_DIR + '/webcritech')

# csse
shutil.rmtree(DATA_DIR + '/csse', ignore_errors=True)
EPIDEMIOLOGY_CSSE_CONFIRMED = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv'
EPIDEMIOLOGY_CSSE_DEATHS = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv'
EPIDEMIOLOGY_CSSE_RECOVERED = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv'
EPIDEMIOLOGY_CSSE_LOOKUP = 'https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv'

dload.save_multi(
    [EPIDEMIOLOGY_CSSE_CONFIRMED, EPIDEMIOLOGY_CSSE_DEATHS, EPIDEMIOLOGY_CSSE_RECOVERED, EPIDEMIOLOGY_CSSE_LOOKUP],
    DATA_DIR + '/csse')

# [1] Demographic data

# manually imported from https://ec.europa.eu/eurostat/web/products-datasets/product?code=demo_r_pjangrp3
