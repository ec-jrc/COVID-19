'''
Definition of parameters commonly used across the project.
'''
from os.path import dirname
import string
import datetime,os

PWD = dirname(__file__)

# general definitions
VALID_FILENAME_CHARS = "-_.() %s%s" % (string.ascii_letters, string.digits)

# base directories for I/O
DIR_DATA = PWD + "/_data"
DIR_GEOSPATIAL = PWD + "/_geospatial"

DIR_DATA_DAY=DIR_DATA+"/"+datetime.datetime.utcnow().strftime('%Y-%m-%d')
if not os.path.isdir(DIR_DATA_DAY+"/webcritech"):
    os.makedirs(DIR_DATA_DAY+"/webcritech")

DIR_REPORTING = PWD+ "/_report/"+datetime.datetime.utcnow().strftime('%Y-%m-%d')
if not os.path.isdir(DIR_REPORTING):
    os.makedirs(DIR_REPORTING)

# epidemiology
FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS = DIR_DATA_DAY + "/webcritech/jrc-covid-19-all-days-by-regions.csv"
FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES = DIR_DATA_DAY + "/webcritech/jrc-covid-19-all-days-by-country.csv"
FILE_EPIDEMIOLOGY_CSSE_CONFIRMED = DIR_DATA + "/csse/time_series_covid19_confirmed_global.csv"
FILE_EPIDEMIOLOGY_CSSE_DEATHS = DIR_DATA + "/csse/time_series_covid19_deaths_global.csv"
FILE_EPIDEMIOLOGY_CSSE_RECOVERED = DIR_DATA + "/csse/time_series_covid19_recovered_global.csv"

# demography
FILE_DEMOGRAPHY = DIR_DATA + '/demo_r_pjangrp3/demo_r_pjangrp3_1_Data.csv'
FILE_DEMOGRAPHY_B=DIR_DATA + '/region_comb.json'

COUNTRIES_SUBSET = "Austria, Bulgaria, Germany, Hungary, Italy, Malta, Portugal, Slovakia".split(", ")
#COUNTRIES_SUBSET = "Italy".split(", ")

eucmParticip = "Iceland, Montenegro, North Macedonia, Norway, Turkey, United Kingdom, Switzerland"
ecList="Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden"


COUNTRIES_ALL = (ecList+', '+eucmParticip).split(", ")



#COUNTRIES_SUBSET = ['Germany']
