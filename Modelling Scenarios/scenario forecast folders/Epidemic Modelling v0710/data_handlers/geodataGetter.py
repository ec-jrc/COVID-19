import time
from geopy.geocoders import Nominatim
import geopandas
import matplotlib.pyplot as plt
import os

from CONF import DIR_GEOSPATIAL


class geodataGetter():
    '''
    This class takes care of retrieving geospatial information useful to the representation of each analysis area
    '''

    def __init__(self, location="", super_location={}, force_recreate=False):
        self.force_recreate = force_recreate

        if location:
            self.df_location = self.get_geo(location)
            if super_location:
                self.df_super_location = self.get_geo(super_location[location])

    def get_geo(self, location):

        try:
            # first attempt to retrieve geospatial information locally, secondly resort to a web api

            filename = DIR_GEOSPATIAL + "/" + location + ".geojson"

            # todo location non-empty string check addition
            if not os.path.isfile(filename) or self.force_recreate:
                # avoid overloading the Nominatim API
                time.sleep(1.5)

                # import geospatial information from Nominatim API
                geolocator = Nominatim(user_agent="my-application")
                location_retrieved = geolocator.geocode(location, geometry="geojson")

                # save temporarily to file the extracted geojson
                geojson_area = str(location_retrieved.raw["geojson"])
                geojson_area = geojson_area.replace("'", '"')

                os.makedirs(DIR_GEOSPATIAL, exist_ok=True)

                # TEMP_FILENAME = "_temp.geojson"
                with open(filename, "w") as text_file:
                    text_file.write(geojson_area)

            # import geospatial information using geopandas utilities
            df_place = geopandas.read_file(filename)

            # os.remove(TEMP_FILENAME)

        except:

            df_place = '{"type": "FeatureCollection", "features": []}'

        return df_place

    def plot(self):

        fig, ax = plt.subplots(1, 1)

        if self.df_super_location is not None and not isinstance(self.df_super_location, str):
            self.df_super_location.plot(ax=ax)

        if self.df_location is not None and not isinstance(self.df_location, str):
            self.df_location.plot(ax=ax, cmap="seismic")

        return fig
