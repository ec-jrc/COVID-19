import pandas as pd
from data_handlers.dhWebcritech import dhWebcritech, get_geofilter_Webcritech
from util.df_matcher import df_matcher_unique


def areanameGetter(REFERENCE_AREA, area_type="Country", country_filter=""):
    '''
    Function useful to get names of areas to consider in the analysis. This includes
    - getting the full name list of areas covered by data (REFERENCE_AREA);
    - for each element of list of areas REFERENCE_AREA, getting the name of the upped administrative division.
    It works on Webcritech datasets (for populating REFERENCE_AREA) and pjangrp3 datasets (for administrative super-divisions).

    Details:
    - if REFERENCE_AREA is None, skip everything
    - if REFERENCE_AREA=[], REFERENCE_AREA is populated with all areas from the FILE_EPIDEMIOLOGY dataset
    - if FILE_EPIDEMIOLOGY refers to regions, REFERENCE_AREA is potentially limited to those regions
        belonging to country COUNTRY_FILTER;
    Returns the (enriched) list REFERENCE_AREA jointly with (same-length) list SUPER_REFERENCE_AREA, which describes
    (where relevant) the higher level area comprising each area in REFERENCE_AREA. SUPER_REFERENCE_AREA is useful, e.g.,
    in order to produce geoplots where regions are plotted over countries, for easier identification.
    :param REFERENCE_AREA:     list of reference areas to be considered
    :param COUNTRY_FILTER:      string, relevant for region-level analysis, useful to filter regions belonging to the corresponding country
    :param FILE_EPIDEMIOLOGY:   filename of the epidemiology dataset
    :return:
    REFERENCE_AREA (list of reference areas, enriched if starting from [])
    SUPER_REFERENCE_AREA (list of super-reference area for each entry in REFERENCE_AREA)
    '''

    # if REFERENCE_AREA is None, skip everything
    if REFERENCE_AREA is None:
        SUPER_REFERENCE_AREA = {}
        return REFERENCE_AREA, SUPER_REFERENCE_AREA

    # REFERENCE_AREA
    geofilter = get_geofilter_Webcritech(area_type)
    #REFERENCE_COUNTRY=[]
    if REFERENCE_AREA == []:
        df_full = dhWebcritech(area_type=area_type).data
        df = df_full.drop_duplicates(geofilter).sort_values(geofilter)
        if area_type == "Region" and country_filter != "":
            # apply the country filter
            df = df[df["CountryName"] == country_filter]
        REFERENCE_AREA = df[geofilter].tolist()
        #REFERENCE_COUNTRY=df["CountryName"].tolist()
    else:
        df = pd.DataFrame()
        for ra in REFERENCE_AREA:
            df_full = dhWebcritech(area_type=area_type, filters={geofilter: ra}).data
            df = df.append(df_full.drop_duplicates(geofilter).sort_values(geofilter))

    #  SUPER_REFERENCE_AREA
    SUPER_REFERENCE_AREA = {area: "" for area in REFERENCE_AREA}
    if area_type == "Region":
        for area in REFERENCE_AREA:
            SUPER_REFERENCE_AREA[area] = df_matcher_unique(df,filters={geofilter: area})["CountryName"].tolist()[0]

    return REFERENCE_AREA, SUPER_REFERENCE_AREA
