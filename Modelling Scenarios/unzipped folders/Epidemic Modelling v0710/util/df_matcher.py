from Levenshtein import distance

def df_matcher(df,filters):

    # apply filters to data
    for key in filters.keys():
        if isinstance(filters[key], str):
            # if string, try perfect matching => Levenshtein-bounded matching => substring matching
            strict_matches = df[df[key] == filters[key]]
            if strict_matches.shape[0] > 0:
                df = strict_matches
            else:
                quasistrict_matches = df[[distance(x, filters[key]) <= 1 for x in df[key]]]
                if quasistrict_matches.shape[0] > 0:
                    df = quasistrict_matches
                else:
                    df = df[df[key].str.contains(filters[key], na=False)]
        else:
            # if not string, search for exact matches
            df = df[df[key] == filters[key]]
            # self.select(key=key, value=filters[key])

    return df

def df_matcher_unique(df,filters):

    df = df_matcher(df,filters)

    # if more than one matches remain after filtering, force strict name matching to get unique entry
    if df.shape[0] > 1:
        df = df.head(1)#df[df.iloc[0]]
        # todo evalutate using different logic

    return df