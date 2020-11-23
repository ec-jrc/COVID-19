from CONF import VALID_FILENAME_CHARS

def filename_cleaner(filename):
    '''
    Removes from a given filename strings non-valid characters, i.e. those not included in VALID_FILENAME_CHARS.
    Supports dicts.
    '''
    if isinstance(filename,dict):
        for key in filename.keys():
            filename[key] = ''.join(c for c in filename[key] if c in VALID_FILENAME_CHARS)
    else:
        filename = ''.join(c for c in filename if c in VALID_FILENAME_CHARS)

    return filename