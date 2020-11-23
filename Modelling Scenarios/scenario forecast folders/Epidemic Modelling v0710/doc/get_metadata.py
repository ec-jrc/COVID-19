from util.filename_cleaner import filename_cleaner

def get_study_filestring(study_type, model_type, analysis_horizon, reference_area):
    return filename_cleaner(study_type + "_" + str(model_type) + "_" + str(analysis_horizon)) # .replace("'", "").replace("/", "")) # + "_" + str(reference_area)


def get_study_metadata(study_type, model_type, analysis_horizon, reference_area, super_reference_area):
    metadata = {}

    metadata["study_type"] = study_type
    metadata["model_type"] = model_type
    metadata["analysis_horizon"] = analysis_horizon
    metadata["reference_area"] = reference_area
    metadata["super_reference_area"] = super_reference_area
    #metadata["reference_country"] = reference_country

    filename = get_study_filestring(study_type, model_type, analysis_horizon, reference_area)
    # define filenames for all the format we might use (actual use depends on the type of study)
    metadata["filenames"] = {}
    metadata["filenames"]["csv"] = filename + ".csv"
    metadata["filenames"]["json"] = filename + ".json"
    metadata["filenames"]["tex"] = filename + ".tex"
    metadata["filenames"]["xlsx"] = filename + ".xlsx"

    return metadata
