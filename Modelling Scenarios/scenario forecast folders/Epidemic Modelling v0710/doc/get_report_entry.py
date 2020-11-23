def get_report_entry(entry_metadata, analysis_outcomes, analysis_output):
    '''
    Generates entry with information relevant to the creation of reports
    '''


    report_entry = {}

    report_entry["metadata"] = entry_metadata
    report_entry["analysis_outcomes"] = analysis_outcomes
    report_entry["analysis_output"] = analysis_output

    return report_entry