import os
import json
import matplotlib.pyplot as plt
import pandas as pd
import pylatex as pl
import shutil

from data_handlers.geodataGetter import geodataGetter
from util.filename_cleaner import filename_cleaner

from CONF import DIR_REPORTING


def safe_plot_intime_save(file_prefix, file_suffix, object_to_plot,dir="./fig/"):
    '''
    safely tries to issue methods plot_intime()/plot() on object_to_plot and save result to file
    :param file_prefix: prefix to use in file naming
    :param file_suffix: suffix to use in file naming
    :param object_to_plot: object over which to issue the method
    :return: filename of the generated file, or "" if unsuccessful
    '''
    # todo improve figure naming
    try:
        FILENAME = dir + filename_cleaner((file_prefix + "_" + file_suffix + ".eps").replace(" ", "_").encode('ascii',
                                                                                                  'replace').decode())
    except:
        FILENAME = ""

    object_to_plot.savefig(FILENAME, format='eps')

    return FILENAME


def append_figure(mp, FILENAME):
    '''
    Appends figure from file FILENAME to Latex division represented by mp.
    '''
    if FILENAME != "":
        mp.append(pl.NoEscape(r'''
        \centering
        \includegraphics[width=\textwidth]{''' + FILENAME + '''}
        '''))


def append_json(mp, json_string):
    '''
    Appends well-formatted representation of a JSON string to Latex division represented by mp.
    '''
    # mp.append(pl.NoEscape("Study context"))
    json_string = "\\begin{lstlisting}[language=json,firstnumber=1]\n" + json_string + "\n\\end{lstlisting}"
    mp.append(pl.NoEscape(json_string))
    return mp


def print_block(section, entry, vspace):
    space = str(vspace) + r"\textwidth"
    with section.create(pl.MiniPage(width=space)) as mp:
        if isinstance(entry, dict):
            append_json(mp, json.dumps(entry, indent=4))
        elif isinstance(entry, str):
            if os.path.isfile(entry):  # todo add file extension check
                append_figure(mp, entry)
            else:
                pass  # todo implement


def get_value(d, target='csv'):
    '''
    Explores dict data structure and returns list of sub-entries containing key specified by target
    :param d:
    :param target:
    :return:
    '''
    # based on https://stackoverflow.com/questions/49662970/how-can-i-search-for-specific-keys-in-this-nested-dictionary-in-python
    val = filter(None,
                 [[b] if a == target else get_value(b, target) if isinstance(b, dict) else None for a, b in d.items()])
    return [i for b in val for i in b]


def stack_data(report, target):
    analysis_output_stacked = {}
    for key in report.keys():
        if key != "metadata":
            if report[key]:
                analysis_output_stacked[key] = get_value(report[key], target)[0]  # keep the first entry
    return analysis_output_stacked


def save_data(filenames, curr_dict):
    # save analysis_output entries to files
    for key in curr_dict.keys():
        if key == "json":  # todo check is dict, also we assume single output per filetype
            with open(filenames[key], "w") as f:
                json.dump(curr_dict[key], f)
        elif key == "csv":  # todo check also here
            df = pd.DataFrame.from_dict(curr_dict[key], orient='index')
            df.to_csv(filenames[key])
        elif key == "xlsx":  # here we assume a dict of df, and report each df on a separate file
            # save to excel files
            # filename_forecasts = DIR_REPORTING + "/forecasts_" + filename_suffix + ".xlsx"
            with pd.ExcelWriter(filenames[key]) as writer:
                for subkey in curr_dict[key].keys():
                    curr_dict[key][subkey].to_excel(writer, sheet_name=subkey)


def docgenerator(DIR_REPORTING, report, detach=False, makereport=False):
    print("")
    print("*" * 5)
    print("Generating documentation.")

    os.makedirs(DIR_REPORTING, exist_ok=True)
    os.chdir(DIR_REPORTING)

    # go to subfolder named based on study type
    try:
        shutil.rmtree(report["metadata"]["study_type"])
    except:
        pass
    os.makedirs(report["metadata"]["study_type"], exist_ok=True)
    os.chdir(report["metadata"]["study_type"])

    print("\tSaving datasets.")

    # go to subfolder named based on study type
    os.makedirs("data", exist_ok=True)
    os.chdir("data")

    # analysis_output => save data to files (approach changes based on data detaching option)
    if detach:
        for key in report.keys():
            if key != "metadata":
                analysis_output_stacked = {}
                analysis_output_stacked["xlsx"] = stack_data(report, target="xlsx")
                # todo same for other formats
                for subkey in analysis_output_stacked:
                    filenames = {}
                    analysis_output = {}
                    analysis_output[subkey] = analysis_output_stacked[subkey][key]
                    filenames[subkey] = filename_cleaner(report["metadata"]["filenames"][subkey].replace(".", "_" + key + "."))
                    #     # .replace(".", "_" + key + ".").replace(
                    #     # "'", "").replace("/", "").replace(":", "")
                    # filenames[subkey] = ''.join(c for c in filenames[subkey] if c in VALID_FILENAME_CHARS)
                save_data(filenames=filenames, curr_dict=analysis_output)

    else:
        analysis_output_stacked = {}
        analysis_output_stacked["csv"] = stack_data(report, target="csv")
        analysis_output_stacked["json"] = stack_data(report, target="json")
        # todo same for other formats
        filenames = filename_cleaner(report["metadata"]["filenames"])
        save_data(filenames=filenames, curr_dict=analysis_output_stacked)

    if not makereport:
        return

    # analysis_outcomes => generate report
    print("\tGenerating report.")

    # go to subfolder named based on study type
    os.chdir("..")
    os.makedirs("report", exist_ok=True)
    os.chdir("report")
    os.makedirs("./fig")
    os.makedirs("./fig_map")

    # generate documentation (Latex-based)
    # https://jeltef.github.io/PyLaTeX/current/examples/complex_report.html

    geometry_options = {
        "head": "40pt",
        "margin": "0.5in",
        "bottom": "0.6in",
        "includeheadfoot": True
    }
    doc = pl.Document(documentclass="book", geometry_options=geometry_options)
    doc.packages.append(pl.Package('graphics'))
    doc.packages.append(pl.Package('booktabs'))
    doc.packages.append(pl.Package('paracol'))

    # Utility to append headings
    doc.preamble.append(pl.NoEscape(r'''

        \usepackage{tikz}
        \usepackage{color}

        \definecolor{pink}{HTML}{FDDBC7}
        \definecolor{red}{HTML}{B2182B}
        \definecolor{blue}{HTML}{0264A2}
        \definecolor{darkblue}{HTML}{074B8A}

        \newcommand{\HE}[2]{
        \begin{tikzpicture}[overlay,remember picture,font=\huge]
        \node[below right,minimum width={\paperwidth-4cm},minimum height={0.1\paperheight},
        draw=darkblue,line width=1mm,fill=blue,text=white] at (current page.north west) {#1};
        \node[below left,minimum width={04cm},minimum height={0.1\paperheight},
        draw=darkblue,line width=1mm,fill=blue,text=blue] at (current page.north east)
        {\includegraphics[height=4cm]{#2}};
        \end{tikzpicture}
        }

    '''))

    # Utility to pretty print json
    # https://tex.stackexchange.com/questions/83085/how-to-improve-listings-display-of-json-files

    doc.preamble.append(pl.NoEscape(r'''

        \usepackage{bera}% optional: just to have a nice mono-spaced font
        \usepackage{listings}
        \usepackage{xcolor}

        \definecolor{eclipseStrings}{RGB}{42,0.0,255}
        \definecolor{eclipseKeywords}{RGB}{127,0,85}
        \colorlet{numb}{magenta!60!black}

        \lstdefinelanguage{json}{
            basicstyle=\normalfont\ttfamily,
            commentstyle=\color{eclipseStrings}, % style of comment
            stringstyle=\color{eclipseKeywords}, % style of strings
            numbers=left,
            numberstyle=\scriptsize,
            stepnumber=1,
            numbersep=8pt,
            showstringspaces=false,
            breaklines=true,
            frame=lines,
            backgroundcolor=\color{white}, %only if you like
            string=[s]{"}{"},
            comment=[l]{:\ "},
            morecomment=[l]{:"},
            literate=
                *{0}{{{\color{numb}0}}}{1}
                 {1}{{{\color{numb}1}}}{1}
                 {2}{{{\color{numb}2}}}{1}
                 {3}{{{\color{numb}3}}}{1}
                 {4}{{{\color{numb}4}}}{1}
                 {5}{{{\color{numb}5}}}{1}
                 {6}{{{\color{numb}6}}}{1}
                 {7}{{{\color{numb}7}}}{1}
                 {8}{{{\color{numb}8}}}{1}
                 {9}{{{\color{numb}9}}}{1}
        }

    '''))

    #  append map descriptors to metadata
    # todo move geodata fetching to reporting modules
    # entry_metadata["map"] = geodataGetter(location=REFERENCE_AREA, super_location=SUPER_REFERENCE_AREA).plot()

    # We generate files for each figure
    print("\tSaving figures.")
    for key in report.keys():
        if key != "metadata":
            if report[key]["analysis_outcomes"]:
                curr_el = report[key]["analysis_outcomes"]

                for subkey in curr_el.keys():
                    if isinstance(curr_el[subkey], plt.Figure):
                        # we redefine curr_el["key"] as the file path
                        curr_el[subkey] = safe_plot_intime_save(file_prefix=key, file_suffix=subkey,
                                                                object_to_plot=curr_el[subkey])

            area = report[key]["metadata"]["REFERENCE_AREA"]
            map = geodataGetter(location=area, super_location=report["metadata"]["super_reference_area"]).plot()
            if isinstance(map, plt.Figure):
                report[key]["metadata"]["map"] = safe_plot_intime_save(file_prefix=key, file_suffix="map",
                                                                       object_to_plot=map,dir = "./fig_map/")

    # We create a separate section for each report_entries entry
    print("\tGenerating report items.")

    for key in report.keys():
        if key != "metadata":

            # todo reinsert headers
            # # append page header
            header_string = r"\HE{" + key + "}{" + report[key]["metadata"]["map"] + "}"
            doc.append(pl.NoEscape(header_string))

            if report[key]["analysis_outcomes"]:
                curr_el = report[key]["analysis_outcomes"]

                # for each key, generate a report sub-item (subsection)
                for subkey in curr_el.keys():

                    with doc.create(pl.Section(subkey)) as section:

                        if isinstance(curr_el[subkey], dict) and any(
                                isinstance(i, dict) for i in curr_el[subkey].values()):
                            # if curr_el[subkey] is a nested dictionary, loop through keys
                            # (if not nested, we parse it out - via the else branch - as a json)
                            n_items = len(curr_el[subkey])
                            vspace = 1 / n_items
                            for entry in curr_el[subkey].keys():
                                print_block(section, curr_el[entry], vspace)
                        else:
                            vspace = 1
                            print_block(section, curr_el[subkey], vspace)

                # new page
                doc.append(pl.NewPage())

    filename_report = report["metadata"]["filenames"]["tex"]
    doc.generate_tex(filepath=filename_report.replace('.tex', ''))
    os.system("pdflatex '" + filename_report + "'  >/dev/null")
    os.system("pdflatex '" + filename_report + "'  >/dev/null")

    # try:
    #     shutil.rmtree("./fig")
    # except:
    #     pass

    try:
        shutil.rmtree("./fig_map")
    except:
        pass
    if os.path.exists(filename_report.replace('.tex', '.aux')):
        os.remove(filename_report.replace('.tex', '.aux'))
    if os.path.exists(filename_report.replace('.tex', '.log')):
        os.remove(filename_report.replace('.tex', '.log'))
    if os.path.exists(filename_report):
        os.remove(filename_report)
