import math
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import numpy as np
import matplotlib.ticker as ticker


def plot_intime_comparison(cou, df_reference, df_comparison, fields_to_plot=[], labelsuffix_df_reference="",
                           labelsuffix_df_comparison=""):
    df_list = [df_comparison, df_reference]
    linestyle = [':', '-']
    label_suffix = [labelsuffix_df_comparison, labelsuffix_df_reference]

    # fp0 = []
    # appendDelta = True
    # for f in fields_to_plot:
    #     if 'delta_' in f:
    #         appendDelta = False
    #     fp0.append(f)
    #
    # if appendDelta:
    #     for f in fp0:
    #         fields_to_plot.append('delta_' + f)

    ncol = len(fields_to_plot)
    nrow = math.ceil(len(fields_to_plot) / ncol)

    fig, axes = plt.subplots(nrow, ncol)
    fig.suptitle(cou,fontweight="bold")

    # locator = MaxNLocator(nbins=3)  # with 3 bins you will have 4 ticks
    # axes.xaxis.set_major_locator(locator)

    # convert axes to 2D, if needed
    if axes.ndim == 1:
        axes = np.reshape(axes, (-1, axes.size))

    for count in range(nrow * ncol):
        c = count % ncol
        r = int((count - c) / ncol)

        if count > len(fields_to_plot) - 1:
            fig.delaxes(axes[r, c])
        else:
            legend_entries = []
            for ii in range(len(df_list)):
                if fields_to_plot[count] in df_list[ii]:
                    try:
                        df_list[ii].plot(ax=axes[r, c], use_index=True, y=[fields_to_plot[count]],
                                         linestyle=linestyle[ii])
                        legend_entries.append(fields_to_plot[count] + " (" + label_suffix[ii] + ")")
                    except:
                        pass

            axes[r, c].legend(legend_entries,fontsize=7)
            axes[r, c].grid()

            if count>0:
                axes[r, c].get_legend().remove()

            axes[r, c].set_title(fields_to_plot[count],pad=20,fontweight="bold",style='italic')
            axes[r, c].get_xaxis().set_major_locator(plt.NullLocator())
            axes[r, c].get_xaxis().set_minor_locator(plt.MaxNLocator(nbins=3))
            axes[r, c].tick_params(axis='both', which='major', labelsize=7)
            axes[r, c].tick_params(axis='both', which='minor', labelsize=7)

            try:
                # format dates when possible and relevant
                formatter = dates.DateFormatter('%d/%m')
                axes[r, c].xaxis.set_minor_formatter(formatter)
            except:
                pass

    # https://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle
    fig.tight_layout(rect=[0, 0.03, 1, 0.9]) #rect=[0, 0.03, 1, 0.95])
    # plt.show()
    # plt.close()

    return fig