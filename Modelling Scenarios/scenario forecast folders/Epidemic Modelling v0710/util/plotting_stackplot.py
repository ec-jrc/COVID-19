# https://stackoverflow.com/questions/50404913/how-to-change-the-color-palette-for-stackplot-matplotlib
# https://python-graph-gallery.com/251-stacked-area-chart-with-seaborn-style/

# library
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data
x = range(1, 6)
y = [[1, 4, 6, 8, 9], [2, 2, 7, 10, 12], [2, 8, 5, 10, 6]]

# Plot
plt.stackplot(x, y, labels=['A', 'B', 'C'])
plt.legend(loc='upper left')
plt.show()

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import norm
# import matplotlib as mpl
# import matplotlib.font_manager as font_manager
#
# def plot_intime_stackplot(df):
#
#
# file = r'E:\FD\Barren_Mudflat\ChinaCoastal\Provinces\0ProvinceStat.csv'
# #set font property of legend
# font1 = {'family' : 'Times New Roman',
# 'weight' : 'normal',
# 'size'   : 16
# }
#
# #set size of figure
# fig, ax = plt.subplots()
# fig.set_size_inches(15, 7.5)
#
# #read values of dataframe
# value = df.values
# #plot stack area
# sp = ax.stackplot(Year, df.values)
# #set legend
# proxy = [mpl.patches.Rectangle((0,0), 0,0, facecolor=pol.get_facecolor()[0])
# for pol in sp]
# ax.legend(proxy, vol,prop = font1, loc='upper left', bbox_to_anchor=
# (0.01,1), ncol = 6)
#
# plt.xlim(1986,2016)
# plt.xticks([1986,1991,1996,2001,2006,2011,2016],fontproperties='Times New
# Roman', size = '16')
# plt.xlabel('Year',fontproperties='Times New Roman', size = '18')
# plt.ylim(0,1400)
# plt.yticks(np.arange(0,1500,200),fontproperties='Times New Roman', size =
# '16')
# plt.ylabel('Mudflat area (thousand ha)',fontproperties='Times New Roman',
# size = '18')
#
# #save fig: run this code before show()
# plt.savefig(r"E:\FD\Barren_Mudflat\ChinaCoastal\Provinces\stackplot.jpg",
# dpi = 600)
# plt.show()
