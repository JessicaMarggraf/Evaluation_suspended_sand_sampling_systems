# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 10:17:46 2026

@author: jmarggra
"""


import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd


#%% Load data

data = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Data EDF broyage\Data_EDF_crushing.csv', sep = ',')
data_av = data.groupby('Grain_size_class_mum').mean(numeric_only=True)
data_std = data.groupby('Grain_size_class_mum').std(numeric_only=True)

#%% Plot 

# Data
labels = [r'63-200 $\rm{\mu}$m', r'200-500 $\rm{\mu}$m', r'500-1000 $\rm{\mu}$m']
cat1 = np.array(data_av['Percentage_not_crushed'])
cat2 = np.array(data_av['Percentage_crushed_sand'])
cat3 = np.array(data_av['Percentage_crushed_fines'])
err1 = np.array(data_std['Percentage_not_crushed'])
err2 = np.array(data_std['Percentage_crushed_sand'])
err3 = np.array(data_std['Percentage_crushed_fines'])

x = np.arange(len(labels))
width = 0.6

fig, ax = plt.subplots(figsize=(11, 8))

ax.bar(x, cat1, width, label='Not crushed', yerr = err1, capsize = 5, color='olive')
ax.bar(x, cat2, width, bottom=cat1, label=r'Crushed, > 63 $\rm{\mu}$m', color='darkseagreen')
ax.bar(x, cat3, width, bottom=cat1 + cat2, label=r'Crushed, < 63 $\rm{\mu}$m', color='#cab2d6')

ax.errorbar(x, cat1 + cat2,  yerr=err2, fmt='none', capsize=5, ecolor='black')
ax.errorbar(x, cat1 + cat2 + cat3,  yerr=err3, fmt='none', capsize=5, ecolor='black')

ax.text(-0.24, 81, '80 %', fontsize = 18)
ax.text(0.74, 83, '82 %', fontsize = 18)
ax.text(1.74, 73, '71 %', fontsize = 18)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize = 18)
ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
ax.set_ylim(0,100)
ax.set_ylabel('Proportion of initial mass (wt%)', weight = 'bold', fontsize = 20)
ax.legend(loc='lower center', fontsize = 18, ncols = 3, bbox_to_anchor = (0.5,-0.15))

plt.tight_layout()
fig.savefig(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Data EDF broyage\Sand_crushing_pump.png', dpi=300, bbox_inches='tight')






























