# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:30:06 2022

@author: jlaible
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.legend_handler import HandlerTuple
from scipy.stats import linregress

#%% Load data

data_BD = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Comparison_BD_gaugings.csv', sep = ';')
data_P6 = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Comparison_P6_gaugings.csv', sep = ';')
data_P6_pump_point = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Comparison_pump_P6_point.csv', sep = ',') #
data_P6_pump_xs= pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Comparison_pump_P6_XS.csv', sep = ';') #
data_P6_pump_GSD = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Comparison_pump_P6_GSD.csv', sep = ';') #
data_BD_pump = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Comparison_BD_pump_point.csv', sep = ',') #

# cross-sectional area GC - for conversion of total cross-sectional flux to flux in kg/m2s to compare with the literature 
area_GC = pd.read_csv(r'C:\Users\jmarggra\Documents\INRAE\Grenoble_Campus\Autres\Area_cross_section.csv', sep=',')

# from literature
data_Dijkman_Milisic82 = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Dijkman_Milisic_1982_data.csv', sep = ',')
data_BW89_all = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Beverage_Williams_1989\Sand_flux_larger_0.062_Colorado_Mississippi.csv', sep=',')
data_BW89_62_125 = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Beverage_Williams_1989\Sand_flux_0.062_0.125_Colorado_Mississippi.csv', sep=',')
data_BW89_125_250 = pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Beverage_Williams_1989\Sand_flux_0.125_0.25_Colorado_Mississippi.csv', sep=',')

# Load data CNR pump testing
data_FBenoit =  pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Data_FBenoit.csv', sep = ',')
data_SPayen_labo =  pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Data_SPayen_labo.csv', sep = ',')
data_SPayen_field =  pd.read_csv(r'C:\Users\jmarggra\Documents\Sampler_article\Data\Data_SPayen_field.csv', sep = ',')

outpath_figures = r'C:\Users\jmarggra\Documents\Sampler_article\Data\Figures'

#%% Prepare data

# calculate absolute uncertainties
err_F_BD = data_BD['U_F']/100*data_BD['Sand_flux_kg_s']
err_C_BD = data_BD['U_C']/100*data_BD['Sand_concentration_g_l']
err_Q_BD = data_BD['U_Q']/100*data_BD['Q_sampling_m3_s']
err_F_P6 = data_P6['U_F']/100*data_P6['Sand_flux_kg_s']
err_C_P6 = data_P6['U_C']/100*data_P6['Sand_concentration_g_l']
err_Q_P6 = data_P6['U_Q']/100*data_P6['Q_sampling_m3_s']

#%% Cross-sectional regressions - y = ax

# x_range = np.linspace(0,1000,10)

# P6 - BD - Csand
# x = np.array(data_P6['Sand_concentration_g_l'])
# y = np.array(data_BD['Sand_concentration_g_l'])
# x = x[:,np.newaxis]
# slope_Csand_P6_BD, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
# lin_model_Csand_P6_BD = [data_P6['Sand_concentration_g_l'][i]*slope_Csand_P6_BD
#                            for i in range(len(data_P6['Sand_concentration_g_l']))]
# R2_Csand_P6_BD = r2_score(data_BD['Sand_concentration_g_l'], lin_model_Csand_P6_BD)
# x_range1 = np.linspace(0,1,10)
# lin_model_Csand_P6_BD_plot = [x_range1[i]*slope_Csand_P6_BD
#                            for i in range(len(x_range1))]

# # P6 - pump - Csand
# idx = np.isfinite(data_P6_pump_xs['Sand_concentration_P6']) & np.isfinite(data_P6_pump_xs['Sand_concentration_pump']) # only for valid data (no nan)
# pump_Csand = [(data_P6_pump_xs['Sand_concentration_pump'][i]) for i in range(len(data_P6_pump_xs['Sand_concentration_P6'])) if idx[i] == True]
# P6_Csand= [(data_P6_pump_xs['Sand_concentration_P6'][i]) for i in range(len(data_P6_pump_xs['Sand_concentration_P6'])) if idx[i] == True]
# y = np.array(pump_Csand)
# x = np.array(P6_Csand)
# x = x[:,np.newaxis]
# slope_Csand_P6_pump, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
# lin_model_Csand_P6_pump = [P6_Csand[i]*slope_Csand_P6_pump
#                            for i in range(len(P6_Csand))]
# R2_Csand_P6_pump = r2_score(pump_Csand, lin_model_Csand_P6_pump)
# x_range2 = np.linspace(0,400,10)
# lin_model_Csand_P6_pump_plot = [x_range2[i]*slope_Csand_P6_pump
#                            for i in range(len(x_range2))]

# # Pump - BD - Csand
# idx = np.isfinite(data_BD_pump['Sand_concentration_BD_ADCP']) & np.isfinite(data_BD_pump['Sand_concentration_pump']) # only for valid data (no nan)
# BD_point = [(data_BD_pump['Sand_concentration_BD_ADCP'][i]) for i in range(len(data_BD_pump['Sand_concentration_BD_ADCP'])) if idx[i] == True]
# pump_point = [(data_BD_pump['Sand_concentration_pump'][i]) for i in range(len(data_BD_pump['Sand_concentration_pump'])) if idx[i] == True]
# x = np.array(pump_point)
# y = np.array(BD_point)
# x = x[:,np.newaxis]
# slope_Csand_pump_BD_point, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
# lin_model_Csand_pump_BD_point = [pump_point[i]*slope_Csand_pump_BD_point
#                            for i in range(len(pump_point))]
# R2_Csand_pump_BD_point = r2_score(BD_point, lin_model_Csand_pump_BD_point)
# x_range4 = np.linspace(0,5,10)
# lin_model_Csand_pump_BD_point_plot = [x_range4[i]*slope_Csand_pump_BD_point
#                            for i in range(len(x_range4))]

# # P6 - pump - Csand 
# x = np.array(data_P6_pump_point['Sand_concentration_P6'])
# y = np.array(data_P6_pump_point['Sand_concentration_pump'])
# x = x[:,np.newaxis]
# slope_Csand_P6_pump_point, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
# lin_model_Csand_P6_pump_point = [data_P6_pump_point['Sand_concentration_P6'][i]*slope_Csand_P6_pump_point
#                             for i in range(len(data_P6_pump_point['Sand_concentration_P6']))]
# R2_Csand_P6_pump_point = r2_score(data_P6_pump_point['Sand_concentration_pump'], lin_model_Csand_P6_pump_point)
# x_range3 = np.linspace(0,5,10)
# lin_model_Csand_P6_pump_point_plot = [x_range3[i]*slope_Csand_P6_pump_point
#                             for i in range(len(x_range3))]

# # GSD pump - P6
# x_range_GSD = np.linspace(0,1000,20)
# number_GSD_P6_BD = data_P6.dropna(subset = ['D90_mum']).shape[0]

# # P6 - pump - D50 
# y = np.array(data_P6_pump_GSD['D50_P6'])
# x = np.array(data_P6_pump_GSD['D50_pump'])
# log_x = np.log10(x)
# log_y = np.log10(y)
# log_x_reshaped = log_x[:, np.newaxis]

# slope_D50_pump_P6_point_log, _, _, _ = np.linalg.lstsq(log_x_reshaped, log_y, rcond=None)
# lin_model_D50_pump_P6_point_log = [log_x[i]*slope_D50_pump_P6_point_log
#                            for i in range(len(log_x))]
# R2_D50_pump_P6_point_log = r2_score(log_y, lin_model_D50_pump_P6_point_log)
# x_range5_log = np.logspace(np.log10(10), np.log10(1000), 100)
# lin_model_D50_pump_P6_point_log_plot = 10**(slope_D50_pump_P6_point_log*np.log10(x_range5_log))

#%% Cross-sectional regressions - y = ax +b

x_range = np.linspace(0,1000,10)

# P6 - BD - Csand
lin_model_Csand_P6_BD = linregress(data_P6['Sand_concentration_g_l'], data_BD['Sand_concentration_g_l'])
x_range1 = np.linspace(0,1,10000)
lin_model_Csand_P6_BD_plot = [x_range1[i]*lin_model_Csand_P6_BD.slope + lin_model_Csand_P6_BD.intercept
                           for i in range(len(x_range1))]

# P6 - pump - Csand
lin_model_Csand_P6_pump_xs = linregress(data_P6_pump_xs['Sand_concentration_P6'], data_P6_pump_xs['Sand_concentration_pump'])
x_range2 = np.linspace(0,4,10)
lin_model_Csand_P6_pump_xs_plot = [x_range2[i]*lin_model_Csand_P6_pump_xs.slope + lin_model_Csand_P6_pump_xs.intercept
                           for i in range(len(x_range2))]

# Pump - BD - Csand
idx = np.isfinite(data_BD_pump['Sand_concentration_BD_ADCP']) & np.isfinite(data_BD_pump['Sand_concentration_pump']) # only for valid data (no nan)
BD_point = [(data_BD_pump['Sand_concentration_BD_ADCP'][i]) for i in range(len(data_BD_pump['Sand_concentration_BD_ADCP'])) if idx[i] == True]
pump_point = [(data_BD_pump['Sand_concentration_pump'][i]) for i in range(len(data_BD_pump['Sand_concentration_pump'])) if idx[i] == True]

lin_model_Csand_pump_BD = linregress(pump_point, BD_point)
x_range4 = np.linspace(0,5,100000)
lin_model_Csand_pump_BD_plot = [x_range4[i]*lin_model_Csand_pump_BD.slope + lin_model_Csand_pump_BD.intercept
                           for i in range(len(x_range4))]

# P6 - pump - Csand 
lin_model_Csand_P6_pump = linregress(data_P6_pump_point['Sand_concentration_P6'], data_P6_pump_point['Sand_concentration_pump'])
x_range2 = np.linspace(0,4,10)
lin_model_Csand_P6_pump_plot = [x_range2[i]*lin_model_Csand_P6_pump.slope + lin_model_Csand_P6_pump.intercept
                           for i in range(len(x_range2))]

# GSD pump - P6
number_GSD_P6_BD = data_P6.dropna(subset = ['D90_mum']).shape[0]

# P6 - pump - D50 
lin_model_Csand_P6_pump_GSD = linregress(data_P6_pump_GSD['D50_P6'], data_P6_pump_GSD['D50_pump'])
x_range_GSD = np.linspace(0,1000,20)
lin_model_Csand_P6_pump_GSD_plot = [x_range_GSD[i]*lin_model_Csand_P6_pump_GSD.slope + lin_model_Csand_P6_pump_GSD.intercept
                           for i in range(len(x_range_GSD))]


#%%################################################################################################

# Fig3

################################################################################################
fig, ax = plt.subplots(3,2, figsize=(12, 18), dpi=200)

# Delft Csand
p1 = ax[0,0].errorbar(data_P6['Sand_concentration_g_l'], data_BD['Sand_concentration_g_l'],          
        marker = 'D', color = 'darkolivegreen', markersize = 8, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2,
        yerr = err_C_BD, xerr = err_C_P6, elinewidth = 1, capsize = 1.5, zorder = 30)

ax[0,0].plot(2*x_range, x_range,  
        ls = ':', lw = 1, color = 'black')
p4, = ax[0,0].plot(x_range, x_range,  
        ls = '--', lw = 1, color = 'black')
p5, = ax[0,0].plot(x_range, 2*x_range,  
        ls = ':', lw = 1, color = 'black', label = 'Error factor 2')

ax[0,0].plot(x_range1, lin_model_Csand_P6_BD_plot, lw = 1.5, ls = '-', color = 'darkolivegreen', zorder = 0)

ax[0,0].text(0.2, 0.95, ('y = ' + str(np.round(lin_model_Csand_P6_BD.slope,2)) + 'x + ' + str(np.round(lin_model_Csand_P6_BD.intercept,2))), 
        color = 'black', fontsize = 14, transform = ax[0,0].transAxes)
ax[0,0].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(lin_model_Csand_P6_BD.rvalue**2,2))) + ', n = ' + str(len(data_P6))), 
        color = 'black', fontsize = 14, transform = ax[0,0].transAxes)

ax[0,0].text(0.9, 0.83, '1:1', transform =  ax[0,0].transAxes,
        color = 'black', fontsize = 12)
ax[0,0].text(0.05, 0.95, 'a)', fontsize = 14, transform = ax[0,0].transAxes)
ax[0,0].set_ylabel(r'$\mathregular{\overline{C_{sand, DB}}}$ (g/l)', fontsize=18, weight = 'bold')
ax[0,0].set_xlabel(r'$\mathregular{\overline{C_{sand, P6}}}$ (g/l)', fontsize=18, weight = 'bold')
ax[0,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[0,0].set_ylim(0.005,1)
ax[0,0].set_xlim(0.005,1)  

# P6 - BD - D50
p2, = ax[1,0].plot(data_P6['D90_mum'], data_BD['D90_mum'],           
        marker = '^', color = 'brown', markersize = 10, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30,
        label = r'$\mathregular{D_{90}}$')
p3, = ax[1,0].plot(data_P6['D50_mum'], data_BD['D50_mum'],           
        marker = 's', color = 'darkgoldenrod', markersize = 8, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30,
        label = r'$\mathregular{D_{50}}$')
p4, = ax[1,0].plot(data_P6['D10_mum'], data_BD['D10_mum'],           
        marker = 'v', color = 'royalblue', markersize = 10, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30,
        label = r'$\mathregular{D_{10}}$')

ax[1,0].text(0.05, 0.95, 'c)', fontsize = 14, transform = ax[1,0].transAxes)
ax[1,0].plot(2*x_range_GSD, x_range_GSD,  
        ls = ':', lw = 1, color = 'black')
p5, = ax[1,0].plot(x_range_GSD, x_range_GSD,  
        ls = '--', lw = 1, color = 'black')
p6, = ax[1,0].plot(x_range_GSD, 2*x_range_GSD,  
        ls = ':', lw = 1, color = 'black', label = 'Error factor 2')
ax[1,0].text(0.9, 0.83, '1:1', transform =  ax[1,0].transAxes,
        color = 'black', fontsize = 12)
ax[1,0].text(0.2, 0.95, ('n = ' + str(number_GSD_P6_BD)), 
        color = 'black', fontsize = 14, transform = ax[1,0].transAxes)

ax[1,0].set_ylabel(r'$\mathregular{\overline{D_{sand, DB}}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
ax[1,0].set_xlabel(r'$\mathregular{\overline{D_{sand, P6}}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
ax[1,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[1,0].set_xlim(10,1000)
ax[1,0].set_ylim(10,1000) 
ax[1,0].legend(fontsize =14, loc = 'lower right', framealpha = 1)

# Csand Delft pump
ax[2,0].plot(data_BD_pump['Sand_concentration_pump'],data_BD_pump['Sand_concentration_BD_ADCP'],        
        marker = 'o', color = 'sandybrown', markersize = 8, ls = '',  alpha = 0.7,
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 20, label = 'ADCP')
ax[2,0].plot(data_BD_pump['Sand_concentration_pump'],data_BD_pump['Sand_concentration_BD_moulinet'],        
        marker = 'P', color = 'darkred', markersize = 9, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 20, label = 'Reel')

ax[2,0].plot(x_range4, lin_model_Csand_pump_BD_plot, lw = 1.5, ls = '-', color = 'peru', zorder = 0)
ax[2,0].text(0.9, 0.83, '1:1', transform =  ax[2,0].transAxes,
        color = 'black', fontsize = 12)

ax[2,0].plot(2*x_range, x_range,  
        ls = ':', lw = 1, color = 'black')
p4, = ax[2,0].plot(x_range, x_range,  
        ls = '--', lw = 1, color = 'black')
p5, = ax[2,0].plot(x_range, 2*x_range,  
        ls = ':', lw = 1, color = 'black')

ax[2,0].text(0.05, 0.95, 'e)', 
        color = 'black', fontsize = 14, transform = ax[2,0].transAxes)
ax[2,0].text(0.2, 0.95, ('y = ' + str(np.round(lin_model_Csand_pump_BD.slope,2)) + 'x + '+ str(np.round(lin_model_Csand_pump_BD.intercept,2))), 
        color = 'black', fontsize = 14, transform = ax[2,0].transAxes)
ax[2,0].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(lin_model_Csand_pump_BD.rvalue**2,2))) + ', n = ' + str(len(data_BD_pump))), 
        color = 'black', fontsize = 14, transform = ax[2,0].transAxes)

ax[2,0].set_ylabel(r'$\mathregular{C_{sand, DB}}$ (g/l)', fontsize=18, weight = 'bold')
ax[2,0].set_xlabel(r'$\mathregular{C_{sand, pump}}$ (g/l)', fontsize=18, weight = 'bold')
ax[2,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[2,0].set_xlim(0.005,5)
ax[2,0].set_ylim(0.005,5)  

ax[2,0].legend(fontsize =14, loc = 'lower right', framealpha = 1)
leg = ax[2,0].legend(loc = 'lower right', fontsize =14, framealpha = 1)
frame = leg.get_frame()
frame.set_facecolor('white')

# pump Csand
ax[0,1].plot(2*x_range, x_range,  
        ls = ':', lw = 1, color = 'black')
p4, = ax[0,1].plot(x_range, x_range,  
        ls = '--', lw = 1, color = 'black')
p5, = ax[0,1].plot(x_range, 2*x_range,  
        ls = ':', lw = 1, color = 'black', label = 'Error factor 2')

ax[0,1].plot(data_P6_pump_xs['Sand_concentration_P6'], data_P6_pump_xs['Sand_concentration_pump'],            
        marker = 'D', color = 'darkolivegreen', markersize = 8, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 15,
        label = 'XS')
ax[0,1].plot(x_range2, lin_model_Csand_P6_pump_xs_plot, lw = 1.5, ls = '-', color = 'darkolivegreen', zorder = 20)

ax[0,1].text(0.2, 0.95, ('y = ' + str(np.round(lin_model_Csand_P6_pump_xs.slope,2)) + 'x + '  + str(np.round(lin_model_Csand_P6_pump_xs.intercept,2))), 
        color = 'black', fontsize = 14, transform = ax[0,1].transAxes)
ax[0,1].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(lin_model_Csand_P6_pump_xs.rvalue**2,2))) + ', n = ' + str(len(data_P6_pump_xs))), 
        color = 'black', fontsize = 14, transform = ax[0,1].transAxes)
ax[0,1].text(0.8, 0.9, '1:1', transform =  ax[0,1].transAxes,
        color = 'black', fontsize = 12)

ax[0,1].text(0.05, 0.95, 'b)', fontsize = 14, transform = ax[0,1].transAxes)
ax[0,1].set_xlabel(r'$\mathregular{\overline{C_{sand, P6}}}$ (g/l)', fontsize=18, weight = 'bold')
ax[0,1].set_ylabel(r'$\mathregular{\overline{C_{sand, pump}}}$ (g/l)', fontsize=18, weight = 'bold')
ax[0,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[0,1].set_xlim(0.005,1)
ax[0,1].set_ylim(0.005,1) 

# pump GSD
ax[1,1].plot(data_P6_pump_GSD['D90_P6'], data_P6_pump_GSD['D90_pump'],            
        marker = '^', color = 'brown', markersize = 10, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30, 
        label = r'$\mathregular{D_{90}}$')
ax[1,1].plot(data_P6_pump_GSD['D50_P6'], data_P6_pump_GSD['D50_pump'],            
        marker = 's', color = 'darkgoldenrod', markersize = 8, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30, 
        label = r'$\mathregular{D_{50}}$')
ax[1,1].plot(data_P6_pump_GSD['D10_P6'], data_P6_pump_GSD['D10_pump'],            
        marker = 'v', color = 'royalblue', markersize = 10, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30, 
        label = r'$\mathregular{D_{10}}$')

ax[1,1].plot(x_range_GSD, lin_model_Csand_P6_pump_GSD_plot, lw = 1.5, ls = '-', color = 'darkgoldenrod', zorder = 20)
ax[1,1].text(0.2, 0.95, ('y = ' + str(np.round(lin_model_Csand_P6_pump_GSD.slope,2)) + 'x + '  + str(np.round(lin_model_Csand_P6_pump_GSD.intercept,2))), 
        color = 'black', fontsize = 14, transform = ax[1,1].transAxes)
ax[1,1].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(lin_model_Csand_P6_pump_GSD.rvalue**2,2))) + ', n = ' + str(len(data_P6_pump_GSD))), 
        color = 'black', fontsize = 14, transform = ax[1,1].transAxes)


ax[1,1].text(0.9, 0.83, '1:1', transform =  ax[1,1].transAxes,
        color = 'black', fontsize = 12)
ax[1,1].plot(2*x_range_GSD, x_range_GSD,  
        ls = ':', lw = 1, color = 'black')
p4, = ax[1,1].plot(x_range_GSD, x_range_GSD,  
        ls = '--', lw = 1, color = 'black')
p5, = ax[1,1].plot(x_range_GSD, 2*x_range_GSD,  
        ls = ':', lw = 1, color = 'black', label = 'Error factor 2')
# ax[1,1].text(0.2, 0.95, ('n = ' + str(len(data_P6_pump_GSD['D90_pump']))), 
#         color = 'black', fontsize = 14, transform = ax[1,1].transAxes)

ax[1,1].text(0.05, 0.95, 'd)', fontsize = 14, transform = ax[1,1].transAxes)
ax[1,1].set_xlabel(r'$\mathregular{D_{sand, P6}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
ax[1,1].set_ylabel(r'$\mathregular{D_{sand, pump}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
ax[1,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[1,1].set_xlim(10,1000)
ax[1,1].set_ylim(10,1000) 
# ax[1,1].legend(fontsize =14, loc = 'lower right', framealpha = 1)

# P6 - pump -points
ax[2,1].plot(data_P6_pump_point['Sand_concentration_P6'], data_P6_pump_point['Sand_concentration_pump'],            
        marker = 'o', color = 'sandybrown', markersize = 8, ls = '', alpha = 0.7,
        markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 10,
        label = 'Point samples')
ax[2,1].plot(x_range2, lin_model_Csand_P6_pump_plot, lw = 1.5, ls = '-', color = 'peru', zorder = 15)

ax[2,1].plot(2*x_range, x_range,  
        ls = ':', lw = 1, color = 'black')
p4, = ax[2,1].plot(x_range, x_range,  
        ls = '--', lw = 1, color = 'black')
p5, = ax[2,1].plot(x_range, 2*x_range,  
        ls = ':', lw = 1, color = 'black', label = 'Error factor 2')

ax[2,1].text(0.2, 0.95, ('y = ' + str(np.round(lin_model_Csand_P6_pump.slope,2)) + 'x + ' + str(np.round(lin_model_Csand_P6_pump.intercept,2))), 
        color = 'black', fontsize = 14, transform = ax[2,1].transAxes)
ax[2,1].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(lin_model_Csand_P6_pump.rvalue**2,2))) + ', n = ' + str(len(data_P6_pump_point))), 
        color = 'black', fontsize = 14, transform = ax[2,1].transAxes)

ax[2,1].text(0.8, 0.9, '1:1', transform =  ax[2,1].transAxes,
        color = 'black', fontsize = 12)

ax[2,1].text(0.05, 0.95, 'f)', fontsize = 14, transform = ax[2,1].transAxes)
ax[2,1].set_xlabel(r'$\mathregular{C_{sand, P6}}$ (g/l)', fontsize=18, weight = 'bold')
ax[2,1].set_ylabel(r'$\mathregular{C_{sand, pump}}$ (g/l)', fontsize=18, weight = 'bold')
ax[2,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[2,1].set_xlim(0.005,5)
ax[2,1].set_ylim(0.005,5) 

ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')
ax[1,1].set_xscale('log')
ax[1,1].set_yscale('log')
ax[2,0].set_xscale('log')
ax[2,0].set_yscale('log')
ax[2,1].set_xscale('log')
ax[2,1].set_yscale('log')

fig.tight_layout()
figname = 'Fig3'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')

#%%################################################################################################

# # Fig3

# ################################################################################################
# fig, ax = plt.subplots(3,2, figsize=(12, 18), dpi=200)

# # Delft Csand
# p1 = ax[0,0].errorbar(data_P6['Sand_concentration_g_l'], data_BD['Sand_concentration_g_l'],          
#         marker = 'D', color = 'darkolivegreen', markersize = 8, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2,
#         yerr = err_C_BD, xerr = err_C_P6, elinewidth = 1, capsize = 1.5, zorder = 30)

# ax[0,0].plot(2*x_range, x_range,  
#         ls = ':', lw = 1, color = 'black')
# p4, = ax[0,0].plot(x_range, x_range,  
#         ls = '--', lw = 1, color = 'black')
# p5, = ax[0,0].plot(x_range, 2*x_range,  
#         ls = ':', lw = 1, color = 'black', label = 'Error factor 2')

# ax[0,0].plot(x_range1, lin_model_Csand_P6_BD_plot, lw = 1.5, ls = '-', color = 'darkolivegreen', zorder = 0)

# ax[0,0].text(0.2, 0.95, ('y = ' + str(np.round(slope_Csand_P6_BD[0],2)) + 'x'), 
#         color = 'black', fontsize = 14, transform = ax[0,0].transAxes)
# ax[0,0].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(R2_Csand_P6_BD,2))) + ', n = ' + str(len(data_P6))), 
#         color = 'black', fontsize = 14, transform = ax[0,0].transAxes)

# ax[0,0].text(0.9, 0.83, '1:1', transform =  ax[0,0].transAxes,
#         color = 'black', fontsize = 12)
# ax[0,0].text(0.05, 0.95, 'a)', fontsize = 14, transform = ax[0,0].transAxes)
# ax[0,0].set_ylabel(r'$\mathregular{\overline{C_{sand, DB}}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[0,0].set_xlabel(r'$\mathregular{\overline{C_{sand, P6}}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[0,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
# ax[0,0].set_ylim(0.005,1)
# ax[0,0].set_xlim(0.005,1)  

# # P6 - BD - D50
# p2, = ax[1,0].plot(data_P6['D90_mum'], data_BD['D90_mum'],           
#         marker = '^', color = 'brown', markersize = 10, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30,
#         label = r'$\mathregular{D_{90}}$')
# p3, = ax[1,0].plot(data_P6['D50_mum'], data_BD['D50_mum'],           
#         marker = 's', color = 'darkgoldenrod', markersize = 8, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30,
#         label = r'$\mathregular{D_{50}}$')
# p4, = ax[1,0].plot(data_P6['D10_mum'], data_BD['D10_mum'],           
#         marker = 'v', color = 'royalblue', markersize = 10, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30,
#         label = r'$\mathregular{D_{10}}$')

# ax[1,0].text(0.05, 0.95, 'c)', fontsize = 14, transform = ax[1,0].transAxes)
# ax[1,0].plot(2*x_range_GSD, x_range_GSD,  
#         ls = ':', lw = 1, color = 'black')
# p5, = ax[1,0].plot(x_range_GSD, x_range_GSD,  
#         ls = '--', lw = 1, color = 'black')
# p6, = ax[1,0].plot(x_range_GSD, 2*x_range_GSD,  
#         ls = ':', lw = 1, color = 'black', label = 'Error factor 2')
# ax[1,0].text(0.9, 0.83, '1:1', transform =  ax[1,0].transAxes,
#         color = 'black', fontsize = 12)
# ax[1,0].text(0.2, 0.95, ('n = ' + str(number_GSD_P6_BD)), 
#         color = 'black', fontsize = 14, transform = ax[1,0].transAxes)

# ax[1,0].set_ylabel(r'$\mathregular{\overline{D_{sand, DB}}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
# ax[1,0].set_xlabel(r'$\mathregular{\overline{D_{sand, P6}}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
# ax[1,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
# ax[1,0].set_xlim(10,1000)
# ax[1,0].set_ylim(10,1000) 
# ax[1,0].legend(fontsize =14, loc = 'lower right', framealpha = 1)

# # Csand Delft pump
# ax[2,0].plot(data_BD_pump['Sand_concentration_pump'],data_BD_pump['Sand_concentration_BD_ADCP'],        
#         marker = 'o', color = 'sandybrown', markersize = 8, ls = '',  alpha = 0.7,
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 20, label = 'ADCP')
# ax[2,0].plot(data_BD_pump['Sand_concentration_pump'],data_BD_pump['Sand_concentration_BD_moulinet'],        
#         marker = 'P', color = 'darkred', markersize = 9, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 20, label = 'Reel')

# ax[2,0].plot(x_range4, lin_model_Csand_pump_BD_point_plot, lw = 1.5, ls = '-', color = 'peru', zorder = 0)
# ax[2,0].text(0.9, 0.83, '1:1', transform =  ax[2,0].transAxes,
#         color = 'black', fontsize = 12)

# ax[2,0].plot(2*x_range, x_range,  
#         ls = ':', lw = 1, color = 'black')
# p4, = ax[2,0].plot(x_range, x_range,  
#         ls = '--', lw = 1, color = 'black')
# p5, = ax[2,0].plot(x_range, 2*x_range,  
#         ls = ':', lw = 1, color = 'black')

# ax[2,0].text(0.05, 0.95, 'e)', 
#         color = 'black', fontsize = 14, transform = ax[2,0].transAxes)
# ax[2,0].text(0.2, 0.95, ('y = ' + str(np.round(slope_Csand_pump_BD_point[0],2)) + 'x'), 
#         color = 'black', fontsize = 14, transform = ax[2,0].transAxes)
# ax[2,0].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(R2_Csand_pump_BD_point,2))) + ', n = ' + str(len(data_BD_pump))), 
#         color = 'black', fontsize = 14, transform = ax[2,0].transAxes)

# ax[2,0].set_ylabel(r'$\mathregular{C_{sand, DB}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[2,0].set_xlabel(r'$\mathregular{C_{sand, pump}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[2,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
# ax[2,0].set_xlim(0.005,5)
# ax[2,0].set_ylim(0.005,5)  

# ax[2,0].legend(fontsize =14, loc = 'lower right', framealpha = 1)
# leg = ax[2,0].legend(loc = 'lower right', fontsize =14, framealpha = 1)
# frame = leg.get_frame()
# frame.set_facecolor('white')

# # pump Csand
# ax[0,1].plot(2*x_range, x_range,  
#         ls = ':', lw = 1, color = 'black')
# p4, = ax[0,1].plot(x_range, x_range,  
#         ls = '--', lw = 1, color = 'black')
# p5, = ax[0,1].plot(x_range, 2*x_range,  
#         ls = ':', lw = 1, color = 'black', label = 'Error factor 2')

# ax[0,1].plot(data_P6_pump_xs['Sand_concentration_P6'], data_P6_pump_xs['Sand_concentration_pump'],            
#         marker = 'D', color = 'darkolivegreen', markersize = 8, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 15,
#         label = 'XS')
# ax[0,1].plot(x_range2, lin_model_Csand_P6_pump_plot, lw = 1.5, ls = '-', color = 'darkolivegreen', zorder = 20)

# ax[0,1].text(0.2, 0.95, ('y = ' + str(np.round(slope_Csand_P6_pump[0],2)) + 'x'), 
#         color = 'black', fontsize = 14, transform = ax[0,1].transAxes)
# ax[0,1].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(R2_Csand_P6_pump,2))) + ', n = ' + str(len(data_P6_pump_xs))), 
#         color = 'black', fontsize = 14, transform = ax[0,1].transAxes)
# ax[0,1].text(0.8, 0.9, '1:1', transform =  ax[0,1].transAxes,
#         color = 'black', fontsize = 12)

# ax[0,1].text(0.05, 0.95, 'b)', fontsize = 14, transform = ax[0,1].transAxes)
# ax[0,1].set_xlabel(r'$\mathregular{\overline{C_{sand, P6}}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[0,1].set_ylabel(r'$\mathregular{\overline{C_{sand, pump}}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[0,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
# ax[0,1].set_xlim(0.005,1)
# ax[0,1].set_ylim(0.005,1) 

# # pump GSD
# ax[1,1].plot(data_P6_pump_GSD['D90_P6'], data_P6_pump_GSD['D90_pump'],            
#         marker = '^', color = 'brown', markersize = 10, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30, 
#         label = r'$\mathregular{D_{90}}$')
# ax[1,1].plot(data_P6_pump_GSD['D50_P6'], data_P6_pump_GSD['D50_pump'],            
#         marker = 's', color = 'darkgoldenrod', markersize = 8, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30, 
#         label = r'$\mathregular{D_{50}}$')
# ax[1,1].plot(data_P6_pump_GSD['D10_P6'], data_P6_pump_GSD['D10_pump'],            
#         marker = 'v', color = 'royalblue', markersize = 10, ls = '', 
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 30, 
#         label = r'$\mathregular{D_{10}}$')

# ax[1,1].plot(x_range_GSD, lin_model_D50_pump_P6_point_log_plot, lw = 1.5, ls = '-', color = 'darkgoldenrod', zorder = 20)


# ax[1,1].text(0.9, 0.83, '1:1', transform =  ax[1,1].transAxes,
#         color = 'black', fontsize = 12)
# ax[1,1].plot(2*x_range_GSD, x_range_GSD,  
#         ls = ':', lw = 1, color = 'black')
# p4, = ax[1,1].plot(x_range_GSD, x_range_GSD,  
#         ls = '--', lw = 1, color = 'black')
# p5, = ax[1,1].plot(x_range_GSD, 2*x_range_GSD,  
#         ls = ':', lw = 1, color = 'black', label = 'Error factor 2')
# ax[1,1].text(0.2, 0.95, ('n = ' + str(len(data_P6_pump_GSD['D90_pump']))), 
#         color = 'black', fontsize = 14, transform = ax[1,1].transAxes)

# ax[1,1].text(0.05, 0.95, 'd)', fontsize = 14, transform = ax[1,1].transAxes)
# ax[1,1].set_xlabel(r'$\mathregular{D_{sand, P6}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
# ax[1,1].set_ylabel(r'$\mathregular{D_{sand, pump}}$ ($\mathregular{\mu}$m)', fontsize=18, weight = 'bold')
# ax[1,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
# ax[1,1].set_xlim(10,1000)
# ax[1,1].set_ylim(10,1000) 
# # ax[1,1].legend(fontsize =14, loc = 'lower right', framealpha = 1)

# # P6 - pump -points
# ax[2,1].plot(data_P6_pump_point['Sand_concentration_P6'], data_P6_pump_point['Sand_concentration_pump'],            
#         marker = 'o', color = 'sandybrown', markersize = 8, ls = '', alpha = 0.7,
#         markeredgecolor = 'black', markeredgewidth = 0.2, zorder = 10,
#         label = 'Point samples')
# ax[2,1].plot(x_range3, lin_model_Csand_P6_pump_point_plot, lw = 1.5, ls = '-', color = 'peru', zorder = 15)

# ax[2,1].plot(2*x_range, x_range,  
#         ls = ':', lw = 1, color = 'black')
# p4, = ax[2,1].plot(x_range, x_range,  
#         ls = '--', lw = 1, color = 'black')
# p5, = ax[2,1].plot(x_range, 2*x_range,  
#         ls = ':', lw = 1, color = 'black', label = 'Error factor 2')

# ax[2,1].text(0.2, 0.95, ('y = ' + str(np.round(slope_Csand_P6_pump_point[0],2)) + 'x'), 
#         color = 'black', fontsize = 14, transform = ax[2,1].transAxes)
# ax[2,1].text(0.2, 0.89, ('R² = ' + str(np.float64(np.round(R2_Csand_P6_pump_point,2))) + ', n = ' + str(len(data_P6_pump_point))), 
#         color = 'black', fontsize = 14, transform = ax[2,1].transAxes)

# ax[2,1].text(0.8, 0.9, '1:1', transform =  ax[2,1].transAxes,
#         color = 'black', fontsize = 12)

# ax[2,1].text(0.05, 0.95, 'f)', fontsize = 14, transform = ax[2,1].transAxes)
# ax[2,1].set_xlabel(r'$\mathregular{C_{sand, P6}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[2,1].set_ylabel(r'$\mathregular{C_{sand, pump}}$ (g/l)', fontsize=18, weight = 'bold')
# ax[2,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
# ax[2,1].set_xlim(0.005,5)
# ax[2,1].set_ylim(0.005,5) 

# ax[0,0].set_xscale('log')
# ax[0,0].set_yscale('log')
# ax[0,1].set_xscale('log')
# ax[0,1].set_yscale('log')
# ax[1,0].set_xscale('log')
# ax[1,0].set_yscale('log')
# ax[1,1].set_xscale('log')
# ax[1,1].set_yscale('log')
# ax[2,0].set_xscale('log')
# ax[2,0].set_yscale('log')
# ax[2,1].set_xscale('log')
# ax[2,1].set_yscale('log')

# fig.tight_layout()
# figname = 'Fig3'
# fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# # fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')

#%%################################################################################################

# Figure 4

#################################################################################################
# Pump lab and field data to determine independence from isocinetics
# a) Ratio C - Ratio v: pompe avec buse et crepine, parallel et perpendiculaire; fig 19 F Benoit - continu et crepine
# types_FBenoit = data_FBenoit['Type'].unique()
# markerss = ['o', 'o', 'H', 'H', 'd', 'd', 'd', 'D', 'D']
# colorss = ['blue', 'blue', 'cornflowerblue', 'cornflowerblue', 'saddlebrown', 'saddlebrown', 'saddlebrown', 'goldenrod', 'goldenrod']
# sizess = [7, 12, 7, 12, 7, 10, 12, 7, 10] 
# # comment differencier entre parallel et perpendiculaire? 

# # b) 
# types_SPayen_labo = data_SPayen_labo['Label'].unique()
# markers_SP_labo = ['d', 'd', 'D', 'D']
# colorss_SP_labo = ['saddlebrown', 'saddlebrown', 'goldenrod', 'goldenrod']
# # color 1ere et 2eme run differemment ou ce n'est pas important?

types_FBenoit = data_FBenoit['Type'].unique()
markerss = ['o', 'o', 's', 's', '^', '^', '^', 'D', 'D']
colorss = ['indigo', 'indigo', 'crimson', 'crimson', 'teal', 'teal', 'teal', 'goldenrod', 'goldenrod']
sizess = [7, 12, 7, 12, 7, 10, 12, 7, 10] 
# comment differencier entre parallel et perpendiculaire? 

# b) 
types_SPayen_labo = data_SPayen_labo['Label'].unique()
markers_SP_labo = ['^', '^', 'D', 'D']
colorss_SP_labo = ['teal', 'teal', 'goldenrod', 'goldenrod']


fig, ax = plt.subplots(1,1, figsize=(10, 8), dpi=200)

for i in range(len(types_FBenoit)):
    ax.plot(data_FBenoit['Ratio_velocity'][data_FBenoit['Type']==types_FBenoit[i]], 
                  data_FBenoit['Ratio_Conc'][data_FBenoit['Type']==types_FBenoit[i]],            
                marker = markerss[i], color = 'None', markersize = sizess[i], ls = '', 
                markeredgecolor = colorss[i], markeredgewidth = 1.5, zorder = 15)

for i in range(len(types_SPayen_labo)):
    ax.plot(data_SPayen_labo['Ratio_velocity'][data_SPayen_labo['Label']==types_SPayen_labo[i]], 
                  data_SPayen_labo['Ratio_conc_norm'][data_SPayen_labo['Label']==types_SPayen_labo[i]],            
                marker = markers_SP_labo[i], color = 'None', markersize = 7, ls = '', 
                markeredgecolor = colorss_SP_labo[i], markeredgewidth = 1.5, zorder = 15)
        
# labels
p1, = ax.plot(data_FBenoit['Ratio_velocity'][0], data_FBenoit['Ratio_Conc'][0],            
              marker = markerss[0], color = 'None', markersize = sizess[0], ls = '', 
              markeredgecolor = colorss[0], markeredgewidth = 1.5, zorder = 0, label = 'Nozzle parallel')
p2, = ax.plot(data_FBenoit['Ratio_velocity'][30], data_FBenoit['Ratio_Conc'][30],            
              marker = markerss[3], color = 'None', markersize = sizess[0], ls = '', 
              markeredgecolor = colorss[3], markeredgewidth = 1.5, zorder = 0, label = 'Nozzle perpendicular')
p3, = ax.plot(data_FBenoit['Ratio_velocity'][36], data_FBenoit['Ratio_Conc'][36],            
              marker = markerss[4], color = 'None', markersize = sizess[0], ls = '', 
              markeredgecolor = colorss[4], markeredgewidth = 1.5, zorder = 0, label = 'Stainer parallel')
p4, = ax.plot(data_FBenoit['Ratio_velocity'][42], data_FBenoit['Ratio_Conc'][42],            
              marker = markerss[7], color = 'None', markersize = sizess[0], ls = '', 
              markeredgecolor = colorss[7], markeredgewidth = 1.5, zorder = 0, label = 'Stainer perpendicular')

p5, = ax.plot(data_FBenoit['Ratio_velocity'][0], data_FBenoit['Ratio_Conc'][0],            
              marker = markerss[0], color = 'None', markersize = sizess[1], ls = '', 
              markeredgecolor = colorss[0], markeredgewidth = 1.5, zorder = 0, label = 'Nozzle parallel')
p6, = ax.plot(data_FBenoit['Ratio_velocity'][20], data_FBenoit['Ratio_Conc'][20],            
              marker = markerss[3], color = 'None', markersize = sizess[3], ls = '', 
              markeredgecolor = colorss[3], markeredgewidth = 1.5, zorder = 0, label = 'Nozzle perpendicular')
p7, = ax.plot(data_FBenoit['Ratio_velocity'][39], data_FBenoit['Ratio_Conc'][39],            
              marker = markerss[4], color = 'None', markersize = sizess[5], ls = '', 
              markeredgecolor = colorss[4], markeredgewidth = 1.5, zorder = 0, label = 'Stainer parallel')
p8, = ax.plot(data_FBenoit['Ratio_velocity'][44], data_FBenoit['Ratio_Conc'][44],            
              marker = markerss[7], color = 'None', markersize = sizess[5], ls = '', 
              markeredgecolor = colorss[7], markeredgewidth = 1.5, zorder = 0, label = 'Stainer perpendicular')
p9, = ax.plot(data_FBenoit['Ratio_velocity'][40], data_FBenoit['Ratio_Conc'][40],            
              marker = markerss[4], color = 'None', markersize = sizess[1], ls = '', 
              markeredgecolor = colorss[4], markeredgewidth = 1.5, zorder = 0, label = 'Stainer parallel')


ax.hlines(1,0,10, color = 'black')
ax.vlines(1,0,10, color = 'grey', ls = '--')
ax.axhspan(0.8, 1.2, 0, 2.5, color = 'grey', alpha = 0.2, zorder = 0)

ax.text(0.12, 0.025, r'sub-isokinetic', fontsize = 14, transform = ax.transAxes)
ax.text(0.4, 0.025, r'super-isokinetic', fontsize = 14, transform = ax.transAxes)
ax.text(0.32, 0.025, r'isokinetic', fontsize = 14, transform = ax.transAxes, rotation='vertical')
# ax.text(0.8, 0.05, r'$\mathregular{u_{flow} = 1.25 m/s}$', fontsize = 14, transform = ax.transAxes)
# ax.set_ylabel(r'$\mathregular{\frac{C_{sand, sampling}}{C_{sand, ref}}}$ (-)', fontsize=18, weight = 'bold')
ax.set_ylabel(r'Concentration ratio $\mathregular{\frac{C_{sand, sampling}}{C_{sand, ref}}}$ (-)', fontsize=18, weight = 'bold')
ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
ax.set_xlim(0.3,2.3)
ax.set_ylim(0.2,2.1) 
ax.legend(handles=[p1, p2, p3, p4, (p1, p2, p3, p4), (p7, p8), (p5, p6,p9)], labels=['Nozzle parallel', 'Nozzle perpendicular', 'Stainer parallel', 'Stainer perpendicular',
                                                      '6 mm', '10 mm', '14 mm'], 
          handler_map={tuple: HandlerTuple(ndivide=None)}, loc = 'upper right', fontsize = 14)
# ax.set_xlabel(r'$\mathregular{\frac{u_{sampling}}{u_{flow}}}$ (-)', fontsize=18, weight = 'bold')
ax.set_xlabel(r'Velocity ratio $\mathregular{\frac{u_{sampling}}{u_{flow}}}$ (-)', fontsize=18, weight = 'bold')

fig.tight_layout()
figname = 'Fig4'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')

#%%################################################################################################

# Figure 5

#################################################################################################


fig, ax = plt.subplots(1,1, figsize=(10, 8), dpi=200)

ax.plot(data_SPayen_field['Ratio_velocity_buse'], data_SPayen_field['Csand_sampling_inf400_buse_norm'],
            marker = 'o', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 1.5, zorder = 15, label = 'Nozzle parallel')

ax.plot(data_SPayen_field['Ratio_velocity_crepine'], data_SPayen_field['Csand_sampling_inf400_crepine_norm'],
            marker = 'D', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 1.5, zorder = 15, label = 'Stainer perpendicular')

ax.hlines(1,0,10, color = 'black')
ax.vlines(1,0,10, color = 'grey', ls = '--')
ax.axhspan(0.8, 1.2, 0, 2.5, color = 'grey', alpha = 0.2, zorder = 0)

ax.legend(fontsize =14, loc = 'lower right', framealpha = 1)
ax.text(0.1, 0.025, r'sub-isokinetic', fontsize = 14, transform = ax.transAxes)
ax.text(0.4, 0.025, r'super-isokinetic', fontsize = 14, transform = ax.transAxes)
ax.text(0.3, 0.025, r'isokinetic', fontsize = 14, transform = ax.transAxes, rotation='vertical')
ax.set_xlabel(r'Velocity ratio $\mathregular{\frac{u_{sampling}}{u_{flow}}}$ (-)', fontsize=18, weight = 'bold')
ax.set_ylabel(r'Concentration ratio $\mathregular{\frac{C_{sand, sampling}}{C_{sand, ref}}}$ (-)', fontsize=18, weight = 'bold')
ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
ax.set_xlim(0,3)
ax.set_ylim(0.25,1.5) 

fig.tight_layout()
figname = 'Fig5a'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')

#%%
fig, ax = plt.subplots(1,2, figsize=(14, 8), dpi=200)
# a) Nozzle
ax[0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Site']=='I_VR'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Site']=='I_VR'],            
            marker = 'D', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'darkcyan', markeredgewidth = 1.5, zorder = 15, label = 'I_VR')

ax[0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Site']=='A_C'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Site']=='A_C'],            
            marker = 's', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'blue', markeredgewidth = 1.5, zorder = 15, label = 'A_C')

ax[0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R1'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R1'],            
            marker = 'v', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'steelblue', markeredgewidth = 1.5, zorder = 15, label = 'R Serie 1')

ax[0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R2'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R2'],            
            marker = '^', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'darkblue', markeredgewidth = 1.5, zorder = 15, label = 'R Serie 2')

ax[0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R3'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R3'],            
            marker = '>', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'dodgerblue', markeredgewidth = 1.5, zorder = 15, label = 'R Serie 3')

ax[0].text(0.02, 0.95, 'a)', fontsize = 14, transform = ax[0].transAxes)
ax[0].vlines(1,0,10, color = 'grey', ls = '--')
ax[0].set_ylabel(r'$\mathregular{C_{sand}}$ (g/l)', fontsize=18, weight = 'bold')
ax[0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[0].set_xlim(0.3, 2.3)
ax[0].set_ylim(0, 1.8) 
# ax[0].legend(loc = 'upper right', fontsize = 14)

# b) Crepine
ax[1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Site']=='I_VR'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Site']=='I_VR'],            
            marker = 'D', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'darkcyan', markeredgewidth = 1.5, zorder = 15, label = 'I_VR')

ax[1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Site']=='A_C'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Site']=='A_C'],            
            marker = 's', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'blue', markeredgewidth = 1.5, zorder = 15, label = 'A_C')

ax[1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R1'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R1'],            
            marker = 'v', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'steelblue', markeredgewidth = 1.5, zorder = 15, label = 'R_PR Serie 1')

ax[1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R2'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R2'],            
            marker = '^', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'darkblue', markeredgewidth = 1.5, zorder = 15, label = 'R_PR Serie 2')

ax[1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R3'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R3'],            
            marker = '>', color = 'None', markersize = 7, ls = '', 
            markeredgecolor = 'dodgerblue', markeredgewidth = 1.5, zorder = 15, label = 'R_PR Serie 3')
        
ax[1].vlines(1,0,10, color = 'grey', ls = '--')

ax[1].text(0.02, 0.95, 'b)', fontsize = 14, transform = ax[1].transAxes)
# ax[1].set_ylabel(r'$\mathregular{\frac{C_{sand, stainer}}{C_{sand, nozzle}}}$ (-)', fontsize=18, weight = 'bold')
ax[1].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[1].set_xlim(0.3, 2.3)
ax[1].set_ylim(0.2, 1.8) 
ax[1].set_yticks([])
ax[1].legend(loc = 'upper right', fontsize = 14)
fig.supxlabel(r'$\mathregular{\frac{u_{sampling}}{u_{flow}}}$ (-)', fontsize=18, weight = 'bold')

fig.tight_layout()
figname = 'Pump_isocinetics_field'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')


#%% Plot field 

fig, ax = plt.subplots(3,2, figsize=(12, 16), dpi=200)
# a) I_VR 
ax[0,0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Site']=='I_VR'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Site']=='I_VR'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 2, zorder = 15, label = 'I_VR')

ax[0,0].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Site']=='I_VR'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Site']=='I_VR'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 2.5, zorder = 30, label = 'I_VR')
ax[0,0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Site']=='I_VR'], 0, 10, color = 'indigo', ls = '-', zorder = 10)
ax[0,0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Site']=='I_VR'], 0, 10, color = 'goldenrod', ls = '-', zorder = 10)
ax[0,0].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='I_VR'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='I_VR'].iloc[0]*1.2, 
                0, 2.5, color = 'indigo', alpha = 0.15, zorder = 0)  
ax[0,0].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='I_VR'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='I_VR'].iloc[0]*1.2, 
                0, 2.5, color = 'goldenrod', alpha = 0.15, zorder = 0)  

# b) A_C
ax[0,1].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Site']=='A_C'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Site']=='A_C'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 2, zorder = 15, label = 'A_C')

ax[0,1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Site']=='A_C'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Site']=='A_C'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 2.5, zorder = 30, label = 'A_C')
ax[0,1].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Site']=='A_C'], 0, 10, color = 'indigo', ls = '-', zorder = 10)
ax[0,1].hlines(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Site']=='A_C'], 0, 10, color = 'goldenrod', ls = '-', zorder = 10)
ax[0,1].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='A_C'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='A_C'].iloc[0]*1.2, 
                0, 2.5, color = 'indigo', alpha = 0.15, zorder = 0)  
ax[0,1].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='A_C'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='A_C'].iloc[0]*1.2, 
                0, 2.5, color = 'goldenrod', alpha = 0.15, zorder = 0)  
# c)
ax[1,0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R1'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R1'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 2, zorder = 15, label = 'R Serie 1')

ax[1,0].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R1'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R1'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 2.5, zorder = 30, label = 'R_PR Serie 1')
ax[1,0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R1'], 0, 10, color = 'indigo', ls = '-', zorder = 10)
ax[1,0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R1'], 0, 10, color = 'goldenrod', ls = '-', zorder = 10)
ax[1,0].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R1'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R1'].iloc[0]*1.2, 
                0, 2.5, color = 'indigo', alpha = 0.15, zorder = 0)  
ax[1,0].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R1'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R1'].iloc[0]*1.2, 
                0, 2.5, color = 'goldenrod', alpha = 0.15, zorder = 0)  

# d)
ax[1,1].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R2'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R2'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 2, zorder = 15, label = 'R Serie 2')

ax[1,1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R2'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R2'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 2.5, zorder = 30, label = 'R_PR Serie 2')
ax[1,1].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R2'], 0, 10, color = 'indigo', ls = '-', zorder = 10)
ax[1,1].hlines(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R2'], 0, 10, color = 'goldenrod', ls = '-', zorder = 10)
ax[1,1].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R2'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R2'].iloc[0]*1.2, 
                0, 2.5, color = 'indigo', alpha = 0.15, zorder = 0)  
ax[1,1].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R2'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R2'].iloc[0]*1.2, 
                0, 2.5, color = 'goldenrod', alpha = 0.15, zorder = 0)  
# e)
ax[2,0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R3'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R3'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 2, zorder = 15, label = 'Nozzle parallel')
ax[2,0].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R3'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R3'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 2.5, zorder = 30, label = 'Stainer perpendicular')
ax[2,0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R3'], 0, 10, color = 'indigo', ls = '-', zorder = 10) #, label = r'$\overline{C_{series, nozzle}}$')
ax[2,0].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R3'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R3'].iloc[0]*1.2, 
                0, 2.5, color = 'indigo', alpha = 0.15, zorder = 0)    
ax[2,0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R3'], 0, 10, color = 'goldenrod', ls = '-', zorder = 10) #, label = r'$\overline{C_{series, nozzle}}$')
ax[2,0].axhspan(data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R3'].iloc[0]*0.8, 
                data_SPayen_field['Mean_Csand_serie_inf400_crepine'][data_SPayen_field['Serie']=='R3'].iloc[0]*1.2, 
                0, 2.5, color = 'goldenrod', alpha = 0.15, zorder = 0)    

ax[0,0].text(0.02, 0.92, 'a)', fontsize = 20, transform = ax[0,0].transAxes)
ax[0,0].vlines(1,0,10, color = 'grey', ls = '--')
ax[0,0].tick_params(axis = 'both', which = 'major', labelsize = 20)
ax[0,0].set_xlim(0, 2.3)
ax[0,0].set_ylim(0, 0.45) 
ax[0,1].vlines(1,0,10, color = 'grey', ls = '--')
ax[0,1].text(0.02, 0.92, 'b)', fontsize = 20, transform = ax[0,1].transAxes)
ax[0,1].tick_params(axis = 'both', which = 'major', labelsize = 20)
ax[0,1].set_xlim(0, 2.3)
ax[0,1].set_ylim(0, 1)
ax[1,0].vlines(1,0,10, color = 'grey', ls = '--')
ax[1,0].text(0.02, 0.92, 'c)', fontsize = 20, transform = ax[1,0].transAxes)
ax[1,0].tick_params(axis = 'both', which = 'major', labelsize = 20)
ax[1,0].set_xlim(0, 2.3)
ax[1,0].set_ylim(0, 1.9) 
ax[1,1].vlines(1,0,10, color = 'grey', ls = '--')
ax[1,1].text(0.02, 0.92, 'd)', fontsize = 20, transform = ax[1,1].transAxes)
ax[1,1].tick_params(axis = 'both', which = 'major', labelsize = 20)
ax[1,1].set_xlim(0, 2.3)
ax[1,1].set_ylim(0, 1) 
ax[2,0].vlines(1,0,10, color = 'grey', ls = '--')
ax[2,0].text(0.02, 0.92, 'e)', fontsize = 20, transform = ax[2,0].transAxes)
ax[2,0].tick_params(axis = 'both', which = 'major', labelsize = 20)
ax[2,0].set_xlim(0, 2.3)
ax[2,0].set_ylim(0, 1.3) 

ax[2,0].text(0.05, 0.025, r'sub-isokinetic', fontsize = 20, transform = ax[2,0].transAxes)
ax[2,0].text(0.6, 0.025, r'super-isokinetic', fontsize = 20, transform = ax[2,0].transAxes)
ax[2,0].text(0.39, 0.025, r'isokinetic', fontsize = 20, transform = ax[2,0].transAxes, rotation='vertical')
ax[2,0].legend(loc = 'lower right', fontsize = 24, framealpha = 1) #, bbox_to_anchor = (0.5, -1), ncol = 2)
fig.supylabel(r'$\mathregular{C_{sand}}$ (g/l)', fontsize=24, weight = 'bold')
fig.supxlabel(r'Velocity ratio $\mathregular{u_{sampling}/u_{flow}}$ (-)', fontsize=24, weight = 'bold')

fig.tight_layout()
figname = 'Fig7'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
fig.savefig(outpath_figures+ '\\' + figname + '.pdf', dpi=300, bbox_inches='tight')
fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')
#################################################################################################

#%% Plot field 

fig, ax = plt.subplots(5,1, figsize=(6, 16), dpi=200)
# a) Serie 1 
ax[0].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Site']=='I_VR'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Site']=='I_VR'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 1.5, zorder = 15, label = 'I_VR')

ax[0].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Site']=='I_VR'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Site']=='I_VR'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 1.5, zorder = 15, label = 'I_VR')
ax[0].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Site']=='I_VR'], 0, 10, color = 'indigo', ls = '-')

# 2) 
ax[1].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Site']=='A_C'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Site']=='A_C'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 1.5, zorder = 15, label = 'A_C')

ax[1].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Site']=='A_C'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Site']=='A_C'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 1.5, zorder = 15, label = 'A_C')
ax[1].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Site']=='A_C'], 0, 10, color = 'indigo', ls = '-')

# 3)
ax[2].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R1'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R1'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 1.5, zorder = 15, label = 'R Serie 1')

ax[2].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R1'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R1'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 1.5, zorder = 15, label = 'R_PR Serie 1')
ax[2].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R1'], 0, 10, color = 'indigo', ls = '-')

# 4)
ax[3].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R2'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R2'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 1.5, zorder = 15, label = 'R Serie 2')

ax[3].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R2'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R2'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 1.5, zorder = 15, label = 'R_PR Serie 2')
ax[3].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R2'], 0, 10, color = 'indigo', ls = '-')

# 5)
ax[4].plot(data_SPayen_field['Ratio_velocity_buse'][data_SPayen_field['Serie']=='R3'], 
        data_SPayen_field['Csand_sampling_inf400_buse'][data_SPayen_field['Serie']=='R3'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'indigo', markeredgewidth = 1.5, zorder = 15, label = 'Nozzle parallel')

ax[4].plot(data_SPayen_field['Ratio_velocity_crepine'][data_SPayen_field['Serie']=='R3'], 
        data_SPayen_field['Csand_sampling_inf400_crepine'][data_SPayen_field['Serie']=='R3'],            
            marker = 'D', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'goldenrod', markeredgewidth = 1.5, zorder = 15, label = 'Stainer perpendicular')
ax[4].hlines(data_SPayen_field['Mean_Csand_serie_inf400_buse'][data_SPayen_field['Serie']=='R3'], 0, 10, color = 'indigo', ls = '-') #, label = r'$\overline{C_{series, nozzle}}$')
       

ax[0].text(0.02, 0.9, 'a)', fontsize = 14, transform = ax[0].transAxes)
ax[0].vlines(1,0,10, color = 'grey', ls = '--')
ax[0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[0].set_xlim(0, 2.3)
ax[0].set_ylim(0, 0.5) 
ax[1].vlines(1,0,10, color = 'grey', ls = '--')
ax[1].text(0.02, 0.9, 'b)', fontsize = 14, transform = ax[1].transAxes)
ax[1].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[1].set_xlim(0, 2.3)
ax[1].set_ylim(0, 1) 
ax[2].vlines(1,0,10, color = 'grey', ls = '--')
ax[2].text(0.02, 0.9, 'c)', fontsize = 14, transform = ax[2].transAxes)
ax[2].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[2].set_xlim(0, 2.3)
ax[2].set_ylim(0, 1.9) 
ax[3].vlines(1,0,10, color = 'grey', ls = '--')
ax[3].text(0.02, 0.9, 'd)', fontsize = 14, transform = ax[3].transAxes)
ax[3].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[3].set_xlim(0, 2.3)
ax[3].set_ylim(0, 1) 
ax[4].vlines(1,0,10, color = 'grey', ls = '--')
ax[4].text(0.02, 0.9, 'e)', fontsize = 14, transform = ax[4].transAxes)
ax[4].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[4].set_xlim(0, 2.3)
ax[4].set_ylim(0, 1.3) 

ax[4].legend(loc = 'lower right', fontsize = 14, framealpha = 1) #, bbox_to_anchor = (0.5, -1), ncol = 2)
fig.supylabel(r'$\mathregular{C_{sand}}$ (g/l)', fontsize=18, weight = 'bold')
fig.supxlabel(r'$\mathregular{\frac{u_{sampling}}{u_{flow}}}$ (-)', fontsize=18, weight = 'bold')

fig.tight_layout()
figname = 'FigS1'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')
#################################################################################################
#%% Alpha graph

ISO_size_classes = pd.read_csv(r'C:\Users\jmarggra\Documents\INRAE\Grenoble_Campus\Autres\ISO_size_classes.csv', sep = ';')
ISO_size_classes = ISO_size_classes.iloc[:,1]
ISO_size_classes = [i for i in ISO_size_classes]

# For classified distribution
size_classes_inf = np.array(ISO_size_classes[0:-1])* 1e-6 / 2
size_classes_sup = np.array(ISO_size_classes[1:])* 1e-6 / 2
size_classes_center = [10**((np.log10(size_classes_inf[i]) + np.log10(size_classes_sup[i]))/2) 
                     for i in range(len(size_classes_inf))]
size_classes_center= np.array(size_classes_center)

colnamess1 = ['Date', 'Start_sampling', 'End_sampling', 'Sampler', 'Q_sampling_m3_s', 'Stage_sampling_m', 'spm_sampling_g_l', 
                                   'Sand_concentration_g_l', 'Sand_flux_kg_s', 'No_samples_used_for_ISO_gsd', 'D10_mum', 'D50_mum', 'D90_mum',
                                   'Fine_concentration_g_l', 'Fine_flux_kg_s','U_C', 'U_Q', 'U_F', 'Method']
colnamess = colnamess1 + ISO_size_classes

#%% BD
import pickle
path_folder = r'C:\Users\jmarggra\Documents\INRAE\Grenoble_Campus'
sampling_dates_BD = pd.read_csv(r'C:\Users\jmarggra\Documents\INRAE\Grenoble_Campus\Autres\Sampling_dates_BD.csv', sep = ';')
datee_BD = sampling_dates_BD['Sampling_dates_2']
samples_BD = pd.DataFrame(columns = colnamess1, index = datee_BD)
samples_BD_ISO = pd.DataFrame(columns = colnamess, index = datee_BD)
analysis_data_all = []
for i in range(len(datee_BD)): # sampling_dates   
    with open(path_folder + '\\' + str(datee_BD.iloc[i]) + '_Campus\\' + str(datee_BD.iloc[i]) + '_Analysis_solid_gauging_BD.txt', "rb") as fp:   
        export_data = pickle.load(fp)
        #ISO_gsd_all_BD.append(export_data.ISO_mean_gsd)
        analysis_data_all.append(export_data.analysis_data)
        start_time = np.min(export_data.analysis_data['Time']) 
        end_time = np.max(export_data.analysis_data['Time'])        
        Q = (np.nanmean(export_data.Q_sampling['Value']))
        stage = (np.nanmean(export_data.stage_sampling['Value']))
        spm = (np.nanmean(export_data.spm_sampling['Value']))
        if export_data.summary_SDC is not None:
            sand_C = export_data.summary_SDC['Conc_mean_SDC_g_l'].iloc[0]        
            sand_flux = export_data.summary_SDC['total_sand_flux_SDC_kg_s'].iloc[0]
            methodi = 'SDC'
        else:
            sand_C = export_data.summary_NN['NN_mean_sand_conc_g_l'].iloc[0]
            sand_flux = export_data.summary_NN['NN_total_sand_flux_kg_s'].iloc[0]
            methodi = 'NN'
        if export_data.uncertainty is not None:
            U_C = export_data.uncertainty['U\'_C'].iloc[0]
            U_Q = export_data.uncertainty['U\'_Q'].iloc[0]
            U_F = export_data.uncertainty['U\'F'].iloc[0]
        ISO_gsd = export_data.ISO_mean_gsd
        if ISO_gsd is not None:
            D10 = np.interp(10, ISO_gsd.iloc[0,:], ISO_size_classes)
            D50 = np.interp(50, ISO_gsd.iloc[0,:], ISO_size_classes)
            D90 = np.interp(90, ISO_gsd.iloc[0,:], ISO_size_classes)
            ISO_gsdd = export_data.ISO_mean_gsd.iloc[0]
        else: 
            D10 = None
            D50 = None
            D90 = None            
            ISO_gsdd = [None]*100
        number_samples_for_ISO = export_data.ISO_mean_gsd_number
        fine_C = np.nan
        fine_flux = np.nan                          
   
    samples_a = pd.DataFrame([datee_BD.iloc[i], start_time, end_time, 'BD', Q, stage, spm, sand_C, 
                                  sand_flux, number_samples_for_ISO, D10, D50, D90, fine_C, fine_flux, 
                                  U_C, U_Q, U_F, methodi]).transpose()   
    samples_b = pd.concat([samples_a, pd.DataFrame(ISO_gsdd).transpose()], axis = 1)
    samples_BD_ISO.iloc[i,:] = samples_b
    samples_BD.iloc[i,:] = samples_a
   

analysis_data_all = pd.concat(analysis_data_all)
analysis_data_all = analysis_data_all.dropna(subset = ['Small_nozzle'], how = 'any')
analysis_data_all.reset_index(drop = True, inplace = True)

#%% Prepare and load alpha
xy = np.round(np.arange(0.001, 4.001, 0.001), 3)
vel_possiblerange = np.round(np.arange(0.01, 4.01, 0.01), 3)
BD_grain_size_steps = [75, 90, 100, 110, 130, 150, 210]
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
from Sed_flux_small_straight import sed_flux_small_straight #4
from Sed_flux_small_bend import sed_flux_small_bend #3
from Sed_flux_big_straight import sed_flux_big_straight #2
from Sed_flux_big_bend import sed_flux_big_bend #1

# Create Alpha relations
sed_flux_small_straight(outpath_figures)
sed_flux_small_bend(outpath_figures)        
sed_flux_big_straight(outpath_figures)        
sed_flux_big_bend(outpath_figures)
alpha_small_straight = sed_flux_small_straight.alpha_small_straight
alpha_small_bend = sed_flux_small_bend.alpha_small_bend
alpha_big_straight = sed_flux_big_straight.alpha_big_straight
alpha_big_bend = sed_flux_big_bend.alpha_big_bend

colors = ['tab:blue','tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink']

#%% Determine alpha
analysis_data_all_GSD = analysis_data_all.dropna(subset = ['D50'])
analysis_data_all_GSD.reset_index(drop=True, inplace=True)

alpha_all = []
nozzle_all = []
vel_all = []
nearest_grain_size_all= []
vel_ind_all = []
for i in range(len(analysis_data_all_GSD)):
    meas = analysis_data_all_GSD.iloc[i,:]
    velocity_point_2digits = np.round(meas['Velocity_sampling_point'], 2)
    vel_all.append(velocity_point_2digits)
    # Find nearest grain size
    nearest_grain_size = BD_grain_size_steps.index(find_nearest(BD_grain_size_steps, meas['D50']))
    nearest_grain_size_all.append(nearest_grain_size)
    vel_ind = np.where(vel_possiblerange == velocity_point_2digits)
    vel_ind_all.append(vel_ind)
    #alpha = alpha_big_bend[nearest_grain_size][vel_ind]
       
    # Determine the used nozzle (small or big)
    if meas['Small_nozzle'] == 1:
        if meas['Height_field'] <= 0.2:
            alpha = alpha_small_bend[nearest_grain_size][vel_ind]
            nozzle_all.append(3)
            alpha_all.append(float(alpha))
        else:
            alpha = alpha_small_straight[nearest_grain_size][vel_ind]
            nozzle_all.append(4)
            alpha_all.append(float(alpha))
    if meas['Small_nozzle'] == 0:
        if meas['Height_field'] <= 0.2:
            alpha = alpha_big_bend[nearest_grain_size][vel_ind]
            nozzle_all.append(1)
            alpha_all.append(float(alpha))
        else:
            alpha = alpha_big_straight[nearest_grain_size][vel_ind]
            nozzle_all.append(2)
            alpha_all.append(float(alpha))

dataa = pd.DataFrame([nozzle_all, alpha_all, vel_all, analysis_data_all_GSD['D50'], vel_ind_all, nearest_grain_size_all, 
                      analysis_data_all_GSD['Concentration_sand_g_l'], analysis_data_all_GSD['Sand_flux_point_kg_s_m2']]).transpose()
dataa.columns = ['Nozzle', 'alpha', 'Velocity_sampling_point', 'D50', 'Vel_ind', 'nearest_grain_size', 
                 'Concentration_sand_g_l', 'Sand_flux_point_kg_s_m2']
dataa['Flux_before_correction'] = dataa['Sand_flux_point_kg_s_m2']/dataa['alpha']
dataa['Concentration_without_corr'] = dataa['Flux_before_correction']/dataa['Velocity_sampling_point']
dataa['Concentration_corr'] = dataa['Concentration_sand_g_l']*dataa['alpha']

#%% Add our data to alpha graphs
linestyles = [(0, (1, 5)), '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 10)), '-']

fig, ax = plt.subplots(2, 2, figsize=(12, 8), dpi=200)

# Big bend
for i in range(len(alpha_big_bend)):
    ax[0,0].plot(vel_possiblerange, alpha_big_bend[i], lw =1.5, linestyle=linestyles[i], 
                 color=colors[i], label= str(BD_grain_size_steps[i]) + ' $\mathregular{\mu}$m')

ax[0,0].xaxis.set_ticklabels([])
ax[0,0].set_xlim(0.5, 2.3)
ax[0,0].set_ylim(0.6, 2.5)
ax[0,0].text(0.05, 0.9, 'a)', fontsize = 14, transform = ax[0,0].transAxes)

# Big straight
for i in range(len(alpha_big_straight)):
    ax[0,1].plot(vel_possiblerange, alpha_big_straight[i], lw =1.5, linestyle=linestyles[i], 
                 color=colors[i],)

ax[0,1].xaxis.set_ticklabels([])
ax[0,1].yaxis.tick_right()
ax[0,1].yaxis.set_label_position("right")
ax[0,1].set_xlim(0.5, 2.3)
ax[0,1].set_ylim(0.6, 2.5)
ax[0,1].text(0.05, 0.9, 'b)', fontsize = 14, transform = ax[0,1].transAxes)

# small bend
for i in range(len(alpha_small_bend)):
    ax[1,0].plot(vel_possiblerange, alpha_small_bend[i], lw =1.5, linestyle=linestyles[i], 
                 color=colors[i])
ax[1,0].set_xlim(0.5, 2.3)
ax[1,0].set_ylim(0.6, 2.5)
ax[1,0].text(0.05, 0.9, 'c)', fontsize = 14, transform = ax[1,0].transAxes)

# small straight
for i in range(len(alpha_small_straight)):
    ax[1,1].plot(vel_possiblerange, alpha_small_straight[i], lw =1.5, linestyle=linestyles[i], 
                 color=colors[i])
ax[1,1].yaxis.tick_right()
ax[1,1].yaxis.set_label_position("right")
ax[1,1].set_xlim(0.5, 2.3)
ax[1,1].set_ylim(0.6, 2.5)
ax[1,1].text(0.05, 0.9, 'd)', fontsize = 14, transform = ax[1,1].transAxes)

# add our data
for i in range(len(analysis_data_all_GSD)):
    if nozzle_all[i] == 1:
        ax[0,0].plot(vel_all[i], alpha_all[i],
                     marker = 'D', color = 'black', alpha = 0.5, markersize = 6, markeredgewidth = 0)
    if nozzle_all[i] == 2:
        ax[0,1].plot(vel_all[i], alpha_all[i],
                     marker = 'D', color = 'black', alpha = 0.5, markersize = 6, markeredgewidth = 0)
    if nozzle_all[i] == 3:
        ax[1,0].plot(vel_all[i], alpha_all[i],
                     marker = 'D', color = 'black', alpha = 0.5, markersize = 6, markeredgewidth = 0)
    if nozzle_all[i] == 4:
        ax[1,1].plot(vel_all[i], alpha_all[i],
                     marker = 'D', color = 'black', alpha = 0.5, markersize = 6, markeredgewidth = 0)


ax[0,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[0,1].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[1,0].tick_params(axis = 'both', which = 'major', labelsize = 14)
ax[1,1].tick_params(axis = 'both', which = 'major', labelsize = 14)

fig.supxlabel('Velocity (m/s)', fontsize = 18, weight = 'bold')
fig.supylabel(r'$\mathregular{\alpha}$ (-)', fontsize = 18, weight = 'bold', x=0.01, y=0.5)
fig.legend(loc = 'lower center', fontsize = 14, bbox_to_anchor = (0.5, -0.12), ncol = 4)

fig.tight_layout()
figname = 'Alpha_graph_our_data'
plt.savefig(outpath_figures + '\\' + figname + '.png', dpi=300, bbox_inches='tight')
plt.savefig(outpath_figures + '\\' + figname + '.eps', dpi=300, bbox_inches='tight')


#%% Plot histogram alpha

plt.hist(alpha_all, bins = 40)

alpha_gre_1_2 = sum(i > 1.2 for i in alpha_all)
alpha_gre_1_2_perc = alpha_gre_1_2/len(alpha_all)

alpha_low_0_8 = sum(i < 0.8 for i in alpha_all)
alpha_low_0_8_perc = alpha_low_0_8/len(alpha_all)

#%% 

plt.plot(dataa['Concentration_sand_g_l'], dataa['Concentration_without_corr'], '.')
plt.plot(dataa['Concentration_sand_g_l'], dataa['Concentration_sand_g_l']/dataa['Concentration_corr'], '.')


#%% Compare with literature data

data_Dijkman_Milisic82 = data_Dijkman_Milisic82.dropna(subset = ['USP'], how = 'any')


# Plot
fig, ax = plt.subplots(1,1, figsize=(10, 10), dpi=200)


p1 = ax.errorbar(data_P6['Sand_flux_kg_s']/area_GC['Area_m2'], data_BD['Sand_flux_kg_s']/area_GC['Area_m2'],          
        marker = 'D', color = 'orange', markersize = 8, ls = '', 
        markeredgecolor = 'black', markeredgewidth = 0.4,
        yerr = err_F_BD/area_GC['Area_m2'], xerr = err_F_P6/area_GC['Area_m2'], elinewidth = 1, 
        capsize = 1.5, zorder = 50, label = 'this study > 62 $\mathregular{\mu}$m')


# in kg/m2s
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='3_10'], 
        data_Dijkman_Milisic82['DBS1'][data_Dijkman_Milisic82['Fraction']=='3_10'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'salmon', markeredgewidth = 1, zorder = 15, label = 'Dijkman & Milisic (1982) 72 - 145 $\mathregular{\mu}$m')
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='3_10'], 
        data_Dijkman_Milisic82['DBB1'][data_Dijkman_Milisic82['Fraction']=='3_10'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'salmon', markeredgewidth = 1, zorder = 15)
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='3_10'], 
        data_Dijkman_Milisic82['DBS2'][data_Dijkman_Milisic82['Fraction']=='3_10'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'salmon', markeredgewidth = 1, zorder = 15)
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='3_10'], 
        data_Dijkman_Milisic82['DBB2'][data_Dijkman_Milisic82['Fraction']=='3_10'],            
            marker = 'o', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'salmon', markeredgewidth = 1, zorder = 15)

ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='10_20'], 
        data_Dijkman_Milisic82['DBS1'][data_Dijkman_Milisic82['Fraction']=='10_20'],            
            marker = 's', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'darkred', markeredgewidth = 1, zorder = 15, label = 'Dijkman & Milisic (1982) 145 - 230 $\mathregular{\mu}$m')
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='10_20'], 
        data_Dijkman_Milisic82['DBB1'][data_Dijkman_Milisic82['Fraction']=='10_20'],            
            marker = 's', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'darkred', markeredgewidth = 1, zorder = 15)
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='10_20'], 
        data_Dijkman_Milisic82['DBS2'][data_Dijkman_Milisic82['Fraction']=='10_20'],            
            marker = 's', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'darkred', markeredgewidth = 1, zorder = 15)
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='10_20'], 
        data_Dijkman_Milisic82['DBB2'][data_Dijkman_Milisic82['Fraction']=='10_20'],            
            marker = 's', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'darkred', markeredgewidth = 1, zorder = 15)

ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='larger20'], 
        data_Dijkman_Milisic82['DBS1'][data_Dijkman_Milisic82['Fraction']=='larger20'],            
            marker = '^', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'tan', markeredgewidth = 1, zorder = 15, label = r'Dijkman & Milisic (1982) $>$230 $\mathregular{\mu}$m')
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='larger20'], 
        data_Dijkman_Milisic82['DBB1'][data_Dijkman_Milisic82['Fraction']=='larger20'],            
            marker = '^', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'tan', markeredgewidth = 1, zorder = 15)
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='larger20'], 
        data_Dijkman_Milisic82['DBS2'][data_Dijkman_Milisic82['Fraction']=='larger20'],            
            marker = '^', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'tan', markeredgewidth = 1, zorder = 15)
ax.plot(data_Dijkman_Milisic82['USP'][data_Dijkman_Milisic82['Fraction']=='larger20'], 
        data_Dijkman_Milisic82['DBB2'][data_Dijkman_Milisic82['Fraction']=='larger20'],            
            marker = '^', color = 'None', markersize = 8, ls = '', 
            markeredgecolor = 'tan', markeredgewidth = 1, zorder = 15)

ax.plot(data_BW89_all['Sand_flux_P61_g_m2s']/1000, data_BW89_all['Sand_flux_DB_g_m2s']/1000,            
            marker = 'P', color = 'teal', markersize = 10, ls = '', 
            markeredgecolor = 'black', markeredgewidth = 0.5, zorder = 20 ,label = 'Beverage & Williams (1989) > 62 $\mathregular{\mu}$m')
ax.plot(data_BW89_62_125['Sand_flux_P61_g_m2s']/1000, data_BW89_62_125['Sand_flux_DB_g_m2s']/1000,            
            marker = 'X', color = 'dodgerblue', markersize = 10, ls = '', 
            markeredgecolor = 'black', markeredgewidth = 0.5, zorder = 20 ,label = 'Beverage & Williams (1989) 62 - 125 $\mathregular{\mu}$m')
ax.plot(data_BW89_125_250['Sand_flux_P61_g_m2s']/1000, data_BW89_125_250['Sand_flux_DB_g_m2s']/1000,            
            marker = '*', color = 'skyblue', markersize = 12, ls = '', 
            markeredgecolor = 'black', markeredgewidth = 0.5, zorder = 20 ,label = 'Beverage & Williams (1989) 125 - 250 $\mathregular{\mu}$m')

ax.plot(2*x_range, x_range,  
        ls = ':', lw = 1, color = 'black')
ax.plot(x_range, x_range,  
        ls = '-', lw = 1, color = 'black')
ax.plot(x_range, 2*x_range,  
        ls = ':', lw = 1, color = 'black',)
ax.text(2,1.6,'1:1', fontsize = 16)
ax.text(0.3,1.8,'Error of a \nfactor of 2', fontsize = 16)

# handles, labels = ax.get_legend_handles_labels()
# order = [0,1,2,3,4,5,8,6,7]
# legend = ax.legend([handles[i] for i in order], [labels[i] for i in order],
#                    framealpha = 1, facecolor = 'white', fontsize =14, loc = 'upper left')
legend = ax.legend( framealpha = 1, facecolor = 'white', fontsize =14, loc = 'upper left') #, bbox_to_anchor = (1.01, 0.9), ncol = 1)
legend.get_frame().set_alpha(1.0)
legend.get_frame().set_facecolor('white')
ax.set_xlabel(r'$\mathregular{\Phi_{sand, P6/P61} \; (kg/m^2s)}$', fontsize=20, weight = 'bold')
ax.set_ylabel(r'$\mathregular{\Phi_{sand, DB} \;(kg/m^2s)}$', fontsize=20, weight = 'bold')
ax.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.001, 3)
ax.set_ylim(0.001,3) 

fig.tight_layout()
figname = 'Fig_comp_literature'
fig.savefig(outpath_figures+ '\\' + figname + '.png', dpi=300, bbox_inches='tight')
# fig.savefig(outpath_figures+ '\\' + figname + '.eps', dpi=300, bbox_inches='tight')































