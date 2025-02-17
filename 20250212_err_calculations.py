# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 11:38:31 2024

@author: Jaka Bizjak, IJS CEU
Predloga za grafični prikaz rezultatov in meritev pri projektu Care4Climate

Predloga podaja usmeritve glede grafičnega prikaza meritev in rezultatov. Pri 
tem se uporabi prosto dostopen repozitorij za oblikovanje grafov za objavo v 
znanstvenih revijah, dostopno na:
    https://github.com/garrettj403/SciencePlots
    https://doi.org/10.5281/zenodo.10206719
Uporabi se stil 'science' in barvno shemo 'muted', pri čemer se barvo #DDCC77 
v shemi zamenja z #DDAA33. Grafi imajo gridline ('grid'). Velikost pisave je 
10 pt, fizikalne veličine se pišejo poševno (med $$), enote pokončno znotraj 
oglatih oklepajev ([]). V primeru točkovnega grafa se uporabi 'scatter' (primer 
2. graf). Decimalno ločilo je , in ne . Območje na x in y osi se prilagaja 
prikazanim podatko. Legenda se nahaja v levem zgornjem kotu, razen če to 
moti prikaz podatkov.
Velikost slike se prilagaja. Širina slike, ki je v dokumentu če celotno širino,
je 7 inchev (privzeto naj bo 7x4). Priporočljivo je uporabiti razmerje širina:višina 4:3 ali 16:9, 
širino slike pa med 4 in 7 inchi. V primeru dveh slik v eni vrsti, se uporabi
širina 3,5 incha.
Za prikaz istega tipa podatkov, naj se pri vseh grafičnih prikazih uporabi iste 
barve. Npr. zemeljski plin je vedno prikazan z rumeno, elektrika z zeleno itd. 
Celoten seznam barv in pripadajočih veličin:
    List of colors
    i     code    colour  what do we use it for
    0     #CC6677 red     stroški
    1     #332288 blue    seštevki spremenljivk
    2     #DDAA33 yellow  zemeljski plin
    3     #117733 green   elektrika
    4     #88CCEE cyan    TGP
    5     #882255 wine
    6     #44AA99 teal    temperatura
    7     #999933 olive   energija
    8     #AA4499 purple
    9     #DDDDDD grey
V nadaljevanju sta prikazana dva primer, ki sta bila uporabljena v PN7_1 in
sicer za črtni in točkovni graf.

"""
#pip install SciencePlots
import numpy as np
import pandas as pd
from pandas import DataFrame
import locale
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import scienceplots
import datetime
import os
import seaborn as sns
import matplotlib.colors as mcolors

# oblikovanje grafa
plt.style.use(['science', 'muted', 'grid']) # glej https://github.com/garrettj403/SciencePlots
plt.rcParams.update({
    'font.size': 10,        # Font size for text
    'axes.titlesize': 10,   # Font size for axes titles
    'axes.labelsize': 10,   # Font size for x and y labels
    'xtick.labelsize': 10,  # Font size for x-axis tick labels
    'ytick.labelsize': 10,  # Font size for y-axis tick labels
    'legend.fontsize': 10,  # Font size for legend
    'lines.linewidth': 2,    # Default line width
})
colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] # list of default colours
colors[2] = '#DDAA33'# it changes #DDCC77
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# List of colors
# i     code    colour  what do we use it for
# 0     #CC6677 red     stroški
# 1     #332288 blue    seštevki spremenljivk
# 2     #DDAA33 yellow  zemeljski plin
# 3     #117733 green   elektrika
# 4     #88CCEE cyan    TGP
# 5     #882255 wine
# 6     #44AA99 teal    temperatura
# 7     #999933 olive   energija
# 8     #AA4499 purple
# 9     #DDDDDD grey
locale.setlocale(locale.LC_NUMERIC, "de_DE.UTF-8") # zamenjava . z , v grafih
plt.rcParams['axes.formatter.use_locale'] = True

os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "0" # tole je pomembno, ker matplotlib prireja velikost slike gleden na ekran!!!
os.environ["MATPLOTLIB_FORCE_DPI"] = "300"

start = datetime.datetime.now() # stopwatch the time of loop calculation

# Load the Excel file
str_excel_1 = 'MN_model_results.xlsx'
str_excel_2 = 'MN_model_results_n_C_r.xlsx'
 
take_from = np.arange(0.1, 0.51, 0.01)
take_from_str = ["%.2f" % number for number in take_from]

frac_O2 = {'phi_O2': np.arange(21, 96, 1)}
df_CO_err, df_H2_err, df_CH4_err, df_CO2_err, df_N2_err, df_T_err, df_CGE_err, df_CGE_tot_err, df_GY_err, df_LHV_err = [pd.DataFrame(frac_O2) for i in range(10)]

# empty dfs
df_phi_21 = pd.DataFrame() # NS model, phi_O = 21 %
df_phi_60 = pd.DataFrame() # NS model, phi_O = 60 %

# dfs of results for S model, phi_O = 21 % and 60 % 
df_4 = pd.read_excel(str_excel_1, sheet_name='stoichiometric_model_ERs')
df_5 = pd.read_excel(str_excel_1, sheet_name='stoichiometric_model_ERs_phi_60')

for value in take_from_str:
    str_sheet = str(value)
    df_1 = pd.read_excel(str_excel_1, sheet_name=str_sheet)
    df_2 = pd.read_excel(str_excel_2, sheet_name=str_sheet)
    df_phi_21 = pd.concat([df_phi_21, df_1[df_1['phi_O'] == 21]], ignore_index=True)
    df_phi_60 = pd.concat([df_phi_60, df_1[df_1['phi_O'] == 60]], ignore_index=True)
    df_CO_err[str(value)] = 1 - df_1['CO']/df_2['CO']
    df_H2_err[str(value)] = 1 - df_1['H2']/df_2['H2']
    df_CH4_err[str(value)] = 1 - df_1['CH4']/df_2['CH4']
    df_CO2_err[str(value)] = 1 - df_1['CO2']/df_2['CO2']
    df_N2_err[str(value)] = 1 - df_1['N2']/df_2['N2']
    df_T_err[str(value)] = 1 - df_1['T [K]']/df_2['T [K]']
    df_CGE_err[str(value)] = 1 - df_1['CGE']/df_2['CGE']
    df_CGE_tot_err[str(value)] = 1 - df_1['CGE_tot']/df_2['CGE_tot']
    df_GY_err[str(value)] = 1 - df_1['GY']/df_2['GY']
    df_LHV_err[str(value)] = 1 - df_1['LHV_g [MJ/kg]']/df_2['LHV_g [MJ/kg]']


''' 
# Sample data for the loop (replace with actual dfs)
dfs = ([df_CO_err, df_H2_err, df_CH4_err, df_CO2_err, df_N2_err, df_T_err, df_CGE_err, df_CGE_tot_err, df_GY_err, df_LHV_err])  # Your dataframes
save_as = (['CO', 'H$_2$', 'CH$_4$', 'CO$_2$', 'N$_2$', 'T', 'CGE', 'CGE_tot', 'GY', 'LHV'])  # Labels

# make every fig separately
for i, df in enumerate(dfs):
    array = dfs[i].to_numpy()[:, 1:]
    fig1, ax1 = plt.subplots(figsize=(3.5, 3.5))
    contour = ax1.contourf(np.arange(0.21,0.96,0.01), take_from, np.transpose(array*100), levels=30, cmap='bwr', norm=mcolors.TwoSlopeNorm(vmin=-20, vmax=20, vcenter=0))  # Contour fill plot
    fig1.colorbar(contour, ax=ax1)  # Add colorbar
    # Add labels and title
    ax1.set_xlabel('$\\varphi_{O_2}$ [\%]')
    ax1.set_ylabel('$\lambda$ [/]')
    if save_as[i] == 'T':
        ax1.set_title('relativno odstopanje $T$ [\%]')
    elif save_as[i] == 'CGE':
        ax1.set_title('relativno odstopanje $CGE$ [\%]')
    elif save_as[i] == 'CGE_tot':
        ax1.set_title('relativno odstopanje $CGE_{tot}$ [\%]')
    elif save_as[i] == 'GY':
        ax1.set_title('relativno odstopanje plina $GY$ [\%]')
    elif save_as[i] == 'LHV':
        ax1.set_title('relativno odstopanje $LHV_{pp}$ [\%]')
    else:
        ax1.set_title('relativno odstopanje volumskega deleža ' + save_as[i] + ' [\%]')
    fig1.tight_layout()
    save = save_as[i].replace('$', '')
    save = save.replace('_', '')
    fig1.savefig('results_err_' + save +  '.png', format='png', dpi=300, bbox_inches='tight')
    plt.close('all')

# make all six figs as subfigs in one fig
fig2, ax2 = plt.subplots(3, 2, figsize=(7, 9.89))  # 2 rows, 3 columns

# Loop through the data and plot each one in a subplot
for i, df in enumerate(dfs[:-4]):
    # Convert dataframe to numpy array and process it
    array = df.to_numpy()[:, 1:]
    # Determine the row and column for the subplot (2 rows, 3 columns)
    row = i // 2  # Integer division for row index (0 or 1)
    col = i % 2   # Modulo for column index (0, 1, or 2)
    # Create contour plot in the correct subplot (ax2[row, col])
    contour = ax2[row, col].contourf(np.arange(0.21, 0.96, 0.01)*100, take_from, np.transpose(array*100), levels=30, cmap='bwr', norm=mcolors.TwoSlopeNorm(vmin=-20, vmax=20, vcenter=0))  # Contour plot
    fig2.colorbar(contour, ax=ax2[row, col])  # Add colorbar
    ax2[row, col].set_xlabel('$\\varphi_{O_2}$ [\%]')  # X-axis label
    ax2[row, col].set_ylabel('$\lambda$ [/]')  # Y-axis label
    # Set the title based on the 'save_as' list
    if save_as[i] == 'T':
        ax2[row, col].set_title('relativno odstopanje $T$ [\%]')  # Title for temperature
    else:
        ax2[row, col].set_title(f'relativno odstopanje volumskega deleža {save_as[i]} [\%]')  # Title for other variables

fig2.tight_layout()  # Adjust layout for better spacing
fig2.savefig('results_err_phi_i.png', format='png', dpi=300, bbox_inches='tight')  # Save the figure as a PNG

plt.close('all')  # Close the plot after saving and showing

# make all four figs as subfigs in one fig
fig5, ax6 = plt.subplots(2, 2, figsize=(7, 6.6))  # 2 rows, 3 columns
# Loop through the data and plot each one in a subplot
for i, df in enumerate(dfs[-4:]):
    # Convert dataframe to numpy array and process it
    array = df.to_numpy()[:, 1:]
    array = np.array(array, dtype=np.float64)  # Convert to float
    # Determine the row and column for the subplot (2 rows, 3 columns)
    row = i // 2  # Integer division for row index (0 or 1)
    col = i % 2   # Modulo for column index (0, 1, or 2)
    i = i + 6
    # Create contour plot in the correct subplot (ax6[row, col])
    contour = ax6[row, col].contourf(np.arange(0.21, 0.96, 0.01)*100, take_from, np.transpose(array*100), levels=30, cmap='bwr', norm=mcolors.TwoSlopeNorm(vmin=-20, vmax=20, vcenter=0))  # Contour plot
    fig5.colorbar(contour, ax=ax6[row, col])  # Add colorbar
    ax6[row, col].set_xlabel('$\\varphi_{O_2}$ [\%]')  # X-axis label
    ax6[row, col].set_ylabel('$\lambda$ [/]')  # Y-axis label
    # Set the title based on the 'save_as' list
    if save_as[i] == 'CGE':
        ax6[row, col].set_title('relativno odstopanje $CGE$ [\%]')
    elif save_as[i] == 'CGE_tot':
        ax6[row, col].set_title('relativno odstopanje $CGE_{tot}$ [\%]')
    elif save_as[i] == 'GY':
        ax6[row, col].set_title('relativno odstopanje $GY$ [\%]')
    elif save_as[i] == 'LHV':
        ax6[row, col].set_title('relativno odstopanje $LHV_{pp}$ [\%]')
fig5.tight_layout()  # Adjust layout for better spacing
fig5.savefig('results_err_phi_i_2.png', format='png', dpi=300, bbox_inches='tight')  # Save the figure as a PNG
'''
# calculation of RMS betwen NS and S models for phi = 21 % and 60 % 
RMS_21_res = []
RMS_60_res = []
RMS_21_ER_02_res = []
RMS_60_ER_02_res = []
x = ['CO', 'H2', 'CH4', 'CO2', 'N2', 'T [K]'] # variable to consider
for species in x:
    RMS_21 = np.sqrt(np.sum((df_phi_21[species].to_numpy() - df_4[species].to_numpy())**2)/np.size(np.arange(0.1, 0.51, 0.01)))
    RMS_60 = np.sqrt(np.sum((df_phi_60[species].to_numpy() - df_5[species].to_numpy())**2)/np.size(np.arange(0.1, 0.51, 0.01)))
    filter_1 = df_phi_21.loc[df_phi_21["ER"] >= 0.2, species].to_numpy()
    filter_2 = df_4.loc[df_4["ER"] >= 0.2, species].to_numpy()
    filter_3 = df_phi_60.loc[df_phi_60["ER"] >= 0.2, species].to_numpy()
    filter_4 = df_5.loc[df_5["ER"] >= 0.2, species].to_numpy()
    RMS_21_ER_02 = np.sqrt(np.sum((filter_1 - filter_2)**2)/np.size(np.arange(0.2, 0.51, 0.01)))
    RMS_60_ER_02 = np.sqrt(np.sum((filter_3 - filter_4)**2)/np.size(np.arange(0.2, 0.51, 0.01)))

    RMS_21_res.append(RMS_21)  # Append results to the list
    RMS_60_res.append(RMS_60)  # Append results to the list
    RMS_21_ER_02_res.append(RMS_21_ER_02)  # Append results to the list
    RMS_60_ER_02_res.append(RMS_60_ER_02)  # Append results to the list
    
# Create a DataFrame with species as the index
df_RMS = pd.DataFrame({
    'Species': x,
    'RMS_21': RMS_21_res,
    'RMS_60': RMS_60_res,    
    'RMS_21 ER > 0,2': RMS_21_ER_02_res,
    'RMS_60 ER > 0,2': RMS_60_ER_02_res
})

# Display the DataFrame
print(df_RMS)
    
plt.close('all')  # Close the plot after saving and showing

end = datetime.datetime.now()
print('\nCalculation time:', end - start) 