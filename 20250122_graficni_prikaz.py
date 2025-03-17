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
str_excel_3 = 'MN_model_results_contour_data.xlsx' 

 
take_from = np.arange(0.1, 0.51, 0.01)
take_from_str = ["%.2f" % number for number in take_from]

frac_O2 = {'phi_O2': np.arange(21, 96, 1)}
df_CO, df_H2, df_CH4, df_CO2, df_N2, df_T, df_CGE, df_CGE_tot, df_GY, df_LHV = [pd.DataFrame(frac_O2) for i in range(10)]

df_phi_21 = pd.DataFrame() # NS model, phi_O = 21 %
df_phi_60 = pd.DataFrame() # NS model, phi_O = 60 %

for value in take_from_str:
    str_sheet = str(value)
    #print(str_sheet)
    df_1 = pd.read_excel(str_excel_1, sheet_name=str_sheet)
    df_phi_21 = pd.concat([df_phi_21, df_1[df_1['phi_O'] == 21]], ignore_index=True)
    df_phi_60 = pd.concat([df_phi_60, df_1[df_1['phi_O'] == 60]], ignore_index=True)
    #print(df_phi_21)
    df_CO[str(value)] = df_1['CO']
    df_H2[str(value)] = df_1['H2']
    df_CH4[str(value)] = df_1['CH4']
    df_CO2[str(value)] = df_1['CO2']
    df_N2[str(value)] = df_1['N2']
    df_T[str(value)] = df_1['T [K]']
    df_CGE[str(value)] = df_1['CGE']
    df_CGE_tot[str(value)] = df_1['CGE_tot']
    df_GY[str(value)] = df_1['GY']
    df_LHV[str(value)] = df_1['LHV_g [MJ/kg]']
    

# Sample data for the loop (replace with actual dfs)
dfs = ([df_CO, df_H2, df_CH4, df_CO2, df_N2, df_T, df_CGE, df_CGE_tot, df_GY, df_LHV])  # Your dataframes
save_as = (['CO', 'H$_2$', 'CH$_4$', 'CO$_2$', 'N$_2$', 'T', 'CGE', 'CGE_tot', 'GY', 'LHV'])  # Labels

df_2 = pd.read_excel(str_excel_1, sheet_name='reaktoro_NS_model')
df_3 = pd.read_excel(str_excel_1, sheet_name='stoichiometric_model')
df_4 = pd.read_excel(str_excel_1, sheet_name='stoichiometric_model_ERs')
df_5 = pd.read_excel(str_excel_1, sheet_name='stoichiometric_model_ERs_phi_60')

# make table of max and min values of i and locations of max and min values
data = []
for i, df in enumerate(dfs):
    max_value = df.iloc[:, 1:].max().max()
    max_loc = df.iloc[:, 1:].stack().idxmax()
    min_value = df.iloc[:, 1:].min().min()
    min_loc = df.iloc[:, 1:].stack().idxmin()
    data.append({
            "": save_as[i].replace('$', ''),   
            "Max Value": max_value,    
            "Max O2": max_loc[0] + 21, 
            "Max ER": max_loc[-1],     
            "Min Value": min_value,    
            "Min O2": min_loc[0] + 21, 
            "Min ER": min_loc[-1]      
        })
# Convert list of dictionaries to a DataFrame
min_max = pd.DataFrame(data)
# Display the DataFrame
print(min_max)

with pd.ExcelWriter(str_excel_1, engine='openpyxl', mode='a',if_sheet_exists='overlay') as writer:
    min_max.to_excel(writer, sheet_name='min_max', index=False, startrow=0, startcol=0)

# get min and max of CGE
CGE_max = np.max([df_CGE.to_numpy()[:, 1:], df_CGE_tot.to_numpy()[:, 1:]])
CGE_min = np.min([df_CGE.to_numpy()[:, 1:], df_CGE_tot.to_numpy()[:, 1:]])


# make every fig separately
for i, df in enumerate(dfs):
    array = dfs[i].to_numpy()[:, 1:]
    array = np.array(array, dtype=np.float64)  # Convert to float
    fig1, ax1 = plt.subplots(figsize=(3.5, 3.5))
    # get min and max of CGE
    if save_as[i] not in ['CGE', 'CGE_tot']:
        contour = ax1.contourf(np.arange(0.21,0.96,0.01), take_from, np.transpose(array), levels=30, cmap='viridis')  # Contour fill plot
    else:
        contour = ax1.contourf(np.arange(0.21,0.96,0.01), take_from, np.transpose(array), levels=30, cmap='viridis', vmin=CGE_min, vmax=CGE_max)
    # Add labels and title
    ax1.set_xlabel('$\\varphi_{O_2}$ [\%]')
    ax1.set_ylabel('$\lambda$ [/]')
    if save_as[i] == 'T':
        ax1.set_title('temperatura uplinjanja [K]')
    elif save_as[i] == 'CGE':
        ax1.set_title('izkoristek uplinjanja $CGE$ [/]')
    elif save_as[i] == 'CGE_tot':
        ax1.set_title('celotni izkoristek uplinjanja $CGE_{tot}$ [/]')
    elif save_as[i] == 'GY':
        ax1.set_title('količina produktnega plina $GY$ [m$^3$/kg]')
    elif save_as[i] == 'LHV':
        ax1.set_title('kurilnost plina $LHV_{pp}$ [MJ/kg]')
    else:
        ax1.set_title('volumski delež ' + save_as[i] + ' [\%]')
    cbar = plt.colorbar(contour, ax=ax1)  # Attach the colorbar to the correct axis
    fig1.tight_layout()
    save = save_as[i].replace('$_', '')
    save = save.replace('$', '')
    fig1.savefig('results_' + save +  '.png', format='png', dpi=300, bbox_inches='tight')
    plt.close('all')
    with pd.ExcelWriter(str_excel_3, engine='openpyxl', mode='a',if_sheet_exists='overlay') as writer:
        dfs[i].to_excel(writer, sheet_name=save, index=False, startrow=0, startcol=0)

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
    contour = ax2[row, col].contourf(np.arange(0.21, 0.96, 0.01)*100, take_from, np.transpose(array), levels=30, cmap='viridis')  # Contour plot
    fig2.colorbar(contour, ax=ax2[row, col])  # Add colorbar
    ax2[row, col].set_xlabel('$\\varphi_{O2}$ [\%]')  # X-axis label
    ax2[row, col].set_ylabel('$\lambda$ [/]')  # Y-axis label
    # Set the title based on the 'save_as' list
    if save_as[i] == 'T':
        ax2[row, col].set_title('temperatura uplinjanja $T$ [K]')  # Title for temperature
    else:
        ax2[row, col].set_title(f'volumski delež {save_as[i]} [\%]')  # Title for other variables

fig2.tight_layout()  # Adjust layout for better spacing
fig2.savefig('results_phi_i.png', format='png', dpi=300, bbox_inches='tight')  # Save the figure as a PNG
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
    if save_as[i] == 'GY' or 'LHV':
        contour = ax6[row, col].contourf(np.arange(0.21, 0.96, 0.01)*100, take_from, np.transpose(array), levels=30, cmap='viridis')  # Contour plot
    else:
        contour = ax6[row, col].contourf(np.arange(0.21, 0.96, 0.01)*100, take_from, np.transpose(array), levels=30, cmap='viridis', vmin=CGE_min, vmax=CGE_max)  # Contour plot
    fig5.colorbar(contour, ax=ax6[row, col])  # Add colorbar
    ax6[row, col].set_xlabel('$\\varphi_{O_2}$ [\%]')  # X-axis label
    ax6[row, col].set_ylabel('$\lambda$ [/]')  # Y-axis label
    # Set the title based on the 'save_as' list
    if save_as[i] == 'CGE':
        ax6[row, col].set_title('izkoristek uplinjanja $CGE$ [/]')
    elif save_as[i] == 'CGE_tot':
        ax6[row, col].set_title('celotni izkoristek uplinjanja $CGE_{tot}$ [/]')
    elif save_as[i] == 'GY':
        ax6[row, col].set_title('količina produktnega plina $GY$ [m$^3$/kg]')
    elif save_as[i] == 'LHV':
        ax6[row, col].set_title('kurilnost plina $LHV_{pp}$ [MJ/kg]')
fig5.tight_layout()  # Adjust layout for better spacing
fig5.savefig('results_phi_i_2.png', format='png', dpi=300, bbox_inches='tight')  # Save the figure as a PNG


plt.close('all')  # Close the plot after saving and showing

# plot of gas composition at ER = 0.30
fig3, ax3 = plt.subplots(figsize=(7,4))
ax4 = ax3.twinx()
p1, = ax3.plot(df_CO['phi_O2'], df_CO['0.30'], label='CO')
p2, = ax3.plot(df_H2['phi_O2'], df_H2['0.30'], label='H$_2$')
p3, = ax3.plot(df_CH4['phi_O2'], df_CH4['0.30'], label='CH$_4$')
p4, = ax3.plot(df_CO2['phi_O2'], df_CO2['0.30'], label='CO$_2$')
p5, = ax3.plot(df_N2['phi_O2'], df_N2['0.30'], label='N$_2$')
p6, = ax4.plot(df_T['phi_O2'], df_T['0.30'], color=colors[7], linestyle='--', label='temperatura')
ax3.set_xlabel('volumski delež kisika $\\varphi_{O_2}$ [\%]')
ax4.set_ylabel('temperatura $T$ [K]')
ax3.set_ylabel('volumski delež $\\varphi_i$ [\%]')
ax3.set_title('$\lambda$ = 0,30')
ax3.set_ylim(top=60)
ax4.grid(None)
ax3.legend(handles=[p1, p2, p3, p4, p5, p6], loc='upper left', ncol=2)
fig3.tight_layout()  # Adjust layout for better spacing
fig3.savefig('results_ER_0_30.png', format='png', dpi=300, bbox_inches='tight')
plt.close('all')  # Close the plot after saving and showing


with plt.style.context(['science', 'grid', 'scatter']):
    fig4, ax4 = plt.subplots(figsize=(7,4))
    ax5 = ax4.twinx()
    p1, = ax4.plot(df_2['case'], df_2['CO'], label='CO', color=colors[0], marker='o')
    p2, = ax4.plot(df_2['case'], df_2['H2'], label='H$_2$', color=colors[1], marker='o')
    p3, = ax4.plot(df_2['case'], df_2['CH4'], label='CH$_4$', color=colors[2], marker='o')
    p4, = ax4.plot(df_2['case'], df_2['CO2'], label='CO$_2$', color=colors[3], marker='o')
    p5, = ax4.plot(df_2['case'], df_2['N2'], label='N$_2$', color=colors[4], marker='o')
    p6, = ax5.plot(df_2['case'], df_2['T [K]'], color=colors[7], label='temperatura', marker='o')
    p11, = ax4.plot(df_3['case'], df_3['CO'], label='CO', color=colors[0], marker='x')
    p12, = ax4.plot(df_3['case'], df_3['H2'], label='H$_2$', color=colors[1], marker='x')
    p13, = ax4.plot(df_3['case'], df_3['CH4'], label='CH$_4$', color=colors[2], marker='x')
    p14, = ax4.plot(df_3['case'], df_3['CO2'], label='CO$_2$', color=colors[3], marker='x')
    p15, = ax4.plot(df_3['case'], df_3['N2'], label='N$_2$', color=colors[4], marker='x')
    p16, = ax5.plot(df_3['case'], df_3['T [K]'], color=colors[7], label='temperatura', marker='x')
    p17, = ax4.plot([], [], 'ko', label='nestehiometrični')
    p18, = ax4.plot([], [], 'kx', label='stehiometrični')
    ax4.set_xlabel('primer')
    ax5.set_ylabel('temperatura $T$ [K]')
    ax4.set_ylabel('volumski delež $\\varphi_i$ [\%]')
    ax4.set_ylim(top=100)
    ax5.set_ylim(bottom=0, top=1600)
    ax5.grid(None)
    ax4.legend(handles=[p1, p2, p3, p4, p5, p6, p17, p18], loc='upper left', ncol=2)
    fig4.tight_layout()  # Adjust layout for better spacing
    fig4.savefig('results_NS_vs_S.png', format='png', dpi=300, bbox_inches='tight')

plt.close('all')  # Close the plot after saving and showing

# plot of gas composition at phi_O = 0.21
fig6, ax7 = plt.subplots(figsize=(7,4))
ax4 = ax7.twinx()
p1, = ax7.plot(df_phi_21['ER'], df_phi_21['CO'], label='CO', color=colors[0])
p2, = ax7.plot(df_phi_21['ER'], df_phi_21['H2'], label='H$_2$', color=colors[1])
p3, = ax7.plot(df_phi_21['ER'], df_phi_21['CH4'], label='CH$_4$', color=colors[2])
p4, = ax7.plot(df_phi_21['ER'], df_phi_21['CO2'], label='CO$_2$', color=colors[3])
p5, = ax7.plot(df_phi_21['ER'], df_phi_21['N2'], label='N$_2$', color=colors[4])
p6, = ax4.plot(df_phi_21['ER'], df_phi_21['T [K]'], color=colors[7], label='temperatura')
p21, = ax7.plot(df_4['ER'], df_4['CO'], ls='--', color=colors[0])
p22, = ax7.plot(df_4['ER'], df_4['H2'], ls='--', color=colors[1])
p23, = ax7.plot(df_4['ER'], df_4['CH4'], ls='--', color=colors[2])
p24, = ax7.plot(df_4['ER'], df_4['CO2'], ls='--', color=colors[3])
p25, = ax7.plot(df_4['ER'], df_4['N2'], ls='--', color=colors[4])
p26, = ax4.plot(df_4['ER'], df_4['T [K]'], color=colors[7], linestyle='--')
p27, = ax4.plot([], [], 'k', label='ne-stehiometrični')
p28, = ax4.plot([], [], 'k', ls='--', label='stehiometrični')
ax7.set_xlabel('razmernik zraka $\lambda$ [/]')
ax4.set_ylabel('temperatura $T$ [K]')
ax7.set_ylabel('volumski delež $\\varphi_i$ [\%]')
ax7.set_title('volumski delež kisika $\\varphi_{O_2}$ = 21 \%')
ax7.set_ylim(top=65)
ax4.grid(None)
ax7.legend(handles=[p1, p2, p3, p4, p5, p6, p27, p28], loc='upper left', ncol=2)
fig6.tight_layout()  # Adjust layout for better spacing
fig6.savefig('results_phi_O_21_S_NS.png', format='png', dpi=300, bbox_inches='tight')
plt.close('all')  # Close the plot after saving and showing

# plot of gas composition at phi_O = 0.21
fig6, ax7 = plt.subplots(figsize=(7,4))
ax4 = ax7.twinx()
p1, = ax7.plot(df_phi_60['ER'], df_phi_60['CO'], label='CO', color=colors[0])
p2, = ax7.plot(df_phi_60['ER'], df_phi_60['H2'], label='H$_2$', color=colors[1])
p3, = ax7.plot(df_phi_60['ER'], df_phi_60['CH4'], label='CH$_4$', color=colors[2])
p4, = ax7.plot(df_phi_60['ER'], df_phi_60['CO2'], label='CO$_2$', color=colors[3])
p5, = ax7.plot(df_phi_60['ER'], df_phi_60['N2'], label='N$_2$', color=colors[4])
p6, = ax4.plot(df_phi_60['ER'], df_phi_60['T [K]'], color=colors[7], label='temperatura')
p21, = ax7.plot(df_5['ER'], df_5['CO'], ls='--', color=colors[0])
p22, = ax7.plot(df_5['ER'], df_5['H2'], ls='--', color=colors[1])
p23, = ax7.plot(df_5['ER'], df_5['CH4'], ls='--', color=colors[2])
p24, = ax7.plot(df_5['ER'], df_5['CO2'], ls='--', color=colors[3])
p25, = ax7.plot(df_5['ER'], df_5['N2'], ls='--', color=colors[4])
p26, = ax4.plot(df_5['ER'], df_5['T [K]'], color=colors[7], linestyle='--')
p27, = ax4.plot([], [], 'k', label='ne-stehiometrični')
p28, = ax4.plot([], [], 'k', ls='--', label='stehiometrični')
ax7.set_xlabel('razmernik zraka $\lambda$ [/]')
ax4.set_ylabel('temperatura $T$ [K]')
ax7.set_ylabel('volumski delež $\\varphi_i$ [\%]')
ax7.set_title('volumski delež kisika $\\varphi_{O_2}$ = 60 \%')
ax7.set_ylim(top=65)
ax4.grid(None)
ax7.legend(handles=[p1, p2, p3, p4, p5, p6, p27, p28], loc='upper left', ncol=2)
fig6.tight_layout()  # Adjust layout for better spacing
fig6.savefig('results_phi_O_60_S_NS.png', format='png', dpi=300, bbox_inches='tight')
plt.close('all')  # Close the plot after saving and showing

end = datetime.datetime.now()
print('\nCalculation time:', end - start) 
