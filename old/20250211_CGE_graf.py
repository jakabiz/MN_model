
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
 
take_from = np.arange(0.1, 0.51, 0.01)
take_from_str = ["%.2f" % number for number in take_from]

frac_O2 = {'phi_O2': np.arange(21, 96, 1)}
df_CO, df_H2, df_CH4, df_CO2, df_N2, df_H2O, df_M_b, df_LHV_b, df_LHV = [pd.DataFrame(frac_O2) for i in range(9)]

df_phi_21 = pd.DataFrame()

for value in take_from_str:
    str_sheet = str(value)
    #print(str_sheet)
    df_1 = pd.read_excel(str_excel_1, sheet_name=str_sheet)
    df_phi_21 = pd.concat([df_phi_21, df_1[df_1['phi_O'] == '21,0']], ignore_index=True)
    #print(df_phi_21)
    df_CO[str(value)] = df_1['n_CO [kmol]']
    df_H2[str(value)] = df_1['n_H2 [kmol]']
    df_CH4[str(value)] = df_1['n_CH4 [kmol]']
    df_CO2[str(value)] = df_1['n_CO2 [kmol]']
    df_N2[str(value)] = df_1['n_N2 [kmol]']
    df_H2O[str(value)] = df_1['n_H2O [kmol]']
    df_M_b[str(value)] = df_1['M_b [kg/kmol]']
    df_LHV_b[str(value)] = df_1['LHV_b [MJ/kg]']
    df_LHV[str(value)] = df_1['LHV_g [MJ/kg]']

dfs = ([df_CO, df_H2, df_CH4, df_CO2, df_N2, df_H2O, df_M_b.replace(',', '.', regex=True).astype(float), df_LHV_b.replace(',', '.', regex=True).astype(float), df_LHV])

arrays = [df.to_numpy()[:, 1:] for df in dfs]

n_CO = arrays[0]
n_H2 = arrays[1]
n_CH4 = arrays[2]
n_CO2 = arrays[3]
n_N2 = arrays[4]
M_b = arrays[6]
LHV_b = arrays[7]
LHV_g = arrays[8]
n_tot = arrays[0] + arrays[1] + arrays[2] + arrays[3] + arrays[4]

CGE = (n_CO*28 + n_H2*2 + n_CH4*20 + n_CO2*44 + n_N2*28)*LHV_g/(M_b*LHV_b)

alpha = 1.41799
beta = 0.59367
ER = np.ones([np.size(np.arange(21,96,1)), np.size(np.arange(10,51,1))])
phi_O = np.ones([np.size(np.arange(21,96,1)), np.size(np.arange(10,51,1))])
ERs = np.arange(0.1, 0.51, 0.01)
phi_Os = np.arange(0.21, 0.96, 0.01)
for i in range(0, np.size(np.arange(10,51,1))):
    ER[:,i] = ER[:,i]*ERs[i]
for i in range(0, np.size(np.arange(21,96,1))):
    phi_O[i,:] = phi_O[i,:]*phi_Os[i]
    
z = ER*(1 + alpha/4 - beta/2)
CGE_tot = (n_CO*28 + n_H2*2 + n_CH4*20 + n_CO2*44 + n_N2*28)*LHV_g/(M_b*LHV_b + z*2*16*0.000559*(phi_O*100 - 21)**2)

fig1, ax1 = plt.subplots(figsize=(3.5, 3.5))
contour = ax1.contourf(np.arange(0.21,0.96,0.01), take_from, np.transpose(CGE), levels=30, cmap='viridis')  # Contour fill plot
fig1.colorbar(contour, ax=ax1)  # Add colorbar
# Add labels and title
ax1.set_xlabel('$\\varphi_{O2}$ [\%]')
ax1.set_ylabel('$ER$ [/]')
ax1.set_title('izkoristek uplinjanja $CGE$ [/]')
#ax1.set_title('celotni izkoristek uplinjanja $CGE_{tot}$ [/]')

fig1.tight_layout()
fig1.savefig('results_CGE_v0.png', format='png', dpi=300, bbox_inches='tight')
plt.close('all')

# fig 2 - cge tot
fig2, ax2 = plt.subplots(figsize=(3.5, 3.5))
contour = ax2.contourf(np.arange(0.21,0.96,0.01), take_from, np.transpose(CGE_tot), levels=30, cmap='viridis')  # Contour fill plot
fig2.colorbar(contour, ax=ax2)  # Add colorbar
# Add labels and title
ax2.set_xlabel('$\\varphi_{O2}$ [\%]')
ax2.set_ylabel('$ER$ [/]')
ax2.set_title('izkoristek uplinjanja $CGE$ [/]')
#ax2.set_title('celotni izkoristek uplinjanja $CGE_{tot}$ [/]')

fig2.tight_layout()
fig2.savefig('results_CGE_tot_v0.png', format='png', dpi=300, bbox_inches='tight')
plt.close('all')
