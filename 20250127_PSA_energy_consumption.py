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
from scipy.optimize import curve_fit


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
str_excel = 'MN_model_results.xlsx'
str_sheet = 'PSA_energy_consumption'

def constrained_quadratic(x, a):
    x_0 = 21  # Point where y = 0
    return a * (x - x_0)**2 



df = pd.read_excel(str_excel, sheet_name=str_sheet) 

# Fit a linear trendline (degree=1 for a straight line)
x = np.linspace(21, 95, 100)
coef_2 = np.polyfit(df['phi_O2 [%]'], df['E [MJ/kg O2]'], deg=2)  # Returns slope and intercept
fit_2 = np.polyval(coef_2, x)  # Compute y-values of the trendline
coef_3 = np.polyfit(df['phi_O2 [%]'], df['E [MJ/kg O2]'], deg=3)  # Returns slope and intercept
fit_3 = np.polyval(coef_3, x)  # Compute y-values of the trendline

popt, _ = curve_fit(constrained_quadratic, df['phi_O2 [%]'], df['E [MJ/kg O2]'])
a_fit = popt[0]  # Extract the scalar value of 'a'

equation_label = f'$E = {a_fit:.6f}'+'(\\varphi_{O2} - 21)^2$'


fig2, ax3 = plt.subplots(figsize=(5, 2.9))

ax3.scatter(df.loc[df['source'] == 'Banaszkiewicz', 'phi_O2 [%]'], df.loc[df['source'] == 'Banaszkiewicz', 'E [MJ/kg O2]'], label='Banaszkiewicz et al.', s=5)
ax3.scatter(df.loc[df['source'] == 'Šulc', 'phi_O2 [%]'], df.loc[df['source'] == 'Šulc', 'E [MJ/kg O2]'], label='Šulc et al.', s=5)
#ax3.plot(x, fit_2, ls='--', label='polinom 2. stopnje')
#ax3.plot(x, fit_3, ls='--', label='polinom 3. stopnje')
ax3.plot(x, constrained_quadratic(x, a_fit), ls='--', label=equation_label.replace('.', ','), color=colors[2])


# Title and grid
#ax3.set_title(opis)
ax3.set_xlabel('volumski delež kisika v zraku $\\varphi_{O2}$ [\%]')
ax3.set_ylabel('raba energije $E$ [MJ/kg O$_2$]')

# Add legend
ax3.legend(loc='upper left', frameon=1)

# Apply tight layout
plt.tight_layout()

# Optionally save the figure
fig2.savefig('PSA_energy_consumption.png', format='png', dpi=300, bbox_inches='tight')




end = datetime.datetime.now()
print('\nCalculation time:', end - start) 
