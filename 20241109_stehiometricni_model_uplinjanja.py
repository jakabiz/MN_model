# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 13:05:33 2024

@author: Jaka Bizjak, IJS CEU

Stehiometrični ravnotežni model uplinjanja biomase v sotočnem uplinjevalniku -
- prirejeno po:
    Pikš K. Modeliranje uplinjanja biomase. Published online 2022. 
     Accessed February 23, 2024. https://repozitorij.uni-lj.si/IzpisGradiva.php?id=137012
"""
import numpy as np
import pandas as pd
from pandas import DataFrame
import locale
import matplotlib
matplotlib.use('Qt5Agg')  # Ensure Qt backend is being used
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import scienceplots
import datetime
import os
import openpyxl
from openpyxl import load_workbook
from scipy.linalg import solve
import reaktoro
from reaktoro import *
from scipy.optimize import fsolve

# oblikovanje grafa
plt.style.use(['science', 'muted', 'grid']) # glej https://github.com/garrettj403/SciencePlots
plt.rcParams.update({
    'font.size': 10,        # Font size for text
    'figure.titlesize': 10, # Font size for figure suptitle
    'axes.titlesize': 10,   # Font size for axes titles
    'axes.labelsize': 10,   # Font size for x and y labels
    'xtick.labelsize': 10,  # Font size for x-axis tick labels
    'ytick.labelsize': 10,  # Font size for y-axis tick labels
    'legend.fontsize': 10,  # Font size for legend
    'lines.linewidth': 2    # Default line width
})
colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] # list of default colours
colors[2] = '#DDAA33'# it changes #DDCC77
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# List of colors
# i     code    colour  what do we use it for
# 0     #CC6677 red     stroški, moč
# 1     #332288 blue    seštevki spremenljivk, energija
# 2     #DDAA33 yellow  zemeljski plin, tok
# 3     #117733 green   elektrika, napetost
# 4     #88CCEE cyan    TGP
# 5     #882255 wine
# 6     #44AA99 teal    temperatura
# 7     #999933 olive   kapaciteta
# 8     #AA4499 purple
# 9     #DDDDDD grey
locale.setlocale(locale.LC_NUMERIC, "de_DE.UTF-8") # zamenjava . z , v grafih
plt.rcParams['axes.formatter.use_locale'] = True

os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "0" # tole je pomembno, ker matplotlib prireja velikost slike gleden na ekran!!!
os.environ["MATPLOTLIB_FORCE_DPI"] = "300"

start = datetime.datetime.now()  # Stopwatch the time of loop calculation

output_file = 'gas_composition_stoichiometric_model.txt'
str_excel = 'MN_model_results.xlsx'
str_sheet ='stoichiometric_model_ERs'

# Check if the file exists; if so, remove it
if os.path.exists(output_file):
    os.remove(output_file)

# Basic stehiometric model

# constant values
log_K_H2 = 0.
log_K_C = 0.
log_K_O2 = 0.
M_C = 12.011 # kg/kmol
M_H = 1.008 # kg/kmol
M_O = 15.999 # kg/kmol
M_N = 14.007 # kg/kmol
R = 8.314 # kJ/kmolK - gas constant
p_0 = 1.013e5
T_0 = 25 + 273.15 # K - initail temperature

# enthalpy of formation
dH_0_CH4 = -74873 # kJ/kmol
dH_0_H2 = 0 # kJ/kmol
dH_0_CO = -110527 # kJ/kmol
dH_0_CO2 = -393522 # kJ/kmol
dH_0_N2 = 0 # kJ/kmol
dH_0_H2O_l = -285830 # kJ/kmol
dH_0_H2O_g = -241826 # kJ/kmol
dH_0_C = 0 # kJ/kmol
dH_vap = 40.7e3 # kJ/kmol - enthalpy of water vaporization

# Defining reaction constants as function of temperature
# Load data form NIST-JANAF tables
data_CH4 = pd.read_csv('NIST_JANAF_CH4.txt', sep=",", header=1)
data_CO = pd.read_csv('NIST_JANAF_CO.txt', sep=",", header=1)
data_CO2 = pd.read_csv('NIST_JANAF_CO2.txt', sep=",", header=1)
data_H2O = pd.read_csv('NIST_JANAF_H2O.txt', sep=",", header=1)

# Define a function to clean the data
def clean_data(df):
    # Drop NaN or Inf values
    df = df.dropna()
    df = df[np.isfinite(df['T/K']) & np.isfinite(df['logKf'])]
    
    # Convert to float if necessary
    df['T/K'] = df['T/K'].astype(float)
    df['logKf'] = df['logKf'].astype(float)
    
    return df

# Clean each dataset
data_CH4 = clean_data(data_CH4)
data_CO = clean_data(data_CO)
data_CO2 = clean_data(data_CO2)
data_H2O = clean_data(data_H2O)

# Calculate linear fits for logKf as a function of T
coef_CH4 = np.polyfit(1/data_CH4['T/K'].values, data_CH4['logKf'].values, 1)
coef_CO = np.polyfit(1/data_CO['T/K'].values, data_CO['logKf'].values, 1)
coef_CO2 = np.polyfit(1/data_CO2['T/K'].values, data_CO2['logKf'].values, 1)
coef_H2O = np.polyfit(1/data_H2O['T/K'].values, data_H2O['logKf'].values, 1)

# Plot of data from NIST-JANAF tables and coresponding linear fits
fig1, ax1 = plt.subplots(figsize=(7,4))
ax1.plot(1/data_CH4['T/K'], data_CH4['logKf'], 'x', color=colors[0], label='CH4')
ax1.plot(1/data_CH4['T/K'], np.polyval(coef_CH4, 1/data_CH4['T/K']), color=colors[0])
ax1.plot(1/data_CO['T/K'], data_CO['logKf'], 'x', color=colors[1], label='CO')
ax1.plot(1/data_CO['T/K'], np.polyval(coef_CO, 1/data_CO['T/K']), color=colors[1])
ax1.plot(1/data_CO2['T/K'], data_CO2['logKf'], 'x', color=colors[2], label='CO2')
ax1.plot(1/data_CO2['T/K'], np.polyval(coef_CO2, 1/data_CO2['T/K']), color=colors[2])
ax1.plot(1/data_H2O['T/K'], data_H2O['logKf'], 'x', color=colors[3], label='H2O')
ax1.plot(1/data_H2O['T/K'], np.polyval(coef_H2O, 1/data_H2O['T/K']), color=colors[3])
ax1.set_xlim([0,0.0051])
ax1.set_ylim([-20,110])
ax1.legend()
ax1.set_xlabel('$1/T$ [1/K]')
ax1.set_ylabel('$\log K_f$ [/]')
fig1.tight_layout()  # Adjust layout for better spacing
fig1.savefig('results_log_K_T.png', format='png', dpi=300, bbox_inches='tight')

#ax1.set_title('Equilibrium Constant vs. Temperature (Arrhenius Plot)')

# define functions for calculating log_K for given substance depending on temperature
def log_K_CH4(T):
    return np.polyval(coef_CH4, 1/T)
def log_K_CO(T):
    return np.polyval(coef_CO, 1/T)
def log_K_CO2(T):
    return np.polyval(coef_CO2, 1/T)
def log_K_H2O(T):
    return np.polyval(coef_H2O, 1/T)

# logK plot of oxidation and reduction reactions
def log_K_1(T): # C + O2 -> CO2 - 1st oxydation reaction
    return log_K_CO2(T) - (log_K_C + log_K_O2) 
def log_K_2(T): # C + 1/2 O2 -> CO - 2nd oxydation reaction
    return log_K_CO(T) - (log_K_C + 0.5*log_K_O2) 
def log_K_3(T): # H2 + 1/2 O2 -> H2O - 3rd oxydation reaction
    return log_K_H2O(T) - (log_K_H2 + 0.5*log_K_O2) 
def log_K_4(T): # C + CO2 <-> 2CO - Boudouard reaction (reduction)
    return 2*log_K_CO(T) - (log_K_C + log_K_CO2(T)) 
def log_K_5(T): # C + H2O <-> CO + H2 - char and steam reduction reaction
    return log_K_CO(T) + log_K_H2 - (log_K_C + log_K_H2O(T))
def log_K_6(T): # CO + H2O <-> CO2 + H2 - Water Gas Shift reaction (reduction)
    return log_K_CO2(T) + log_K_H2 - (log_K_CO(T) + log_K_H2O(T))
def log_K_7(T): # C + 2H2 <-> CH4 - methanization reaction (reduction)
    return log_K_CH4(T) - (log_K_C + 2*log_K_H2)
def log_K_8(T): # CH4 + H2O <-> CO2 + 3H2 - steam reforming reaction (reduction)
    return log_K_CO2(T) + 3*log_K_H2 - (log_K_CH4(T) + log_K_H2O(T))
T = np.linspace(200, 2500, 200)
fig1, ax1 = plt.subplots(figsize=(7,4))
ax1.plot(1/T, log_K_1(T), label='popolno zgorevanje ogljika')
ax1.plot(1/T, log_K_2(T), label='nepopolno zgorevanje ogljika')
ax1.plot(1/T, log_K_3(T), label='popolno zgorevanje vodika')
ax1.plot(1/T, log_K_4(T), label='Boudouardova reakcija')
ax1.plot(1/T, log_K_5(T), label='reakcija oglja z vodo')
ax1.plot(1/T, log_K_6(T), label='reakcija CO z vodno paro (WGS)')
ax1.plot(1/T, log_K_7(T), label='reakcija metanacije')
ax1.plot(1/T, log_K_8(T), label='parni reforming metana')

# Add secondary x-axis (Temperature T)
ax2 = ax1.secondary_xaxis('top')
ax2.set_xlabel('$T$ [K]')

# Get primary x-axis ticks
primary_ticks = ax1.get_xticks()

# Filter valid values (avoid division by zero)
valid_ticks = primary_ticks[primary_ticks > 0]

# Set secondary x-axis ticks at the same locations
ax2.set_xticks(valid_ticks)
ax2.set_xticklabels([f"{int(1/tick)}" for tick in valid_ticks])  # Convert 1/T to T and format



#ax1.set_xlim([0,0.0051])
ax1.set_ylim([-50,150])
ax1.legend(loc='upper left', ncol=2)
ax1.set_xlabel('$1/T$ [1/K]')
ax1.set_ylabel('$\log K$ [/]')
fig1.tight_layout()  # Adjust layout for better spacing
fig1.savefig('results_log_K_T_rections.png', format='png', dpi=300, bbox_inches='tight')

'''
# get proximate and ultimate analysis data from literature saved in files
biomass = pd.read_csv('biomass.txt', sep=",")
model_Piks = pd.read_csv('gas_composition_basic_model_Piks.txt', sep=",", header=0)

all_results = []
ERs = np.arange(0.1, 0.51, 0.01)
for ER in ERs:
    '''
'''
cases = biomass.shape[0]
for i in np.arange(0, cases, 1):
    
    # define biomass
    w_C = biomass['w_C'][i]/100
    w_H = biomass['w_H'][i]/100
    w_O = biomass['w_O'][i]/100
    w_N = biomass['w_N'][i]/100
    w_S = 0
    MC = biomass['MC'][i]/100
    ER = biomass['ER'][i]
    
    phi_N = biomass['phi_N'][i]
    phi_O = biomass['phi_O'][i]
    '''
'''
    phi_N = 0.79
    phi_O = 0.21
    
    w_C = 0.521
    w_H = 0.062
    w_O = 0.412
    w_N = 0.005
    w_S = 0
    MC = 0.10  # Moisture content
    
    alpha = w_H*M_C/(w_C*M_H)
    beta = w_O*M_C/(w_C*M_O)
    gamma = w_N*M_C/(w_C*M_N)
    x = MC*(M_C + alpha*M_H + beta*M_O + gamma*M_N)/(2*M_H + M_O)
    z = ER*(1 + alpha/4 - beta/2)
    
    M_b = (M_C + alpha*M_H + beta*M_O + gamma*M_N) # kg/kmol
    HHV_b = (0.3491*w_C + 1.1783*w_H + 0.1005*w_S - 0.1034*w_O - 0.0151*w_N)*100*1000*M_b # kJ/kmol - higher heating value of biomass (multipled by 100, so mass fractions are in % - empirical equation)
    
    n_C = 0.05 # we assume 5 % of carbon is left as fixed carbon
    n_N2 = 0.5*gamma + phi_N/phi_O*z
    
    
    dH=1000 # initial guess of energy balance so the while loop starts
    temp = 273 # initial temperature
    while np.absolute(dH) > 200:
        K_CH4 = log_K_CH4(temp)
        K_CO = log_K_CO(temp)
        K_CO2 = log_K_CO2(temp)
        K_H2O = log_K_H2O(temp)
        
        # calculation of gas composition
        # define system of equations
        def equations(vars):
            n_CH4, n_CO, n_CO2, n_H2O, n_H2 = vars
            n_tot = n_CH4 + n_H2 + n_H2O + n_CO + n_CO2 + n_N2
            eq_C = n_CO + n_CO2 + n_CH4 + n_C - 1 
            eq_H = 2*n_H2 + 4*n_CH4 + 2*n_H2O - alpha -2*x
            eq_O = n_CO + 2*n_CO2 + n_H2O - beta - x - 2*z
            eq_WGS = n_CO*n_H2O*10**(K_CO2 + log_K_H2 - (K_CO + K_H2O)) - n_CO2*n_H2
            eq_MET = n_H2**2*10**(K_CH4 - (log_K_C + 2*log_K_H2)) - n_CH4
            return [eq_C, eq_H, eq_O, eq_WGS, eq_MET]
        
        # Initial guess for the solution
        initial_guess = [0.01, 0.1, 0.1, 0.5, 0.1]  # Initial guess for n_CH4, n_CO, n_CO2, n_H2O, n_H2
        
        # Solve the system using fsolve
        n_CH4, n_CO, n_CO2, n_H2O, n_H2 = fsolve(equations, initial_guess)
        
        n_tot = n_CH4 + n_CO + n_CO2 + n_H2O + n_H2 + n_N2 - n_H2O
        
        # volume fractions in produced gas
        phi_CH4 = n_CH4/n_tot
        phi_CO = n_CO/n_tot
        phi_CO2 = n_CO2/n_tot
        phi_H2 = n_H2/n_tot
        phi_H2O = n_H2O/n_tot
        phi_N2 = n_N2/n_tot
        
        # c_p for all gases depending on temperature
        cp_CH4 = R*(1.702 + 9.081e-3*((T_0 + temp)/2) + (-2.164e-6)*(4*((T_0 - temp)/2)**2 - temp*T_0)/3) # J/kmolK
        cp_H2 = R*(3.249 + 0.422e-3*((T_0 + temp)/2) + 0.083e5/(temp*T_0)) # kJ/kmolK
        cp_CO = R*(3.376 + 0.557e-3*((T_0 + temp)/2) + (-0.031e5)/(temp*T_0)) # kJ/kmolK
        cp_CO2 = R*(5.457 + 1.047e-3*((T_0 + temp)/2) + (-1.157e5)/(temp*T_0)) # kJ/kmolK
        cp_N2 = R*(3.28 + 0.593e-3*((T_0 + temp)/2) + (0.04e5)/(temp*T_0)) # kJ/kmolK
        cp_H2O = R*(3.47 + 1.45e-3*((T_0 + temp)/2) + (0.121e5)/(temp*T_0)) # kJ/kmolK
        cp_C = R*(1.771 + 0.771e-3*((T_0 + temp)/2) + (-0.867e5)/(temp*T_0)) # kJ/kmolK
        
        dH_b = HHV_b + dH_0_CO2 + alpha/2*dH_0_H2O_l # kJ/kmol - formation enthalpy of biomass
        dH_H2O_g = cp_H2O*(temp - T_0) + dH_0_H2O_g # kJ/kmol
        dH_H2 = cp_H2*(temp - T_0) + dH_0_H2 # kJ/kmol
        dH_CO = cp_CO*(temp - T_0) + dH_0_CO # kJ/kmol
        dH_CO2 = cp_CO2*(temp - T_0) + dH_0_CO2 # kJ/kmol
        dH_N2 = cp_N2*(temp - T_0) + dH_0_N2 # kJ/kmol
        dH_CH4 = cp_CH4*(temp - T_0) + dH_0_CH4 # kJ/kmol
        dH_C = cp_C*(temp - T_0) + dH_0_C # kJ/kmol
        
        # enthalpy of reactants and products
        dH_reaktanti = dH_b + x*(dH_0_H2O_g + dH_vap)
        dH_produkti = n_H2*dH_H2 + n_CO*dH_CO + n_CO2*dH_CO2 + n_CH4*dH_CH4 + n_H2O*dH_H2O_g + n_N2*dH_N2 + n_C*dH_C 
        dH = dH_produkti - dH_reaktanti
        #print('dH: ', dH)
        # if dH is greater than 100, then increase temperature
        temp = temp + 1
    print(dH)
    M_g = (M_C + M_O)*phi_CO + 2*M_H*phi_H2 + (M_C + 4*M_H)*phi_CH4 + (M_C + 2*M_O)*phi_CO2 + 2*M_N*phi_N2
    
    LHV_b = HHV_b - (9*w_H)*dH_vap # kJ/kmol - LHV of biomass
    LHV_b_m = LHV_b/(M_b*1000) # MJ/kg
    LHV_g = (phi_CO*12.622 + phi_H2*10.788 + phi_CH4*35.814)*1000*1000*R*T_0/p_0
    LHV_g_V = (phi_CO*12.622 + phi_H2*10.788 + phi_CH4*35.814)*1000
    LHV_g_m = LHV_g/(M_g*1000) # MJ/kg
    CGE = (n_CO*(M_C+M_O) + n_H2*2*M_H + n_CH4*(M_C + 4*M_H))*LHV_g_m/(M_b*LHV_b_m)
    CGE_tot = (n_CO*(M_C+M_O) + n_H2*2*M_H + n_CH4*(M_C + 4*M_H))*LHV_g_m/(M_b*LHV_b_m + z*2*M_O*(0.000559*(phi_O*100-21)**2)) # including energy consumption for PSA
    # E_PSA is based on fit on experimental data - see PSA_energy_consumption.png and update if necessery
    GY = n_tot*R*T_0*1000/(M_b*p_0)

    result = [MC*100, w_C*100, w_H*100, w_O*100, w_N*100, 
              M_b, LHV_b_m, np.round(ER, 2), np.round(phi_N*100, 2), np.round(phi_O*100, 2), 
              temp, phi_CO*100, phi_H2*100, phi_CH4*100, phi_CO2*100, phi_H2O*100, phi_N2*100, n_C, M_g, LHV_g_V/1000, LHV_g_m, CGE, CGE_tot, GY,
              n_CO, n_H2, n_CH4, n_CO2, n_H2O, n_N2]
    all_results.append(result)
    print('ER = ', np.round(ER,2), ' ; phi_O = ', np.round(phi_O, 2))
    # Create DataFrame from results
df_results = pd.DataFrame(all_results, columns=['MC', 'w_C', 'w_H', 'w_O', ' w_N', 
                                                    'M_b [kg/kmol]', 'LHV_b [MJ/kg]', 'ER', 'phi_N', 'phi_O', 
                                                    'T [K]', 'CO', 'H2', 'CH4', 'CO2', 'H2O', 'N2', 'n_C_r [kmol]', 'M_g [kg/kmol]', 'LHV_g [MJ/m3]', 'LHV_g [MJ/kg]', 'CGE', 'CGE_tot', 'GY',
                                                    'n_CO [kmol]', 'n_H2 [kmol]', 'n_CH4 [kmol]', 'n_CO2 [kmol]', 'n_H2O [kmol]', 'n_N2 [kmol]'])
        


# Write the DataFrame to the sheet, starting from the next row
with pd.ExcelWriter(str_excel, engine='openpyxl', mode='a',if_sheet_exists='overlay') as writer:
    df_results.to_excel(writer, sheet_name=str_sheet, index=False, startrow=0, startcol=0)

'''
plt.close('all')

end = datetime.datetime.now()
print('\nCalculation time:', end - start)
