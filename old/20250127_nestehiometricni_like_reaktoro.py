# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 11:48:32 2024

@author: Jaka Bizjak, IJS CEU
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
from scipy.linalg import solve
import reaktoro
from reaktoro import *

plt.close('all')
'''
# oblikovanje grafa
plt.style.use(['science', 'muted', 'grid']) # glej https://github.com/garrettj403/SciencePlots
plt.rcParams.update({
    'font.size': 10,        # Font size for text
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
'''
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "0" # tole je pomembno, ker matplotlib prireja velikost slike gleden na ekran!!!
os.environ["MATPLOTLIB_FORCE_DPI"] = "300"

start = datetime.datetime.now() # stopwatch the time of loop calculation

str_excel = 'MN_model_results.xlsx'
str_sheet ='reaktoro_NS_model'

#*****************************************************************************

# Constants
M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007  # kg/kmol
R = 8.314  # kJ/kmolK
T_0 = 25 + 273.15 # STPC - K
p_0 = 1.013e5

# enthalpy of formation
dH_0_CO = -110527 # kJ/kmol
dH_0_H2 = 0 # kJ/kmol
dH_0_CH4 = -74873 # kJ/kmol
dH_0_CO2 = -393522 # kJ/kmol
dH_0_N2 = 0 # kJ/kmol
dH_0_H2O_l = -285830 # kJ/kmol
dH_0_H2O_g = -241826 # kJ/kmol
dH_0_N2 = 0 # kJ/kmol
dH_0_C = 0 # kJ/kmol
dH_vap = 40.7e3 # kJ/kmol - enthalpy of water vaporization
# Standard Gibbs free energies of formation (kJ/kmol)
dG_f = np.array([
    -137163,   # CO
    0,         # H2
    -50768,    # CH4
    -394389,   # CO2
    -228582   # H2O
])

# Atomic composition matrix (C, H, O, N)
Sigma = np.array([
    [1, 0, 1],  # CO
    [0, 2, 0],  # H2
    [1, 4, 0],  # CH4
    [1, 0, 2],  # CO2
    [0, 2, 1]  # H2O
])

# Read biomass data from file
biomass = pd.read_csv('biomass.txt', sep=",")

all_results = []

cases = biomass.shape[0]
abc = np.arange(15, 16, 1)
for i in abc:
    w_C = biomass['w_C'][i] / 100
    w_H = biomass['w_H'][i] / 100
    w_O = biomass['w_O'][i] / 100
    w_N = biomass['w_N'][i] / 100
    w_S = 0
    MC = biomass['MC'][i] / 100  # Moisture content
    ER = biomass['ER'][i]    # Equivalence ratio
    phi_N = biomass['phi_N'][i] # volume fraction of N2 in gasifying agent
    phi_O = biomass['phi_O'][i] # volume fraction of O2 in gasifying agent
    
    # Calculate atomic ratios
    alpha = w_H * M_C / (w_C * M_H)
    beta = w_O * M_C / (w_C * M_O)
    gamma = w_N * M_C / (w_C * M_N)
    y = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)
    z = ER * (1 + alpha / 4 - beta / 2)
    
    HHV_b = (0.3491*w_C + 1.1783*w_H + 0.1005*w_S - 0.1034*w_O - 0.0151*w_N)*100*1000*(M_C + alpha*M_H + beta*M_O + gamma*M_N) # kJ/kmol - higher heating value of biomass (multipled by 100, so mass fractions are in % - empirical equation)
    M_b = (M_C + alpha*M_H + beta*M_O + gamma*M_N) # kg/kmol
    
    n_N2 = 0.5*gamma + z*phi_N/phi_O # 
    n_C = 0.05
    
    # initial guesses of temp and enthalpy 
    temp = 900
    dH = 1000
    
    while np.absolute(dH) > 100:
    
        #RAND algoritem
        '''
        # initial guesses
        n_CO = 0.6
        n_H2 = 0.6
        n_CH4 = 0.03
        n_CO2 = 0.2
        n_H2O = 0.03
        
        psi = np.array([
            1e-6, # psi_C
            1e-6, # psi_H
            1e-6  # psi_O
            ])
        
        n_array = np.array([
            n_CO,
            n_H2,
            n_CH4,
            n_CO2,
            n_H2O])
        
        epsilon = 1
        count = 0
        
        # Lists to store count, epsilon and other values for ploting
        count_values = []
        epsilon_values = []
        n1_values, n2_values, n3_values, n4_values, n5_values, omega_val = [], [], [], [], [], []
    
        fig1, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.set_xlabel('Iteration Count')
        ax1.set_ylabel('$\epsilon$')
        ax2.set_ylabel('$n_i$/$\omega$')
        ax1.set_title('Convergence Over Iterations')
        ax2.grid()
        ax1.grid()
        #ax1.legend(handels=[line, line1, line2, line3, line4, line5])
        ax2.set_ylim(0, 1.01)  # Example: Set n_array y-limits
    
        # Set up the initial plot with empty data
        line, = ax1.plot([], [], 'k', label='epsilon')  # 'bo' stands for blue circle markers
        line1, = ax2.plot([], [], label='n_CO')
        line2, = ax2.plot([], [], label='n_H2')
        line3, = ax2.plot([], [], label='n_CH4')
        line4, = ax2.plot([], [], label='n_CO2')
        line5, = ax2.plot([], [], label='n_H2O')
        line6, = ax2.plot([], [], label='omega')
        ax1.legend(handles=[line, line1, line2, line3, line4, line5, line6], ncol=2)
        # to run GUI event loop
        plt.ion()
        
        while epsilon > 0.01:
            # Matrix construction
            Gamma = np.diag(np.sum(n_array)/n_array)
            Jacobian = np.block([
                [Gamma, -Sigma],
                [Sigma.T, np.zeros((np.shape(Sigma)[-1], np.shape(Sigma)[-1]))]
                ])
            B = np.zeros(np.shape(Jacobian)[-1])
            for i in range(np.size(n_array)):
                B[i] = -(-dG_f[i]/(R * T) + np.log(n_array[i]/np.sum(n_array)) - np.dot(Sigma[i, :], psi))
            B[5:] = -(np.dot(Sigma.T, n_array) - np.dot(np.array([1 - n_C, alpha + 2*y, beta + y +2*z]).T, 1))
            x = solve(Jacobian, B)
            delta_n = x[:5]
            delta_psi = x[5:]
            
            # calculation of omega
            pos = np.where(delta_n < 0)[0]
            a = delta_n[pos]
            value = np.abs((n_array[pos])/delta_n[pos])
            omega = np.min(np.append(value, 1))
            epsilon = np.max(np.abs(delta_n/n_array))
            n_array = n_array + omega*delta_n + 1e-4
            psi = psi + omega*delta_psi + 1e-4
    
            count +=1
            # Store count and epsilon values
            count_values.append(count)
            epsilon_values.append(epsilon)
            n1_values.append(n_array[0])
            n2_values.append(n_array[1])
            n3_values.append(n_array[2])
            n4_values.append(n_array[3])
            n5_values.append(n_array[4])
            omega_val.append(omega)
            
            # Update plot with new data
            line.set_xdata(count_values)
            line.set_ydata(epsilon_values)
            line1.set_xdata(count_values)
            line1.set_ydata(n1_values)
            line2.set_xdata(count_values)
            line2.set_ydata(n2_values)
            line3.set_xdata(count_values)
            line3.set_ydata(n3_values)
            line4.set_xdata(count_values)
            line4.set_ydata(n4_values)
            line5.set_xdata(count_values)
            line5.set_ydata(n5_values)
            line6.set_xdata(count_values)
            line6.set_ydata(omega_val)
            ax1.relim()  # Recalculate limits
            ax1.autoscale_view()  # Autoscale view
            fig1.canvas.draw()  # Redraw the plot
            fig1.canvas.flush_events()  # Ensure the plot updates in real-time
            plt.pause(0.01)  # Pause for a small interval to allow GUI to update
        
            if count >= 100:
                break
        # Disable interactive mode after the loop is done
        plt.ioff()
        #plt.show()  # Make sure the plot window stays open after the loop
        
        print('count = ', count)
        
        '''

        # Load the thermodynamic database and define the system
        db = NasaDatabase("nasa-cea")  # NASA CEA database
        
        formula = ChemicalFormula(f"CH{alpha}O{beta}N{gamma}")
        
        state = ChemicalState(system)
        
        conditions = EquilibriumConditions(system)
        conditions.temperature(temp, "kelvin")
        conditions.pressure(1.013, "bar")
        conditions.setInitialComponentAmounts(b)
        solver = EquilibriumSolver(system)
        results = solver.solve(state, conditions)
        succeeded =  results.succeeded()

        # Extract species amounts as float variables
        n_CO = state.speciesAmount("CO")
        n_H2 = state.speciesAmount("H2")
        n_CH4 = state.speciesAmount("CH4")
        n_CO2 = state.speciesAmount("CO2")
        n_H2O = state.speciesAmount("H2O")
        n_N2 = state.speciesAmount("N2")
        
        n_tot = n_CH4 + n_CO + n_CO2 + n_H2O + n_H2 + n_N2 - n_H2O
        
        # volume fractions in produced gas
        phi_CH4 = n_CH4/n_tot
        phi_CO = n_CO/n_tot
        phi_CO2 = n_CO2/n_tot
        phi_H2 = n_H2/n_tot
        phi_H2O = n_H2O/n_tot
        phi_N2 = n_N2/n_tot
        
        # energy balance        
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
        dH_reaktanti = dH_b + y*(dH_0_H2O_g)# + dH_vap)
        dH_produkti = n_H2*dH_H2 + n_CO*dH_CO + n_CO2*dH_CO2 + n_CH4*dH_CH4 + n_H2O*dH_H2O_g + n_N2*dH_N2 + n_C*dH_C 
        dH = dH_produkti - dH_reaktanti
    
        # if dH is greater than 100, then increase temperature
        temp = temp + 1
        
    M_g = (M_C + M_O)*phi_CO + 2*M_H*phi_H2 + (M_C + 4*M_H)*phi_CH4 + (M_C + 2*M_O)*phi_CO2 + 2*M_N*phi_N2
    
    LHV_b = HHV_b - (9*w_H)*dH_vap # kJ/kmol - LHV of biomass
    LHV_b_m = LHV_b/(M_b*1000) # MJ/kg
    LHV_g = (phi_CO*12.622 + phi_H2*10.788 + phi_CH4*35.814)*1000*1000*R*T_0/p_0
    LHV_g_V = (phi_CO*12.622 + phi_H2*10.788 + phi_CH4*35.814)*1000
    LHV_g_m = LHV_g/(M_g*1000) # MJ/kg
    CGE = n_tot*M_g*LHV_g_m/(M_b*LHV_b_m)
    GY = n_tot*R*T_0*1000/(M_b*p_0)

    print('dH: ', dH)
    result = [i+1, MC*100, w_C*100, w_H*100, w_O*100, w_N*100, 
              M_b, LHV_b_m, ER, np.round(phi_N*100, 2), np.round(phi_O*100, 2), 
              temp, phi_CO*100, phi_H2*100, phi_CH4*100, phi_CO2*100, phi_H2O*100, phi_N2*100, M_g, LHV_g_V/1000, LHV_g_m, CGE, GY, succeeded]
    all_results.append(result)

# Create DataFrame from results
df_results = pd.DataFrame(all_results, columns=['case', 'MC', 'w_C', 'w_H', 'w_O', ' w_N', 
                                                'M_b [kg/kmol]', 'LHV_b [MJ/kg]', 'ER', 'phi_N', 'phi_O', 
                                                'T [K]', 'CO', 'H2', 'CH4', 'CO2', 'H2O', 'N2', 'M_g [kg/kmol]', 'LHV_g [MJ/m3]', 'LHV_g [MJ/kg]', 'CGE', 'GY', 'succeeded'])

# Replace all numbers (ints and floats) with comma-separated decimals
def format_numbers(value):
    if isinstance(value, (int, float)):
        return str(value).replace('.', ',')  # Replace '.' with ','
    return value  # Leave non-numeric values untouched

# Apply formatting to each column of the DataFrame
df_results = df_results.apply(lambda col: col.map(format_numbers))

# Write the DataFrame to the sheet, starting from the next row
with pd.ExcelWriter(str_excel, engine='openpyxl', mode='a',if_sheet_exists='overlay') as writer:
    df_results.to_excel(writer, sheet_name=str_sheet, index=False, header =False, startrow=25, startcol=0)



end = datetime.datetime.now()
print('\nCalculation time:', end - start) 
