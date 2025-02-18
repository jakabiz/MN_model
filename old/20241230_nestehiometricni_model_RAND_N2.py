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

#*****************************************************************************

# Constants
M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007  # kg/kmol
R = 8.314  # kJ/kmolK
T_0 = 25 + 273.15 # STPC - K

# enthalpy of formation
dH_f_CO = -110527 # kJ/kmol
dH_f_H2 = 0 # kJ/kmol
dH_f_CH4 = -74873 # kJ/kmol
dH_f_CO2 = -393522 # kJ/kmol
dH_f_N2 = 0 # kJ/kmol
dH_f_H2O_l = -285830 # kJ/kmol
dH_f_H2O_g = -241826 # kJ/kmol
dH_f_N2 = 0 # kJ/kmol
dH_f_C = 0 # kJ/kmol
dH_vap = 40.7e3 # kJ/kmol - enthalpy of water vaporization
# Standard Gibbs free energies of formation (kJ/kmol)
dG_f = np.array([
    -137163,   # CO
    0,         # H2
    -50768,    # CH4
    -394389,   # CO2
    -228582,   # H2O
    0          # N2 
])

# Atomic composition matrix (C, H, O, N)
Sigma = np.array([
    [1, 0, 1, 0],  # CO
    [0, 2, 0, 0],  # H2
    [1, 4, 0, 0],  # CH4
    [1, 0, 2, 0],  # CO2
    [0, 2, 1, 0],  # H2O
    [0, 0, 0, 2]   # N2
])

E = np.shape(Sigma)[-1]
S = np.shape(Sigma)[0]

# Read biomass data from file
biomass = pd.read_csv('biomass.txt', sep=",")
model_Piks = pd.read_csv('gas_composition_basic_model_Piks.txt', sep=",", header=0)
stehiometric_model = pd.read_csv('gas_composition_stehiometric_model.txt', sep=",", header=0)

all_results = []

cases = 1#biomass.shape[0]

for i in np.arange(0, cases, 1):
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
    
    n_N2 = 0.5*gamma + z*phi_N/phi_O # 
    n_C = 0.05
    
    # add looping over temperatures
    T = 964
    
    #RAND algoritem
    
    # initial guesses
    n_CO = 0.6
    n_H2 = 0.6
    n_CH4 = 0.03
    n_CO2 = 0.2
    n_H2O = 0.03
    
    psi = np.array([
        1e-6, # psi_C
        1e-6, # psi_H
        1e-6, # psi_O
        1e-6  # psi_N
        ])
    
    n_array = np.array([
        n_CO,
        n_H2,
        n_CH4,
        n_CO2,
        n_H2O,
        n_N2
        ])
    
    epsilon = 1
    count = 0
    
    # Lists to store count, epsilon and other values for ploting
    count_values = []
    epsilon_values = []
    n1_values, n2_values, n3_values, n4_values, n5_values, n6_values, omega_val = [], [], [], [], [], [], []

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
    line6, = ax2.plot([], [], label='n_N2')
    line7, = ax2.plot([], [], label='omega')
    ax1.legend(handles=[line, line1, line2, line3, line4, line5, line6, line7], ncol=2)
    # to run GUI event loop
    plt.ion()
    
    while epsilon > 0.01:
        # Matrix construction
        Gamma = np.diag(np.sum(n_array)/n_array)
        Jacobian = np.block([
            [Gamma, -Sigma],
            [Sigma.T, np.zeros((np.shape(Sigma)[-1], np.shape(Sigma)[-1]))]
            ])
        B = np.zeros(E+S)
        for i in range(np.size(n_array)):
            B[i] = -(-dG_f[i]/(R * T) + np.log(n_array[i]/np.sum(n_array)) - np.dot(Sigma[i, :], psi))
        B[S:] = -(np.dot(Sigma.T, n_array) - np.dot(np.array([1 - n_C, alpha + 2*y, beta + y +2*z, n_N2]).T, 1))
        x = solve(Jacobian, B)
        delta_n = x[:S]
        delta_psi = x[S:]
        
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
        n6_values.append(n_array[5])
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
        line6.set_ydata(n6_values)
        line7.set_xdata(count_values)
        line7.set_ydata(omega_val)
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
    
    
        
    db = NasaDatabase("nasa-cea") # database for thermodynamic properties
    gases = GaseousPhase("CH4 O2 CO2 CO H2O H2") # product gases
    system = ChemicalSystem(db, gases)
    solver = EquilibriumSolver(system)



end = datetime.datetime.now()
print('\nCalculation time:', end - start) 
