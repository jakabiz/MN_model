# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:49:01 2024

@author: Jaka Bizjak, IJS CEU

Ne-stehiometrični ravnotežni model uplinjanja biomase v sotočnem uplinjevalniku
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
from scipy.optimize import minimize
import datetime
import os

start = datetime.datetime.now() # stopwatch the time of loop calculation

output_file = 'gas_composition_non_stehiometric_model.txt'

# Check if the file exists; if so, remove it
if os.path.exists(output_file):
    os.remove(output_file)

# Create an empty DataFrame with the required columns
df = pd.DataFrame(columns=['T [K]', 'ER', 'CO', 'H2', 'CH4', 'CO2', 'H2O', 'N2'])
df.to_csv(output_file, index=False)  # Write the header

# Constants
M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007  # kg/kmol
R = 8.314  # kJ/kmolK
T_0 = 25 + 273.15 # STPC - K

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

# vector of standard Gibbs free energies of formation (kJ/kmol)
delta_Gf = np.array([
    -137163,   # CO
    0,         # H2
    -50768,    # CH4
    -394389,   # CO2
    -228582,   # H2O
    0          # N2
])

# Matrix of elements (atomic matrix)
#         C, H, O, N
sigma = [[1, 0, 1, 0],  # CO
         [0, 2, 0, 0],  # H2
         [1, 4, 0, 0],  # CH4
         [1, 0, 2, 0],  # CO2
         [0, 2, 1, 0],  # H2O
         [0, 0, 0, 2]]  # N2


# Read biomass data from file
biomass = pd.read_csv('biomass.txt', sep=",")
model_Piks = pd.read_csv('gas_composition_basic_model_Piks.txt', sep=",", header=0)
stehiometric_model = pd.read_csv('gas_composition_stehiometric_model.txt', sep=",", header=0)

all_results = []

cases = 1#biomass.shape[0]

# calculation for every case
for i in np.arange(0, cases, 1):
    w_C = biomass['w_C'][i] / 100
    w_H = biomass['w_H'][i] / 100
    w_O = biomass['w_O'][i] / 100
    w_N = biomass['w_N'][i] / 100
    w_S = 0
    MC = biomass['MC'][i] / 100  # Moisture content
    ER = biomass['lambda'][i]    # Equivalence ratio
    phi_N = biomass['phi_N'][i]  # volume fraction of N in gasifing agent
    phi_O = biomass['phi_O'][i]  # volume fraction of O in gasifing agent
    
    # Calculate atomic ratios
    alpha = w_H * M_C / (w_C * M_H)
    beta = w_O * M_C / (w_C * M_O)
    gamma = w_N * M_C / (w_C * M_N)
    y = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)
    z = ER * (1 + alpha / 4 - beta / 2)
    
    # total number of atoms 
    N_C = 1,                                  # C
    N_H = alpha + 2 * y,                      # H
    N_O = beta + 2 * z + y,                   # O
    N_N = (gamma + 2 * z * phi_N/phi_O) / 2]] # N
