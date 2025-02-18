# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 12:16:02 2024

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


plt.rcParams.update({
    'font.size': 10,        # Font size for text
    'axes.titlesize': 10,   # Font size for axes titles
    'axes.labelsize': 10,   # Font size for x and y labels
    'xtick.labelsize': 10,  # Font size for x-axis tick labels
    'ytick.labelsize': 10,  # Font size for y-axis tick labels
    'legend.fontsize': 10,  # Font size for legend
    'lines.linewidth': 2,   # Default line width
    'lines.markersize': 5    # Default marker size for all plots
})

start = datetime.datetime.now() # stopwatch the time of loop calculation

output_file = 'gas_composition_non_stehiometric_model.txt'

# Check if the file exists; if so, remove it
if os.path.exists(output_file):
    os.remove(output_file)

# Create an empty DataFrame with the required columns
df = pd.DataFrame(columns=['T [K]', 'ER', 'H2', 'CO', 'CO2', 'CH4', 'N2'])
df.to_csv(output_file, index=False)  # Write the header

# Constants
M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007  # kg/kmol
R = 8.314  # kJ/kmolK
T_0 = 25 + 273.15 # STPC - K

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

# Standard Gibbs free energies of formation (kJ/kmol)
delta_Gf = np.array([
    -137163,   # CO
    0,         # H2
    -50768,    # CH4
    -394389,   # CO2
    -228582,   # H2O
    0,         # N2
    0          # Solid carbon residual
])

# Atomic composition matrix (C, H, O, N)
a_matrix = np.array([
    [1, 0, 1, 0],  # CO
    [0, 2, 0, 0],  # H2
    [1, 4, 0, 0],  # CH4
    [1, 0, 2, 0],  # CO2
    [0, 2, 1, 0],  # H2O
    [0, 0, 0, 2],  # N2
    [1, 0, 0, 0]   # Solid carbon residual
])

# Read biomass data from file
biomass = pd.read_csv('biomass.txt', sep=",")
model_Piks = pd.read_csv('gas_composition_basic_model_Piks.txt', sep=",", header=0)
stehiometric_model = pd.read_csv('gas_composition_stehiometric_model.txt', sep=",", header=0)

all_results = []

cases = 3#biomass.shape[0]

for i in np.arange(0, cases, 1):
    w_C = biomass['w_C'][i] / 100
    w_H = biomass['w_H'][i] / 100
    w_O = biomass['w_O'][i] / 100
    w_N = biomass['w_N'][i] / 100
    w_S = 0
    MC = biomass['MC'][i] / 100  # Moisture content
    ER = biomass['lambda'][i]    # Equivalence ratio
    
    # Calculate atomic ratios
    alpha = w_H * M_C / (w_C * M_H)
    beta = w_O * M_C / (w_C * M_O)
    gamma = w_N * M_C / (w_C * M_N)
    x = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)
    z = ER * (1 + alpha / 4 - beta / 2)
    lambda_air = 3.76
    
    # Total atoms
    def total_atoms(z, x):
        A_C = 1  # Carbon from biomass
        A_H = alpha + 2 * x
        A_O = beta + 2 * z + x
        A_N = (gamma + 2 * z * lambda_air) / 2
        return A_C, A_H, A_O, A_N
    
    # Gibbs free energy function
    def gibbs_energy(n):
        n_total = np.sum(n) + 1e-10
        y = n / n_total + 1e-10
        G = np.sum(n * delta_Gf) + np.sum(n * R * temp * np.log(y))
        return G
    
    # Constraints
    def create_constraints(z, x):
        A_C, A_H, A_O, A_N = total_atoms(z, x)
        return [
            {'type': 'eq', 'fun': lambda n: np.dot(a_matrix[:, 0], n) - A_C},  # Carbon balance
            {'type': 'eq', 'fun': lambda n: np.dot(a_matrix[:, 1], n) - A_H},  # Hydrogen balance
            {'type': 'eq', 'fun': lambda n: np.dot(a_matrix[:, 2], n) - A_O},  # Oxygen balance
            {'type': 'eq', 'fun': lambda n: np.dot(a_matrix[:, 3], n) - A_N},  # Nitrogen balance
            {'type': 'eq', 'fun': lambda n: n[6] - 0.05 * A_C}                 # Solid carbon residual
        ]
    
    # Initial guesses and bounds
    #initial_guess = [0.2, 0.2, 0.1, 0.3, 0.3, 0.2, 0.05]
    A_C, A_H, A_O, A_N = total_atoms(z, x)
    initial_guess = [
        0.5 * A_C,  # CO
        0.2 * A_H,  # H2
        0.05 * A_C, # CH4
        0.3 * A_C,  # CO2
        0.2 * A_H,  # H2O
        0.8 * A_N,  # N2
        0.05 * A_C  # Solid carbon residual
    ]
    bounds = [(1e-10, None)] * 7
    
    # Perform optimization
    dH=1000 # initial guess of energy balance so the while loop starts
    temp = 873 # initial temperature
    while np.absolute(dH) > 100:
        constraints = create_constraints(z, x)

        result = minimize(
            gibbs_energy,
            initial_guess,
            bounds=bounds,
            constraints=constraints,
            method='SLSQP'
        )
        
        # Output results
        if result.success:
            n_opt = result.x  # Optimal mole numbers for each component
            total_moles = np.sum(n_opt[:6])  # Total moles for gas components, excluding solid carbon
        
            # Calculate mole fractions for each gas component
            n_CO, n_H2, n_CH4, n_CO2, n_H2O, n_N2, n_C = n_opt  
        
        #else:
           # print("Optimization failed:", result.message)
            
        n_tot = n_CH4 + n_CO + n_CO2 + n_H2O + n_H2 + n_N2 - n_H2O
        
        # volume fractions in produced gas
        phi_CH4 = n_CH4/n_tot
        phi_CO = n_CO/n_tot
        phi_CO2 = n_CO2/n_tot
        phi_H2 = n_H2/n_tot
        phi_H2O = n_H2O/n_tot
        phi_N2 = n_N2/n_tot
        
        # energy balance
        HHV_b = (0.3491*w_C + 1.1783*w_H + 0.1005*w_S - 0.1034*w_O - 0.0151*w_N)*100*1000*(M_C + alpha*M_H + beta*M_O + gamma*M_N) # kJ/kmol - higher heating value of biomass (multipled by 100, so mass fractions are in % - empirical equation)
        
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
        dH_reaktanti = dH_b + x*(dH_0_H2O_g)# + dH_vap)
        dH_produkti = n_H2*dH_H2 + n_CO*dH_CO + n_CO2*dH_CO2 + n_CH4*dH_CH4 + n_H2O*dH_H2O_g + n_N2*dH_N2 + n_C*dH_C 
        dH = dH_produkti - dH_reaktanti
        # if dH is greater than 100, then increase temperature
        temp = temp + 1
        
    tabelca = {
        '': ['T [K]', 'ER [/]', 'H2 [%]', 'CO [%]', 'CO2 [%]', 'CH4 [%]', 'N2 [%]'],
        'non-stehiometric model': [temp, ER, phi_H2*100, phi_CO*100, phi_CO2*100, phi_CH4*100, phi_N2*100],
        'stehiometric model': [stehiometric_model['T [K]'][i], stehiometric_model['ER'][i], stehiometric_model['H2'][i], stehiometric_model['CO'][i], stehiometric_model['CO2'][i], stehiometric_model['CH4'][i], stehiometric_model['N2'][i]],
        'Pikš basic model': [np.nan, ER, model_Piks['H2'][i], model_Piks['CO'][i], model_Piks['CO2'][i], model_Piks['CH4'][i], model_Piks['N2'][i]]
        }
    print('case ' + str(i+1), '\n', pd.DataFrame(tabelca).round(2).T.to_string(header=False), '\n')
    
    # Collect results for this case
    result = [temp, ER, phi_H2 * 100, phi_CO * 100, phi_CO2 * 100, phi_CH4 * 100, phi_N2 * 100]
    all_results.append(result) 
    
# Create DataFrame from results
df_results = pd.DataFrame(all_results, columns=['T [K]', 'ER', 'H2', 'CO', 'CO2', 'CH4', 'N2'])

# Append results to file
df_results.to_csv(output_file, mode='a', header=False, index=False)
    
end = datetime.datetime.now()
print('\nCalculation time:', end - start)
