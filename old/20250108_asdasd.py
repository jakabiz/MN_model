# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:33:31 2025

@author: jakab
"""

import numpy as np
import pandas as pd
from pandas import DataFrame
import cantera as ct
from scipy.optimize import root, minimize
from pint import UnitRegistry
import datetime
import os

ureg = UnitRegistry()
Q_ = ureg.Quantity

def to_si(quant):
    """Converts a Pint Quantity to magnitude at base SI units."""
    return quant.to_base_units().magnitude

start = datetime.datetime.now()  # Stopwatch the time of loop calculation

output_file = 'gas_composition_nonstehiometric_model.txt'

# Check if the file exists; if so, remove it
if os.path.exists(output_file):
    os.remove(output_file)

# Create an empty DataFrame with the required columns
df = pd.DataFrame(columns=['T [K]', 'ER', 'H2', 'CO', 'CO2', 'CH4', 'N2'])
df.to_csv(output_file, index=False)  # Write the header

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

biomass = pd.read_csv('biomass.txt', sep=",")
all_results = []

cases = biomass.shape[0]
#for i in np.arange(0, cases, 1):
i=0
w_C = biomass['w_C'][i] / 100
w_H = biomass['w_H'][i] / 100
w_O = biomass['w_O'][i] / 100
w_N = biomass['w_N'][i] / 100
w_S = 0
MC = biomass['MC'][i] / 100  # Moisture content
ER = biomass['ER'][i]    # Equivalence ratio
phi_N = biomass['phi_N'][i] # volume fraction of N2 in gasifying agent
phi_O = biomass['phi_O'][i] # volume fraction of O2 in gasifying agent

M_C, M_H, M_O, M_N = 12.01, 1.008, 16.00, 14.01  # g/mol
# Calculate atomic ratios
alpha = w_H * M_C / (w_C * M_H)
beta = w_O * M_C / (w_C * M_O)
gamma = w_N * M_C / (w_C * M_N)

y = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)
z = ER * (1 + alpha / 4 - beta / 2)

# Calculate initial moles of N2 (inert gas)
n_N2 = 0.5 * gamma + z * (phi_N / phi_O)
n_C = 0.05

# Products: CO, H2, CH4, CO2, H2O
components = ['CO', 'H2', 'CH4', 'CO2', 'H2O']
elemental_comp = np.array([
    [1, 0, 1, 1, 0],  # Carbon
    [0, 2, 4, 0, 2],  # Hydrogen
    [1, 0, 0, 2, 1],  # Oxygen
])

# Initial molar amounts of each element
initial_moles_elements = np.array([
    1.0 - n_C,  # Carbon
    alpha + 2 * y,  # Hydrogen
    beta + y + 2*z,  # Oxygen (including O2 from air and moisture)
])

T = 900
dH=10000 # initial guess of energy balance so the while loop starts
while np.absolute(dH) > 100:
    temperature = Q_(T, 'K')  # Set temperature
    pressure = Q_(1, 'atm')  # Set pressure
    # ambiental conditions
    T_0 = 25 + 273 # K - initail temperature
    p_tot = 1.013e5 # Pa - gasification pressure
    R = 8.314 # kJ/kmolK - gas constant
    temp = temperature.magnitude
    
    
    def lagrange_system_minimize(x, temperature, pressure, components, gas, elemental_comp, initial_moles_elements):
        """System of equations for Lagrange multiplier approach for use with 'minimize'."""
        moles = x[:len(components)]
        multipliers = x[len(components):]
        
        mole_fractions = moles / np.sum(moles)
        
        # Standard-state Gibbs free energy
        gibbs = np.zeros(len(components))
        for idx, comp in enumerate(components):
            gas.TPX = to_si(temperature), to_si(Q_(1, 'atm')), f'{comp}:1.0'
            gibbs[idx] = gas.gibbs_mole
        gibbs *= Q_('J/kmol')
        
        gas_constant = Q_(ct.gas_constant, 'J/(kmol*K)')
        chemical_potentials = gibbs + gas_constant * temperature * np.log(mole_fractions * pressure / Q_(1.0, 'atm'))
        
        # Elemental balance equations
        moles_elements = np.dot(elemental_comp, moles)
        element_equations = moles_elements - initial_moles_elements
        
        # Lagrange multiplier equations
        multiplier_equations = to_si(chemical_potentials + np.dot(multipliers, elemental_comp) * Q_('J/kmol'))
        
        return np.concatenate((element_equations, multiplier_equations))

    # Solve for equilibrium
    gas = ct.Solution('gri30.yaml')
    
    # Initial guess
    x0 = [0.7, 0.7, 0.05, 0.2, 0.1] + [1e3] * elemental_comp.shape[0]
    
    # Bounds for each mole fraction (positive)
    bounds = [(0, None)] * len(x0)  # Ensures mole fractions are non-negative
    
    # Solve with 'minimize'
    sol = minimize(
        lambda x: np.sum(lagrange_system_minimize(x, temperature, pressure, components, gas, elemental_comp, initial_moles_elements)**2), 
        x0, 
        bounds=bounds,
        options={'maxiter': 500}
    )
    print('case ', i+1)
    print('T: ', T, ' K')
    print('Root-finding success:', sol.success)
    print('Equilibrium mole fractions:')
    if sol.success:
        moles = sol.x[:len(components)]
        n_CO, n_H2, n_CH4, n_CO2, n_H2O = sol.x[:len(components)]
        n_tot = n_CH4 + n_CO + n_CO2 + n_H2 + n_N2
        mole_fractions = moles / n_tot
        phi_CO, phi_H2, phi_CH4, phi_CO2, phi_H2O = mole_fractions
        phi_N2 = n_N2/n_tot
        for idx, comp in enumerate(components):
            print(f'{comp}: {mole_fractions[idx]:.3f}')
        print('N2: ', np.round(n_N2/n_tot, 3))
        
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
        dH_reaktanti = dH_b + y*(dH_0_H2O_g)# + dH_vap)
        dH_produkti = n_H2*dH_H2 + n_CO*dH_CO + n_CO2*dH_CO2 + n_CH4*dH_CH4 + n_H2O*dH_H2O_g + n_N2*dH_N2 + n_C*dH_C 
        dH = dH_produkti - dH_reaktanti
        # if dH is greater than 100, then increase temperature
        print('dH: ', dH)
    T = T + 1
# Collect results for this case
result = [T, ER, phi_H2 * 100, phi_CO * 100, phi_CO2 * 100, phi_CH4 * 100, phi_N2 * 100]
all_results.append(result)
    
# Create DataFrame from results
df_results = pd.DataFrame(all_results, columns=['T [K]', 'ER', 'H2', 'CO', 'CO2', 'CH4', 'N2'])

# Append results to file
df_results.to_csv(output_file, mode='a', header=False, index=False)

end = datetime.datetime.now()
print('\nCalculation time:', end - start)