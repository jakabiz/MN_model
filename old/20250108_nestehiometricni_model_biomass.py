# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:33:31 2025

@author: jakab
"""

# Import necessary libraries
import numpy as np
import pandas as pd
from pandas import DataFrame
import cantera as ct
from scipy.optimize import root, minimize
from pint import UnitRegistry
import datetime
import os

# Initialize unit registry for unit conversions using Pint
ureg = UnitRegistry()
Q_ = ureg.Quantity

# Function to convert Pint Quantity to base SI units
def to_si(quant):
    """Converts a Pint Quantity to magnitude at base SI units."""
    return quant.to_base_units().magnitude

# Start a stopwatch to track the calculation time
start = datetime.datetime.now()

# Output file for storing gas composition results
output_file = 'gas_composition_nonstehiometric_model.txt'

# Check if the output file already exists; if so, remove it to start fresh
if os.path.exists(output_file):
    os.remove(output_file)

# Create an empty DataFrame with specified columns to store results
df = pd.DataFrame(columns=['T [K]', 'ER', 'H2', 'CO', 'CO2', 'CH4', 'N2'])
df.to_csv(output_file, index=False)  # Write the header to the output file

# Define standard enthalpy of formation (in kJ/kmol) for different species
dH_0_CH4 = -74873  # CH4
dH_0_H2 = 0  # H2
dH_0_CO = -110527  # CO
dH_0_CO2 = -393522  # CO2
dH_0_N2 = 0  # N2 (inert)
dH_0_H2O_l = -285830  # H2O (liquid)
dH_0_H2O_g = -241826  # H2O (gas)
dH_0_C = 0  # Carbon (elemental)
dH_vap = 40.7e3  # Enthalpy of water vaporization (kJ/kmol)

# Read biomass data from the file 'biomass.txt' into a DataFrame
biomass = pd.read_csv('biomass.txt', sep=",")
all_results = []  # List to collect results for each biomass case

# Number of biomass cases to iterate over
cases = biomass.shape[0]

# Example: iterate over biomass case i (currently set to 3)
i = 0
# Extract the respective biomass properties for this case
w_C = biomass['w_C'][i] / 100  # Weight fraction of Carbon
w_H = biomass['w_H'][i] / 100  # Weight fraction of Hydrogen
w_O = biomass['w_O'][i] / 100  # Weight fraction of Oxygen
w_N = biomass['w_N'][i] / 100  # Weight fraction of Nitrogen
w_S = 0  # Sulfur is assumed to be 0 (not given in biomass data)
MC = biomass['MC'][i] / 100  # Moisture content (fraction)
ER = biomass['ER'][i]  # Equivalence ratio
phi_N = biomass['phi_N'][i]  # Volume fraction of Nitrogen in the gasifying agent
phi_O = biomass['phi_O'][i]  # Volume fraction of Oxygen in the gasifying agent

# Molar masses for Carbon, Hydrogen, Oxygen, and Nitrogen (in g/mol)
M_C, M_H, M_O, M_N = 12.01, 1.008, 16.00, 14.01

# Calculate atomic ratios for H/C, O/C, and N/C based on biomass composition
alpha = w_H * M_C / (w_C * M_H)
beta = w_O * M_C / (w_C * M_O)
gamma = w_N * M_C / (w_C * M_N)

# Calculate moisture contribution (y) and equivalence ratio-related adjustment (z)
y = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)
z = ER * (1 + alpha / 4 - beta / 2)

# Calculate initial moles of Nitrogen (n_N2) and Carbon (n_C)
n_N2 = 0.5 * gamma + z * (phi_N / phi_O)
n_C = 0.05  # Assumed initial mole fraction for Carbon

# List of products formed during biomass gasification
components = ['CO', 'H2', 'CH4', 'CO2', 'H2O']
# Elemental composition of the products (C, H, O) for each component
elemental_comp = np.array([
    [1, 0, 1, 1, 0],  # Carbon
    [0, 2, 4, 0, 2],  # Hydrogen
    [1, 0, 0, 2, 1],  # Oxygen
])

# Initial moles of each element based on the atomic ratios and biomass composition
initial_moles_elements = np.array([
    1.0 - n_C,  # Carbon
    alpha + 2 * y,  # Hydrogen
    beta + y + 2 * z,  # Oxygen (including O2 from air and moisture)
])

# Initial guess for temperature and energy balance
T = 900  # Starting temperature in K
dH = 10000  # Initial guess for energy balance (large enough to start the while loop)
# While loop to adjust temperature until energy balance is within a threshold
while np.absolute(dH) > 100:
    # Define temperature and pressure using Pint Quantity objects
    temperature = Q_(T, 'K')  # Temperature in K
    pressure = Q_(1, 'atm')  # Pressure in atm

    # Environmental conditions and constants
    T_0 = 25 + 273  # Initial temperature in K (25°C converted to K)
    p_tot = 1.013e5  # Total pressure in Pa (1 atm)
    R = 8.314  # Gas constant in kJ/kmol·K (universal)
    temp = temperature.magnitude  # Get temperature value as a float

    # Define the Lagrange system of equations for solving equilibrium
    def lagrange_system(x, temperature, pressure, components, gas, elemental_comp, initial_moles_elements):
        """System of equations for Lagrange multiplier approach."""
        moles = x[:len(components)]  # Extract mole amounts for each component
        multipliers = x[len(components):]  # Extract Lagrange multipliers

        # Calculate mole fractions for each component
        mole_fractions = moles / np.sum(moles)
        
        # Initialize an array for Gibbs free energy for each component
        gibbs = np.zeros(len(components))
        for idx, comp in enumerate(components):
            gas.TPX = to_si(temperature), to_si(Q_(1, 'atm')), f'{comp}:1.0'
            gibbs[idx] = gas.gibbs_mole  # Get Gibbs free energy for component
        gibbs *= Q_('J/kmol')  # Convert Gibbs energy to J/kmol
        
        # Calculate chemical potentials for each component (Gibbs + RT ln(P))
        gas_constant = Q_(ct.gas_constant, 'J/(kmol*K)')  # Gas constant in J/(kmol·K)
        chemical_potentials = gibbs + gas_constant * temperature * np.log(mole_fractions * pressure / Q_(1.0, 'atm'))
        
        # Calculate the elemental balance equations
        moles_elements = np.dot(elemental_comp, moles)  # Elemental balances for C, H, O
        element_equations = moles_elements - initial_moles_elements  # Difference from initial amounts
        
        # Define the Lagrange multiplier equations
        multiplier_equations = to_si(chemical_potentials + np.dot(multipliers, elemental_comp) * Q_('J/kmol'))
        
        # Return all the equations (elemental balances + multiplier equations)
        return np.concatenate((element_equations, multiplier_equations)) 

    # Set up the Cantera gas solution
    gas = ct.Solution('gri30.yaml')

    # Initial guess for mole fractions of each component (arbitrary values to start)
    x0 = [0.7, 0.7, 0.05, 0.2, 0.1] + [1e3] * elemental_comp.shape[0]
    
    # Define bounds for the solution variables (moles must be >= 0)
    bounds = [(0, None)]  # This means x must be >= 0
    
    # Solve the system of equations using the 'root' function from scipy
    sol = root(
        lagrange_system, x0, method='lm',  # Use Levenberg-Marquardt method
        args=(temperature, pressure, components, gas, elemental_comp, initial_moles_elements), 
        options={'maxiter': 500}  # Set maximum iterations
    )
    
    # Output the case results
    print('case ', i+1)
    print('T: ', T, ' K')  # Print the temperature for this iteration
    print('Root-finding success:', sol.success)  # Check if root-finding was successful
    print('Equilibrium mole fractions:')
    
    # If the solution was successful, extract and print the mole fractions
    if sol.success:
        moles = sol.x[:len(components)]
        n_CO, n_H2, n_CH4, n_CO2, n_H2O = sol.x[:len(components)]
        n_tot = n_CH4 + n_CO + n_CO2 + n_H2 + n_N2  # Total moles of all components
        mole_fractions = moles / n_tot  # Calculate mole fractions
        
        # Print the mole fractions of each component
        phi_CO, phi_H2, phi_CH4, phi_CO2, phi_H2O = mole_fractions
        phi_N2 = n_N2 / n_tot
        for idx, comp in enumerate(components):
            print(f'{comp}: {mole_fractions[idx]:.3f}')
        print('N2: ', np.round(n_N2/n_tot, 3))  # Print N2 mole fraction
        
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
        print('dH: ', dH)
        # if dH is greater than 100, then increase temperature
        temp = temp + 1
        
    T += 1  # Increase temperature to try again if energy balance is not within tolerance

# Collect results for this case
result = [T, ER, phi_H2 * 100, phi_CO * 100, phi_CO2 * 100, phi_CH4 * 100, phi_N2 * 100]
all_results.append(result)
    
# Create DataFrame from results
df_results = pd.DataFrame(all_results, columns=['T [K]', 'ER', 'H2', 'CO', 'CO2', 'CH4', 'N2'])

# Append results to file
df_results.to_csv(output_file, mode='a', header=False, index=False)

end = datetime.datetime.now()
print('\nCalculation time:', end - start)