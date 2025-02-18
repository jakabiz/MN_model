# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 11:31:48 2025

@author: Jaka Bizjak, IJS CEU

equilibrium non-stoichiometric model of biomass gasification, summarized by 
Computational Thermodynamics by Kyle Niemeyer 
https://kyleniemeyer.github.io/computational-thermo/content/intro.html
"""

import numpy as np
import cantera as ct
import pandas as pd
from scipy.optimize import root, minimize

from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

# for convenience:
def to_si(quant):
    '''Converts a Pint Quantity to magnitude at base SI units.
    '''
    return quant.to_base_units().magnitude

# Known information

components = ['CO', 'O2', 'CO2']
moles_initial = np.array([1.0, 0.5, 0.0])

# Elemental makeup of components
elemental_comp = np.array([
    [1, 0, 1], # carbon
    [1, 2, 2], # oxygen
    ])

temperature = Q_(2500, 'K')
pressures = [1, 10] * Q_('atm')

def lagrange_system(x, temperature, pressure, components, 
                    gas, elemental_comp, moles_initial):
    '''System of equations for Lagrange multiplier approach.
    '''
    moles = np.array([x[0], x[1], x[2]])
    multipliers = np.array([x[3], x[4]])
    
    mole_fractions = moles / np.sum(moles)
    
    # get standard-state Gibbs free energy of each component
    gibbs = np.zeros(len(components))
    for idx, comp in enumerate(components):
        gas.TPX = (
            to_si(temperature), to_si(Q_(1, 'atm')),
            f'{comp}:1.0'
            )
        gibbs[idx] = gas.gibbs_mole
        
    gibbs *= Q_('J/kmol')
    
    gas_constant = Q_(ct.gas_constant, 'J/(kmol*K)')
    chemical_potentials = (
        gibbs + gas_constant * temperature * np.log(
            mole_fractions * pressure / Q_(1.0, 'atm')
            )
        )

    # initial molar amounts of each element
    initial_moles_elements = np.dot(elemental_comp, moles_initial)
    moles_elements = np.dot(elemental_comp, moles)
    
    # We can take advantage of element-wise operations with these arrays,
    # and concisely evaluate all the equations
    element_equations = moles_elements - initial_moles_elements
    multiplier_equations = to_si(
        chemical_potentials + 
        np.dot(multipliers, elemental_comp) * Q_('J/kmol')
        )
    
    # Return the set of equations joined together
    return np.concatenate((element_equations, multiplier_equations))

# Solve at first pressure

pressure = pressures[0]
gas = ct.Solution('gri30.yaml')

# initial guesses
x0 = [1.0, 1.0, 1.0, 1e6, 1e6]

sol = root(
    lagrange_system, x0, method='lm',
    args=(temperature, pressure, components, gas, elemental_comp, moles_initial)
    )

print('Root-finding algorithm success: ', sol.success)
print(f'Function evaluation (should be small): {sol.fun}')
print('Number of function evaluations: ', sol.nfev)
print()

moles = sol.x[0:3]
mole_fractions = moles / np.sum(moles)
print(f'Mole fractions at {pressure: .1f}:')
for idx, comp in enumerate(components):
    print(f'{comp}: {mole_fractions[idx]: .3f}')
    
gas = ct.Solution("gri30.yaml")  # GRI-Mech 3.0 database   
gas.TPX = 2500, ct.one_atm, "CO2:1"  # Set temperature, pressure, and composition

# Get the Gibbs free energy in J/kmol
gibbs_free_energy = gas.gibbs_mole

# Output Gibbs free energy for CO2
print(gibbs_free_energy)
