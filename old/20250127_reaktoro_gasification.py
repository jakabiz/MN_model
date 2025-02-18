# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:37:47 2025

@author: jakab
"""

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
from reaktoro import *
from reaktplot import *

start = datetime.datetime.now() # stopwatch the time of loop calculation

db = NasaDatabase("nasa-cea")
'''
# Inputs -----------------------------------------------------------------------------
massC              = 49.30  # g/mol
massH              = 5.5    # g/mol
massO              = 45.2   # g/mol
HHV                = 18.933 # kJ/g
Q                  = 0.05   # heat loss (%) with respect to HHV
fAirMin            = 0.70   # minimum value for fAir (fuel to air mass ratio)
fAirDelta          = 0.1    # step-size for fAir values
fAirNum            = 30     # number of steps for fAir
output_all_species = False   # if false, output only a few species
MC = 0.1
# ------------------------------------------------------------------------------------

nC = massC / 12.011
nH = massH / 1.00797
nO = massO / 15.9994

a = nH / nC
b = nO / nC
'''
# Constants
M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007  # kg/kmol
R = 8.314  # kJ/kmolK
T_0 = 25 + 273.15 # STPC - K
p_0 = 1.013e5 # Pa
# equilibrium model for syngas composition calculation based on free Gibbs energy minimization
w_C = 0.521
w_H = 0.062
w_O = 0.412
w_N = 0.005
w_S = 0
MC = 0.10  # Moisture content

# Calculate atomic ratios
alpha = w_H * M_C / (w_C * M_H)
beta = w_O * M_C / (w_C * M_O)
gamma = w_N * M_C / (w_C * M_N)
y = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)
M_b = (M_C + alpha*M_H + beta*M_O + gamma*M_N) # kg/kmol
HHV_b = (0.3491*w_C + 1.1783*w_H + 0.1005*w_S - 0.1034*w_O - 0.0151*w_N)*100*1000*M_b # kJ/kmol - higher heating value of biomass (multipled by 100, so mass fractions are in % - empirical equation)


h0CO2 = -393.522  # kJ/g
h0H2O = -285.830   # kJ/g


formula = ChemicalFormula(f"C{0.95}H{alpha}O{beta}")

h0fuel = HHV_b/1000 + h0CO2 + 0.5*h0H2O*alpha  # in kJ/kmol
heatDuty = HHV_b/1000 + h0CO2 + 0.5*h0H2O*alpha  # in kJ/kmol

stmodelparams = StandardThermoModelParamsConstant()
stmodelparams.G0 = 1.0e+3
stmodel = StandardThermoModelConstant(stmodelparams)

species = Species()
species = species.withName("Fuel")
species = species.withElements(formula.elements())
species = species.withAggregateState(AggregateState.CondensedPhase)
species = species.withStandardThermoModel(stmodel)

db.addSpecies(species)

gases = GaseousPhase(speciate("C H O N"))
condensedphases = CondensedPhases("Fuel H2O(l)")

system = ChemicalSystem(db, condensedphases, gases)

specs = EquilibriumSpecs(system)
specs.pressure()
specs.enthalpy()

solver = EquilibriumSolver(specs)

table = Table()

ERs = np.arange(0.15,0.41, 0.01)
phi_Os = np.arange(0.21, 0.96, 0.01)

for ER in ERs:
    for phi_O in phi_Os:
        phi_N = 1 - phi_O
    
        state = ChemicalState(system)
        state.temperature(T_0, "K")
        state.pressure(p_0, "Pa")
        state.set("Fuel", 1.0, "mol")
        state.set('H2O', y, 'mol')
        state.set("O2", ER*(1 + alpha/4 - beta/2), "mol")
        state.set("N2", ER*(1 + alpha/4 - beta/2)*phi_N/phi_O, "mol")
    
        props = ChemicalProps(state)
        conditions = EquilibriumConditions(specs)
    
        conditions.pressure(props.pressure())
        conditions.enthalpy(heatDuty, "kJ")
    
        conditions.setLowerBoundTemperature(500, "K")
        conditions.setUpperBoundTemperature(3500.0, "K")
    
        result = solver.solve(state, conditions)
    
        assert result.succeeded()
    
        table.column("ER") << ER
        table.column("phi_O") << phi_O
        table.column("Temperature") << state.temperature()
        table.column(f"n(Gases)") << state.props().phaseProps("GaseousPhase").amount() - state.speciesAmount("H2O")
        
        #table.column("n(C(gr))")    << state.speciesAmount("C(gr)")
        table.column("n(CH4)")      << state.speciesAmount("CH4")
        table.column("n(CO)")       << state.speciesAmount("CO")
        table.column("n(CO2)")      << state.speciesAmount("CO2")
        table.column("n(H2)")       << state.speciesAmount("H2")
        table.column("n(H2O)")      << state.speciesAmount("H2O")
        table.column("n(N2)")       << state.speciesAmount("N2")
    
table.save("reaktoro_oxygen_frac.txt")


end = datetime.datetime.now()
print('\nCalculation time:', end - start) 