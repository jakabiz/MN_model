# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 11:49:54 2025
updates:
    220250129 - like reaktoro model from https://reaktoro.org/applications/biomass-gasification/biomass-gasification.html

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
import openpyxl
from openpyxl import load_workbook
from scipy.linalg import solve
import reaktoro
from reaktoro import *

plt.close('all')

start = datetime.datetime.now() # stopwatch the time of loop calculation

str_excel = 'MN_model_results_n_C_r_moist.xlsx'

db = NasaDatabase("nasa-cea")  # NASA CEA database

# Constants
M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007  # kg/kmol
R = 8.314  # kJ/kmolK
T_0 = 25 + 273.15 # STPC - K
p_0 = 1.013e5 # Pa

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
dH_vap = -40.7e3 # kJ/kmol - enthalpy of water vaporization

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

# equilibrium model for syngas composition calculation based on free Gibbs energy minimization
w_C = 0.521
w_H = 0.062
w_O = 0.412
w_N = 0.005
w_S = 0
MC = 0.40  # Moisture content

# Calculate atomic ratios
alpha = w_H * M_C / (w_C * M_H)
beta = w_O * M_C / (w_C * M_O)
gamma = w_N * M_C / (w_C * M_N)
y = MC * (M_C + alpha * M_H + beta * M_O + gamma * M_N) / (2 * M_H + M_O)

M_b = (M_C + alpha*M_H + beta*M_O + gamma*M_N) # kg/kmol
#M_b = formula.molarMass() * 1000  # from kg/mol to g/mol - also possible

#n_N2 = 0.5*gamma + z*phi_N/phi_O #do not need it
#n_C = 0.05

HHV_b = (0.3491*w_C + 1.1783*w_H + 0.1005*w_S - 0.1034*w_O - 0.0151*w_N)*100*1000*M_b # kJ/kmol - higher heating value of biomass (multipled by 100, so mass fractions are in % - empirical equation)

formula = ChemicalFormula(f"CH{alpha}O{beta}N{gamma}")

dH_0_b = HHV_b + dH_0_CO2 + 0.5*dH_0_H2O_l*alpha  # in kJ/kmol
heat_duty = dH_0_b + y*(dH_0_H2O_g + dH_vap) # enthalpy of reactants in kJ/kmol

stmodelparams = StandardThermoModelParamsConstant()
stmodelparams.G0 = 1.0e+3
stmodel = StandardThermoModelConstant(stmodelparams)

species = Species()
species = species.withName('biomass')
species = species.withElements(formula.elements())
species = species.withAggregateState(AggregateState.CondensedPhase)
species = species.withStandardThermoModel(stmodel)
db.addSpecies(species)

gases = GaseousPhase(speciate("C H O N")) # defines elements in the system
condensedphases = CondensedPhases("biomass H2O(l) C(gr)") # defines condensed species in the system
system = ChemicalSystem(db, condensedphases, gases)
specs = EquilibriumSpecs(system)
specs.pressure()
specs.enthalpy()

solver = EquilibriumSolver(specs)

ERs = np.arange(0.1, 0.51, 0.01)
#ERs = np.array([0.1]) # use it just to chech calculation time of one ER
phi_Os = np.arange(0.21, 0.96, 0.01)


for ER in ERs:
    all_results = []
    z = ER * (1 + alpha/4 - beta/2)
    for phi_O in phi_Os:
        phi_N = 1 - phi_O # volume fraction of N2 in gasifying agent
        
        state = ChemicalState(system)
        state.temperature(T_0, "K")
        state.pressure(p_0, "Pa")
        state.set("biomass", 1.0, "mol")
        state.set('H2O', y, 'mol')
        state.set("O2", z, "mol")
        state.set("N2", z*phi_N/phi_O, "mol")
    
        props = ChemicalProps(state)
        conditions = EquilibriumConditions(specs)
    
        conditions.pressure(props.pressure())
        conditions.enthalpy(heat_duty/1000, "kJ") # devided by 1000 because amount of biomass, water and air is in mol
    
        conditions.setLowerBoundTemperature(600, "K")
        conditions.setUpperBoundTemperature(3200, "K")
    
        results = solver.solve(state, conditions)
        succeeded =  results.succeeded()
        
        temp = state.temperature()
        # Extract species amounts as float variables
        n_CO = state.speciesAmount("CO")
        n_H2 = state.speciesAmount("H2")
        n_CH4 = state.speciesAmount("CH4")
        n_CO2 = state.speciesAmount("CO2")
        n_H2O = state.speciesAmount("H2O")
        n_N2 = state.speciesAmount("N2")
        n_C_r = state.speciesAmount("C(gr)") # only if you define solid residual as condesed species

        n_tot = state.props().phaseProps("GaseousPhase").amount() - state.speciesAmount("H2O") # total amount of produced gas on dry basis
        
        # volume fractions in produced gas
        phi_CH4 = n_CH4/n_tot
        phi_CO = n_CO/n_tot
        phi_CO2 = n_CO2/n_tot
        phi_H2 = n_H2/n_tot
        phi_H2O = n_H2O/n_tot
        phi_N2 = n_N2/n_tot    
        
        M_g = (M_C + M_O)*phi_CO + 2*M_H*phi_H2 + (M_C + 4*M_H)*phi_CH4 + (M_C + 2*M_O)*phi_CO2 + 2*M_N*phi_N2 # molar mass of produces gas - kg/kmol
        
        LHV_b = HHV_b - (9*w_H)*dH_vap # kJ/kmol - LHV of biomass
        LHV_b_m = LHV_b/(M_b*1000) # MJ/kg
        LHV_g = (phi_CO*12.622 + phi_H2*10.788 + phi_CH4*35.814)*1000*1000*R*T_0/p_0
        LHV_g_V = (phi_CO*12.622 + phi_H2*10.788 + phi_CH4*35.814)*1000
        LHV_g_m = LHV_g/(M_g*1000) # MJ/kg
        CGE = n_tot*M_g*LHV_g_m/(M_b*LHV_b_m)
        CGE_tot = n_tot*M_g*LHV_g_m/(M_b*LHV_b_m + z*2*M_O*(0.000559*(phi_O*100-21)**2)) # including energy consumption for PSA
        # E_PSA is based on fit on experimental data - see PSA_energy_consumption.png and update if necessery
        GY = n_tot*R*T_0*1000/(M_b*p_0)
    
        result = [MC*100, w_C*100, w_H*100, w_O*100, w_N*100, 
                  M_b, LHV_b_m, np.round(ER, 2), np.round(phi_N*100, 2), np.round(phi_O*100, 2), 
                  temp, phi_CO*100, phi_H2*100, phi_CH4*100, phi_CO2*100, phi_H2O*100, phi_N2*100, n_C_r, M_g, LHV_g_V/1000, LHV_g_m, CGE, CGE_tot, GY, succeeded,
                  n_CO, n_H2, n_CH4, n_CO2, n_H2O, n_N2]
        all_results.append(result)
        print('ER = ', np.round(ER,2), ' ; phi_O = ', np.round(phi_O, 2))
        # Create DataFrame from results
    df_results = pd.DataFrame(all_results, columns=['MC', 'w_C', 'w_H', 'w_O', ' w_N', 
                                                        'M_b [kg/kmol]', 'LHV_b [MJ/kg]', 'ER', 'phi_N', 'phi_O', 
                                                        'T [K]', 'CO', 'H2', 'CH4', 'CO2', 'H2O', 'N2', 'n_C_r [kmol]', 'M_g [kg/kmol]', 'LHV_g [MJ/m3]', 'LHV_g [MJ/kg]', 'CGE', 'CGE_tot', 'GY', 'succeeded',
                                                        'n_CO [kmol]', 'n_H2 [kmol]', 'n_CH4 [kmol]', 'n_CO2 [kmol]', 'n_H2O [kmol]', 'n_N2 [kmol]'])
    
    # Replace all numbers (ints and floats) with comma-separated decimals
    #def format_numbers(value):
    #    if isinstance(value, (int, float)):
    #        return str(value).replace('.', ',')  # Replace '.' with ','
    #    return value  # Leave non-numeric values untouched
    #    # Apply formatting to each column of the DataFrame
    #df_results = df_results.apply(lambda col: col.map(format_numbers))
    
    # Write the DataFrame to the sheet, starting from the next row
    with pd.ExcelWriter(str_excel, engine='openpyxl', mode='a',if_sheet_exists='overlay') as writer:
        df_results.to_excel(writer, sheet_name=f"{ER:.2f}", index=False, startrow=0, startcol=0)



end = datetime.datetime.now()
print('\nCalculation time:', end - start) 