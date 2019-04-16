"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

from: CANTERA
"""

import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Simulation parameters

mechanism = 'gri30.cti' #'utils/lu19.cti' #'gri30.cti'

p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]

yF = 0.055  # fuel
yO2 = 0.232  # oxidizer

for yF in [0.055,0.0570,0.06,0.07]:
    yOx = (1 - yF) * 0.232
    yN2 = (1 - yF) * 0.768

    # reactants in mass fraction
    reactants = ('CH4:%f O2:%f N2:%f' % (yF, yOx, yN2))

    # Calculate Bilger Mixture fraction
    yC = yF * 0.748753
    yH = yF * 0.251247
    yO = yOx
    mixtureFractionBilger = (2 * yC / 12.01 + yH / (2 * 1.01) - (yO - 0.232) / 16) / \
                            (2 * (0.75 / 12.01) + 0.25 / (2 * 1.01) + 0.232 / 16)

    print('Bilger Mixture Fraction is: %f' % mixtureFractionBilger)

    width = 0.01  # m
    initial_grid = np.linspace(0, width, 150)

    loglevel = 1  # amount of diagnostic output (0 to 8)

    # IdealGasMix object used to compute mixture properties, set to the state of the
    # upstream fuel-air mixture
    gas = ct.Solution(mechanism)
    gas.TPY = Tin, p, reactants

    # Set up flame object
    f = ct.FreeFlame(gas, initial_grid)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.06)
    f.show_solution()

    # Solve with mixture-averaged transport model
    f.transport_model = 'Mix'
    f.solve(loglevel=loglevel, auto=True)

    # Solve with the energy equation enabled
    f.save('lu19_adiabatic.xml', 'mix', 'solution with mixture-averaged transport')
    f.show_solution()
    print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))

    # Solve with multi-component transport properties
    f.transport_model = 'Multi'
    f.solve(loglevel)  # don't use 'auto' on subsequent solves
    f.show_solution()
    print('multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))
    # f.save('lu19_adiabatic.xml','multi', 'solution with multicomponent transport')

    # write the velocity, temperature, density, and mole fractions to a CSV file
    # f.write_csv('lu19_adiabatic.csv', quiet=False)

    # plt.plot(f.grid,f.Y[gas.species_index('CH4')])
    plt.plot(f.grid, f.T)
    plt.show(block=False)