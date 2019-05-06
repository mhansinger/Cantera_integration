# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:44:34 2015
Write out all flame quantities of interest
@author: julian
"""
import csv

def write_full_csv(flame, filename, species='Y', quiet=True):
        """
        Write the velocity, temperature, density, and species profiles
        to a CSV file.

        :param filename:
            Output file name
        :param species:
            Attribute to use obtaining species profiles, e.g. ``X`` for
            mole fractions or ``Y`` for mass fractions.
        """

        z = flame.grid
        T = flame.T
        u = flame.u
        HeatRelease = flame.heat_release_rate
        ProdRate = flame.gas.species_names

        for i, val in enumerate(flame.gas.species_names):
            ProdRate[i] = 'ProdRate-%s' %flame.gas.species_names[i]

 

        csvfile = open(filename, 'w')
        writer = csv.writer(csvfile)
        writer.writerow(['z', 'u',
                         'T', 'rho', 'viscosity', 'thermalConductivity', 'cp', 'HeatRelease'] + flame.gas.species_names + ProdRate)
        for n in range(flame.flame.n_points):
            flame.set_gas_state(n)

            writer.writerow([z[n], u[n], T[n], flame.gas.density, flame.gas.viscosity, flame.gas.thermal_conductivity, flame.gas.cp,  HeatRelease[n] ] +
                            list(getattr(flame.gas, species)) + list(flame.gas.net_production_rates))
        csvfile.close()
        if not quiet:
            print("Solution saved to '{0}'.".format(filename))