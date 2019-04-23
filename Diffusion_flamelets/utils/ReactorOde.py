'''
This demonstrates an approach for solving problems where Cantera's reactor
network model cannot be configured to describe the system in question. Here,
Cantera is used for evaluating thermodynamic properties and kinetic rates while
an external ODE solver is used to integrate the resulting equations. In this
case, the SciPy wrapper for VODE is used, which uses the same variable-order BDF
methods as the Sundials CVODES solver used by Cantera.
'''

import numpy as np


class ReactorOde(object):
    def __init__(self, gas):
        # Parameters of the ODE system and auxiliary data are stored in the
        # ReactorOde object.
        self.gas = gas
        self.P = gas.P

    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """

        # State vector is [T, Y_1, Y_2, ... Y_K]
        self.gas.set_unnormalized_mass_fractions(y[1:])
        self.gas.TP = y[0], self.P
        rho = self.gas.density

        wdot = self.gas.net_production_rates
        dTdt = - (np.dot(self.gas.partial_molar_enthalpies, wdot) /
                  (rho * self.gas.cp))
        dYdt = wdot * self.gas.molecular_weights / rho

        return np.hstack((dTdt, dYdt))