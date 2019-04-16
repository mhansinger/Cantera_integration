###############################################################
#
# Constant pressure or Constant-volume reactor,
#             adiabatic kinetics simulation.
#
###############################################################

#import :

import sys
from cantera import *
import numpy as np
import matplotlib.pyplot as plt
import csv


#################################################################
# Prepare your run
#################################################################
#Import the gas with the selected mechanism
gas = Solution('Lu19.xml')

# Create the reservoir that represents the atmosphere and fill it with air at 1 bar
air = Solution('air.xml')

#Number of species in the.cti file.
m=gas.n_species

#Find fuel, nitrogen, and oxygen indices
fuel_species = 'CH4'
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

if ifuel < 0:
    raise "fuel species "+fuel_species+" not present!"

if gas.n_atoms(fuel_species,'O') > 0 or  gas.n_atoms(fuel_species,'N') > 0:
    raise "Error: only hydrocarbon fuels are supported."

#################
#Enter general parameters

	#Stoechiometry
print ""
print "-------------------------------------------------------- "
print "    THERMO PROPERTIES: "
print "-------------------------------------------------------- "
print ""
phi        = input('Enter Stoichiometric ratio phi : ')
phi        = float(phi)
print ""

		#Air composition
air_N2_O2_ratio = 3.76
stoich_O2 = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')

		#Mass fraction vector
x = np.zeros(m,'d')
x[ifuel] = phi
x[io2] = stoich_O2
#x[in2] = stoich_O2*air_N2_O2_ratio
x[in2] = 0      # pure O2

	# Specify intial pressures and temperature of the reactor
Ti = input('Enter temperature (in kelvin) : ')
Ti = float(Ti)       # Kelvin

Pi = input('Enter pressure (in bar) : ')
Pi = float(Pi)*1e5         # Pascal


	#Set initial thermodynamic state of the gas
gas.TPX = Ti, Pi, x

T_ref = 298.15
gas0 = Solution('Lu19.xml')
gas0.TP = T_ref, Pi
h0 = gas0.partial_molar_enthalpies    # J/kmol

#################################################################
# Program starts here
#################################################################
#Create the batch reactor, constant Volume
r   = IdealGasReactor(gas)  	# IdealGasConstPressureReactor for constant P


#Specify the conditions: Pression or Volume constant
print "--------------------------------------------------- "
print "    Equilibirum conditions: "
print "--------------------------------------------------- "
print ""
print "For a constant volume equilibrium, enter :      UV "
print "For a constant pressure equilibrium, enter :    HP "
print ""
cond  = raw_input('Specify the equilibrium condition : ')
cond  = str(cond)
print ""
while cond != 'HP' and cond != 'UV':
     print "You must choose between UV and HP !  "
     cond  = raw_input('Specify the equilibrium condition : ')
     cond  = str(cond)

	#Particular case of a constant-pressure reactor
if cond == 'HP':
	# Define a wall between the reactor and the environment, and
	# make it flexible, so that the pressure in the reactor is held
	# at the environment pressure.
     env 		    = Reservoir(air)
     w 		            = Wall(r,env)
     w.expansion_rate_coeff = 1.0e6   # set expansion parameter. dV/dt = KA(P_1 - P_2)
     w.area                 = 1.0       # set wall area


# Now create a reactor network consisting of the single batch reactor
# Reason: the only way to advance reactors in time is through a network
sim = ReactorNet([r])
#sim.atol = 1e-18  # abs tolerances for chemistry integrator
#sim.rtol = 1e-11 # relative tolerance

#################
#Computational properties: we're going to advance the network in time

print ""
print "-------------------------------------------------------- "
print "    COMPUTATIONAL PROPERTIES: "
print "-------------------------------------------------------- "
print ""

	# Initial simulation time
time = 0.05

	# Specify the number of time steps
nt        = input('Enter maximum time: ')
nt        = float(nt)

	# Specify the time step
dtms	  = input('Enter the time step (in ms): ')
dtms	  = float(dtms)

dt = dtms * 1.0e-3 #s

#################
#Run the simulation

	#parameters
tim = [] #np.zeros(nt,'d')
temp = [] #np.zeros(nt,'d')
press = [] #np.zeros(nt,'d')
#mfrac = [] #np.zeros([nt,m],'d')
Y_CH4 = []
Y_O2 = []
Y_H2 = []
Y_O = []
Y_OH = []
Y_H2O = []
Y_CO = []
Y_CO2 = []
Y_H = []
RR_CH4 = []
RR_O2 = []
RR_H2 = []
RR_O = []
RR_OH =[]
RR_H2O = []
RR_CO = []
RR_CO2 = []
RR_H = []        # -> to check effects on transport
Heat_release = []
t_ign=[]    # capture only the delta t 0.24192 < t < 0.24194 s where the reaction occurs
t_min = 0.22
t_max = 0.24
t_space = 0.0001   # used in the plot

    # create index for interested species
indexCH4 = gas.species_index('CH4')
indexO2 = gas.species_index('O2')
indexH2 = gas.species_index('H2')
indexO = gas.species_index('O')
indexOH = gas.species_index('OH')
indexH2O = gas.species_index('H2O')
indexCO = gas.species_index('CO')
indexCO2 = gas.species_index('CO2')
indexH = gas.species_index('H')

# Initialize plot for run-time check on temperature
#fig = plt.gcf()
#fig.show()
#fig.canvas.draw()

	#Loop for nt time steps of dt seconds.
print 'time [s] ,   T [K] ,   p [Pa] ,   h [J/kg]'

#for n in range(nt):
while time < nt:
    time += dt
    sim.advance(time)
    tim.append(time)
    temp.append(r.T)
    press.append(r.thermo.P)
    Y_CH4.append(r.thermo['CH4'].Y)
    Y_O2.append(r.thermo['O2'].Y)
    Y_H2.append(r.thermo['H2'].Y)
    Y_O.append(r.thermo['O'].Y)
    Y_OH.append(r.thermo['OH'].Y)
    Y_H2O.append(r.thermo['H2O'].Y)
    Y_CO.append(r.thermo['CO'].Y)
    Y_CO2.append(r.thermo['CO2'].Y)
    Y_H.append(r.thermo['H'].Y)

    if (time > t_min) and (time < t_max):
        t_ign.append(time)
        # kmol/m3/s *kg/kmol = kg/m3/s
        RR_CH4.append(gas.net_production_rates[indexCH4]*gas.molecular_weights[indexCH4])
        RR_O2.append(gas.net_production_rates[indexO2]*gas.molecular_weights[indexO2])
        RR_H2.append(gas.net_production_rates[indexH2]*gas.molecular_weights[indexH2])
        RR_O.append(gas.net_production_rates[indexO]*gas.molecular_weights[indexO])
        RR_OH.append(gas.net_production_rates[indexOH]*gas.molecular_weights[indexOH])
        RR_H2O.append(gas.net_production_rates[indexH2O]*gas.molecular_weights[indexH2O])
        RR_CO.append(gas.net_production_rates[indexCO]*gas.molecular_weights[indexCO])
        RR_CO2.append(gas.net_production_rates[indexCO2]*gas.molecular_weights[indexCO2])

        Heat = 0.0
        for i in range(m):
            Heat -= gas.net_production_rates[i] * h0[i]     # kmol/m3/s * J/kmol = W / m3
        Heat_release.append(Heat)

        # following rates only as check for Lu13
        RR_H.append(gas.net_production_rates[indexH] * gas.molecular_weights[indexH])
    #for i in range(m):
     #   mfrac[n,i]=r.thermo[gas.species_name(i)].Y
    print '%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T,
                                           r.thermo.P, r.thermo.h)
    # update plot for T
    #plt.plot(tim, temp)
    #fig.canvas.draw()

#################################################################
# Save your results if needed
#################################################################
l1 =len(Y_CH4)
l2 =len(t_ign)

# write output CSV file for importing into Excel
if cond == 'HP':
     csvfile = 'Reactor_HP_Lu19.csv'
elif cond == 'UV':
     csvfile = 'Reactor_UV_Lu19.csv'

csv_file = csvfile
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile,delimiter=',')
    writer.writerow(['t','T','Y_CH4', 'Y_O2', 'Y_H2','Y_O','Y_OH', 'Y_H2O', 'Y_CO', 'Y_CO2','Y_H'])
    for i in range(0, l1, 100):
        writer.writerow([tim[i],temp[i],str(Y_CH4[i])[1:-1],str(Y_O2[i])[1:-1],str(Y_H2[i])[1:-1],str(Y_O[i])[1:-1],str(Y_OH[i])[1:-1],str(Y_H2O[i])[1:-1],str(Y_CO[i])[1:-1],str(Y_CO2[i])[1:-1],str(Y_H[i])[1:-1]])
print 'output written to '+csvfile

# write output CSV file for importing into Excel
# Save file #1
if cond == 'HP':
     csvfile = 'Reactor_HP_Lu19_RRi.csv'
elif cond == 'UV':
     csvfile = 'Reactor_UV_Lu19_RRi.csv'

# Save file #2
csv_file = csvfile
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile,delimiter=',')
    writer.writerow(['t', 'RR.CH4', 'RR.O2', 'RR.H2', 'RR.O', 'RR.OH', 'RR.H2O', 'RR.CO', 'RR.CO2', 'HeatRelease', 'RR.H'])
    for i in range(0, l2):
        writer.writerow([t_ign[i], RR_CH4[i], RR_O2[i], RR_H2[i], RR_O[i], RR_OH[i], RR_H2O[i], RR_CO[i], RR_CO2[i],
                         Heat_release[i], RR_H[i]])
print 'output written to '+csvfile

'''
#################################################################
# Plot your results
#################################################################
# plot the results if matplotlib is installed.
# see http://matplotlib.sourceforge.net to get it
args = sys.argv
if len(args) > 1 and args[1] == '-plot':
	import matplotlib.pyplot as plt
        plt.clf()
        plt.subplot(2,2,1)
        plt.plot(tim,temp[:])
        plt.xlabel('Time (s)');
        plt.ylabel('Temperature (K)');
        plt.subplot(2,2,2)
        plt.plot(tim,Y_OH)
        plt.xlabel('Time (s)');
        plt.ylabel('OH Mass Fraction');
        plt.subplot(2,2,3)
        plt.plot(tim,Y_H);
        plt.xlabel('Time (s)');
        plt.ylabel('H Mass Fraction');
        plt.subplot(2,2,4)
        plt.plot(tim,Y_H2);
        plt.xlabel('Time (s)');
        plt.ylabel('H2 Mass Fraction');
        plt.show()
else:
    print """To view a plot of these results, run this script with the option -plot"""
'''