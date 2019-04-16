'''
This is makes use of PyFoam to loop over different reactor simulations

last change: 4.3.2019
'''

import numpy as np
import pandas as pd
from os.path import join
import sys, re
import os
from shutil import copyfile, copytree
import fileinput

from PyFoam.Execution.ConvergenceRunner import ConvergenceRunner
from PyFoam.Execution.UtilityRunner import UtilityRunner
from PyFoam.LogAnalysis.BoundingLogAnalyzer import BoundingLogAnalyzer
from PyFoam.RunDictionary.SolutionFile import SolutionFile
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
import Ofpp


# parent directory
main_dir = '/home/max/HDD2_Data/OF4_Simulations/CanteraIgnition/DataGeneration/IgniteLu19/Cases'

source_dir = '/home/max/HDD2_Data/OF4_Simulations/CanteraIgnition/DataGeneration/IgniteLu19/source_files'

base_dir_name = 'Lu19_f'

solver = 'canteraOnlyReactingFoam'

lu_19_species = ['C2H2', 'C2H4', 'C2H6', 'CH2CO', 'CH2O', 'CH3', 'CH3OH', 'CH4','CO', 'CO2', 'H', 'H2', 'H2O', 'H2O2', 'HO2', 'N2', 'O', 'O2', 'OH']

T_start = 1100

#range of scalar dissipation rates
f_range=np.linspace(0.005,0.3,150)

os.chdir(main_dir)

# now loop over different chi values
for f in f_range:

    os.chdir(main_dir)

    print('Case: ', f)

    this_dir_name = base_dir_name+str(f)
    if os.path.isdir(this_dir_name):
        print('%s directory exists' % this_dir_name)
    else:
        os.mkdir(this_dir_name)

    # copy: system, constant, 0 to current directory
    try:
        copytree(join(source_dir,'0'),join(this_dir_name,'0'))
        copytree(join(source_dir,'constant'), join(this_dir_name,'constant'))
        copytree(join(source_dir,'system'), join(this_dir_name,'system'))
    except:
        print('Could not copy files from sourceFiles')

    # change into the case directory
    os.chdir(this_dir_name)
    print()

    # correct the flamelet properties
    for line in fileinput.input('0/CH4',inplace=True):
        print(line.replace('uniform', 'uniform  %s; //' % str(f))),

    y_O2 = (1 - f)*0.2315
    y_N2 = 1-y_O2 -f
    for line in fileinput.input('0/O2',inplace=True):
        print(line.replace('uniform', 'uniform  %s; //' % str(y_O2))),

    for line in fileinput.input('0/N2',inplace=True):
        print(line.replace('uniform', 'uniform  %s; //' % str(y_N2))),

    # for line in fileinput.input('0/H',inplace=True):
    #     print(line.replace('uniform', 'uniform  %s; //' % str(0))),

    run=ConvergenceRunner(BoundingLogAnalyzer(),argv=[solver],silent=True)
    print('Starting: ',solver)
    run.start()
    print('Done!\n')

    # remove all the log files
    print('Cleaning up!')
    allFiles = os.listdir('.')

    rmFiles = [f for f in allFiles if f.startswith('PyFoam')]

    for file in rmFiles:
        try:
            os.remove(file)
        except:
            print('Could not remove: ', file)

