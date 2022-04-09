# %% Script description
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Sep 22 17:29:01 2020
@author: oliviakosterlitz
"""

script_description = """
Code will take the data from one treatement with various populations,
and it will analyze different measures for transfer.
"""

# %% Load packages
import argparse
import sys
import os
import numpy
import time
import pandas as pd
from gillespy2 import ODESolver
import matplotlib.pyplot as pyplot
import seaborn
from SimFunctions import *

# %% System arguments
parser = argparse.ArgumentParser(
          description = script_description,
          formatter_class = argparse.RawDescriptionHelpFormatter)
         
parser.add_argument("-t", "--treatment", required = True, help = "example, T0")
parser.add_argument("-c", "--parameter_csv", required = True, help = "Treatments_master.csv")
parser.add_argument("-f", "--csv_folder", required = True)
parser.add_argument("-o", "--output_folder", required = True, help = "Output folder")

args = parser.parse_args()
                
treatment = str(args.treatment)        
parameter_csv = args.parameter_csv
csvfilefolder = args.csv_folder
output_folder = args.output_folder

# %% Troubleshooting aruguments
#from GitHub.Programs.SimFunctions import *
#treatment = "T1"
#parameter_csv = "/Users/oliviakosterlitz/Dropbox/Pop3/GitHub/Simulations/BasicRun/Treatments/SimSetup_inputCSV_example_with_incubation.csv"
#csvfilefolder = "/Users/oliviakosterlitz/Dropbox/Pop3/GitHub/Simulations/BasicRun/SimData/T1"
#output_folder = "/Users/oliviakosterlitz/Dropbox/Pop3/GitHub/Simulations/BasicRun"

#treatment = "T16"
#parameter_csv = "/Users/oliviakosterlitz/Dropbox/Pop3/Simulations/Submission/Output/Treatments/Treatments_master_with_incubation.csv"
#csvfilefolder = "/Users/oliviakosterlitz/Dropbox/Pop3/Simulations/Submission/Output/SimData/T59"
#output_folder = "/Users/oliviakosterlitz/Dropbox/Pop3/Simulations/Submission/Output"

# %% Load treatment parameters
parameters_all = pd.read_csv(parameter_csv)
parameters = parameters_all[parameters_all['Treatment_ID'] == treatment]
num_experiments = int(parameters['Num_of_experiments'])
num_I2_populations = int(parameters['Num_sims_I2'])
num_I1_populations = int(parameters['Num_sims_I1'])

# %% Load customizations
time_series_intervals = 0.25 
ALPHA = 0.25
Num_of_trans_incubation_time = 50

# %% Load Gillespie simulation data
parent_dir = os.getcwd()
csvfiles = os.listdir(csvfilefolder)

# Build treatment dictionary
data = {}
data['parameters'] = {}
data['experiments'] = {}
data['LDM'] = {}
data['SIM'] = {}
data['ASM'] = {}
data['TDR'] = {}

# Load the data for the treatment (i.e. various experiments with various number of replicate  populations with different incubation types; I1, I2, I3)
for file in csvfiles:
    if str(treatment+"_E") in file:
        sim = file.split(".")[0].split("_")
        experiment = int(sim[1][1:])
        incubation = str(sim[2])
        population = int(sim[3][1:])        
        if experiment not in data['experiments'].keys():
            data['experiments'][experiment] = {}
            data['experiments'][experiment]['I2'] = [None] * num_I2_populations
            data['experiments'][experiment]['I1'] = [None] * num_I1_populations
        if incubation == 'IL' or incubation == "I1":
            data['experiments'][experiment]["I1"][population]=pd.read_csv(csvfilefolder + '/' + file)
        if incubation == 'IS' or incubation == "I2":
            data['experiments'][experiment]["I2"][population]=pd.read_csv(csvfilefolder + '/' + file)
# Continue loading treatment parameters
data['parameters']['psiD'] = float(parameters['psiD'])
data['parameters']['psiR'] = float(parameters['psiR'])
data['parameters']['psiT'] = float(parameters['psiT'])
data['parameters']['psiF'] = float(parameters['psiF'])
data['parameters']['gammaDR'] = float(parameters['gammaDR'])
data['parameters']['gammaTR'] = float(parameters['gammaTR'])
data['parameters']['gammaDF'] = float(parameters['gammaDF'])
data['parameters']['gammaTF'] = float(parameters['gammaTF'])
data['parameters']['tauD'] = float(parameters['tauD'])
data['parameters']['tauT'] = float(parameters['tauT'])
data['parameters']['D_0'] = float(parameters['D_0'])
data['parameters']['R_0'] = float(parameters['R_0'])
data['parameters']['T_0'] = float(parameters['T_0'])
data['parameters']['F_0'] = float(parameters['F_0'])
data['parameters']['I1_t'] = float(parameters['Cutoff_I1_t'])
data['parameters']['I2_t'] = float(parameters['Cutoff_I2_t'])
if data['parameters']['I1_t'] < data['parameters']['I2_t']:
    data['parameters']['I2_t'] = data['parameters']['I1_t']

# %% Numerical simulation
psiD_sec=rate_hourtosecond(data['parameters']['psiD'])
psiR_sec=rate_hourtosecond(data['parameters']['psiR'])
psiT_sec=rate_hourtosecond(data['parameters']['psiT'])
psiF_sec=rate_hourtosecond(data['parameters']['psiF'])

gammaDR_sec=rate_hourtosecond(data['parameters']['gammaDR'])
gammaTR_sec=rate_hourtosecond(data['parameters']['gammaTR'])
gammaDF_sec=rate_hourtosecond(data['parameters']['gammaDF'])
gammaTF_sec=rate_hourtosecond(data['parameters']['gammaTF'])

tauDsec=rate_hourtosecond(data['parameters']['tauD'])
tauT_sec=rate_hourtosecond(data['parameters']['tauT'])

D_init=data['parameters']['D_0']
R_init=data['parameters']['R_0']
T_init=data['parameters']['T_0']
F_init=data['parameters']['F_0']

ODE_runtime_sec=time_hourtosecond(data['parameters']['I1_t'])

# Initiate the model 
ODEmodel = PlasmidDynamics(runtime = ODE_runtime_sec,
                            psiD = psiD_sec,
                            psiR = psiR_sec,
                            psiT = psiT_sec,
                            psiF = psiF_sec,
                            gammaDR = gammaDR_sec,
                            gammaTR = gammaTR_sec,
                            gammaDF = gammaDF_sec,
                            gammaTF = gammaTF_sec,
                            tauD = tauDsec,
                            tauT = tauT_sec,
                            D_init = D_init,
                            R_init = R_init,
                            T_init = T_init,
                            F_init = F_init)

# Run the model
resultODE = ODEmodel.run(solver=ODESolver)

# %% Incubation time for non-LDM methods

if resultODE[0]['Transconjugant'][-1] < Num_of_trans_incubation_time:
    tcrit = ((resultODE[0]['time'][-1])/3600) #If the number of transconjugants doesn't pass the threshold set, pick the last time point 
else: 
    tcrit = ((resultODE[0]['time'][resultODE[0]['Transconjugant'] > Num_of_trans_incubation_time][0])/3600)

# %% LDM calculation
# Time series data
data['LDM']['time_series'] = numpy.arange(time_series_intervals, data['parameters']['I2_t'], time_series_intervals)  # make all of the time points at certain time intervals
data['LDM']['TimeSeries'] = {}

# Loop through the experiment and creat a sub list for each (i.e., rep = experiment)
for rep in range(num_experiments):  # create the empty placeholders to calculate the LDM at each time interval for each experiment
    data['LDM']['TimeSeries'][rep] = [None] * len(data['LDM']['time_series'])

# store the first time a transconjugant appears for each population in each experiment
for rep in range(num_experiments):
    data['LDM'][rep] = {}
    data['LDM'][rep]['firstTimes'] = [None] * num_I2_populations
    for pop in range(num_I2_populations):
        if data['experiments'][rep]['I2'][pop]['Transconjugant'][len(data['experiments'][rep]['I2'][pop]['Transconjugant'])-1] > 0:
            data['LDM'][rep]['firstTimes'][pop] = data['experiments'][rep]['I2'][pop]['time'][data['experiments'][rep]['I2'][pop]['Transconjugant'] >= 1].iloc[0]
        else:
            data['LDM'][rep]['firstTimes'][pop] = len(data['experiments'][rep]['I2'][pop]['Transconjugant'])

for t in range(len(data['LDM']['time_series'])):  # loops through all of the time points
    t_hours = data['LDM']['time_series'][t]  # current time point
    t_seconds = data['LDM']['time_series'][t]*3600
    for rep in range(num_experiments):
        data['LDM'][rep][t_hours] = {} 
        data['LDM'][rep][t_hours]['count_nonturbid'] = len([i for i in data['LDM'][rep]['firstTimes'] if i > t_seconds])        
        data['LDM'][rep][t_hours]['p0'] = (data['LDM'][rep][t_hours]['count_nonturbid']/num_I2_populations)
        data['LDM'][rep][t_hours]['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Donor'])
        data['LDM'][rep][t_hours]['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Recipient'])
        data['LDM']['TimeSeries'][rep][t] = LDMmetric(t_hours,
                                                      data['LDM'][rep][t_hours]['p0'],
                                                      data['parameters']['D_0'],
                                                      data['parameters']['R_0'], 
                                                      data['LDM'][rep][t_hours]['Dt'], 
                                                      data['LDM'][rep][t_hours]['Rt'])
        data['LDM'][rep][t_hours]['LDM'] = data['LDM']['TimeSeries'][rep][t]

# LDM estimate at the incubation time corresponding to on average 1 transconjugant determined using the numerical solution
data['LDM']['time_secs'] = (resultODE[0]['time'][resultODE[0]['Transconjugant']>1][0])
data['LDM']['time_hours'] = (data['LDM']['time_secs'])/3600

for rep in range(num_experiments):
    data['LDM'][rep]['criteria_t'] = {}
    data['LDM'][rep]['criteria_t']['count_nonturbid'] = len([i for i in data['LDM'][rep]['firstTimes'] if i > data['LDM']['time_secs']])   
    data['LDM'][rep]['criteria_t']['p0'] = (data['LDM'][rep]['criteria_t']['count_nonturbid']/num_I2_populations)
    data['LDM'][rep]['criteria_t']['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['LDM']['time_secs']]['Donor'])
    data['LDM'][rep]['criteria_t']['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['LDM']['time_secs']]['Recipient'])
    data['LDM'][rep]['criteria_t']['LDM'] = LDMmetric(data['LDM']['time_hours'],
                                                      data['LDM'][rep]['criteria_t']['p0'],
                                                      data['parameters']['D_0'],
                                                      data['parameters']['R_0'], 
                                                      data['LDM'][rep]['criteria_t']['Dt'], 
                                                      data['LDM'][rep]['criteria_t']['Rt'])
    
# %% SIM calculation
data['SIM']['time_series'] = numpy.arange(time_series_intervals, data['parameters']['I1_t'], time_series_intervals)
data['SIM']['time_series_I2'] = numpy.arange(time_series_intervals, data['parameters']['I2_t'], time_series_intervals)
data['SIM']['TimeSeries'] = {}
data['SIM']['TimeSeries_I2'] = {}
for rep in range(num_experiments):
    data['SIM']['TimeSeries'][rep] = [None] * len(data['SIM']['time_series'])
    data['SIM']['TimeSeries_I2'][rep] = {}
    data['SIM'][rep] = {}
    data['SIM'][rep]['I2'] = {}
    for pop in range(num_I2_populations):
            data['SIM'][rep]['I2'][pop] = {}
            data['SIM']['TimeSeries_I2'][rep][pop] = [None] * len(data['SIM']['time_series_I2'])
for t in range(len(data['SIM']['time_series'])):
    t_hours = data['SIM']['time_series'][t]
    t_seconds = data['SIM']['time_series'][t]*3600
    for rep in range(num_experiments):
        data['SIM'][rep][t_hours] = {}
        data['SIM'][rep][t_hours]['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Donor'])
        data['SIM'][rep][t_hours]['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Recipient'])
        data['SIM'][rep][t_hours]['Tt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Transconjugant'])
        data['SIM'][rep][t_hours]['N0'] = data['parameters']['D_0'] + data['parameters']['R_0']
        data['SIM'][rep][t_hours]['Nt'] = data['SIM'][rep][t_hours]['Dt'] + data['SIM'][rep][t_hours]['Rt'] + data['SIM'][rep][t_hours]['Tt']
        if (data['SIM'][rep][t_hours]['Dt'] > 0 and data['SIM'][rep][t_hours]['Rt'] > 0):
            data['SIM']['TimeSeries'][rep][t] = SIMmetric(t_hours,
                                                          data['SIM'][rep][t_hours]['N0'], 
                                                          data['SIM'][rep][t_hours]['Nt'], 
                                                          data['SIM'][rep][t_hours]['Dt'], 
                                                          data['SIM'][rep][t_hours]['Rt'], 
                                                          data['SIM'][rep][t_hours]['Tt'])
            data['SIM'][rep][t_hours]['SIM'] = data['SIM']['TimeSeries'][rep][t]  
        if t_hours < data['parameters']['I2_t']:
            for pop in range(num_I2_populations):
                data['SIM'][rep]['I2'][pop][t_hours] = {}
                data['SIM'][rep]['I2'][pop][t_hours]['Dt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Donor'])
                data['SIM'][rep]['I2'][pop][t_hours]['Rt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Recipient'])
                data['SIM'][rep]['I2'][pop][t_hours]['Tt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Transconjugant'])
                data['SIM'][rep]['I2'][pop][t_hours]['N0'] = data['parameters']['D_0'] + data['parameters']['R_0']
                data['SIM'][rep]['I2'][pop][t_hours]['Nt'] = data['SIM'][rep]['I2'][pop][t_hours]['Dt'] + data['SIM'][rep]['I2'][pop][t_hours]['Rt'] + data['SIM'][rep]['I2'][pop][t_hours]['Tt']
                if data['SIM'][rep]['I2'][pop][t_hours]['Dt'] > 0 and data['SIM'][rep]['I2'][pop][t_hours]['Rt'] > 0:
                    data['SIM']['TimeSeries_I2'][rep][pop][t] = SIMmetric(t_hours,
                                                                  data['SIM'][rep]['I2'][pop][t_hours]['N0'], 
                                                                  data['SIM'][rep]['I2'][pop][t_hours]['Nt'], 
                                                                  data['SIM'][rep]['I2'][pop][t_hours]['Dt'], 
                                                                  data['SIM'][rep]['I2'][pop][t_hours]['Rt'], 
                                                                  data['SIM'][rep]['I2'][pop][t_hours]['Tt'])
                    data['SIM'][rep]['I2'][pop][t_hours]['SIM'] = data['SIM']['TimeSeries_I2'][rep][pop][t]

data['SIM']['time_hours'] = tcrit
data['SIM']['time_secs'] = int(time_hourtosecond(data['SIM']['time_hours']))
for rep in range(num_experiments):
    data['SIM'][rep]['criteria_t'] = {}
    data['SIM'][rep]['criteria_t']['N0'] = data['parameters']['D_0'] + data['parameters']['R_0']
    data['SIM'][rep]['criteria_t']['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['SIM']['time_secs']]['Donor'])
    data['SIM'][rep]['criteria_t']['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['SIM']['time_secs']]['Recipient'])
    data['SIM'][rep]['criteria_t']['Tt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['SIM']['time_secs']]['Transconjugant'])
    data['SIM'][rep]['criteria_t']['Nt'] = data['SIM'][rep]['criteria_t']['Dt'] + data['SIM'][rep]['criteria_t']['Rt'] + data['SIM'][rep]['criteria_t']['Tt']
    if data['SIM'][rep]['criteria_t']['Dt'] > 0 and data['SIM'][rep]['criteria_t']['Rt'] > 0:
        data['SIM'][rep]['criteria_t']['SIM'] = SIMmetric(data['SIM']['time_hours'],
                                                          data['SIM'][rep]['criteria_t']['N0'], 
                                                          data['SIM'][rep]['criteria_t']['Nt'], 
                                                          data['SIM'][rep]['criteria_t']['Dt'], 
                                                          data['SIM'][rep]['criteria_t']['Rt'], 
                                                          data['SIM'][rep]['criteria_t']['Tt'])
    for pop in range(num_I2_populations):
        data['SIM'][rep]['criteria_t'][pop] = {}
        data['SIM'][rep]['criteria_t'][pop]['N0'] = data['parameters']['D_0'] + data['parameters']['R_0']
        data['SIM'][rep]['criteria_t'][pop]['Dt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['SIM']['time_secs']]['Donor'])
        data['SIM'][rep]['criteria_t'][pop]['Rt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['SIM']['time_secs']]['Recipient'])
        data['SIM'][rep]['criteria_t'][pop]['Tt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['SIM']['time_secs']]['Transconjugant'])
        data['SIM'][rep]['criteria_t'][pop]['Nt'] = data['SIM'][rep]['criteria_t'][pop]['Dt'] + data['SIM'][rep]['criteria_t'][pop]['Rt'] + data['SIM'][rep]['criteria_t'][pop]['Tt']
        if data['SIM'][rep]['criteria_t'][pop]['Dt'] > 0 and data['SIM'][rep]['criteria_t'][pop]['Rt'] > 0:
            data['SIM'][rep]['criteria_t'][pop]['SIM'] = SIMmetric(data['SIM']['time_hours'],
                                                              data['SIM'][rep]['criteria_t'][pop]['N0'], 
                                                              data['SIM'][rep]['criteria_t'][pop]['Nt'], 
                                                              data['SIM'][rep]['criteria_t'][pop]['Dt'], 
                                                              data['SIM'][rep]['criteria_t'][pop]['Rt'], 
                                                              data['SIM'][rep]['criteria_t'][pop]['Tt'])
            
# %% ASM calculation
data['ASM']['time_series'] = numpy.arange(time_series_intervals, data['parameters']['I1_t'], time_series_intervals)
data['ASM']['time_series_I2'] = numpy.arange(time_series_intervals, data['parameters']['I2_t'], time_series_intervals)
data['ASM']['TimeSeries'] = {}
data['ASM']['TimeSeries_I2'] = {}
for rep in range(num_experiments):
    data['ASM']['TimeSeries'][rep] = [None] * len(data['ASM']['time_series'])
    data['ASM']['TimeSeries_I2'][rep] = {}
    data['ASM'][rep] = {}
    data['ASM'][rep]['I2'] = {}
    for pop in range(num_I2_populations):
            data['ASM'][rep]['I2'][pop] = {}
            data['ASM']['TimeSeries_I2'][rep][pop] = [None] * len(data['ASM']['time_series_I2'])
for t in range(len(data['ASM']['time_series'])):
    t_hours = data['ASM']['time_series'][t]
    t_seconds = data['ASM']['time_series'][t]*3600
    for rep in range(num_experiments):
        data['ASM'][rep][t_hours] = {}
        data['ASM'][rep][t_hours]['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Donor'])
        data['ASM'][rep][t_hours]['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Recipient'])
        data['ASM'][rep][t_hours]['Tt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Transconjugant'])
        if data['ASM'][rep][t_hours]['Dt'] > 0 and data['ASM'][rep][t_hours]['Rt'] > 0:
            data['ASM'][rep][t_hours]['psiD'] = growthrate(t_hours, data['parameters']['D_0'], data['ASM'][rep][t_hours]['Dt'])
            data['ASM'][rep][t_hours]['psiR'] = growthrate(t_hours, data['parameters']['R_0'], data['ASM'][rep][t_hours]['Rt'])
            data['ASM']['TimeSeries'][rep][t] = ASMmetric(t_hours,
                                                          data['ASM'][rep][t_hours]['psiD'],
                                                          data['ASM'][rep][t_hours]['psiR'],
                                                          data['parameters']['psiT'],
                                                          data['parameters']['D_0'],
                                                          data['parameters']['R_0'],
                                                          data['ASM'][rep][t_hours]['Dt'], 
                                                          data['ASM'][rep][t_hours]['Rt'], 
                                                          data['ASM'][rep][t_hours]['Tt'])   
            data['ASM'][rep][t_hours]['ASM'] = data['ASM']['TimeSeries'][rep][t] 
        if t_hours < data['parameters']['I2_t']:
            for pop in range(num_I2_populations):
                data['ASM'][rep]['I2'][pop][t_hours] = {}
                data['ASM'][rep]['I2'][pop][t_hours]['Dt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Donor'])
                data['ASM'][rep]['I2'][pop][t_hours]['Rt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Recipient'])
                data['ASM'][rep]['I2'][pop][t_hours]['Tt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Transconjugant'])
                if data['ASM'][rep]['I2'][pop][t_hours]['Dt'] > 0 and data['ASM'][rep]['I2'][pop][t_hours]['Rt'] > 0:
                    data['ASM'][rep]['I2'][pop][t_hours]['psiD'] = growthrate(t_hours, data['parameters']['D_0'], data['ASM'][rep]['I2'][pop][t_hours]['Dt'])
                    data['ASM'][rep]['I2'][pop][t_hours]['psiR'] = growthrate(t_hours, data['parameters']['R_0'], data['ASM'][rep]['I2'][pop][t_hours]['Rt'])
                    data['ASM']['TimeSeries_I2'][rep][pop][t] = ASMmetric(t_hours,
                                                                  data['ASM'][rep]['I2'][pop][t_hours]['psiD'],
                                                                  data['ASM'][rep]['I2'][pop][t_hours]['psiR'],
                                                                  data['parameters']['psiT'],
                                                                  data['parameters']['D_0'],
                                                                  data['parameters']['R_0'],
                                                                  data['ASM'][rep]['I2'][pop][t_hours]['Dt'], 
                                                                  data['ASM'][rep]['I2'][pop][t_hours]['Rt'], 
                                                                  data['ASM'][rep]['I2'][pop][t_hours]['Tt'])   
                    data['ASM'][rep]['I2'][pop][t_hours]['ASM'] = data['ASM']['TimeSeries_I2'][rep][pop][t]        
                
data['ASM']['time_hours'] = tcrit
data['ASM']['time_secs'] = int(time_hourtosecond(data['ASM']['time_hours']))
for rep in range(num_experiments):
    data['ASM'][rep]['criteria_t'] = {}
    data['ASM'][rep]['criteria_t']['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['ASM']['time_secs']]['Donor'])
    data['ASM'][rep]['criteria_t']['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['ASM']['time_secs']]['Recipient'])
    data['ASM'][rep]['criteria_t']['Tt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['ASM']['time_secs']]['Transconjugant'])
    if data['ASM'][rep]['criteria_t']['Dt'] and data['ASM'][rep]['criteria_t']['Rt']:
        data['ASM'][rep]['criteria_t']['psiD'] = growthrate(data['ASM']['time_hours'], data['parameters']['D_0'], data['ASM'][rep]['criteria_t']['Dt'])
        data['ASM'][rep]['criteria_t']['psiR'] = growthrate(data['ASM']['time_hours'], data['parameters']['R_0'], data['ASM'][rep]['criteria_t']['Rt'])
        data['ASM'][rep]['criteria_t']['ASM'] = ASMmetric(data['ASM']['time_hours'],
                                                          data['ASM'][rep]['criteria_t']['psiD'],
                                                          data['ASM'][rep]['criteria_t']['psiR'],
                                                          data['parameters']['psiT'],
                                                          data['parameters']['D_0'],
                                                          data['parameters']['R_0'],
                                                          data['ASM'][rep]['criteria_t']['Dt'], 
                                                          data['ASM'][rep]['criteria_t']['Rt'], 
                                                          data['ASM'][rep]['criteria_t']['Tt'])  
    for pop in range(num_I2_populations):
        data['ASM'][rep]['criteria_t'][pop] = {}
        data['ASM'][rep]['criteria_t'][pop]['Dt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['ASM']['time_secs']]['Donor'])
        data['ASM'][rep]['criteria_t'][pop]['Rt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['ASM']['time_secs']]['Recipient'])
        data['ASM'][rep]['criteria_t'][pop]['Tt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['ASM']['time_secs']]['Transconjugant'])
        if data['ASM'][rep]['criteria_t'][pop]['Dt'] and data['ASM'][rep]['criteria_t'][pop]['Rt']:
            data['ASM'][rep]['criteria_t'][pop]['psiD'] = growthrate(data['ASM']['time_hours'], data['parameters']['D_0'], data['ASM'][rep]['criteria_t'][pop]['Dt'])
            data['ASM'][rep]['criteria_t'][pop]['psiR'] = growthrate(data['ASM']['time_hours'], data['parameters']['R_0'], data['ASM'][rep]['criteria_t'][pop]['Rt'])
            data['ASM'][rep]['criteria_t'][pop]['ASM'] = ASMmetric(data['ASM']['time_hours'],
                                                              data['ASM'][rep]['criteria_t'][pop]['psiD'],
                                                              data['ASM'][rep]['criteria_t'][pop]['psiR'],
                                                              data['parameters']['psiT'],
                                                              data['parameters']['D_0'],
                                                              data['parameters']['R_0'],
                                                              data['ASM'][rep]['criteria_t'][pop]['Dt'], 
                                                              data['ASM'][rep]['criteria_t'][pop]['Rt'], 
                                                              data['ASM'][rep]['criteria_t'][pop]['Tt'])  
# %% TDR calculation
data['TDR']['time_series'] = numpy.arange(time_series_intervals, data['parameters']['I1_t'], time_series_intervals)
data['TDR']['time_series_I2'] = numpy.arange(time_series_intervals, data['parameters']['I2_t'], time_series_intervals)
data['TDR']['TimeSeries'] = {}
data['TDR']['TimeSeries_I2'] = {}
for rep in range(num_experiments):
    data['TDR']['TimeSeries'][rep] = [None] * len(data['TDR']['time_series'])
    data['TDR']['TimeSeries_I2'][rep] = {}
    data['TDR'][rep] = {}
    data['TDR'][rep]['I2'] = {}
    for pop in range(num_I2_populations):
            data['TDR'][rep]['I2'][pop] = {}
            data['TDR']['TimeSeries_I2'][rep][pop] = [None] * len(data['TDR']['time_series_I2'])
for t in range(len(data['TDR']['time_series'])):
    t_hours = data['TDR']['time_series'][t]
    t_seconds = data['TDR']['time_series'][t]*3600
    for rep in range(num_experiments):
        data['TDR'][rep][t_hours] = {}
        data['TDR'][rep][t_hours]['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Donor'])
        data['TDR'][rep][t_hours]['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Recipient'])
        data['TDR'][rep][t_hours]['Tt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == t_seconds]['Transconjugant'])
        if data['TDR'][rep][t_hours]['Dt'] > 0 and data['TDR'][rep][t_hours]['Rt'] > 0:
            data['TDR']['TimeSeries'][rep][t] = TDRmetric(t_hours,
                                                          data['TDR'][rep][t_hours]['Dt'], 
                                                          data['TDR'][rep][t_hours]['Rt'], 
                                                          data['TDR'][rep][t_hours]['Tt'])
            data['TDR'][rep][t_hours]['TDR'] = data['TDR']['TimeSeries'][rep][t] 
        if t_hours < data['parameters']['I2_t']:
            for pop in range(num_I2_populations):
                data['TDR'][rep]['I2'][pop][t_hours] = {}
                data['TDR'][rep]['I2'][pop][t_hours]['Dt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Donor'])
                data['TDR'][rep]['I2'][pop][t_hours]['Rt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Recipient'])
                data['TDR'][rep]['I2'][pop][t_hours]['Tt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == t_seconds]['Transconjugant'])
                if data['TDR'][rep]['I2'][pop][t_hours]['Dt'] > 0 and data['TDR'][rep]['I2'][pop][t_hours]['Rt'] > 0:
                    data['TDR']['TimeSeries_I2'][rep][pop][t] = TDRmetric(t_hours,
                                                                          data['TDR'][rep]['I2'][pop][t_hours]['Dt'], 
                                                                          data['TDR'][rep]['I2'][pop][t_hours]['Rt'], 
                                                                          data['TDR'][rep]['I2'][pop][t_hours]['Tt'])
                    data['TDR'][rep]['I2'][pop][t_hours]['TDR'] = data['TDR']['TimeSeries_I2'][rep][pop][t]
        
data['TDR']['time_hours'] = tcrit
data['TDR']['time_secs'] = int(time_hourtosecond(data['TDR']['time_hours']))
for rep in range(num_experiments):
    data['TDR'][rep]['criteria_t'] = {}
    data['TDR'][rep]['criteria_t']['Dt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['TDR']['time_secs']]['Donor'])
    data['TDR'][rep]['criteria_t']['Rt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['TDR']['time_secs']]['Recipient'])
    data['TDR'][rep]['criteria_t']['Tt'] = int(data['experiments'][rep]['I1'][0][data['experiments'][rep]['I1'][0]['time'] == data['TDR']['time_secs']]['Transconjugant'])
    if data['TDR'][rep]['criteria_t']['Dt'] and data['TDR'][rep]['criteria_t']['Rt']:
        data['TDR'][rep]['criteria_t']['TDR'] = TDRmetric(data['TDR']['time_hours'],
                                                          data['TDR'][rep]['criteria_t']['Dt'], 
                                                          data['TDR'][rep]['criteria_t']['Rt'], 
                                                          data['TDR'][rep]['criteria_t']['Tt'])
    for pop in range(num_I2_populations):
        data['TDR'][rep]['criteria_t'][pop] = {}
        data['TDR'][rep]['criteria_t'][pop]['Dt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['TDR']['time_secs']]['Donor'])
        data['TDR'][rep]['criteria_t'][pop]['Rt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['TDR']['time_secs']]['Recipient'])
        data['TDR'][rep]['criteria_t'][pop]['Tt'] = int(data['experiments'][rep]['I2'][pop][data['experiments'][rep]['I2'][pop]['time'] == data['TDR']['time_secs']]['Transconjugant'])
        if data['TDR'][rep]['criteria_t'][pop]['Dt'] and data['TDR'][rep]['criteria_t'][pop]['Rt']:
            data['TDR'][rep]['criteria_t'][pop]['TDR'] = TDRmetric(data['TDR']['time_hours'],
                                                                   data['TDR'][rep]['criteria_t'][pop]['Dt'], 
                                                                   data['TDR'][rep]['criteria_t'][pop]['Rt'], 
                                                                   data['TDR'][rep]['criteria_t'][pop]['Tt'])  
        
# %% Density figures
if not os.path.exists(output_folder + '/' + "Analysis_figures"):
    os.makedirs(output_folder + '/' + "Analysis_figures")
if not os.path.exists(output_folder + '/' + "Analysis_figures" + '/' + "Densityonly"):
    os.makedirs(output_folder + '/' + "Analysis_figures" + '/' + "Densityonly")


density_filename_I1 = output_folder + '/' + "Analysis_figures" + '/' + "Densityonly" + "/" + str(treatment) + "_I1_simulation_density.pdf"
fig, ax = pyplot.subplots(1,1, figsize=(8,5)) # pyplot.subplots(rows, columns, figsize=(width,length))
donor_color = "red"
recipient_color = "blue"
transconjugant_color = "purple"
plasmid_free_donor_color = "grey"
for rep in range(4): 
    ax.plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['Donor'], linestyle='-', alpha = ALPHA, color=donor_color)
    ax.plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['Recipient'], linestyle='-', alpha = ALPHA, color=recipient_color)
    ax.plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['Transconjugant'], linestyle='-', alpha = ALPHA, color=transconjugant_color)
    ax.plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['PlasmidFreeDonor'], linestyle='-', alpha = ALPHA, color=plasmid_free_donor_color)
ax.plot([], [], ' ', label=r'$\psi_D=$'+str(data['parameters']['psiD']))
ax.plot([], [], ' ', label=r'$\psi_R=$'+str(data['parameters']['psiR']))
ax.plot([], [], ' ', label=r'$\psi_T=$'+str(data['parameters']['psiT']))
ax.plot([], [], ' ', label=r'$\psi_F=$'+str(data['parameters']['psiF']))
ax.plot([], [], ' ', label=r'$\gamma_{D}=$'+str(data['parameters']['gammaDR']))
ax.plot([], [], ' ', label=r'$\gamma_{T}=$'+str(data['parameters']['gammaTR']))
ax.plot([], [], ' ', label=r'$\tau_{D}=$'+str(data['parameters']['tauD']))
ax.plot([], [], ' ', label=r'$\tau_{T}=$'+str(data['parameters']['tauT']))
ax.plot([], [], ' ', label='selected incubation times:')
ax.axvline(x=data['parameters']['I2_t'], linestyle='--', color = 'darkgrey', label=('I2 =' + str(round(data['parameters']['I2_t'],2))))
ax.axvline(x=data['parameters']['I1_t'], linestyle='--', color = 'grey', label=('I1 =' + str(round(data['parameters']['I1_t'],2))))
ax.axvline(x=data['LDM']['time_hours'], linestyle='-.', color = 'brown', label=r'$t_{LDM}=$' + str(round(data['LDM']['time_hours'],2)))
ax.axvline(x=data['SIM']['time_hours'], linestyle='-.', color = 'orange', label=r'$t_{x}=$' + str(round(data['SIM']['time_hours'],2)))
ax.set_ylabel("population density (cells/mL)")
ax.set_xlabel("time (hours)")
ax.set_yscale('log')
ax.set_ylim(bottom = 0.9)
ax.legend(bbox_to_anchor=(1.01,1), facecolor='white', edgecolor='none', framealpha=0.9, prop={'size': 8})
fig.savefig(density_filename_I1, bbox_inches='tight')    


density_filename_I2 = output_folder + '/' + "Analysis_figures" + '/' + "Densityonly" + "/" + str(treatment) + "_I2_simulation_density.pdf"
fig, ax = pyplot.subplots(1,1, figsize=(8,5)) # pyplot.subplots(rows, columns, figsize=(width,length))
donor_color = "red"
recipient_color = "blue"
transconjugant_color = "purple"
plasmid_free_donor_color = "grey"
for rep in range(4): 
    ax.plot(data['experiments'][rep]['I2'][0]['time']/3600, data['experiments'][rep]['I2'][0]['Donor'], linestyle='-', alpha = ALPHA, color=donor_color)
    ax.plot(data['experiments'][rep]['I2'][0]['time']/3600, data['experiments'][rep]['I2'][0]['Recipient'], linestyle='-', alpha = ALPHA, color=recipient_color)
    ax.plot(data['experiments'][rep]['I2'][0]['time']/3600, data['experiments'][rep]['I2'][0]['Transconjugant'], linestyle='-', alpha = ALPHA, color=transconjugant_color)
    ax.plot(data['experiments'][rep]['I2'][0]['time']/3600, data['experiments'][rep]['I2'][0]['PlasmidFreeDonor'], linestyle='-', alpha = ALPHA, color=plasmid_free_donor_color)
ax.plot([], [], ' ', label=r'$\psi_D=$'+str(data['parameters']['psiD']))
ax.plot([], [], ' ', label=r'$\psi_R=$'+str(data['parameters']['psiR']))
ax.plot([], [], ' ', label=r'$\psi_T=$'+str(data['parameters']['psiT']))
ax.plot([], [], ' ', label=r'$\psi_F=$'+str(data['parameters']['psiF']))
ax.plot([], [], ' ', label=r'$\gamma_{D}=$'+str(data['parameters']['gammaDR']))
ax.plot([], [], ' ', label=r'$\gamma_{T}=$'+str(data['parameters']['gammaTR']))
ax.plot([], [], ' ', label=r'$\tau_{D}=$'+str(data['parameters']['tauD']))
ax.plot([], [], ' ', label=r'$\tau_{T}=$'+str(data['parameters']['tauT']))
ax.plot([], [], ' ', label='selected incubation times:')
ax.axvline(x=data['parameters']['I2_t'], linestyle='--', color = 'darkgrey', label=('I2 =' + str(round(data['parameters']['I2_t'],2))))
ax.axvline(x=data['parameters']['I1_t'], linestyle='--', color = 'grey', label=('I1 & I3 =' + str(round(data['parameters']['I1_t'],2))))
ax.axvline(x=data['LDM']['time_hours'], linestyle='-.', color = 'brown', label=r'$t_{LDM}=$' + str(round(data['LDM']['time_hours'],2)))
ax.axvline(x=data['SIM']['time_hours'], linestyle='-.', color = 'orange', label=r'$t_{x}=$' + str(round(data['SIM']['time_hours'],2)))
ax.set_ylabel("population density (cells/mL)")
ax.set_xlabel("time (hours)")
ax.set_yscale('log')
ax.set_ylim(bottom = 0.9)
ax.legend(bbox_to_anchor=(1.01,1), facecolor='white', edgecolor='none', framealpha=0.9, prop={'size': 8})
fig.savefig(density_filename_I2, bbox_inches='tight')    
    
# %% Time series figure
if not os.path.exists(output_folder + '/' + "Analysis_figures/Timeseriesonly"):
    os.makedirs(output_folder + '/' + "Analysis_figures/Timeseriesonly")

LDM_timeseries_filename = output_folder + '/' + "Analysis_figures/Timeseriesonly" + "/" + str(treatment) + "__LDM_timeseries.pdf"
TimeSeriesfigure(data, 'LDM', num_experiments, filename = LDM_timeseries_filename)

SIM_timeseries_filename = output_folder + '/' + "Analysis_figures/Timeseriesonly" + "/" + str(treatment) + "__SIM_timeseriesI1.pdf"
TimeSeriesfigure(data, 'SIM', num_experiments, SIM_timeseries_filename)

ASM_timeseries_filename = output_folder + '/' + "Analysis_figures/Timeseriesonly" + "/" + str(treatment) + "__ASM_timeseriesI1.pdf"
TimeSeriesfigure(data, 'ASM', num_experiments, ASM_timeseries_filename)

TDR_timeseries_filename = output_folder + '/' + "Analysis_figures/Timeseriesonly" + "/" + str(treatment) + "__TDR_timeseriesI2.pdf"
TimeSeriesfigure(data, 'TDR', num_experiments, TDR_timeseries_filename)

# %% Combination figure
combined_filename_I1 = output_folder + '/' + "Analysis_figures" + "/" + str(treatment) + "_I1_combined.pdf"
fig, ax = pyplot.subplots(2,1, figsize=(6,8)) # pyplot.subplots(rows, columns, figsize=(width,length))
donor_color = "red"
recipient_color = "blue"
transconjugant_color = "purple"
plasmid_free_donor_color = "grey"
for rep in range(4): 
    ax[0].plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['Donor'], linestyle='--', alpha = ALPHA, color=donor_color)
    ax[0].plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['Recipient'], linestyle='--', alpha = ALPHA, color=recipient_color)
    ax[0].plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['Transconjugant'], linestyle='--', alpha = ALPHA, color=transconjugant_color)
    ax[0].plot(data['experiments'][rep]['I1'][0]['time']/3600, data['experiments'][rep]['I1'][0]['PlasmidFreeDonor'], linestyle='--', alpha = ALPHA, color=plasmid_free_donor_color)
ax[0].plot(resultODE['time']/3600, resultODE['Donor'], linestyle='-', color=donor_color, label='numerical D')
ax[0].plot(resultODE['time']/3600, resultODE['Recipient'], linestyle='-', color=recipient_color, label='numerical R')
ax[0].plot(resultODE['time']/3600, resultODE['Transconjugant'], linestyle='-', color=transconjugant_color, label='numerical T')
ax[0].plot(resultODE['time']/3600, resultODE['PlasmidFreeDonor'], linestyle='-', color=plasmid_free_donor_color, label='numerical F')
ax[0].plot([], [], ' ', label=r'$\psi_D=$'+str(data['parameters']['psiD']))
ax[0].plot([], [], ' ', label=r'$\psi_R=$'+str(data['parameters']['psiR']))
ax[0].plot([], [], ' ', label=r'$\psi_T=$'+str(data['parameters']['psiT']))
ax[0].plot([], [], ' ', label=r'$\psi_F=$'+str(data['parameters']['psiF']))
ax[0].plot([], [], ' ', label=r'$\gamma_{D}=$'+str(data['parameters']['gammaDR']))
ax[0].plot([], [], ' ', label=r'$\gamma_{T}=$'+str(data['parameters']['gammaTR']))
ax[0].plot([], [], ' ', label=r'$\tau_{D}=$'+str(data['parameters']['tauD']))
ax[0].plot([], [], ' ', label=r'$\tau_{T}=$'+str(data['parameters']['tauT']))
ax[0].plot([], [], ' ', label='selected incubation times:')
ax[0].axvline(x=data['parameters']['I2_t'], linestyle='--', color = 'darkgrey', label=('I2 =' + str(round(data['parameters']['I2_t'],2))))
ax[0].axvline(x=data['parameters']['I1_t'], linestyle='--', color = 'grey', label=('I1 =' + str(round(data['parameters']['I1_t'],2))))
ax[0].axvline(x=data['LDM']['time_hours'], linestyle='-.', color = 'brown', label=r'$t_{LDM}=$' + str(round(data['LDM']['time_hours'],2)))
ax[0].axvline(x=data['SIM']['time_hours'], linestyle='-.', color = 'orange', label=r'$t_{x}=$' + str(round(data['SIM']['time_hours'],2)))
ax[0].set_ylabel("population density (cells/mL)")
ax[0].set_yscale('log')
ax[0].set_ylim(bottom = 0.9)
ax[0].legend(bbox_to_anchor=(1.01,1), facecolor='white', edgecolor='none', framealpha=0.9, prop={'size': 8})

ax[1].plot(data['SIM']['time_series'], ([data['parameters']['gammaDR']]*len(data['SIM']['time_series'])), alpha = 0)
for rep in range(4):
    ax[1].plot(data['TDR']['time_series'], data['TDR']['TimeSeries'][rep], linestyle='-', alpha = ALPHA, color='cyan')
    ax[1].plot(data['ASM']['time_series'], data['ASM']['TimeSeries'][rep], linestyle='-', alpha = ALPHA, color='green')
    ax[1].plot(data['SIM']['time_series'], data['SIM']['TimeSeries'][rep], linestyle='-', alpha = ALPHA, color='orange')
    ax[1].plot(data['LDM']['time_series'], data['LDM']['TimeSeries'][rep], linestyle='-', alpha = ALPHA, color='brown')
ax[1].set_yscale('log')
ax[1].plot([], [], linestyle='-', color='brown', label = 'LDM')
ax[1].plot([], [], linestyle='-', color='orange', label = 'SIM')
ax[1].plot([], [], linestyle='-', color='cyan', label = 'TDR')
ax[1].plot([], [], linestyle='-', color='green', label = 'ASM')
ax[1].axhline(y=data['parameters']['gammaDR'], linestyle='--', color = 'grey', label=r'$\gamma_{DR}=$'+str(data['parameters']['gammaDR']))
ax[1].axvline(x=data['parameters']['I2_t'], linestyle='--', color = 'darkgrey')
ax[1].axvline(x=data['parameters']['I1_t'], linestyle='--', color = 'grey')
ax[1].axvline(x=data['LDM']['time_hours'], linestyle='-.', color = 'brown')
ax[1].axvline(x=data['SIM']['time_hours'], linestyle='-.', color = 'orange')
ax[1].set_ylabel("Conjugation rate estimate")
ax[1].set_xlabel("sampling time (hours)")
ax[1].legend(bbox_to_anchor=(1.01,1), facecolor='white', edgecolor='none', framealpha=0.9, prop={'size': 8})
fig.savefig(combined_filename_I1, bbox_inches='tight')    


# %% Output criteria_t estimates
Exp = [None] * (num_experiments+1)
LDM = [None] * (num_experiments+1)
ASM = [None] * (num_experiments+1)
SIM = [None] * (num_experiments+1)
TDR = [None] * (num_experiments+1)
LDM_p0 = [None] * (num_experiments+1)
sampled_Tt = [None] * (num_experiments+1)

Exp_pop = [None] * (num_I2_populations * num_experiments)
LDM_pop = [None] * (num_I2_populations * num_experiments)
ASM_pop = [None] * (num_I2_populations * num_experiments)
SIM_pop = [None] * (num_I2_populations * num_experiments)
TDR_pop = [None] * (num_I2_populations * num_experiments)
p0_pop = [None] * (num_I2_populations * num_experiments)
Tt_pop = [None] * (num_I2_populations * num_experiments)

Exp[0] = 'time'
LDM[0] = data['LDM']['time_hours']
SIM[0] = data['SIM']['time_hours']
ASM[0] = data['ASM']['time_hours']
TDR[0] = data['TDR']['time_hours']
LDM_p0[0] = "n/a"
sampled_Tt[0] = "n/a"

for rep in range(num_experiments):
    Exp[rep+1] = 'E' + str(rep) + '_I1'
    Exp_I2_prefix = 'E' + str(rep) + '_I2_P'
    LDM[rep+1] = data['LDM'][rep]['criteria_t']['LDM']
    ASM[rep+1] = data['ASM'][rep]['criteria_t']['ASM']
    SIM[rep+1] = data['SIM'][rep]['criteria_t']['SIM']
    TDR[rep+1] = data['TDR'][rep]['criteria_t']['TDR']
    LDM_p0[rep+1] = data['LDM'][rep]['criteria_t']['p0']
    sampled_Tt[rep+1] = data['ASM'][rep]['criteria_t']['Tt']
    for pop in range(num_I2_populations):
        Exp_pop[((rep*num_I2_populations)+pop)] = Exp_I2_prefix + str(pop)
        Tt_pop[((rep*num_I2_populations)+pop)] = data['ASM'][rep]['criteria_t'][pop]['Tt']
        if data['SIM'][rep]['criteria_t'][pop]['Dt'] > 0 and data['SIM'][rep]['criteria_t'][pop]['Rt'] > 0: 
            SIM_pop[((rep*num_I2_populations)+pop)] = data['SIM'][rep]['criteria_t'][pop]['SIM']
            ASM_pop[((rep*num_I2_populations)+pop)] = data['ASM'][rep]['criteria_t'][pop]['ASM']
            TDR_pop[((rep*num_I2_populations)+pop)] = data['TDR'][rep]['criteria_t'][pop]['TDR']
            
        
Exp.extend(Exp_pop)
LDM.extend(LDM_pop)
ASM.extend(ASM_pop)
SIM.extend(SIM_pop)
TDR.extend(TDR_pop)
LDM_p0.extend(p0_pop)
sampled_Tt.extend(Tt_pop)

zippedList = list(zip(Exp, LDM, SIM, ASM, TDR, LDM_p0, sampled_Tt))
AllEstimates = pd.DataFrame(zippedList, columns = ['Experiments', 'LDM', 'SIM', 'ASM', 'TDR', 'p0', 'Tt'])

if not os.path.exists(output_folder + '/' + 'Analysis_estimates'):
    os.makedirs(output_folder + '/' + 'Analysis_estimates')
    
outfile = output_folder + '/' + 'Analysis_estimates' + "/" + treatment + "_estimates.csv"
AllEstimates.to_csv(outfile)
# %% Output time series estimates
for t in range(len(data['SIM']['time_series'])):    
    Exp = [None] * (num_experiments+1)
    LDM = [None] * (num_experiments+1)
    ASM = [None] * (num_experiments+1)
    SIM = [None] * (num_experiments+1)
    TDR = [None] * (num_experiments+1)
    LDM_p0 = [None] * (num_experiments+1)
    sampled_Tt = [None] * (num_experiments+1)
    
    Exp_pop = [None] * (num_I2_populations * num_experiments)
    LDM_pop = [None] * (num_I2_populations * num_experiments)
    ASM_pop = [None] * (num_I2_populations * num_experiments)
    SIM_pop = [None] * (num_I2_populations * num_experiments)
    TDR_pop = [None] * (num_I2_populations * num_experiments)
    p0_pop = [None] * (num_I2_populations * num_experiments)
    Tt_pop = [None] * (num_I2_populations * num_experiments)
    
    Exp[0] = 'time'
    LDM[0] = data['LDM']['time_hours']
    SIM[0] = data['SIM']['time_hours']
    ASM[0] = data['ASM']['time_hours']
    TDR[0] = data['TDR']['time_hours']
    LDM_p0[0] = "n/a"
    sampled_Tt[0] = "n/a"
    
    if t >= len(data['LDM']['time_series']):
        for rep in range(num_experiments):
            Exp[rep+1] = 'E' + str(rep) + '_I1'
            Exp_I2_prefix = 'E' + str(rep) + '_I2_P'
            LDM[rep+1] = None
            ASM[rep+1] = data['ASM']['TimeSeries'][rep][t]
            SIM[rep+1] = data['SIM']['TimeSeries'][rep][t]
            TDR[rep+1] = data['TDR']['TimeSeries'][rep][t]
            LDM_p0[rep+1] = None
            sampled_Tt[rep+1] = data['ASM'][rep][data['ASM']['time_series'][t]]['Tt']
    else:
        for rep in range(num_experiments):
            Exp[rep+1] = 'E' + str(rep) + '_I1'
            Exp_I2_prefix = 'E' + str(rep) + '_I2_P'
            LDM[rep+1] = data['LDM']['TimeSeries'][rep][t]
            ASM[rep+1] = data['ASM']['TimeSeries'][rep][t]
            SIM[rep+1] = data['SIM']['TimeSeries'][rep][t]
            TDR[rep+1] = data['TDR']['TimeSeries'][rep][t]
            LDM_p0[rep+1] = data['LDM'][rep][data['LDM']['time_series'][t]]['p0']
            sampled_Tt[rep+1] = data['ASM'][rep][data['ASM']['time_series'][t]]['Tt']
            for pop in range(num_I2_populations):
                Exp_pop[((rep*num_I2_populations)+pop)] = Exp_I2_prefix + str(pop)
                Tt_pop[((rep*num_I2_populations)+pop)] = data['ASM'][rep]['I2'][pop][data['ASM']['time_series'][t]]['Tt']
                if data['SIM'][rep]['I2'][pop][data['SIM']['time_series'][t]]['Dt'] > 0 and data['SIM'][rep]['I2'][pop][data['SIM']['time_series'][t]]['Rt'] > 0: 
                    SIM_pop[((rep*num_I2_populations)+pop)] = data['SIM'][rep]['I2'][pop][data['SIM']['time_series'][t]]['SIM']
                    ASM_pop[((rep*num_I2_populations)+pop)] = data['ASM'][rep]['I2'][pop][data['ASM']['time_series'][t]]['ASM']
                    TDR_pop[((rep*num_I2_populations)+pop)] = data['TDR'][rep]['I2'][pop][data['TDR']['time_series'][t]]['TDR']
                    
     #data['SIM'][2]['I2'][0][3.0]['Dt']
    
    Exp.extend(Exp_pop)
    LDM.extend(LDM_pop)
    ASM.extend(ASM_pop)
    SIM.extend(SIM_pop)
    TDR.extend(TDR_pop)
    LDM_p0.extend(p0_pop)
    sampled_Tt.extend(Tt_pop)
    
    zippedList = list(zip(Exp, LDM, SIM, ASM, TDR, LDM_p0, sampled_Tt))
    AllEstimates = pd.DataFrame(zippedList, columns = ['Experiments', 'LDM', 'SIM', 'ASM', 'TDR', 'p0', 'Tt'])
    
    if not os.path.exists(output_folder + '/' + 'Analysis_estimates'):
        os.makedirs(output_folder + '/' + 'Analysis_estimates')
    if not os.path.exists(output_folder + '/' + 'Analysis_estimates' + '/' + 'Time_series'):
        os.makedirs(output_folder + '/' + 'Analysis_estimates' + '/' + 'Time_series')
        
    outfile = output_folder + '/' + 'Analysis_estimates' + "/" + 'Time_series' + "/" + treatment + "_" + str(data['TDR']['time_series'][t]) + "_estimates.csv"
    AllEstimates.to_csv(outfile)
    










