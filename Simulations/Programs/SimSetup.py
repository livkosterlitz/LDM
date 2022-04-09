#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 2 13:35:44 2020

@author: oliviakosterlitz
"""

script_description = """
Code will use a deterministic simulation to determine incubation 
times given user specified conditions, which is provided by the 
input CSV file with treatement information:
    Treatment_ID : T# where the # is a unique identifier 	
    psiD	
    psiR	
    psiT	
    psiF	
    gammaDR	
    gammaTR	
    gammaDF	
    gammaTF	
    tauD	
    tauT
    D_0	
    R_0	
    T_0	
    F_0	
    Cutoff_I1	
    Cutoff_I2
    Num_sims_I1
    Num_sims_I2
    Num_of_experiments
"""

import argparse
import os
import sys
import pandas
import numpy as np
import gillespy2
from gillespy2 import ODESolver
from SimFunctions import *

#troubleshooting 
# from GitHub.Programs.SimFunctions import *
# SimRun_path = "SimRun.py"
# input_CSV = "/Users/oliviakosterlitz/Dropbox/Pop3/Simulations/Submission/GitHub/BasicRun/Treatments/SimSetup_inputCSV_example.csv"           
# output_folder = "/Users/oliviakosterlitz/Dropbox/Pop3/Simulations/Submission/GitHub/BasicRun/"

parser = argparse.ArgumentParser(
         description = script_description,
         formatter_class = argparse.RawDescriptionHelpFormatter)
         
parser.add_argument("-s", "--SimRun_path", required = True, help = "Path to the SimRun.py program")
parser.add_argument("-c", "--csv_input", required = True, help = "Input CSV treatment information as column headers")
parser.add_argument("-o", "--output_folder", required = True, help = "Output folder")

args = parser.parse_args()

SimRun_path = args.SimRun_path
input_CSV = args.csv_input
output_folder = args.output_folder

df = pandas.read_csv(input_CSV)

ODE_test_time = 24
ODE_runtime_sec=time_hourtosecond(ODE_test_time)

if os.path.isdir(output_folder + '/' + "ODE_figures") == False:
    os.mkdir(output_folder + '/' + "ODE_figures")
if os.path.isdir(output_folder + '/' + "ODE_data") == False:
    os.mkdir(output_folder + '/' + "ODE_data")

T_cutoff_t_list = [None]*len(df['Treatment_ID'])
Pop_size_cutoff_t_list = [None]*len(df['Treatment_ID'])

for i in range(len(df['Treatment_ID'])):
    
    # growth rates (psi)
    psiD=df['psiD'][i]
    psiR=df['psiR'][i]
    psiT=df['psiT'][i]
    psiF=df['psiF'][i]
    psiD_sec=rate_hourtosecond(psiD)
    psiR_sec=rate_hourtosecond(psiR)
    psiT_sec=rate_hourtosecond(psiT)
    psiF_sec=rate_hourtosecond(psiF)
    
    # transfer rates (gamma)
    gammaDR=df['gammaDR'][i]
    gammaTR=df['gammaTR'][i]
    gammaDF=df['gammaDF'][i]
    gammaTF=df['gammaTF'][i]
    gammaDR_sec=rate_hourtosecond(gammaDR)
    gammaTR_sec=rate_hourtosecond(gammaTR)
    gammaDF_sec=rate_hourtosecond(gammaDF)
    gammaTF_sec=rate_hourtosecond(gammaTF)    
    
    # segregation rates (tau)
    tauD=df['tauD'][i]
    tauT=df['tauT'][i]
    tauD_sec=rate_hourtosecond(tauD)
    tauT_sec=rate_hourtosecond(tauT)
    
    # initial densities    
    D_0=df['D_0'][i]
    R_0=df['R_0'][i]
    T_0=df['T_0'][i]
    F_0=df['F_0'][i]
    
    # the ODE population density that determines the incubation time decisions
    Pop_size_cutoff = df['Cutoff_I1'][i]
    T_cutoff = df['Cutoff_I2'][i]
    
    # initiate the plasmid population model 
    ODEmodel = PlasmidDynamics(runtime = ODE_runtime_sec,
                               psiD = psiD_sec,
                               psiR = psiR_sec,
                               psiT = psiT_sec,
                               psiF = psiF_sec,
                               gammaDR = gammaDR_sec,
                               gammaTR = gammaTR_sec,
                               gammaDF = gammaDF_sec,
                               gammaTF = gammaTF_sec,
                               tauD = tauD_sec,
                               tauT = tauT_sec,
                               D_init = D_0,
                               R_init = R_0,
                               T_init = T_0,
                               F_init = F_0)

    # run the ODE model
    resultODE = ODEmodel.run(solver=ODESolver)
    
    # Output the numerical solution
    to_csv_custom(resultODE, nametag = output_folder + '/' + "ODE_data" + "/" + str(df['Treatment_ID'][i]))
    
    # The long incubation is at the time when the population density is equal to the density provided (i.e., long_cutoff)
    total_population_array = resultODE[0]['Recipient']+resultODE[0]['Donor']+resultODE[0]['Transconjugant']+resultODE[0]['PlasmidFreeDonor']
    Pop_size_incubation = (resultODE[0]['time'][total_population_array >= Pop_size_cutoff][0])/3600
    Pop_size_cutoff_t_list[i] = Pop_size_incubation
    
    # The short incubation is at the time when the transconjugant density is equal to the density provided (i.e., short_cutoff) 
    T_cutoff_incubation = (resultODE[0]['time'][resultODE[0]['Transconjugant'] >= T_cutoff][0])/3600
    T_cutoff_t_list[i] = T_cutoff_incubation

    # make a figure for the treatment
    makefigure = ODEfigure([], 
                            psiD=psiD, 
                            psiR=psiR, 
                            psiT=psiT,
                            psiF=psiF,
                            gammaDR=gammaDR, 
                            gammaTR=gammaTR,
                            gammaDF=gammaDF,
                            gammaTF=gammaTF,
                            tauD=tauD,
                            tauT=tauT,
                            inc_T = T_cutoff_incubation,
                            inc_Pop_size = Pop_size_incubation,
                            ODE_solution = resultODE, 
                            filename = output_folder + '/' + "ODE_figures/" + "ODE_" + str(df['Treatment_ID'][i]) + ".pdf")
    
## export incubation times for each treatment ##
df['Cutoff_I1_t'] = Pop_size_cutoff_t_list
df['Cutoff_I2_t'] = T_cutoff_t_list

# make a folder to hold treatment + incubation info if the folder doesn't already exist
if os.path.isdir(output_folder + '/' + "Treatments") == False:
    os.mkdir(output_folder + '/' + 'Treatments')
df.to_csv(output_folder + '/' + "Treatments/" + str(input_CSV.split('/')[-1].split(".")[0]) + '_with_incubation.csv')

# make a folder to hold shell script if the folder doesn't already exist
if os.path.isdir(output_folder + '/' + "SimData") == False:
    os.mkdir(output_folder + '/' + 'SimData')

shellfolder = output_folder + "/SimData"

for i in range(len(df['Treatment_ID'])):

    # make a shell script for the treatment to run the simulations
    treatment = str(df['Treatment_ID'][i])
    
    treatment_folder = "SimData/" + treatment
    
    if os.path.isdir(output_folder + '/' + treatment_folder) == False:
        os.mkdir(output_folder + '/' + treatment_folder)
    
    f = open(output_folder + '/' + treatment_folder + '/' + treatment +'.sh', 'w')

    # write lines for each experiment
    for e in range(df['Num_of_experiments'][i]):

        # write lines for simulations with the population size cutoff
        for p in range(df['Num_sims_I1'][i]):
            incubation = "I1" # used to be IL for long incubation
            shell_line (SimRun_path, df, f, i, e, incubation, p)
    
        # write lines for simulations with the transconjugant size cutoff
        for t in range(df['Num_sims_I2'][i]):
            incubation = "I2" # used to be IS for short incubation
            shell_line (SimRun_path, df, f, i, e, incubation, t)

    f.close()
    







