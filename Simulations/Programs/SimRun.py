#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 13:37:59 2020

@author: oliviakosterlitz
"""

import argparse
import sys
import pandas
import gillespy2
from gillespy2 import NumPySSASolver
from SimFunctions import *
import numpy as np

script_description = """
Code will use run stochastic simulation by running the GillesPy packages
NumPySSASolver to determine incubation times given user specified conditions. 
Input for the run is as follows:
    output_filename
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
    t
"""

parser = argparse.ArgumentParser(
         description = script_description,
         formatter_class = argparse.RawDescriptionHelpFormatter)
         
parser.add_argument("-o", "--output_filename", required = True, help = "The output filename")
parser.add_argument("-pD", "--psiD", required = True, help = "The growth rate of the donor")
parser.add_argument("-pR", "--psiR", required = True, help = "The growth rate of the recipient")
parser.add_argument("-pT", "--psiT", required = True, help = "The growth rate of the transconjugant")
parser.add_argument("-pF", "--psiF", required = True, help = "The growth rate of the plasmid-free recipient")
parser.add_argument("-gDR", "--gammaDR", required = True, help = "The transfer rate from donor to recipient")
parser.add_argument("-gTR", "--gammaTR", required = True, help = "The transfer rate from transconjugant to recipient")
parser.add_argument("-gDF", "--gammaDF", required = True, help = "The transfer rate from donor to plasmid-free recipient")
parser.add_argument("-gTF", "--gammaTF", required = True, help = "The transfer rate from transconjugant to plasmid-free recipient")
parser.add_argument("-tD", "--tauD", required = True, help = "The segregation rate from donor to plasmid-free recipient")
parser.add_argument("-tT", "--tauT", required = True, help = "The segregation rate from transconjugant to recipient")
parser.add_argument("-D0", "--D_0", required = True, help = "The initial density of the donor")
parser.add_argument("-R0", "--R_0", required = True, help = "The initial density of the recipient")
parser.add_argument("-T0", "--T_0", required = True, help = "The initial density of the transconjugant")
parser.add_argument("-F0", "--F_0", required = True, help = "The initial density of the plasmid-free donor")
parser.add_argument("-t", "--t", required = True, help = "incubation time")

args = parser.parse_args()
                
output_filename = args.output_filename           
psiD = float(args.psiD)
psiR = float(args.psiR)
psiT = float(args.psiT)
psiF = float(args.psiF)
gammaDR = float(args.gammaDR)
gammaTR = float(args.gammaTR)
gammaDF = float(args.gammaDF)
gammaTF = float(args.gammaTF)
tauD = float(args.tauD)
tauT = float(args.tauT)
D_0 = float(args.D_0)
R_0 = float(args.R_0)
T_0 = float(args.T_0)
F_0 = float(args.F_0)
t = float(args.t)

# Convert the rates to "per second". This is necessary for using gillespy2.
runtime_sec=time_hourtosecond(t)

# growth rates
psiD_sec=rate_hourtosecond(psiD)
psiR_sec=rate_hourtosecond(psiR)
psiT_sec=rate_hourtosecond(psiT)
psiF_sec=rate_hourtosecond(psiF)

# transfer rates
gammaDR_sec=rate_hourtosecond(gammaDR)
gammaTR_sec=rate_hourtosecond(gammaTR)
gammaDF_sec=rate_hourtosecond(gammaDF)
gammaTF_sec=rate_hourtosecond(gammaTF)

# segregation rates
tauD_sec=rate_hourtosecond(tauD)
tauT_sec=rate_hourtosecond(tauT)

# Initiate the model 
SSAmodel = PlasmidDynamics( runtime=runtime_sec,
                            psiD=psiD_sec,
                            psiR=psiR_sec,
                            psiT=psiT_sec,
                            psiF=psiF_sec,
                            gammaDR=gammaDR_sec,
                            gammaTR=gammaTR_sec,
                            gammaDF=gammaDF_sec,
                            gammaTF=gammaTF_sec,
                            tauD=tauD_sec,
                            tauT=tauT_sec,
                            D_init=D_0,
                            R_init=R_0,
                            T_init=T_0,
                            F_init=F_0)


resultSSA = SSAmodel.run(solver=NumPySSASolver)

# Save the file
to_csv_custom(resultSSA, nametag = output_filename)








