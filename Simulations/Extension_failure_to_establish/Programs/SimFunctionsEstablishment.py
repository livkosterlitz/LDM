#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 2 14:31:32 2020

@author: oliviakosterlitz
"""
import gillespy2
import numpy as np

# %% #### Functions for SimSetup.py ####

def rate_hourtosecond (rate):
    secondrate = rate/3600
    return (secondrate)

def time_hourtosecond (time):
    seconds = int(time*3600)
    return (seconds)

class PlasmidDynamics(gillespy2.Model):
     def __init__(self, runtime, 
                          psiD, 
                          psiR, 
                          psiT,
                          psiF,
                          gammaDR, 
                          gammaTR,
                          gammaDF,
                          gammaTF,
                          tauD,
                          tauT,
                          D_init, 
                          R_init, 
                          T_init,
                          F_init,
                          parameter_values=None):
            # initialize Model
            gillespy2.Model.__init__(self, name="PlasmidDynamics")
            
            # parameters
            rate1 = gillespy2.Parameter(name='psiD', expression= psiD)
            rate2 = gillespy2.Parameter(name='psiR', expression= psiR)
            rate3 = gillespy2.Parameter(name='psiT', expression = psiT)
            rate4 = gillespy2.Parameter(name='psiF', expression = psiF)
            rate5 = gillespy2.Parameter(name='gammaDR', expression=gammaDR)
            rate6 = gillespy2.Parameter(name='gammaTR', expression=gammaTR)
            rate7 = gillespy2.Parameter(name='gammaDF', expression=gammaDF)
            rate8 = gillespy2.Parameter(name='gammaTF', expression=gammaTF)
            rate9 = gillespy2.Parameter(name='tauD', expression=tauD)
            rate10 = gillespy2.Parameter(name='tauT', expression=tauT)
            self.add_parameter([rate1,rate2,rate3,rate4,rate5,rate6,rate7,rate8,rate9,rate10])
            
            # species
            D = gillespy2.Species(name='Donor', initial_value= D_init)
            R = gillespy2.Species(name='Recipient', initial_value= R_init)
            T = gillespy2.Species(name='Transconjugant', initial_value= T_init)
            F = gillespy2.Species(name='PlasmidFreeDonor', initial_value= F_init)
            self.add_species([D,R,T,F])
            
            ## reactions
            
            # growth of donors
            r1 = gillespy2.Reaction(name="r1",reactants={D:1}, products={D:2},
                   rate=rate1)
            
            # growth of recipients
            r2 = gillespy2.Reaction(name="r2",reactants={R:1}, products={R:2},
                    rate=rate2)
            
            # growth of transconjugants
            r3 = gillespy2.Reaction(name="r3",reactants={T:1}, products={T:2},
                    rate=rate3)

            # growth of plasmid free donors
            r4 = gillespy2.Reaction(name="r4",reactants={F:1}, products={F:2},
                    rate=rate4)
            
            # conjugation between donors and recipients
            r5 = gillespy2.Reaction(name="r5",reactants={D:1, R:1}, products={D:1,T:1},
                    rate=rate5)
            
            # conjugation between transconjugants and recipients
            r6 = gillespy2.Reaction(name="r6",reactants={T:1, R:1}, products={T:2},
                    rate=rate6)
            
            # conjugation between donors and plasmid free donors
            r7 = gillespy2.Reaction(name="r7",reactants={D:1, F:1}, products={D:2},
                    rate=rate7)

            # conjugation between transconjugants and plasmid free donors
            r8 = gillespy2.Reaction(name="r8",reactants={T:1, F:1}, products={T:1,D:1},
                    rate=rate8)
            
            # segregation event from donors
            r9 = gillespy2.Reaction(name="r9",reactants={D:1}, products={F:1},
                    rate=rate9)

            # segregation event from transconjugants
            r10 = gillespy2.Reaction(name="r10",reactants={T:1}, products={R:1},
                    rate=rate10)
            
            self.add_reaction([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10])
            self.timespan(np.linspace(0,runtime,runtime+1))
        
def ODEfigure(population, 
           psiD=None, 
           psiR=None, 
           psiT=None,
           psiF=None,
           gammaDR=None, 
           gammaTR=None,
           gammaDF=None,
           gammaTF=None,
           tauD=None,
           tauT=None,
           inc_T=None,
           inc_Pop_size=None,
           ODE_solution=False, 
           filename=None, 
           show=True, 
           legends=True):
    import numpy as np
    import matplotlib.pyplot as pyplot
    import seaborn
    fig, ax = pyplot.subplots(1,1, figsize=(8,5)) # pyplot.subplots(rows, columns, figsize=(width,length))
    donor_color = "red"
    recipient_color = "blue"
    transconjugant_color = "purple"
    plasmid_free_donor_color = "grey"
    for pop in population:
        ax.plot(pop['time']/3600, pop['Donor'], linestyle='--', color=donor_color)
        ax.plot(pop['time']/3600, pop['Recipient'], linestyle='--', color=recipient_color)
        ax.plot(pop['time']/3600, pop['Transconjugant'], linestyle='--', color=transconjugant_color)
    if ODE_solution != False:
        ax.plot(ODE_solution['time']/3600, ODE_solution['Donor'], linestyle='-', color=donor_color, label='D')
        ax.plot(ODE_solution['time']/3600, ODE_solution['Recipient'], linestyle='-', color=recipient_color, label='R')
        ax.plot(ODE_solution['time']/3600, ODE_solution['Transconjugant'], linestyle='-', color=transconjugant_color, label='T')
        ax.plot(ODE_solution['time']/3600, ODE_solution['PlasmidFreeDonor'], linestyle='-', color=plasmid_free_donor_color, label='F')
    if psiD != None:
        ax.plot([], [], ' ', label=r'$\psi_D=$'+str(psiD))
        ax.plot([], [], ' ', label=r'$\psi_R=$'+str(psiR))
        ax.plot([], [], ' ', label=r'$\psi_T=$'+str(psiT))
        ax.plot([], [], ' ', label=r'$\psi_F=$'+str(psiF))
        ax.plot([], [], ' ', label=r'$\gamma_{D}=$'+str(gammaDR))
        ax.plot([], [], ' ', label=r'$\gamma_{T}=$'+str(gammaTR))
        ax.plot([], [], ' ', label=r'$\tau_{D}=$'+str(tauD))
        ax.plot([], [], ' ', label=r'$\tau_{T}=$'+str(tauT))
        ax.plot([], [], ' ', label='selected incubation times:')
        ax.axvline(x=inc_T, linestyle='--', color = 'darkgrey', label=('I2 =' + str(round(inc_T,2))))
        ax.axvline(x=inc_Pop_size, linestyle='--', color = 'grey', label=('I1 =' + str(round(inc_Pop_size,2))))
    ax.set_ylabel("population density (cells/mL)")
    ax.set_xlabel("time (hours)")
    ax.set_yscale('log')
    ax.set_ylim(bottom = 0.9)
    ax.legend(bbox_to_anchor=(1.01,1), facecolor='white', edgecolor='none', framealpha=0.9, prop={'size': 8})
    if(filename):
        fig.savefig(filename, bbox_inches='tight')

def output_name (treatment, experiment, incubation, population):
    experiment = 'E' + str(experiment)
    population = 'P' + str(population)
    output_name = treatment + '_' + experiment + '_' + incubation + '_' + population   
    return (output_name)
    
def shell_line (SimRun_path, df, f, i, e, incubation, replicate):
    treatment = str(df['Treatment_ID'][i])
    f.write('time' + ' ')
    f.write('python3' + ' ')
    f.write(SimRun_path + ' ')
    f.write('--output_filename ' + output_name(treatment, e, incubation, replicate) + ' ')
    f.write('--psiD ' + str(df['psiD'][i]) + ' ')
    f.write('--psiR ' + str(df['psiR'][i]) + ' ')
    f.write('--psiT ' + str(df['psiT'][i]) + ' ')
    f.write('--psiF ' + str(df['psiF'][i]) + ' ')
    f.write('--gammaDR ' + str(df['gammaDR'][i]) + ' ')
    f.write('--gammaTR ' + str(df['gammaTR'][i]) + ' ')
    f.write('--gammaDF ' + str(df['gammaDF'][i]) + ' ')
    f.write('--gammaTF ' + str(df['gammaTF'][i]) + ' ')
    f.write('--tauD ' + str(df['tauD'][i]) + ' ')
    f.write('--tauT ' + str(df['tauT'][i]) + ' ')
    if incubation == 'I3': 
        f.write('--D_0 ' + str(0) + ' ')
        f.write('--R_0 ' + str(0) + ' ')
        f.write('--T_0 ' + str(df['D_0'][i]) + ' ')
        f.write('--F_0 ' + str(0) + ' ')
    else:
        f.write('--D_0 ' + str(df['D_0'][i]) + ' ')
        f.write('--R_0 ' + str(df['R_0'][i]) + ' ')
        f.write('--T_0 ' + str(df['T_0'][i]) + ' ')
        f.write('--F_0 ' + str(df['F_0'][i]) + ' ')
    if incubation == 'I2':     
        f.write('--t ' + str(df['Cutoff_I2_t'][i]) + ' ')
    else:
        f.write('--t ' + str(df['Cutoff_I1_t'][i]) + ' ')
    f.write('\n')        
# %% #### Functions for SimAnalysis.py ####
# Calculates LDM conjugation rate
def LDMmetric (runtime, p0, D_0, R_0, D_t, R_t):
    num = (-np.log(p0)) * (np.log(D_t * R_t) - np.log(D_0 * R_0))
    den = runtime * ((D_t * R_t) - (D_0 * R_0))
    LDM_gamma_estimate=num/den
    return(LDM_gamma_estimate)   

def correctedLDMmetric(D_0, R_0, psiD, psiR, psiT, runtime, pi, Pnt):
    import scipy.stats
    import numpy as np
    factor1 = (-np.log(Pnt)*(psiD+psiR))/(D_0*R_0)
    factor1_den = np.exp((psiD+psiR)*runtime)
    factor2_den = scipy.special.hyp2f1(1, (psiD+psiR)/psiT, 1 + ((psiD+psiR)/psiT), pi/(pi-1))
    factor3_den = scipy.special.hyp2f1(1, (psiD+psiR)/psiT, 1 + ((psiD+psiR)/psiT), (pi/(pi-1))*np.exp(-psiT*runtime))
    LDM = factor1 * (1/(factor1_den * factor2_den - factor3_den))
    return(LDM)

# Calculates SIM conjugation rate
def SIMmetric (runtime, N_0, N_t, D_t, R_t, T_t):
    SIM_gamma_estimates=((np.log(N_t/N_0)/runtime)*(np.log(1+(N_t*T_t)/(R_t * D_t))) * (1/(N_t - N_0)))
    return(SIM_gamma_estimates)

# Calculates ASM conjugation rate
def ASMmetric (runtime, psiD, psiR, psiT, D_0, R_0, D_t, R_t, T_t):
    ASM_gamma_estimates=(psiD+psiR-psiT)*(T_t)*(1/((D_t*R_t)-(D_0*R_0*np.exp(psiT*runtime))))
    return(ASM_gamma_estimates)

def growthrate (time, initial, final):
    growthrate = (1/time)*(np.log(final) - np.log(initial))
    return(growthrate)

def TDRmetric (runtime, D_t, R_t, T_t):
    TDR_tao_estimates=((1/runtime)*(T_t/(D_t*R_t)))
    return (TDR_tao_estimates)

# Time series figure
def TimeSeriesfigure(data, estimate, num_experiments, filename, timeseries = None):
    import numpy as np
    import matplotlib.pyplot as pyplot
    import seaborn
    import matplotlib.ticker as mticker
    if estimate == 'LDM':
        color_select = 'brown'
    elif estimate == 'SIM':
        color_select = 'orange'
    elif estimate == 'ASM':
        color_select = 'green'
    else:
        color_select = 'cyan'
    fig, ax = pyplot.subplots(1,1, figsize=(6,3.5)) # pyplot.subplots(rows, columns, figsize=(width,length))
    ax.plot(data['SIM']['time_series'], ([data['parameters']['gammaDR']]*len(data['SIM']['time_series'])), alpha = 0)
    for rep in range(num_experiments):
        ax.plot(data[estimate]['time_series'], data[estimate]['TimeSeries'][rep], linestyle='-', alpha = 0.25, color=color_select)
    ax.set_yscale('log')
    ax.plot([], [], linestyle='-', color=color_select, label = estimate)
    ax.axhline(y=data['parameters']['gammaDR'], linestyle='--', color = 'grey', label=r'$\gamma_{DR}=$'+str(data['parameters']['gammaDR']))
    ax.axvline(x=data['parameters']['I2_t'], linestyle='--', color = 'darkgrey', label=('I2 =' + str(round(data['parameters']['I2_t'],2))))
    ax.axvline(x=data['parameters']['I1_t'], linestyle='--', color = 'grey', label=('I1 =' + str(round(data['parameters']['I1_t'],2))))
    ax.axvline(x=data['LDM']['time_hours'], linestyle='-.', color = 'brown', label=r'$t_{LDM}=$' + str(round(data['LDM']['time_hours'],2)))
    ax.axvline(x=data['SIM']['time_hours'], linestyle='-.', color = 'orange', label=r'$t_{x}=$' + str(round(data['SIM']['time_hours'],2)))
    ax.set_ylabel(estimate + " Estimate")
    ax.set_xlabel("sampling time (hours)")
    ax.legend(bbox_to_anchor=(1.01,1), facecolor='white', edgecolor='none', framealpha=0.9, prop={'size': 8})
    fig.savefig(filename, bbox_inches='tight')
       
def to_csv_custom(data, nametag=None):
      """
      outputs the Results to one or more .csv files in a new directory.

      :param nametag: allows the user to optionally "tag" the directory and included files. Defaults to the model
      name.
      :type nametag: str
      :param path: the location for the new directory and included files. Defaults to model location.
      :type path: str
      """
      import csv
      import os

      identifier = nametag
      directory = os.getcwd()
      
      # multiple trajectories
      if isinstance(data.data, list):
          for i, trajectory in enumerate(data.data):  # write each CSV file
              filename = os.path.join(directory, str(identifier) +".csv")
              field_names = []
              for species in trajectory:  # build the header
                  field_names.append(species)
              with open(filename, 'w', newline='') as csv_file:
                  csv_writer = csv.writer(csv_file)
                  csv_writer.writerow(field_names)  # write the header
                  for n,time in enumerate(trajectory['time']):  # write all lines of the CSV file
                      this_line=[]
                      for species in trajectory:  # build one line of the CSV file
                          this_line.append(trajectory[species][n])
                      csv_writer.writerow(this_line)  
    
def dilution_function(number):
    import scipy.stats
    dilution_list = []
    dilution_list.append(number)
    while number > 0:
        number = number * 0.1
        number = (scipy.stats.poisson.rvs(mu = number, size = 1))[0]
        dilution_list.append(number)
    return(dilution_list)

def plating_function(dilution_series):
    import numpy as np
    import scipy.stats
    dilution = np.array(dilution_series) * 0.1
    plating_list = [(scipy.stats.poisson.rvs(mu = xi, size = 1))[0] for xi in dilution]
    return(plating_list)

def extinction_function(plating_list, extinction_probability):
    import scipy.stats
    after_extinction = [scipy.stats.binom.rvs(xi, 1-extinction_probability, size=1)[0] for xi in plating_list]
    return(after_extinction)

def density_calculation_function(plating_series):
    x = [e for e in plating_series if 10 <= e <= 300]
    if len(x) > 0 :
        counts = min(x)
    elif max(plating_series) == 0:
       return(0)
    elif max(plating_series) < 10:
        counts = max(plating_series)
    else: 
        counts = min(plating_series)
    dilution = plating_series.index(counts)
    density = counts * (10**dilution) * 10
    return(density)

def extinction_correction_function(density, extinction_probability):
    correction = round(density * (1/(1-extinction_probability)))
    return(correction)

def dilution_plating_ext (density, extinction_probability):
    dil = dilution_function(density)
    plated = plating_function(dil)
    ext = extinction_function(plated, extinction_probability)
    cal = density_calculation_function(ext)
    correct = extinction_correction_function(cal, extinction_probability)
    return(correct)

def density_init (density, extinction_probability):
    correct = dilution_plating_ext(density, extinction_probability)
    if correct == 0: 
        while correct == 0:
            correct = dilution_plating_ext(density, extinction_probability)
    return(correct)


# import pandas as pd
# import numpy as np
# import statistics
# density_list = [1e9, 1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1]
# extinction_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# mean_frame = {'Density': density_list}
# var_frame = {'Density': density_list}
# for i in range(len(mean_frame['Density'])):
#     for l in extinction_list:
#         mean_frame[l] = [0] * len(mean_frame['Density'])
#         var_frame[l] = [0] * len(var_frame['Density'])

# for i in range(len(mean_frame['Density'])):
#     for l in extinction_list:
#         dlist = []
#         number = mean_frame['Density'][i]
#         extinction_probability = l
#         #df[l][i] = [number, l]
#         for r in range(10000):
#             ds = dilution_function(number)
#             plating_list = plating_function(ds)
#             ext = extinction_function(plating_list, extinction_probability) 
#             densitycal = density_calculation_function(ext)
#             densitycorrected = extinction_correction_function(densitycal, extinction_probability)
#             dlist.append(densitycorrected)
#         mean_frame[l][i] = np.average(dlist)
#         var_frame[l][i] = np.var(dlist)
#         print('new treatment')
#         print(number)
#         print(l)
#         print(np.average(dlist))
#         print(np.var(dlist))

# mean_df = pd.DataFrame(mean_frame)
# var_df = pd.DataFrame(var_frame)
# mean_df.to_csv('mean_r10000.csv')
# var_df.to_csv('variance_r10000.csv')
        





