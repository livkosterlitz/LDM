These programs create a stochastic simulation framework for modeling bacterial population dynamics to calculate conjugation rates. 

# The plasmid population model 

These programs implement a mathematical model describing the change in density of donors, recipients, transconjugants, and plasmid-free donors over time (given by dynamic variables *D<sub>t</sub>*, *R<sub>t</sub>*, *T<sub>t</sub>*, and *F<sub>t</sub>*, respectively). In the model, each population type grows exponentially at a population-specific growth rate (*&psi;<sub>D<sub>*, *&psi;<sub>R<sub>*, *&psi;<sub>T<sub>*, and *&psi;<sub>F<sub>*, respectively). In addition, each plasmid-containing population (donors and transconjugants) decreases as a result of plasmid loss due to segregation events. Transconjugants are transformed into plasmid-free recipients at the transconjugant segregation rate (*&tau;<sub>T<sub>*). Similarly, donors are transformed into plasmid-free donors at the donor segregation rate (*&tau;<sub>D<sub>*). Therefore, each plasmid-free population (plasmid-free donors and recipients) increases due to these segregation events at the population-specific segregation rates (*&tau;<sub>D<sub>*, *&tau;<sub>T<sub>*, respectively). Furthermore, the transconjugant density increases as a result of conjugation events from donors to recipients at the donor-to-recipient conjugation rate (*&gamma;<sub>DR<sub>*) and from existing transconjugants to recipients at the transconjugant-to-recipient conjugation rate (*&gamma;<sub>TR<sub>*). The recipient density decreases due to these conjugation events. Similarly, the donor density increases as a result of conjugation events at the transconjugant-to-plasmid-free-donor conjugation rate (*&gamma;<sub>TF<sub>*) and donors-to-plasmid-free-donor conjugation rate (*&gamma;<sub>DF<sub>*). The plasmid-free donor density decreases due to these conjugation events. 
<p align="center">
  <img width="700" src="https://github.com/livkosterlitz/LDMsimulations/blob/main/Images/plasmid_population_model.png">
</p>  

<p align="center">
  <img width="350" src="https://github.com/livkosterlitz/LDMsimulations/blob/main/Images/model_equations.png">
</p>

# Conjugation rate estimates
Currently, the programs estimate conjugation rate with 4 metrics: 
| Method | Acronym | Year  | Estimate |
| :--- | :--- | :---  | :--- |
| Levin _et. al._ method | [TDR](https://doi.org/10.1016/0147-619X(79)90043-X) | 1979  | <img width="100" src="https://github.com/livkosterlitz/LDMsimulations/blob/main/Images/TDR.png"> |
| Simonsen _et. al._ method <br>(aka. end-point method) | [SIM](https://doi.org/10.1099/00221287-136-11-2319) | 1990 | <img width="300" src="https://github.com/livkosterlitz/LDMsimulations/blob/main/Images/SIM.png"> |
| Approximate extended Simonsen method | [ASM](https://doi.org/10.1016/j.plasmid.2022.102627) | 2021 |  <img  width="400" src="https://github.com/livkosterlitz/LDMsimulations/blob/main/Images/ASM.png"> |
| Luria-Delbrück method | [LDM](https://doi.org/10.1101/2021.01.06.425583) | 2021 | <img  width="375" src="https://github.com/livkosterlitz/LDMsimulations/blob/main/Images/LDM.png"> |

# Program Workflow

|Steps| Script |Step description| 
| :--- | :--- | :--- | 
| [Step 1](#Step-1) | `SimSetup.py` | Pass simulation parameters to determine the necessary incubation time to run the stochastic simulation |
| [Step 2](#Step-2) | `SimRun.py` | Run the stochastic simulation with the calculated incubation time |
| [Step 3](#Step-3) | `SimAnalysis.py` | Combine the stochastic simulation results from many populations  |

A basic example is used to illustrate the workflow and provide details for all of the steps. 

# Program Requirements
The programs (`SimSetup.py`, `SimRun.py`, `SimAnalysis.py`) were written in [Python](https://www.python.org/download/releases/2.7/) and depend on several python packages:
  * [Pandas](https://pypi.org/project/pandas/)
  * [Numpy](https://pypi.org/project/numpy/)
  * [Matplotlib](https://pypi.org/project/matplotlib/)
  * [Seaborn](https://pypi.org/project/seaborn/)
  * [Gillespy2](https://pypi.org/project/gillespy2/)

## <a name="Step-1"></a> **Step 1** `SimSetup.py`

The user specifies the treatment(s) by providing modeling parameter values and initial density for each dynamic variable (*D<sub>t</sub>*, *R<sub>t</sub>*, *T<sub>t</sub>*, and *F<sub>t</sub>*). The `SimSetup.py` program runs a deterministic simulation on the treatment(s) using the plasmid population model. This is a quick computation capturing the average population dynamics of each treatment used to determine the proper incubation times for calculating the various metrics of interest (TDR, SIM, ASM, and LDM). The incubation time is a required input for running the stochastic simulation. Since each conjugation metric has different requirements, two different incubation options are created.
 
* I1: The total population density increases until a user-specified population density is reached. The incubation time needed to reach the specified density is calculated using a numerical simulation. 
* I2: The total population increases in density until a user-specified transconjugant density is reached. The incubation time needed to reach the specified transconjugant population density is calculated using a numerical simulation. This incubation time is typically shorter than I1 and therefore is used to run the additional populations that are necessary for calculating the LDM metric.

#### Inputs
The `SimSetup.py` program uses the custom functions provided by `SimFunction.py` to create shell scripts for running simulations based on the following custom inputs:
* -s/--SimRun_path: The **full** computer path to the SimRun.py program 
* -o/--output_folder: The folder for the output files
* -c/--csv_input: A user-created input CSV containing plasmid population model parameters and the cutoff criteria to determine the incubation time ([template input CSV](https://github.com/livkosterlitz/LDMsimulations/blob/main/BasicRun/SimSetup_inputCSV_example.csv), more details in table below)

| Column name | Example | Description | Units |
| :--- | :--- | :--- | :--- |
| Treatment_ID | T0 | user-specified treatment identifier | |
| psiD | 0.5 | donor growth rate | h<sup>-1 |
| psiR | 0.5 | recipient growth rate | h<sup>-1 |
| psiT | 0.5 | transconjugant growth rate | h<sup>-1 |
| psiF | 0.5 | plasmid-free donor growth rate | h<sup>-1 |
| gammaDR | 1 x 10<sup>-9 | donor-to-recipient conjugation rate | ml cfu<sup>-1</sup> h<sup>-1 |
| gammaTR | 1 x 10<sup>-9 | transconjugant-to-recipient conjugation rate | ml cfu<sup>-1</sup> h<sup>-1 |
| gammaDF | 1 x 10<sup>-9 | donor-to-plasmid-free-donor conjugation rate | ml cfu<sup>-1</sup> h<sup>-1 |
| gammaTF | 1 x 10<sup>-9 | transconjugant-to-plasmid-free-donor conjugation rate | ml cfu<sup>-1</sup> h<sup>-1 |
| tauD | 0 | donor segregation rate | h<sup>-1 |
| tauT | 0 | transconjugant segregation rate | h<sup>-1 |
| D_0 | 1 x 10<sup>2 | initial donor density | cfu ml<sup>-1 |
| R_0 | 1 x 10<sup>2 | initial recipient density | cfu ml<sup>-1 |
| T_0 | 0 | initial transconjugant density | cfu ml<sup>-1 |
| F_0 | 0 | initial plasmid-free donor density | cfu ml<sup>-1 |
| Cutoff_I1 | 1 x 10<sup>7 | average population density to terminate simulation | cfu ml<sup>-1 |
| Cutoff_I2 | 1 x 10<sup>2 | average transconjugant density to terminate simulation | cfu ml<sup>-1 |
| Num_sims_I1 | 1 | number of simulations with I1 incubation condition |  |
| Num_sims_I2 | 99 | number of simulations with I2 incubation condition |  |
| Num_of_experiments| 100 | An integer for the number of experiments. Note: The number of populations simulated within each experiment is specified in columns Num_sims_I1 and Num_sims_I2. |  |

#### Terminal command
```bash
python Programs/SimSetup.py -s /Users/Fullpath/Programs/SimRun.py -c Treatments/SimSetup_inputCSV_example.csv -o Anyuserfolder/
```
#### Outputs
The program has the following outputs which are placed in separate folders:
1. Treatments: [example CSV output](https://github.com/livkosterlitz/LDMsimulations/blob/main/BasicRun/Treatments/SimSetup_inputCSV_example_with_incubation.csv)
> For each treatment (i.e., row in the input CSV), the program determines incubation times based on cutoff criteria provided by the user: (I1) maximum population size or (I2) maximum population size of transconjugants. The output CSV file is identical to the input CSV except for the two columns which return the incubation times (hours) that will be used for running stochastic simulations with the `SimRun.py` program. 
2. ODE_figures: [example ODE plot](https://github.com/livkosterlitz/LDMsimulations/blob/main/BasicRun/ODE_figures/ODE_T1.pdf)
> For each treatment, a plot is generated for the deterministic simulation. 
3. ODE_data
> For each treatment, a csv file is generated alongside the ODE figure.  
4. SimData: [example Shell script](https://github.com/livkosterlitz/LDMsimulations/blob/main/BasicRun/SimData/T1/T1.sh)
> For each treatment, one folder is generated containing one shell script (e.g., SimData/T1/T1.sh). Each row of the shell script is a stochastic simulation with the `SimRun.py` program. The number of rows will depend on the number of simulations requested by the user. 

## <a name="Step-2"></a> **Step 2**  Shell script for `SimRun.py`

Within each treatment's shell script, one row corresponds to one stochastic simulation with the `SimRun.py` program. 

### Execute shell script in terminal to simulate mating cultures
```bash
cd BasicRun/SimData/T1/

chmod a+x T1.sh

./T1.sh
```
#### `SimRun.py` Input
* Various user-specified inputs passed initially through `SimSetup.py`
* -t/--t: incubation time calculated in `SimSetup.py`
* -o/--output_filename: A unique filename to store the stochastic simulation results. There are four fields in the output filename all separated by an underscore (e.g., T1_E0_I1_P1.csv)

| Field | Example | Description |
| :--- | :--- | :--- | 
| Treatment | T1 | User specified treatment identifier in the input CSV connected to the input parameters (*&psi;<sub>D<sub>*, *&psi;<sub>R<sub>*, *&psi;<sub>T<sub>*, *&psi;<sub>F<sub>*, *&tau;<sub>T<sub>*, *&tau;<sub>D<sub>*, *&gamma;<sub>DR<sub>*, *&gamma;<sub>TR<sub>*, *&gamma;<sub>DF<sub>*, *&gamma;<sub>TF<sub>*) and initial values (*D<sub>t</sub>*, *R<sub>t</sub>*, *T<sub>t</sub>*, *F<sub>t</sub>*)
| Experiment | E0 | A numeric identifier for each experiment| 
| Incubation condition | I1 | Indicates the incubation condition: I1 or I2 | 
| Population replicate | P1 | A numeric identifer for a replicate population. The number of population replicates for each incubation condition is user-specified in the input CSV | 

#### Example command in the shell script to generate one simulated population with the `SimRun.py` program
```bash
python3 /Users/Fullpath/Programs/SimRun.py --output_filename T1_E0_I1_P0 --psiD 0.5 --psiR 0.5 --psiT 0.5 --psiF 0.5 --gammaDR 1e-09 --gammaTR 1e-06 --gammaDF 1e-06 --gammaTF 1e-09 --tauD 0 --tauT 0 --D_0 100.0 --R_0 100.0 --T_0 0 --F_0 0 --t 21.6
```
#### `SimRun.py` output
The program outputs a single CSV file with the specified output filename consisting of various columns:
| Column | Description |
| :--- | :--- | 
| Time| A time counter in seconds | 
| Donor| Donor population size at a given second in the incubation | 
| Plasmid-free-donor| Plasmid-free-donor population size at a given second in the incubation | 
| Recipient| Recipient population size at a given second in the incubation | 
| Transconjugant| Transconjugant population size at a given second in the incubation | 

## <a name="Step-3"></a> **Step 3**  `SimAnalysis.py`

This program takes the stochastically simulated mating assays and calculates the various conjugation rate estimates at various incubation times. 

#### Inputs

The `SimAnalysis.py` program uses the simulated mating assays generated in [Step 2](#Step-2) from the `SimRun.py` program to estimate the conjugation rate using various metrics and provide various diagnostic plots describing the simulations based on the following custom inputs:
* -t/--treatment: The name of the treatment (e.g., T0)
* -c/--parameter_csv: The **full** computer path to the treatments CSV output generated in [Step 1](#Step-1) ([example CSV output](https://github.com/livkosterlitz/LDMsimulations/blob/main/BasicRun/Treatments/SimSetup_inputCSV_example_with_incubation.csv))
* -f/--csv_folder: The **full** computer path to the data folder containing the output files from [Step 2](#Step-2) (e.g., SimData/T1/)
* -o/--output_folder: The **full** computer path to deposit the output files

#### Terminal command
```bash
python /Users/Fullpath/Programs/SimAnalysis.py -t T1 -c Treatments/SimSetup_inputCSV_example_with_incubation.csv -f SimData/T1 -o Anyuserfolder/
```
#### Outputs
The program has two output folders:
1. Analysis estimates
> Put explanation and examples here....
2. Analysis figures
> Put explanation and examples here......
