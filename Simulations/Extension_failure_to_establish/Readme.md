These programs alter the main simulation framework to add in two additional sources of stochasticity: plating and extinction probabilities. Below highlights within the workflow what is modified and the changes needed by the user. The output files from [Step 3](#Step-3) are in the [Step3_output folder](https://github.com/livkosterlitz/LDM/tree/main/Simulations/Extension_failure_to_establish/Step3_output) and a CSV file (i.e., [Treatments_establishment.csv](https://github.com/livkosterlitz/LDM/blob/main/Simulations/Extension_failure_to_establish/Treatments_establishment.csv)) describing each treatments parameter settings and the relevant figure in the manuscript. 

# Conjugation rate estimates
These programs estimate conjugation rate with 2 metrics: 
| Method | Acronym | Year  | Estimate |
| :--- | :--- | :---  | :--- |
| Simonsen _et. al._ method <br>(aka. end-point method) | [SIM](https://doi.org/10.1099/00221287-136-11-2319) | 1990 | <img width="300" src="https://github.com/livkosterlitz/LDM/blob/main/Simulations/Images/SIM/SIM.png"> |
| Luria-Delbrück method <br>(general formula) | [LDM](https://doi.org/10.1101/2021.01.06.425583) | 2021 | <img  width="375" src="https://github.com/livkosterlitz/LDM/blob/main/Simulations/Images/LDM_general/LDM_general.png"> |

# Program Workflow

|Steps| Script |Step description|
| :--- | :--- | :--- |
| [Step 1](#Step-1) | `SimSetupEstablishment.py` | Adds a transconjugant monoculture (i.e., Num_sims_I3) to each experiment used to estimate the transconjugant growth rate. This is an extra variable that is needed for the general formula of the LDM estimate. |
| Step 2 | `SimRun.py` | No change. |
| [Step 3](#Step-3) | `SimAnalysisEstablishment.py` | New functions incorporated to add in the two sources of stochasticity.  |

## <a name="Step-1"></a> **Step 1** `SimSetupEstablishment.py`

Since the general formula of the LDM estimate that incorporates non-zero extinction probabilities has different requirements, there is a new incubation option added.
 
* I3: The transconjugant density increases until a user-specified population density is reached. The same incubation time is used for I3 that is selected for I1. 

#### Inputs
* -c/--csv_input: Two columns are added to the user-created input CSV.

| Column name | Example | Description | 
| :--- | :--- | :--- |
| Num_sims_I3 | 1 | number of simulations with I1 incubation condition. By default, this should match the number of I1 simulations. |
| Extinction_probability| 0.25 | Extinction probability for the donor, recipient and transconjugant cell types.|

## <a name="Step-3"></a> **Step 3**  `SimAnalysisEstablishment.py`

This script adds in the two sources of stochasticity before calculating the conjugation rate.
* Each estimate of population density goes through dilution plating. This uses the scipy package in python and the poisson distribution by setting the mean to the 'actual' population density generate in Step 2 by the Gillespie algorithm. 
* Each estimate of the population density goes through sampling using the extinction probability set by the user. This also uses the scipy package in python and the binomial distribution where the probability of a success is set to 1 - extinction probability. 
* The density estimates used to calculate the conjugation rates are are first 'corrected' by the extinction probability. Given that this correction isn't perfect given the random number generation of the last step, this effectively adds noise to the calculations.


