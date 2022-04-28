These programs alter the main simulation framework to add in two additional sources of stochasticity: plating and extinction probabilities. Below highlights within the workflow what is modified and the changes needed by the user. 

# Conjugation rate estimates
Currently, the programs estimate conjugation rate with 4 metrics: 
| Method | Acronym | Year  | Estimate |
| :--- | :--- | :---  | :--- |
| Simonsen _et. al._ method <br>(aka. end-point method) | [SIM](https://doi.org/10.1099/00221287-136-11-2319) | 1990 | <img width="300" src="https://github.com/livkosterlitz/LDM/blob/main/Simulations/Images/SIM/SIM.png"> |
| Luria-Delbrück method <br>(general formula) | [LDM](https://doi.org/10.1101/2021.01.06.425583) | 2021 | <img  width="375" src="https://github.com/livkosterlitz/LDM/blob/main/Simulations/Images/LDM_general/LDM_general.png"> |



|Steps| Script |Step description|
| :--- | :--- | :--- |
| [Step 1](#Step-1) | `SimSetupEstablishment.py` | Altered SimSetup script that adds in a third type of incubation condition (i.e., Num_sims_I3), which incubates a mono-culture of transconjugants used to estimate the transconjugant growth rate. This is an extra variable that is needed for the more general form of the LDM estimate, which is the necessary form when there is a positive probability that a transconjugant cell failed to establish. |
| [Step 2](#Step-2) | `SimRun.py` | No change. |
| [Step 3](#Step-3) | `SimAnalysisEstablishment.py` | New functions are added to calculate the LDM and SIM.  |
