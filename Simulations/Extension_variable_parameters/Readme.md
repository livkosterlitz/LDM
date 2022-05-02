These programs alter the main simulation framework to add in one additional source of stochasticity: variation in the user-specified parameter values for growth and transfer rates. Below highlights within the workflow what is modified and the changes needed by the user. The output files from [Step 3](#Step-3) are in the [Step3_output folder](https://github.com/livkosterlitz/LDM/tree/main/Simulations/Extension_variable_parameters/Step3_output). 

# Program Workflow

|Steps| Script |Step description|
| :--- | :--- | :--- |
| [Step 1](#Step-1) | `SimSetupStochasticParams.py` | Using a user specified value for the standard deviation, the growth rate and transfer rates are drawn from a normal distribution on a natural logarithm and base 10 logarithm scale, respectively. |
| Step 2 | `SimRun.py` | No change. |
| Step 3 | `SimAnalysis.py` | No change.  |

## <a name="Step-1"></a> **Step 1** `SimSetupStochasticParams.py`

To be able to add variation to the user-specified parameters, the user needs to specify the standard deviation for the distribution which the parameters will be sampled from. 

#### Inputs
* -c/--csv_input: One column is added to the user-created input CSV.

| Column name | Example | Description | 
| :--- | :--- | :--- |
| std_dev | 0.2 | the standard deviation for the normal distribution used to sample the growth rate and transfer rate parameters. |

