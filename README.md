# WHAT-IF: Supporting water infrastructure investment planning with hydroeconomic optimization models 

The WHAT-IF (Water, Hydropower, Agriculture Tool for Investment and Financing) model is an open source decision support tool 
distributed under the GPLv3 license and is described in the [HESS publication](https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-167/):
> WHAT-IF: an open-source decision support tool for water infrastructure investment planning within the Water-Energy-Food-Climate Nexus


The decision-support tool can be used in single-run mode or in scenario mode:

* **WHATIF_main.py** 
Single run of the model, collects the parameters from the data module, creates an hydroeconomic optimization model, solves it and exports the results 
* **WHATIF_scenario.py** 
Multiple runs of the model, collects the parameters from the data module, creates an hydroeconomic optimization model, 
solves it for all specified scenarios and exports the results for each scenario as well as a synthesis of all scenarios
* **WHATIF_scenario_mpc.py** 
Same with the possibility of using the Model Predictive Control (MPC) framework which avoids the perfect foresight assumption

The **bin** folder contains the main libraries for the decision support tool:

* **data_collection.py**
Collects the data from the excel files (MainFile, WaterModule, AgricultureModule, CropMarketModule, EnergyModule) into a python object
* **hydroeconomic_optimization.py** 
Formulates the optimization problem of the water-energy-food management (and investments at a later stage) 
* **result_analysis.py**
Performs the analysis of the results and exports them towards the excel/python files 
* **model_predictive_control.py**
Routines for the Model Predictive Control (MPC) framework

The **Data** folder contains the requiered data for the different modules:

* **MainFile.xlsx** Is the configuration file, specifying the time steps, options and scenarios that are used by the model.
* **WaterModule.xlsx** Data supporting the Water module, including the hydrology (runoff, evapo-transpiration, precipitation, groundwater recharge, catchments …), the reservoirs and the environmental requirements
* **AgricultureModule.xlsx** Data supporting the Agriculture module, including farming zones, farm types, crops and cultures characteristics.
* **CropMarketModule.xlsx** Data supporting the crop market module, including crop markets, demands, value of crops, transport routes and food security constraints.
* **EnergyModule.xlsx** Data supporting the energy production and energy markets modules, including hydropower plants, other power plants, power technologies, fuels, energy markets, demands, value of energy and transmission lines.

The **Documents** folder contains additional information:
* **How_to_Install_WHATIF.docx** is a step by step guide on how to install and run WHAT-IF
* **How_to_Compare_Scenarios** is a step by step guide on how to create and compare scenarios in WHAT-IF
* **How_to_Use_WHATIF_visualize.docx** is a step by step guide on how to use an enhance result viewer
* **How_to_Spot_Common_Errors.docx** lists most common programing errors when using/developing the model
* **Data_Organization.docx** summarizes the organization of the data for the WHAT-IF model and the HESS publication
* **WHATIF_py37.yml** is the conda environment with all the useful packages to run WHAT-IF 

## Usage

The **Documents** folder contains a step by step guide **How_to_Install_WHATIF** describing the process, in brief:
* Install the [Anaconda navigator](https://anaconda.org/anaconda/anaconda-navigator), 
the simpliest way to manage your pyhton packages and versions

* Install the WHATIF_py37 environment from the **Documents** folder

To do so, in the anaconda prompt, run:

`conda env create -f WHATIF_py37.yml` 

* Recommended solvers are cplex (free for academics) or glpk (open-source)
* WHAT-IF is now only released with Python 3.7


## Contributors 
Payet-Burin Raphaël (DTU/COWI - rapp@env.dtu.dk)
Mikkel Kromann (dansk energi)
Peter Bauer-Gottwein (DTU)
Silvio Pereira-Cardenal (COWI)
Kenneth Strzepek (MIT)
