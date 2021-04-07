# WHAT-IF: Supporting water infrastructure investment planning within the water-energy-food nexus 
**Master branch** Comes with example case data - perfect to understand the model and start your own WHAT-IF study case

## The decision support tool

WHAT-IF (Water, Hydropower, Agriculture Tool for Investment and Financing) is an open source decision support tool 
distributed under the GPLv3 license and is described in the [HESS publication](https://doi.org/10.5194/hess-23-4129-2019):
> Payet-Burin, R., Kromann, M., Pereira-Cardenal, S., Strzepek, K. and Bauer-Gottwein, P.: WHAT-IF: An open-source decision support tool for water infrastructure investment planning within the water-energy-food-climate nexus, Hydrol. Earth Syst. Sci., 23(10), 4129–4152, [doi:10.5194/hess-23-4129-2019](https://doi.org/10.5194/hess-23-4129-2019), 2019  

WHAT-IF was developed to analyze water-related investments in an economic context. WHAT-IF has a holistic and bottom-up approach, the management of infrastructure (e.g. reservoirs, irrigation, water supply, energy production) is solved simultaneously in order to maximize welfare benefits considering crop and power markets and trade. 
The hydro-economic optimization framework enables the tool to solve synergies and trade-offs between the water, 
energy and agricultural sector and explore a large range of scenarios considering exogenous climate change and socio-economic drivers. 
The main feedback loops in the model are summarized in the figure below.  

![WHAT-IF model](https://github.com/RaphaelPB/WHAT-IF/blob/master/Documents/images/WHATIF_model.PNG)  

The water resource availability is typically represented by an exogenous rainfall-runoff model, while the natural (e.g. river network, lakes) and engineered (e.g. reservoirs, transfer schemes) flow network is represented internally at the sub-basin and monthly scale. Water users (except hydropower and irrigation) are represented through their demand/value, while ecosystems are represented through environmental flow constraints.

The agriculture sector is represented by rainfed and irrigated agriculture, that produce within crop markets (typically at the national scale), while trade occurs between markets (including e.g. a world market). Crop demand is represented per market considering own-price elasticity. The main links with the water resources are rainfall, surface and groundwater supply to agriculture, using mainly the FAO 33 (Doorenbos and Kassam, 1979) and FAO 56 (Allen et al., 1998) for irrigation requirement and yield water response functions. 

Power plants (e.g. hydropower, thermal, renewables), produce within power markets (typically at the national scale), that trade through transmission lines, while greenhouse gas emissions are traced. A capacity expansion model represents the development of generic power technologies. Power demand is represented as inelastic, but different load segment (sometimes called "time slices") that sub-divide the monthly demand can be defined (e.g. peak and base demand). Capacity factors for power plants can be defined at the monthly and load segment scale and represent limited availability of power plants. For example, capacity factors can be used to represent seasonal and intraday variability of renewable energies. Hydropower production is dependent on river flow and reservoir releases.

While the previously described processes are predefined in the model (with a flexible implementation), a general activity module represents any other process that consumes and/or produces one or several commodities (land, water, power, and crops) connecting sub-basins, crop markets, power markets and agriculture land. This can represent various processes such as: desalinization (consumes energy and produces water), livestock (consumes land, water and crops and produces another food commodity), food processing (consumes energy, water and crops and produces another commodity), and bioenergy (consumes crops or crop residues and produces energy).
The main exogenous drivers are demand for commodities (water, crops, electricity), technology development (e.g. power technologies, yields, efficiencies), external markets (import/export prices), policies (e.g. environmental constraints, food security policies, carbon taxes), and climate change.

Operation of infrastructure (e.g. releases, water supply, energy production, trade, cropping) is solved simultaneously considering physical and policy constraints to maximize welfare benefits. Welfare benefits are the sum of consumer and producer surplus (see Krugman and Wells, 2005). The optimization framework simulates how operation and management responds to new infrastructure and exogenous climate change and socio-economic drivers. As the model is solved in a single iteration for the entire planning horizon, this results in assuming "perfect foresight" in represented infrastructure operation.

## Wiki

For all detailed information check the [wiki page](https://github.com/RaphaelPB/WHAT-IF/wiki):
* [Installing and running WHAT IF](https://github.com/RaphaelPB/WHAT-IF/wiki/Installing-and-running-WHAT-IF)
* [Code structure](https://github.com/RaphaelPB/WHAT-IF/wiki/Code-structure)
* [Creating, running and comparing scenarios](https://github.com/RaphaelPB/WHAT-IF/wiki/Creating,-running-and-comparing-scenarios) 
* and many more...

## Model versions (/branches)
If you look for a specific version of the model (eventually corresponding to a specific paper/study case), change the branch
* **[Master](https://github.com/RaphaelPB/WHAT-IF)** is the main branch containing the most recent public model developments and a synthetic data case to understand the model
* **[HESS_zambezi](https://github.com/RaphaelPB/WHAT-IF/tree/HESS_Zambezi)** is the model as published in [HESS](https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-167/) with the Zambezi dataset
* **[FrontiersInWater_Zambezi](https://github.com/RaphaelPB/WHAT-IF/tree/FrontiersInWater_Zambezi)** is the model as published in (under review) for the investment planning study in the Zambezi River Basin as in the article 
* **[FrontiersInWater_SyntheticCase](https://github.com/RaphaelPB/WHAT-IF/tree/FrontiersInWater_SyntheticCase)** is the model as published in (under review) for the investment planning study in the synthetic case (which is based on the Master branch example) 
* **[Nile_SyntheticCase](https://github.com/RaphaelPB/WHAT-IF/tree/Nile_SyntheticCase)** model set-up described in [this publication](https://www.essoar.org/doi/10.1002/essoar.10504115.1) synthetic case on the nile comparing SDP, MPC, and perfect foresight (we do not recommend to use except to reproduce the experiments from the publication)

The data and the model version are relatively disconnected - but using data from an older set-up might require a few adaptations.
We recommend to always use the latest version of the model.

## Install

The [Installing and running WHAT IF](https://github.com/RaphaelPB/WHAT-IF/wiki/Installing-and-running-WHAT-IF) contains a step by step guide describing the process, in brief:
* Install the [Anaconda navigator](https://anaconda.org/anaconda/anaconda-navigator), 
the simpliest way to manage your pyhton packages and versions
* In the anaconda prompt, run:  
`conda config --add channels conda-forge`  
`conda create -n WHATIF_py38 python=3.8.8`  
`conda install -n WHATIF_py38 openpyxl xlsxwriter xlrd pyomo=5.7.3 pandas=1.2.3 numpy multiprocess ipopt=3.11.1 glpk`  
`conda activate WHATIF_py38`  

* Recommended solvers are glpk (open-source, linear - slow), ipopt (open-source, non-linear), cplex (free for academics, linear - fast) see [Installing extra solvers](https://github.com/RaphaelPB/WHAT-IF/wiki/Installing-extra-solvers)
* WHAT-IF is now only released with Python 3.x (above 3.6, preferably 3.8)

## Contributors 
> Innovation Fund Denmark (grant no. 7038-00015B), COWIfonden (grant no. C-137.02), the Technical University of Denmark (DTU), the Massachusetts Institute of Technology (MIT), and COWI A/S funded the industrial PhD project in which this research was carried out  

Payet-Burin Raphaël (DTU/COWI - rapy@cowi.com)  
Mikkel Kromann (dansk energi)  
[Peter Bauer-Gottwein](https://orbit.dtu.dk/en/persons/peter-bauer-gottwein) (DTU)  
Silvio Pereira-Cardenal (COWI)  
Kenneth Strzepek (MIT)  
