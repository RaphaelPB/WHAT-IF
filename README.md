# WHAT-IF: Supporting water infrastructure investment planning with hydroeconomic optimization models 
**Master branch**

## The decision support tool

The WHAT-IF (Water, Hydropower, Agriculture Tool for Investment and Financing) model is an open source decision support tool 
distributed under the GPLv3 license and is described in the [HESS publication](https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-167/):
> WHAT-IF: an open-source decision support tool for water infrastructure investment planning within the Water-Energy-Food-Climate Nexus  

WHAT-IF has a holistic and bottom-up approach, the management of local infrastructure (e.g. reservoirs, irrigation, water supply, energy production) 
is solved simultaneously in order to maximize welfare benefits considering crop and power markets and trade. 
The hydro-economic optimization framework enables the tool to solve synergies and trade-offs between the water, 
energy and agricultural sector and explore a large range of scenarios considering exogenous climate change and socio-economic drivers. 
The main feedback loops in the model are summarized in the figure below.  

[WHAT-IF model](https://github.com/RaphaelPB/WHAT-IF/blob/master/Documents/images/WHATIF_model.PNG)  

The model operates at a user-defined (monthly is default) timestep and at a catchment scale. 
The hydrological module is represented through timeseries of hydrological cycle water balance variables (rainfall, runoff, evapotranspiration, groundwater recharge). 
Reservoirs can store and release water, while domestic and industrial users withdraw and return water considering supply costs and treatment technologies. 
Individual crops are represented and their water requirements are internally computed in the model using the FAO-56 method, 
while actual yields are calculated using potential yields and the FAO-33 method; livestock water demands can be user-specified. 
Crop production is calculated at the catchment level. Markets and trade are represented at national or subnational level, 
depending on the scope and data availability. 
Hydropower plants are individually represented, while other energy technologies are represented through aggregated power units 
(at the national or subnational level). A capacity expansion model considers investments in additional power capacity.
Power is traded on markets considering the transmission capacity within and between countries.

## Wiki

Check the [wiki page](https://github.com/RaphaelPB/WHAT-IF/wiki) for all detailed information:
* [Code structure](https://github.com/RaphaelPB/WHAT-IF/wiki/Code-structure)
* [Installing and running WHAT IF](https://github.com/RaphaelPB/WHAT-IF/wiki/Installing-and-running-WHAT-IF)
* [Creating, running and comparing scenarios](https://github.com/RaphaelPB/WHAT-IF/wiki/Creating,-running-and-comparing-scenarios) 
* and many more...


## Install

The [Installing and running WHAT IF](https://github.com/RaphaelPB/WHAT-IF/wiki/Installing-and-running-WHAT-IF) contains a step by step guide describing the process, in brief:
* Install the [Anaconda navigator](https://anaconda.org/anaconda/anaconda-navigator), 
the simpliest way to manage your pyhton packages and versions
* In the anaconda prompt, run:  
`conda config --add channels conda-forge`  
`conda create -n WHATIF_py37 python=3.7.3`  
`conda install -n WHATIF_py37 openpyxl xlsxwriter xlrd pyomo=5.6.2 pandas numpy multiprocess ipopt glpk`  
`conda activate WHATIF_py37`  

* Recommended solvers are ipopt (open-source, non-linear), cplex (free for academics, linear) 
* WHAT-IF is now only released with Python 3.7


## Contributors 
> Innovation Fund Denmark (grant no. 7038-00015B), COWIfonden, the Technical University of Denmark (DTU), the Massachusetts Institute of Technology (MIT), and COWI A/S funded the industrial PhD project in which this research was carried out  

Payet-Burin Raphaël (DTU/COWI - rapp@env.dtu.dk)  
Mikkel Kromann (dansk energi)  
Peter Bauer-Gottwein (DTU)  
Silvio Pereira-Cardenal (COWI)  
Kenneth Strzepek (MIT)  
