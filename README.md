# WHAT-IF: Supporting water infrastructure investment planning with hydroeconomic optimization models 

**This is the example branch - perfect start point to create your own WHAT-IF case study !**

The WHAT-IF (Water, Hydropower, Agriculture Tool for Investment and Financing) model is an open source decision support tool 
distributed under the GPLv3 license and is described in the [HESS publication](https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-167/):
> WHAT-IF: an open-source decision support tool for water infrastructure investment planning within the Water-Energy-Food-Climate Nexus  

Check the [wiki page](https://github.com/RaphaelPB/WHAT-IF/wiki) for all detailed information


## Usage

The **Documents** folder contains a step by step guide **How_to_Install_WHATIF** describing the process, in brief:
* Install the [Anaconda navigator](https://anaconda.org/anaconda/anaconda-navigator), 
the simpliest way to manage your pyhton packages and versions

* Install the WHATIF_py37 environment from the **Documents** folder

To do so, in the anaconda prompt, run:

* `conda config --add channels conda-forge`
* `conda create -n WHATIF_py37 python=3.7.3`
* `conda install -n WHATIF_py37 openpyxl xlsxwriter xlrd pyomo=5.6.2 pandas numpy multiprocess ipopt glpk`
* `conda activate WHATIF_py37` 

* Recommended solvers are ipopt (non-linear), cplex (free for academics) or glpk (open-source)
* WHAT-IF is now only released with Python 3.7


## Contributors 
Payet-Burin Raphaël (DTU/COWI - rapp@env.dtu.dk)
Mikkel Kromann (dansk energi)
Peter Bauer-Gottwein (DTU)
Silvio Pereira-Cardenal (COWI)
Kenneth Strzepek (MIT)
