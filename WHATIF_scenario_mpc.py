# -*- coding: utf-8 -*-
"""
    Copyright 2018 Technical University of Denmark
    Copyright 2018 COWI A/S
    
    Authors: 
        Raphaël Payet-Burin
        Mikkel Kromann

    This file is part of WHAT-IF.

    WHAT-IF is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 (GPLv3) as published by the Free Software Foundation,

    WHAT-IF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with WHAT-IF. If not, see http://www.gnu.org/licenses/.

"""
#Import python libraries: pandas, xlswriter, openpyxl, xlrd, pyomo, glpk 
import time
import os
#import pickle
import sys
from pyomo.environ                import SolverFactory      #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from result_analysis              import ResultAnalysis             #Exports results to excel sheets
from model_predictive_control     import SolveScenario
#For Model Predictive Control framework
#from joblib import Parallel, delayed
import multiprocessing as mp

#%% DEFINE PARAMETERS
#OPTIONS
UPDATE = 0 #0 updates all parameters, 1 updates only selected parameters
epll=24 #Maximum number of parallel runs for ensemble forecast 
PARALLEL_scenario = 0 #Run scenarios in parallel (BOTH IS NOT POSSIBLE)
RESULTFOLDER = 'mpc'+'_'+time.strftime("%a%d_%m_%Y_%Hh%M")
EXPORT = 'all' #'all' powerBI files + following, 'xlsx': individual excel files + following, 'txt': selected results + excel of all selected results
SOLVER = 'ipopt' #'cplex'
SOLVERPATH = 0#'~/CoinIpopt/bin/ipopt' # 0 if not needed

#%% Start parallel solving of problems
if __name__ == '__main__':
    #Paths
    paths  = {'data':os.path.join(dirname, 'Data'),
              'output':os.path.join(dirname, 'Results',RESULTFOLDER),
              'export':EXPORT}    
    #Define path to data
    Main            = paths['data'] + os.sep + 'MainFile_ex.xlsx'
    Water           = paths['data'] + os.sep + 'WaterModule_ex.xlsx'
    Agriculture     = paths['data'] + os.sep + 'AgricultureModule_ex.xlsx'
    CropMarket      = paths['data'] + os.sep + 'CropMarketModule_ex.xlsx'
    Energy          = paths['data'] + os.sep + 'EnergyModule_ex.xlsx'
    Parameter       = paths['data'] + os.sep + 'Parameterspy37.txt' #parameters (python dictionnaries) saved as txt    
    #Collect parameters
    t=time.time()
    print('Harvesting parameters...')
    parameters      = Database(update=UPDATE,DataFile=Parameter)
    parameters.harvest_all([Main,Water,Agriculture,CropMarket,Energy]) #
    parameters.save(Parameter) #Save parameters    
    print('*Parameters harvested in '+str(time.time()-t)+' seconds')
    #Create folder
    if not os.path.exists(paths['output']):
        os.makedirs(paths['output'])
    
    #Choose solver
    solver = SolverFactory(SOLVER,executable=SOLVERPATH) if SOLVERPATH != 0 else SolverFactory(SOLVER)
    if solver.name == 'ipopt':
        #solver.options['linear_solver']='ma97' #PARDISO
        #solver.options['mu_strategy']='adaptive' #this option is found to speed up problem resolution but is unstable for MPC
        solver.options['bound_relax_factor']=10**-12
        
    #Run through scenarios
    scenarios = parameters.val['nscenario']
    if PARALLEL_scenario == 1:
        pool = mp.Pool(min(epll,len(scenarios)))    
        parallelresults = pool.starmap(SolveScenario, [(ss,parameters,solver,paths) for ss in scenarios])
        pool.close()
    else:
        parallelresults=[SolveScenario(ss,parameters,solver,paths) for ss in scenarios]
    result_analysis=ResultAnalysis()
    result_analysis.export_scenario_analysis(scenarios,parameters.val['sRefScen'],parallelresults,paths['output'])