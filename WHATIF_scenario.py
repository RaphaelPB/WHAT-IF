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

#Import python packages
from __future__                   import division   #As in python 2.7 division between integers gives an integer - avoid problems
import time
import os
import pickle
import sys
from pyomo.environ                import Suffix, SolverFactory #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Collects results and exports to excel sheets
#Parallel computing
import multiprocess as mp

#%%OPTIONS - MODIFY BY USER
RESULTFOLDER = 'scenarios'+'_'+time.strftime("%d_%m_%Y_%Hh%M") #automatically generates names based on time and date (avoids erasing previous results)
UPDATE = 0 #0 updates all parameters, 1 updates only selected parameters
PARALLEL_scenario = 0 #Run scenarios in parallel 
npll = 3 #Number of paralllel scenario solves
EXPORT = 'all' #'all' powerBI files + following, 'xlsx': individual excel files + following, 'txt': selected results + excel of all selected results
SOLVER = 'ipopt' #'cbc' #'cplex'
SOLVERPATH = 0 #0 is default, precise path only if required by solver

#%%
#Define paths
ResultFolderPath= os.path.join(dirname, 'Results', RESULTFOLDER)
DataFolderPath  = os.path.join(dirname, 'Data')

#Collect parameters
Main            = DataFolderPath + os.sep + 'MainFile_ex.xlsx'
Water           = DataFolderPath + os.sep + 'WaterModule_ex.xlsx'
Agriculture     = DataFolderPath + os.sep + 'AgricultureModule_ex.xlsx'
CropMarket      = DataFolderPath + os.sep + 'CropMarketModule_ex.xlsx'
Energy          = DataFolderPath + os.sep + 'EnergyModule_ex.xlsx'
Investment      = DataFolderPath + os.sep + 'InvestmentModule_ex.xlsx'
Param           = DataFolderPath + os.sep + 'Parameterspy37.txt' #parameters (python dictionnaries) saved as txt

print('Harvesting parameters...')
parameters      = Database(update=UPDATE,DataFile=Param)
parameters.harvest_all([Main,Water,Agriculture,CropMarket,Energy])
investopt=[parameters.val['Options'][k] for k in parameters.val['Options'].keys() if k[0]=='Investment module'] #investment module options
if 1 in investopt or 'continuous' in investopt:
    parameters.harvest_all([Investment])
#Save parameters
parameters.save(Param)
print('*Parameters harvested*')

#%%
#Define single scenario analysis function    
def ScenarioAnalysis(ss,parameters,solver):
    tt=time.time()
    #Collect extra hydrology parameters (hydrology scenario has to be called 'sClimate')
    parameters.load_csv(DataFolderPath,ss)
    #if 'sClimate' in parameters.val.keys() and parameters.val['sClimate'][ss] not in parameters.val['wRunOff'].keys():
    #    parameters.load_hydrology(DataFolderPath,scenarios=[parameters.val['sClimate'][ss]])
    #Create model  
    print('Formulating scenario: ' + ss)
    HOM             = HydroeconomicOptimization(parameters,scenario=ss)
    HOM.model.dual  = Suffix(direction=Suffix.IMPORT)    
    
    #Solve scenario
    print('Solving scenario: ' + ss)
    solverstatus = solver.solve(HOM.model)
    if HOM.model.Options['Investment module'] == 1:
    #Fix selected investment to switch to fully linear model and get dual values
        for ip in HOM.model.ninvphase:
            for inv in HOM.model.ninvest:   
                HOM.model.IbINVEST[ip,inv].fixed = True
        print('...re-solving scenario ' + ss)
        solverstatus = solver.solve(HOM.model)        
    
    #Collect and export results
    print('...saving and exporting results of ' + ss)
    exppath = ResultFolderPath + os.sep + ss + '.xlsx' 
    result_analysis = ResultAnalysis()
    if EXPORT in ['xlsx','all']:
        results=result_analysis.readresults(HOM.model,parameters,solverstatus,PRINT=0,scenario=ss) #   
        result_analysis.export_to_excel(results,exppath,NEWSHEET=1) #export formated results to excel 
    if EXPORT in ['all']:
        result_analysis.export_all_DV(HOM.model,ResultFolderPath,scenario=ss) #export all decision variables
        result_analysis.export_mass_balances(results,ResultFolderPath,scenario=ss) #export mass balances (for powerBI)
        result_analysis.export_index_mapping(HOM.model,os.path.join(ResultFolderPath,'Mapping')) #export index mapping/connections (for powerBI)
    
    #Save selected results for scenario analysis
    output=result_analysis.selectedresults(HOM.model,parameters,scenario=ss)
    outputpath=ResultFolderPath + os.sep + ss + '.txt' 
    print(ss,'saving to text file')
    pickle.dump(output,open(outputpath,"wb"))
    print('scenario '+ ss +' solved in '+str(time.time()-tt)+' seconds')
    return output

#%% Start parallel solving of problems
#Create folders
if not os.path.exists(ResultFolderPath):
    os.makedirs(ResultFolderPath)
    
#Choose solver
solver = SolverFactory(SOLVER,executable=SOLVERPATH) if SOLVERPATH != 0 else SolverFactory(SOLVER)
#if solver.name == 'ipopt':
#    solver.options['linear_solver']='ma97'
    #solver.options['mu_strategy']='adaptive'
    #solver.options['bound_relax_factor']=10**-12
if solver.name=='cplex' and PARALLEL_scenario==1:
    solver.options['threads'] = 1 #limit the number of parallel computing

if __name__ == '__main__':
    scenarios = parameters.val['nscenario']
    if PARALLEL_scenario == 1:
        pool = mp.Pool(min(npll,len(scenarios)))
        parallelresults = pool.starmap(ScenarioAnalysis, [(ss,parameters,solver) for ss in scenarios])
        pool.close()
    else:
        parallelresults=[ScenarioAnalysis(ss,parameters,solver) for ss in scenarios]
#Save results (can be done later)
    result_analysis=ResultAnalysis()
    result_analysis.export_scenario_analysis(scenarios,parameters.val['sRefScen'],parallelresults,ResultFolderPath)
