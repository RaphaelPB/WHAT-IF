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
#Import python libraries: 
import time
import os
import pickle
import sys
from pyomo.environ                import Suffix, SolverFactory      #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory 
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Exports results to excel sheets

#%%OPTIONS - MODIFY BY USER
SCENARIO    = 'DEVsemidry' #Scenario to be run, 'WHATIF_main' is default
RESULTFOLDER= 'WHATIF_main'#+time.strftime("%d_%m_%Y_%Hh%M") #automatically generates names based on time and date (avoids erasing previous results)
NEWSHEET    = 1 #1 creates a new sheet, 0 fills existing sheet
UPDATE      = 0 #0 updates all parameters, 1 updates only selected parameters (through the "Info" sheet in the parameters excel files)
EXPORT      = 'all' #'all' powerBI files + following, 'xlsx': individual excel files + following, 'txt': selected results + excel of all selected results
SOLVER      = 'ipopt' #'cplex'
SOLVERPATH  = 0# 0 is default, precise path only if required by solver #'~/CoinIpopt/bin/ipopt' 

#%% DEFINE PARAMETERS
#Data Folder
DataFolderPath  = os.path.join(dirname, 'Data')
ResultFolderPath= os.path.join(dirname, 'Results', RESULTFOLDER)    
#Excel output file
OutPath         = ResultFolderPath + os.sep + 'modelrun_' + time.strftime("%d_%m_%Y_%Hh%M") + '.xlsx'

#Define path to data
Main            = DataFolderPath + os.sep + 'MainFile.xlsx'
Water           = DataFolderPath + os.sep + 'WaterModule.xlsx'
Agriculture     = DataFolderPath + os.sep + 'AgricultureModule.xlsx'
CropMarket      = DataFolderPath + os.sep + 'CropMarketModule.xlsx'
Energy          = DataFolderPath + os.sep + 'EnergyModule.xlsx'
Investment      = DataFolderPath + os.sep + 'InvestmentModule.xlsx'
Param           = DataFolderPath + os.sep + 'Parameters.txt' #parameters (python dictionnaries) saved as txt

#Collect parameters
t=time.time()
print('Harvesting parameters...')
parameters      = Database(update=UPDATE,DataFile=Param)
parameters.harvest_all([Main,Water,Agriculture,CropMarket,Energy])
#Collect investment module parameters
if parameters.val['Options']['Investment module',parameters.val['sOptions'][SCENARIO]] == 1:
    parameters.harvest_all([Investment]) 
#Save parameters    
parameters.save(Param)
#Collect parameters in .csv files
parameters.load_csv(DataFolderPath,SCENARIO)

print(time.time()-t)
print('*Parameters harvested*')

#%% BUILD HYDROECONOMIC OPTIMIZATION MODEL
t=time.time()
print('Creating Hydroeconomic optimization model...')
HOM             = HydroeconomicOptimization(parameters,scenario=SCENARIO)
#Save duals: shadowprices
HOM.model.dual  = Suffix(direction=Suffix.IMPORT)
print(time.time()-t)
print('*Model created*')

#%% SOLVE MODEL
print('Solving model...')
t=time.time()
solver = SolverFactory(SOLVER,executable=SOLVERPATH) if SOLVERPATH != 0 else SolverFactory(SOLVER)
#if solver.name == 'ipopt':
    #solver.options['linear_solver']='ma97'
    #solver.options['mu_strategy']='adaptive' 
    #solver.options['bound_relax_factor']=0

solverstatus    = solver.solve(HOM.model)#,tee=True,keepfiles=True,logfile=os.path.join(dirname, 'logREAD.log')) #

if HOM.model.Options['Investment module'] == 1:
#Fix selected investment to switch to fully linear model and get dual values
   for ip in HOM.model.ninvphase:
       for inv in HOM.model.ninvest:   
           HOM.model.IbINVEST[ip,inv].fixed = True
   solverstatus = solver.solve(HOM.model)   

print(solverstatus)
print(time.time()-t)

#%%EXPORT RESULTS
print('Exporting results...')
t=time.time()

#Create result folder
if not os.path.exists(ResultFolderPath):
    os.makedirs(ResultFolderPath)

#Read and export results
exppath = ResultFolderPath + os.sep + SCENARIO + '.xlsx' 
result_analysis = ResultAnalysis()
if EXPORT in ['xlsx','all']:
    results=result_analysis.readresults(HOM.model,parameters,solverstatus,scenario=SCENARIO,PRINT=1,VALIDATION=1) #   
    result_analysis.export_to_excel(results,exppath,NEWSHEET=NEWSHEET,VALIDATION=1) #export formated results to excel 
if EXPORT in ['all']:
    result_analysis.export_all_DV(HOM.model,ResultFolderPath,scenario=SCENARIO) #export all decision variables
    result_analysis.export_mass_balances(results,ResultFolderPath,scenario=SCENARIO) #export mass balances (for powerBI)
    result_analysis.export_index_mapping(HOM.model,os.path.join(ResultFolderPath,'Mapping')) #export index mapping/connections (for powerBI)
selectedresults=result_analysis.selectedresults(HOM.model,parameters,scenario=SCENARIO)
pickle.dump(selectedresults,open(os.path.join(ResultFolderPath,SCENARIO + '.txt'),"wb"))
print(time.time()-t)
print('*Result Exported*')
