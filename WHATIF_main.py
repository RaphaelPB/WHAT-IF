# -*- coding: utf-8 -*-
"""
    Copyright 2018 Raphaël Payet-Burin
    Copyright 2018 Mikkel Kromann
    Copyright 2018 Technical University of Denmark
    Copyright 2018 COWI A/S

    This file is part of WHAT-IF.

    WHAT-IF is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 (GPLv3) as published by the Free Software Foundation,

    WHAT-IF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with WHAT-IF. If not, see http://www.gnu.org/licenses/.

"""
#Import python libraries: pandas, xlswriter, openpyxl, xlrd, pyomo, glpk 
from __future__                   import division                   #As in python 2.7 division between integers gives an integer - avoid problems
import time
import os
import pickle
import sys
from pyomo.environ                import Suffix, SolverFactory, Var #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Exports results to excel sheets

#%% DEFINE PARAMETERS
#OPTIONS
    #Data Folder
DataFolderPath  = os.path.join(dirname, 'Data')
ResultFolderPath= os.path.join(dirname, 'Results',time.strftime("%a%d_%m_%Y_%Hh%M"))
    #Excel output file
NEWSHEET        = 1 #1 creates a new sheet, 0 fills existing sheet
OutPath         = ResultFolderPath + os.sep + 'RESULTS_' + time.strftime("%d_%m_%Y_%Hh%M") + '.xlsx'
    #Python obj output file (txt)
ResultExport    = ResultFolderPath + os.sep + 'RESULTS.txt' #results (python objects) saved as txt
ModelExport     = ResultFolderPath + os.sep + 'MODEL.txt' #solved pyomo model (unmodified python objects) saved as txt
    #Load excel files or existing python object
UPDATE          = 0 #0 updates all parameters, 1 updates only selected parameters

#Define path to data
Main            = DataFolderPath + os.sep + 'MainFile.xlsx'
Water           = DataFolderPath + os.sep + 'WaterModule.xlsx'
Agriculture     = DataFolderPath + os.sep + 'AgricultureModule.xlsx'
CropMarket      = DataFolderPath + os.sep + 'CropMarketModule.xlsx'
Energy          = DataFolderPath + os.sep + 'EnergyModule.xlsx'
Param           = DataFolderPath + os.sep + 'Parameters.txt' #parameters (python dictionnaries) saved as txt

#Collect parameters
t=time.time()
print('Harvesting parameters...')
parameters      = Database(update=UPDATE,DataFile=Param)
parameters.harvest_all([Main,Water,Agriculture,CropMarket,Energy]) 

#Save parameters    
parameters.save(Param)
print(time.time()-t)
print('*Parameters harvested*')

#%% BUILD HYDROECONOMIC OPTIMIZATION MODEL
t=time.time()
print('Creating Hydroeconomic optimization model...')
HOM             = HydroeconomicOptimization(parameters,scenario='WHATIF_main')
print(time.time()-t)
print('*Model created*')

#%% SOLVE MODEL
print('Solving model...')
t=time.time()
#Save duals: shadowprices
HOM.model.dual  = Suffix(direction=Suffix.IMPORT)

#Choose solver
solver          = SolverFactory('cplex') 

#Solve
solverstatus    = solver.solve(HOM.model)   

print(solverstatus)
print(time.time()-t)

#%%EXPORT RESULTS
print('Exporting results...')
t=time.time()

#Create result folder
if not os.path.exists(ResultFolderPath):
    os.makedirs(ResultFolderPath)

#Create result python object
results     = ResultAnalysis()
results.readresults(HOM.model,parameters,solverstatus,scenario='WHATIF_main',PRINT=1)
results.exportresults(OutPath,HOM.model,NEWSHEET=NEWSHEET)
print(time.time()-t)
print('*Result Exported*')

#Save results
pickle.dump(results,open(ResultExport,"wb")) 
