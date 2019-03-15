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
#Import python packages
from __future__                   import division   #As in python 2.7 division between integers gives an integer - avoid problems
import time
import os
import pickle
import sys
import pandas as pd
from pyomo.environ                import Suffix, SolverFactory #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Collects results and exports to excel sheets
#Parallel computing
from joblib import Parallel, delayed
import multiprocessing


#Define paths to data
DataFolderPath  = os.path.join(dirname, 'Data')
Main            = DataFolderPath + os.sep + 'MainFile.xlsx'
Water           = DataFolderPath + os.sep + 'WaterModule.xlsx'
Agriculture     = DataFolderPath + os.sep + 'AgricultureModule.xlsx'
CropMarket      = DataFolderPath + os.sep + 'CropMarketModule.xlsx'
Energy          = DataFolderPath + os.sep + 'EnergyModule.xlsx'
Param           = DataFolderPath + os.sep + 'Parameters.txt' #parameters (python dictionnaries) saved as txt

#Define paths to results
ResultFolderPath= os.path.join(dirname, 'Results',time.strftime("%a%d_%m_%Y_%Hh%M"))
ScenarioOutPath = ResultFolderPath + os.sep + 'SCENARIOS_results.xlsx'

#Activate or not investment module
INVEST=1
SAVEALL=1 # 1: save all results, 0: only saves selected results

#Collect parameters
t=time.time()
print('Harvesting parameters...')
parameters      = Database(update=0,DataFile=Param)
parameters.harvest_all([Main,Water,Agriculture,CropMarket,Energy])

#Save parameters
parameters.save(Param)

#Create result folder
if not os.path.exists(ResultFolderPath):
    os.makedirs(ResultFolderPath)

#Define single scenario analysis function    
def ScenarioAnalysis(ss,parameters):
    
    #Create model  
    print('formulating scenario: ' + ss)
    HOM             = HydroeconomicOptimization(parameters,scenario=ss)
    HOM.model.dual  = Suffix(direction=Suffix.IMPORT)    
    
    #Solve scenario
    print('...solving scenario ' + ss)
    solver      = SolverFactory('cplex') #choose solver
    solverstatus= solver.solve(HOM.model)
    
    #Collect and export results
    print('...saving and exporting results of ' + ss)
    results     = ResultAnalysis()
    if SAVEALL==1:
        #Collect and export all results for individual scenario
        if parameters.val['sExportExcel'][ss] != 0:
            results.readresults(HOM.model,parameters,solverstatus,PRINT=0,scenario=ss)
            OutPath = ResultFolderPath + os.sep + ss + '.xlsx'             
            results.exportresults(OutPath,HOM.model,NEWSHEET=1)
        if parameters.val['sExportPython'][ss] != 0:
            Result = ResultFolderPath + os.sep + ss + '.txt'
            pickle.dump(results,open(Result,"wb"))
    
        #Save selected results for scenario analysis
    results.selectedresults(HOM.model,parameters,scenario=ss)
    output={'EconomicBalance':                  results.EconomicBalance,
            'EconomicBalance_co':               results.EconomicBalance_co,
            'WaterEconomicBalance_co':          results.WaterEconomicBalance_co,
            'EnergyEconomicBalance_co':         results.EnergyEconomicBalance_co,
            'AgricultureEconomicBalance_co':    results.AgricultureEconomicBalance_co,
            'AgricultureEconomicBalance2_co':   results.AgricultureEconomicBalance2_co,
            'HydropowerProduction_co':          results.HydropowerProd_co,
            'GrossCultivatedArea_co':           results.GrossCulArea_co,
            'NetIrrArea_co':                    results.NetIrrArea_co,
            'IrrigNetCons_co':                  results.IrrigNetCons_co,
            'EnergyValue_co':                   results.EnergyValue_co,
            'CropValue_co':                     results.CropValue_co,
            'OtherIndicators':                  results.OtherIndicators}
    if HOM.model.Options['Energy market'] == 1:
        output['EnergyExportBenefit_co']    = results.EnergyExportBenefit_co
        output['EnergyExports_co']          = results.EnergyExports_co
        output['EnergyCapacityConst_co']    = results.EnergyCapacityConst_co
        output['SolarCapacityConst_co']     = results.SolarCapacityConst_co
        output['CO2Emission_co']            = results.CO2Emission_co
        output['CropExpBenefit_co']         = results.CropExpBenefit_co
        output['EnergyConsSurplus_co']      = results.EgyConsSurplus_co
        output['EnergyProdSurplus_co']      = results.EgyProdSurplus_co
    if HOM.model.Options['Crop market'] == 1:
        output['AgricultureConsSurplus_co'] = results.AgrConsSurplus_co
        output['AgricultureProdSurplus_co'] = results.AgrProdSurplus_co
    return output

def RelativeResults(data, parameters):
    ref_scen = parameters.val['sRefScen']
    relresult=[]
    for result in data:       
        relresult.append({scen+'_'+ref_scen[scen]:{elem:result[scen][elem]-result[ref_scen[scen]][elem] for elem in result[scen].keys() if elem in result[ref_scen[scen]].keys()} for scen in result.keys() if ref_scen[scen] in result.keys()})
    return relresult
        
#Parallel loops
num_cores = multiprocessing.cpu_count()
scenarios = parameters.val['nscenario']
parallelresults = Parallel(n_jobs=num_cores)(delayed(ScenarioAnalysis)(ss,parameters) for ss in scenarios)

#Save results 
EconomicBalance={scenarios[ss]:parallelresults[ss]['EconomicBalance'] for ss in range(len(scenarios))} #'Total':0,'Water':0,'Agriculture':0,'Energy':0
EconomicBalance_co={scenarios[ss]:parallelresults[ss]['EconomicBalance_co'] for ss in range(len(scenarios))}
WaterEconomicBalance_co={scenarios[ss]:parallelresults[ss]['WaterEconomicBalance_co'] for ss in range(len(scenarios))}
EnergyEconomicBalance_co={scenarios[ss]:parallelresults[ss]['EnergyEconomicBalance_co'] for ss in range(len(scenarios))}
AgricultureEconomicBalance_co={scenarios[ss]:parallelresults[ss]['AgricultureEconomicBalance_co'] for ss in range(len(scenarios))}
AgricultureEconomicBalance2_co={scenarios[ss]:parallelresults[ss]['AgricultureEconomicBalance2_co'] for ss in range(len(scenarios))}
HydropowerProduction_co={scenarios[ss]:parallelresults[ss]['HydropowerProduction_co'] for ss in range(len(scenarios))}
EnergyExportBenefit_co={scenarios[ss]:parallelresults[ss]['EnergyExportBenefit_co'] for ss in range(len(scenarios)) if 'EnergyExportBenefit_co' in parallelresults[ss].keys()}
EnergyExports_co={scenarios[ss]:parallelresults[ss]['EnergyExports_co'] for ss in range(len(scenarios)) if 'EnergyExports_co' in parallelresults[ss].keys()}
EnergyCapacityConst_co={scenarios[ss]:parallelresults[ss]['EnergyCapacityConst_co'] for ss in range(len(scenarios)) if 'EnergyCapacityConst_co' in parallelresults[ss].keys()}
SolarCapacityConst_co={scenarios[ss]:parallelresults[ss]['SolarCapacityConst_co'] for ss in range(len(scenarios)) if 'SolarCapacityConst_co' in parallelresults[ss].keys()}
GrossCultivatedArea_co={scenarios[ss]:parallelresults[ss]['GrossCultivatedArea_co'] for ss in range(len(scenarios))}
CropExpBenefit_co={scenarios[ss]:parallelresults[ss]['CropExpBenefit_co'] for ss in range(len(scenarios)) if 'CropExpBenefit_co' in parallelresults[ss].keys()}
NetIrrArea_co={scenarios[ss]:parallelresults[ss]['NetIrrArea_co'] for ss in range(len(scenarios))}
IrrigNetCons_co={scenarios[ss]:parallelresults[ss]['IrrigNetCons_co'] for ss in range(len(scenarios))}
CO2Emission_co={scenarios[ss]:parallelresults[ss]['CO2Emission_co'] for ss in range(len(scenarios)) if 'CO2Emission_co' in parallelresults[ss].keys()}
EnergyValue_co={scenarios[ss]:parallelresults[ss]['EnergyValue_co'] for ss in range(len(scenarios))}
CropValue_co={scenarios[ss]:parallelresults[ss]['CropValue_co'] for ss in range(len(scenarios))}
AgricultureConsSurplus_co={scenarios[ss]:parallelresults[ss]['AgricultureConsSurplus_co'] for ss in range(len(scenarios)) if 'AgricultureConsSurplus_co' in parallelresults[ss].keys()}
AgricultureProdSurplus_co={scenarios[ss]:parallelresults[ss]['AgricultureProdSurplus_co'] for ss in range(len(scenarios)) if 'AgricultureProdSurplus_co' in parallelresults[ss].keys()}
EnergyConsSurplus_co={scenarios[ss]:parallelresults[ss]['EnergyConsSurplus_co'] for ss in range(len(scenarios)) if 'EnergyConsSurplus_co' in parallelresults[ss].keys()}
EnergyProdSurplus_co={scenarios[ss]:parallelresults[ss]['EnergyProdSurplus_co'] for ss in range(len(scenarios)) if 'EnergyProdSurplus_co' in parallelresults[ss].keys()}
OtherIndicators={scenarios[ss]:parallelresults[ss]['OtherIndicators'] for ss in range(len(scenarios))}

#EXPORT OVERALL RESULTS
results = ResultAnalysis()
writer = pd.ExcelWriter(ScenarioOutPath, engine='openpyxl') #able to save formulas (graphs die however)

#EconomicBalance
data    = [EconomicBalance,EconomicBalance_co,WaterEconomicBalance_co,EnergyEconomicBalance_co,AgricultureEconomicBalance_co,AgricultureEconomicBalance2_co,
           AgricultureConsSurplus_co, AgricultureProdSurplus_co, EnergyConsSurplus_co, EnergyProdSurplus_co]
reldata = RelativeResults(data, parameters)
title   = ['Economic Balance [$]','Economic Balance [$]','Water Economic Balance [$]',' Energy Economic Balance [$]','Agriculture Economic Balance [$]','Agriculture Economic Balance 2[$]',
           'Agr cons surplus [$]','Agr prod surplus [$]','Egy cons surplus [$]', 'Egy prod surplus [$]']
results.exportsheet(writer,'EconomicBalance',data,title)
results.exportsheet(writer,'RelEconomicBalance',reldata,title)

#Other indicators
data    = [OtherIndicators, HydropowerProduction_co, EnergyExportBenefit_co, EnergyExports_co, GrossCultivatedArea_co, CropExpBenefit_co,
           NetIrrArea_co, IrrigNetCons_co, EnergyCapacityConst_co, SolarCapacityConst_co, CO2Emission_co, EnergyValue_co, CropValue_co]
reldata = RelativeResults(data, parameters)
title   = ['Other Indicators','Hydropower Production [kWh/y]','Energy Export benefit [$/y]','Energy Exports [kWh/y]','Gross Cultivated area [ha/y]', 'Crop export benefits [$/y]',
           'Net irrigated area [ha/y]','Agricultural water consumption [m3/y]','Additional Power Invest [kWh/month]','Solar Invest [kWh/month]','CO2 emissions [t/y]','Average energy price [$/kWh]','Average crop price [$/t]']
results.exportsheet(writer,'OtherIndicators',data,title)
results.exportsheet(writer,'RelOtherIndicators',reldata,title)

#Options
results.exportsheet(writer,'Options',[parameters.val['Options']],['OPTIONS'])
writer.save()

#EXPORT TO TXT
AllResultsPath  = ResultFolderPath + os.sep + 'SCENARIOS_results.txt'
AllResults      = {scenarios[ss]:parallelresults[ss] for ss in range(len(scenarios))}
AllResults['ref_scen'] = parameters.val['sRefScen']
AllResults['scenarios'] = parameters.val['nscenario']
pickle.dump(AllResults,open(AllResultsPath,"wb"))