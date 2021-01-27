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
import time
import os
import pickle
import sys
from pyomo.environ                import Suffix, SolverFactory, value, Var #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Collects results and exports to excel sheets
#Parallel computing
import multiprocess as mp
import numpy as np
import itertools
import pandas as pd
import copy

#%%OPTIONS
SPECIALCASE=0#'eg'#'food'#0#
RUNS=300 #'all': does all possible combinations
RANDOMCAPEX=0
RANDOMDISCOUNT=1
LOADSCEN='montecarlo_scenarios_300.txt'#'save'# # 0# : nothing special, 'save' :save montecarlo_scenarios.txt (to result folder), paths : loads scenarios previously saved (placed in Data)
LOADVAR ='modelvars.txt' #'export'# (0: do nothing, 'export': exports solved model variables as txt file, 'whatevername.txt': saves model variables of first run and then uses it for the next runs)
ANALYSIS='invest' # 'reanalysis': NPV with - without all investments 'invest' selected investments for scenario, 'indinv' adds with-without analysis of each individual investment
#Export options
UPDATE = 0 #0 updates all parameters, 1 updates only selected parameters
RESULTFOLDER = 'zam_montecarlo_energyblind300_'
EXPORT =None#'ifpri'# 
#parallel options
PARALLEL_scenario = 1 #Run scenarios in parallel
epll=24
#solver options
THREADS = 1 #LSB_DJOB_NUMPROC
SOLVER = 'cplex' #'cbc' #'cplex'
SOLVERPATH = 0#'~/CoinIpopt/bin/ipopt' #'/home/software/ipopt/3.12/bin/ipopt'#'~/miniconda3/pkgs/coincbc-2.10.5-hab63836_0/bin/cbc'#0 0 is default, precise path only if required by solver 

#%%
#Define single scenario analysis function    
def ScenarioAnalysis(ss,parameters_in,solver):
    tt=time.time()
    np.random.seed()
    #copy parameters
    parameters=copy.deepcopy(parameters_in)
    #Collect extra parameters
    parameters.load_csv(DataFolderPath,ss)
    
    invest={inv:{} for inv in parameters.val['ninvest']}
    #add a Random CAPEX for individual investments and a Discount Rate (general)
    if RANDOMDISCOUNT==1:
        dr=round(np.random.uniform(8,12))
        parameters.val['Options']['Discount rate',parameters.val['sOptions'][ss]]=dr
    if RANDOMDISCOUNT==1 or RANDOMCAPEX==1:
        for inv in parameters.val['ninvest']:            
            if RANDOMCAPEX==1:
                f=np.random.uniform(0.9,1.5)
                parameters.val['iCAPEX'][inv]=round(parameters.val['iCAPEX'][inv]*f)
                invest[inv]['iCAPEX']=parameters.val['iCAPEX'][inv]
            invest[inv]['sDiscount']=dr

    #Create model
    print('Formulating scenario: ' + ss)
    HOM             = HydroeconomicOptimization(parameters,scenario=ss)
    HOM.model.dual  = Suffix(direction=Suffix.IMPORT)
    
    #load variables
    if LOADVAR not in [0,'export']:
        varpath=os.path.join(ResultFolderPath,LOADVAR)
        if os.path.exists(varpath):
            LoadVar=pickle.load(open(varpath,"rb"))
            for var in HOM.model.component_objects(Var):
                if str(var) in LoadVar.keys():
                    for index in var:
                        if index in LoadVar[str(var)].keys():
                            var[index].value=LoadVar[str(var)][index]
                        
    #Solve scenario
    print('...solving scenario: ' + ss)
    solverstatus= solver.solve(HOM.model)
    if HOM.model.Options['Investment module'] in [1,'continuous']:
    #Fix selected investment to switch to fully linear model and get dual values
        for ip in HOM.model.ninvphase:
            for inv in HOM.model.ninvest:
                HOM.model.IbINVEST[ip,inv].fixed = True
                if HOM.model.Options['Investment module'] in ['continuous']:
                    if inv in HOM.model.ninvestb:
                        HOM.model.IbINVESTb[ip,inv].fixed = True
        print('...re-solving scenario ' + ss)
        solverstatus = solver.solve(HOM.model)        
    print('scenario '+ ss +' solved in '+str(time.time()-tt)+' seconds')
    
    #Save variables
    if LOADVAR != 0:
        varpath=os.path.join(ResultFolderPath,LOADVAR)
        if LOADVAR == 'export' or not os.path.exists(varpath):
            LoadVar={str(var):{index:var[index].value for index in var} for var in HOM.model.component_objects(Var)}
            pickle.dump(LoadVar,open(varpath,"wb"))
    
    #Export results
    exppath = ResultFolderPath + os.sep + ss + '.xlsx' 
    result_analysis = ResultAnalysis()
    try:
        if EXPORT in ['txt','all','xlsx','ifpri']:
            output=result_analysis.selectedresults(HOM.model,parameters,scenario=ss)
            outputpath=ResultFolderPath + os.sep + ss + '.txt' 
            pickle.dump(output,open(outputpath,"wb"))
        if EXPORT in ['xlsx','all']:
            results=result_analysis.readresults(HOM.model,parameters,solverstatus,PRINT=0,scenario=ss,VALIDATION=0) #   
            result_analysis.export_to_excel(results,exppath,NEWSHEET=1,VALIDATION=0) #export formated results to excel
        if EXPORT in ['ifpri','all']:          
                result_analysis.export_ifpri_indicators(HOM.model,ResultFolderPath,scenario=ss) 
        if EXPORT in ['all']:
            result_analysis.export_all_DV(HOM.model,ResultFolderPath,scenario=ss) #export all decision variables
            result_analysis.export_mass_balances(results,ResultFolderPath,scenario=ss) #export mass balances (for powerBI)
            result_analysis.export_index_mapping(HOM.model,os.path.join(ResultFolderPath,'Mapping')) #export index mapping/connections (for powerBI)
    except:
        print('results of run '+ss+' could not be exported')
    #Deactivate investment constraints
    if ANALYSIS in ['reanalysis','indinv']:
        HOM.model.invest_after.deactivate()
        HOM.model.invest_year.deactivate()
        HOM.model.invest_binary.deactivate()
    #Save economic balance of optimal invesmtents
    if ANALYSIS in ['reanalysis','indinv']:
        #outputinv=result_analysis.selectedresults(HOM.model,parameters,scenario=ss)
        ebalanceinv=-value(HOM.model.obj)#outputinv['EconomicBalance']['Total']
    
    #Save investment decisions
    if ANALYSIS in ['invest','indinv']:
        #save investment decisions
        for inv in HOM.model.ninvest:
            invest[inv]['INVEST']={ip:round(HOM.model.IbINVEST[ip,inv].value*100) 
                                   for ip in HOM.model.ninvphase}
            invest[inv]['type']=HOM.model.iInvType[inv]
            invest[inv]['INVcap']=HOM.model.iInvCap[inv]
        #delete extra keys (some investments might be in parameters but not specific scenarios)
        for inv in parameters.val['ninvest']:
            if inv not in HOM.model.ninvest:
                invest.pop(inv)
    
    #Save individual NPV of investments for each phase    
    if ANALYSIS in ['indinv']:
        print('...individual scenario analysis ' + ss)
        tt=time.time()
        #NPV with-without analysis for each investment and phase   
        for inv in HOM.model.ninvest: #loop through investments
            invest[inv]['invNPV']={ip:0 for ip in HOM.model.ninvphase}
            invphase=None
            invval=1
            #Find investment (if so) and set to 0
            for ip in HOM.model.ninvphase:               
                if HOM.model.IbINVEST[ip,inv].value>0.05: #ADD threshold value
                    invphase=ip #save phase
                    invval=HOM.model.IbINVEST[ip,inv].value #save value
                    HOM.model.IbINVEST[ip,inv].value=0 #set value to 0
            #Total economic benefits without investment
            if invphase == None: #use directly original output as investment value is 0
                ebalance0=ebalanceinv                
            else: #resolve model without investment
                solverstatus = solver.solve(HOM.model)
                ebalance0=-value(HOM.model.obj)#output0['EconomicBalance']['Total']
            #Investment with-without analysis for each phase
            for ip in HOM.model.ninvphase:
                if invphase == ip: #use directly original output as investment is invested in that phase
                    invest[inv]['invNPV'][ip]=round(ebalanceinv-ebalance0,2)
                else: #resolve model with investment
                    HOM.model.IbINVEST[ip,inv].value=invval
                    solverstatus = solver.solve(HOM.model)
                    ebalance1=-value(HOM.model.obj)
                    invest[inv]['invNPV'][ip]=round(ebalance1-ebalance0,2)
                    #set back investment value to 0
                    HOM.model.IbINVEST[ip,inv].value=0
            #set back investment value to original value
            if invphase != None:
                HOM.model.IbINVEST[invphase,inv].value=invval
    
    #Save NPV of total investment plan (by comparing to without any investment)
    if ANALYSIS in ['reanalysis','indinv']:
        for inv in HOM.model.ninvest:
            for ip in HOM.model.ninvphase:
                HOM.model.IbINVEST[ip,inv].value=0
        solverstatus = solver.solve(HOM.model)
        #output0=result_analysis.selectedresults(HOM.model,parameters,scenario=ss)
        ebalance0=-value(HOM.model.obj)#output0['EconomicBalance']['Total']
        if ANALYSIS in ['indinv']:
            for inv in HOM.model.ninvest: #ADD iteration because of structure of invest
                invest[inv]['NPV']=round(ebalanceinv-ebalance0,2)
        else:
            ebalance=round(ebalanceinv-ebalance0,2)
        print('...additional analysis '+ ss +' solved in '+str(time.time()-tt)+' seconds')
    
    #return agruments    
    if ANALYSIS in ['reanalysis']:
          return ebalance
    if ANALYSIS in ['invest','indinv']:
          return invest

#%% Start parallel solving of problems
if __name__ == '__main__':
    #Define paths
    ResultFolderPath= os.path.join(dirname, 'Results',RESULTFOLDER+time.strftime("%d_%m_%Y_%Hh%M"))
    DataFolderPath  = os.path.join(dirname, 'Data')    
    #Collect parameters
    Main            = DataFolderPath + os.sep + 'MainFile_ex.xlsx'
    Water           = DataFolderPath + os.sep + 'WaterModule_ex.xlsx'
    Agriculture     = DataFolderPath + os.sep + 'AgricultureModule_ex.xlsx'
    CropMarket      = DataFolderPath + os.sep + 'CropMarketModule_ex.xlsx'
    Energy          = DataFolderPath + os.sep + 'EnergyModule_ex.xlsx'
    Investment      = DataFolderPath + os.sep + 'InvestmentModule_ex.xlsx'
    Param           = DataFolderPath + os.sep + 'Parameters.txt' #parameters (python dictionnaries) saved as txt    
    print('Harvesting parameters...')
    parameters      = Database(update=UPDATE,DataFile=Param)
    parameters.harvest_all([Main,Water,Agriculture,CropMarket,Energy])
    investopt=[parameters.val['Options'][k] for k in parameters.val['Options'].keys() if k[0]=='Investment module'] #investment module options
    if 1 in investopt or 'continuous' in investopt:
        parameters.harvest_all([Investment])
    parameters.save(Param) #Save parameters
    print('*Parameters harvested*')
    
    #Create folders
    if not os.path.exists(ResultFolderPath):
        os.makedirs(ResultFolderPath)
        
    #Choose solver
    solver = SolverFactory(SOLVER,executable=SOLVERPATH) if SOLVERPATH != 0 else SolverFactory(SOLVER)
    if solver.name == 'ipopt':
        solver.options['linear_solver']='ma97' 

    if solver.name == 'cplex':
        solver.options['threads']=THREADS
    
    #Generate scenarios
    #Definitions
    nS=parameters.val['nsampling']
    mP=parameters.val['mParameter']
    mV=parameters.val['mValue']
    sn0=parameters.val['nscenario'][0]
    scenarios=[]
    scenparam=[p for p in parameters.val.keys() if str(p)[0]=='s']
    #All scenario combinations
    if RUNS == 'all':
        sP=[[mP[k] for k in mP.keys() if k[0]==pn] for pn in nS]
        pO={nS[k]:k for k in range(len(nS))}
        allcombinations=list(itertools.product(*sP))
        nruns=len(allcombinations)
    else: #Montecarlo
        nruns=RUNS    
    #Create scenarios
    if LOADSCEN in ['save',0]:
        memory={m:{} for m in range(nruns)}
        for m in range(nruns):
            sn='m'+str(m)
            scenarios.append(sn)
            for pn in scenparam:                
                if pn in nS: #parameter with uncertainty
                    if RUNS == 'all': #runs through all possible scenarios
                        pscen=allcombinations[m][pO[pn]]
                        pval=[mV[k] for k in mV.keys() if k[0]==pn and pscen==mP[k]][0]
                    else: #random selection of scenario
                        choices=[k for k in mP.keys() if k[0]==pn]
                        choice=choices[np.random.choice(range(len(choices)))]
                        pscen=mP[choice]
                        pval=mV[choice]
                    parameters.val[pn][sn]=pscen
                    memory[m][pn]=pval              
                else: #parameter without uncertainty
                    parameters.val[pn][sn]=parameters.val[pn][sn0]
        if LOADSCEN == 'save':
            #save the montecarlo scenarios
            savescen={pn:parameters.val[pn] for pn in scenparam}
            savescen.update({'memory':memory,
                             'scenarios':scenarios})
            outpath=os.path.join(ResultFolderPath,'montecarlo_scenarios.txt')
            pickle.dump(savescen,open(outpath,"wb"))
    #Load existing scenarios
    else:
        loadpath=os.path.join(DataFolderPath,LOADSCEN)
        loadscen=pickle.load(open(loadpath,"rb"))
        scenarios=loadscen['scenarios']
        memory=loadscen['memory']
        for pn in scenparam:
            if pn in nS:
                parameters.val[pn]=loadscen[pn]
                #ADD custom special case for synthetic case with silo versus nexus
                if SPECIALCASE in ['noag']:
                    if pn == 'sUserDem':
                        for sn in scenarios:
                            parameters.val[pn][sn]='noag_'+parameters.val[pn][sn]
            else:
                for sn in scenarios:
                    parameters.val[pn][sn]=parameters.val[pn][sn0]
                   
    #Run scenarios
    #scenarios = parameters.val['nscenario']
    if PARALLEL_scenario == 1:
        pool = mp.Pool(min(epll,len(scenarios)))    
        parallelresults = pool.starmap(ScenarioAnalysis, [(ss,parameters,solver) for ss in scenarios])
        pool.close()
    else:
        parallelresults=[ScenarioAnalysis(ss,parameters,solver) for ss in scenarios]
    
    #%%Investment selection results per scenario
    idx=pd.IndexSlice
    if ANALYSIS in ['invest','indinv']:
        #Create dataframe hosting results
        ninvest=[k for k in parallelresults[0].keys()] #investments
        ninvphase=parallelresults[0][ninvest[0]]['INVEST'].keys() #investment phases        
        mindex=pd.MultiIndex.from_product([range(nruns),ninvest,ninvphase], 
                                          names=['nrun','ninvest','ninvphase'])
        columns=np.append(parameters.val['nsampling'],['INVEST','type','INVcap'])
        if RANDOMCAPEX==1:
            columns=np.append(columns,['iCAPEX'])
        if RANDOMDISCOUNT==1:
            columns=np.append(columns,['sDiscount'])
        if ANALYSIS in ['indinv']:
            columns=np.append(columns,['invNPV','NPV'])         
        pinvest=pd.DataFrame(columns=columns,index=mindex,dtype=float)
            
        for k in range(len(scenarios)):
            invest=parallelresults[k]
            sn=scenarios[k]        
            for ninv in invest.keys():
                pinvest.loc[idx[k,ninv,:],'type']=invest[ninv]['type']
                for pn in parameters.val['nsampling']:
                    pinvest.loc[idx[k,ninv,:],pn]=memory[k][pn]
                if RANDOMCAPEX==1:
                    pinvest.loc[idx[k,ninv,:],'iCAPEX']=invest[ninv]['iCAPEX']
                if RANDOMDISCOUNT==1:
                    pinvest.loc[idx[k,ninv,:],'sDiscount']=invest[ninv]['sDiscount']                    
                for ip in ninvphase:
                    pinvest.loc[idx[k,ninv,ip],'INVEST']=invest[ninv]['INVEST'][ip]
                    pinvest.loc[idx[k,ninv,ip],'INVcap']=invest[ninv]['INVEST'][ip]*invest[ninv]['INVcap']/100
                    if ANALYSIS in ['indinv']:
                        pinvest.loc[idx[k,ninv,ip],'NPV']=invest[ninv]['NPV']
                        pinvest.loc[idx[k,ninv,ip],'invNPV']=invest[ninv]['invNPV'][ip]                   
                    
        pinvest.to_csv(os.path.join(ResultFolderPath,'INVEST.csv'),sep=';',decimal=',')
    
    #%%Get NPV of investments for each scenario
    if ANALYSIS in ['reanalysis']:     
        #create result dataframe
        mindex=range(nruns)
        columns=np.append(parameters.val['nsampling'],['NPV'])
        pnpv=pd.DataFrame(columns=columns,index=mindex)
        pnpv.index.name='nrun'
        for k in range(len(scenarios)):
            #NPV=difference of economic balance with and without investments
            NPV=parallelresults[k]
            pnpv.loc[k,'NPV']=round(NPV)
            #save scenario parameters
            for pn in parameters.val['nsampling']:
                pnpv.loc[k,pn]=memory[k][pn]
                
        pnpv.to_csv(os.path.join(ResultFolderPath,'NPV.csv'),sep=';',decimal=',')
        
#%% Define investment plans
