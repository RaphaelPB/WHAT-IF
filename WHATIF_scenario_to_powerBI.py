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

import sys
import os
import pickle
import pandas as pd
dirname = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(dirname, 'bin'))


#%%OPTIONS - MODIFY BY USER
#SCENFILE=os.path.join(dirname,'Data','Scenarios_to_compare.xlsx')

SHEET='investAll'
FOLDERNAME='iALL_mpc3b_naThu07_05_2020_17h08'
DIFFMODE=0 # if=1 exports relative results to POWERBI, instead of absolute (=0)
SCENFILE='Scenarios_to_compare.xlsx'
result_path = os.path.join(dirname,'Results',FOLDERNAME)
#%%#####################################
#               FUNCTIONS

#Export panda dataframe to .csv
def pd_to_csv(folder,file,pdata):
    file = os.path.join(folder,file)
    if isinstance(pdata,pd.DataFrame):
        return pdata.to_csv(file,sep=';',decimal=',',index='False')
#Transfor dictionary in panda dataframe (for specific dic structure)
def dict_to_pd(dic):
    if dic != {}:
        #creates dataframe from nested dic, moves left column to index, and first index to column
        pdata=pd.concat({k: pd.DataFrame(v).T for k, v in dic.items()}, axis=0)
        if not pdata.empty:
            pdata=pdata.stack().unstack(0)
        return pdata

def aggregate_scenarios_to_csv(scenarios,vardic,outpath,keytype='tuple',indexname=0,elist=0,refscen=0,renamecol=1):
    ##scenarios: list of scenarios
    ##vardic: dictionary containing all elements for alls scenarios, in the form:
    #vardic ={scen:{elem:{...} for elem in elist} for scen in scenarios}
    ##outpath: folder to export to
    ##keytype: how to read vardic[scen][elem], 
    #'tuple' = dic in the form {(i1,...,in):result for i1 in I1 ... for in in In}
    #'nested' = dic in the form {i1:...{in:result for in in In}... for i1 in I1}
    ##indexname : names of the indexes [i1_name,...,in_name]
    ##elist : list of all elements to be taken from vardic - if =0, all elements are exported
    ##refscen : Use differential (A-B) values of scenarios according to defined reference scenario: 
    #refscen = {scen: relative_scen for scen in scenarios}, relative_scen has to be in vardic.keys()
    ##renamecol: 1 = renames column by elem, 0 = keeps default name (or leaves empty)
    
    #create outpath
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    #diff function - does panda-pandaref + keeps panda where pandaref does not exist
    def diff_panda(panda,pandaref):
        pandaf=panda.subtract(pandaref)
        pandaf=pandaf.fillna(panda)
        return pandaf
    
    #generate elist if not passed
    if elist==0:
        elist=[]
        for scen in scenarios:
            for key in vardic[scen].keys():
                if key != 'AlCULAREAXXXX': #solve troubles
                    elist.append(key)        
        elist=set(elist)
    
    for elem in elist:
        if refscen==0:
            #assemble different scenarios
            if keytype=='tuple':
                frames=[pd.Series(vardic[scen][elem]) for scen in scenarios if elem in vardic[scen].keys()] 
            elif keytype=='nested':
                frames=[dict_to_pd(vardic[scen][elem]) for scen in scenarios if vardic[scen][elem]!={}]
        else:
            #assemble different scenarios and subtract ref scenario
            if keytype=='tuple':
                frames=[diff_panda(pd.Series(vardic[scen][elem]),pd.Series(vardic[refscen[scen]][elem])) for scen in scenarios if elem in vardic[scen].keys()]
            elif keytype=='nested':
                frames=[diff_panda(dict_to_pd(vardic[scen][elem]),dict_to_pd(vardic[refscen[scen]][elem])) for scen in scenarios if vardic[scen][elem]!={}]
        
        #concatane to single dataframe
        if refscen==0:
            pdata=pd.concat(frames,keys=[scen for scen in scenarios if elem in vardic[scen].keys()])
        else:
            pdata=pd.concat(frames,keys=[scen+'_'+refscen[scen] for scen in scenarios if elem in vardic[scen].keys()])
        if keytype=='tuple':
            pdata=pdata.to_frame() 
        #elif keytype=='nested':
         #   pdata=pd.concat(frames,keys=[scen for scen in scenarios if elem in vardic[scen].keys()])
        
        if not pdata.empty:
            #rename columns and indexes
            if indexname != 0:
                indexname[elem].insert(0,'scenario')
                pdata.index.names=indexname[elem]
                if renamecol==1:
                    pdata.columns=[elem]
            #export to csv
            filename=str(elem)+'.csv'
            pd_to_csv(outpath,filename,pdata=pdata)

        if indexname != 0:
            for elem in indexname.keys():
                if elem not in elist:
                    indexname[elem].insert(0,'scenario')
                    pdata=pd.DataFrame(columns=indexname[elem]+[elem])
                    filename=str(elem)+'.csv'
                    pdata.to_csv(os.path.join(outpath,filename),sep=';',decimal=',',index=False)
                    
#%%#####################################
#        COLLECT INFORMATION
        
sceninfo=pd.read_excel(SCENFILE,sheet_name=SHEET, skiprows=1, index_col=[0], engine='xlrd')
refscen=sceninfo.to_dict()['refscen']
scenarios=[s for s in refscen.keys()]
#scen_to_load=[s for s in refscen.keys()]


#%%#####################################
#               POWERBI

#%%DECISION VARIABLES
#create directory
DVpath=os.path.join(result_path,'DecisionVariables')

#Varaibles index names  #ADD: Variables can have different indexs (eg CULAREA)
VarIndex={
 'AcEXTPROD':['nyear', 'ncmarket', 'ncrop'],
 'AcPROD':['nyear', 'nfzone', 'ncrop'],
 'AcSUPPLY':['nyear', 'ncmarket', 'ncrop','ncdstep'],
 'AcTRANS':['nyear', 'nctrans', 'ncrop'],
 'AlCULAREA':['nyear', 'nfzone', 'nflied', 'nfieldculture'],
 'CULAREA':['nyear','nfzone','nculture'],
 'AwSUPPLY':['ntime', 'nfzone', 'nculture'],
 'DUMYCROP':['nyear', 'nfzone', 'ncrop'],
 'DUMYEFLOW':['ntime', 'neflow'],
 'DUMYSTOR':['nres'],
 'DUMYWATER':['ntime','ncatch'],
 'EeGENCAP':['nyear', 'nptech','npmarket'],
 'EeGENPROD':['ntime', 'npload', 'nptech','npmarket'],
 'EeHPPROD':['ntime', 'npload', 'nhpp'],
 'EeOPPROD':['ntime', 'npload', 'nopp'],
 'EeSUPPLY':['ntime', 'npload', 'npmarket'],
 'EeTRANS':['ntime', 'npload', 'ntransline'],
 'EwHPDISCHARGE':['ntime', 'npload', 'nhpp'],
 'WwGWSTORAGE':['ntime', 'naquifer'],
 'WwOUTFLOW':['ntime', 'ncatch'],
 'WwRSTORAGE':['ntime', 'nres'],
 'WwSUPPLY':['ntime', 'nuser'],
 'WwTRANSFER':['ntime', 'ntransfer'],
 'energy_shadow':['ntime','npload','npmarket'],
 'crop_shadow':['nyear','ncmarket','ncrop','ncdstep'],
 'water_shadow':['ntime','ncatch']}

#load all decision variables (DV) from model runs
scen_to_load=set([s for s in refscen.keys()]+[s for s in refscen.values() if s==s]) 
ScenarioDV={scen:pickle.load(open(os.path.join(result_path,scen+'_DV.txt'),"rb")) for scen in scen_to_load}
#Assemble and Export to csv
REF = 0 if DIFFMODE == 0 else refscen
aggregate_scenarios_to_csv(scenarios,ScenarioDV,DVpath,keytype='tuple',indexname=VarIndex,elist=0,refscen=REF,renamecol=1)


#%%BALANCES (Economic, Energy, Crop, Water)
#create directory
Bpath=os.path.join(result_path,'Balances')
#balances index names
BalIndex={'EconomicBalance':['nyear','ncountry'],
          'EnergyBalance':['ntime','npmarket'],
          'CropBalance':['nyear','ncmarket'],
          'WaterBalance':['ntime','ncatch']}

        
#load all balances from model runs
ScenarioB={scen:pickle.load(open(os.path.join(result_path,scen+'_Balances.txt'),"rb")) for scen in scen_to_load}    
#Assemble and Export to csv
REF = 0 if DIFFMODE == 0 else refscen
aggregate_scenarios_to_csv(scenarios,ScenarioB,Bpath,keytype='nested',indexname=BalIndex,elist=0,refscen=REF,renamecol=0)

#for k in range(4):
#    frames=[dict_to_pd(ScenarioB[scen][k]) for scen in scenarios if ScenarioB[scen][k]!={}] #assemble Balances from different scenarios
#    if frames != []:
#        pdata=pd.concat(frames,keys=scenarios)  #concatane into single dataframe
#        pdata.index.names=BalIndex[k] #rename indexes
#        filename=['EconomicBalance.csv','EnergyBalance.csv','CropBalance.csv','WaterBalance.csv'][k]+'.csv'
#        pd_to_csv(Balancepath,filename,pdata=pdata)
        
