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
import locale
#set decimal and separator (for handling of csv files)
langlocale = locale.getdefaultlocale()[0]
locale.setlocale(locale.LC_ALL, langlocale)
CSVDECIMAL = locale.localeconv()['decimal_point']
CSVSEPARATOR = ','
if CSVDECIMAL == ',':
    CSVSEPARATOR = ';'

#%%OPTIONS - MODIFY BY USER
SHEET='ifpri'
FOLDERNAME='nexus_300_ifpriidx'
DIFFMODE=0 # if=1 exports relative results to POWERBI, instead of absolute (=0)
SCENFILE='Scenarios_to_compare.xlsx'
result_path = os.path.join(dirname,'Results',FOLDERNAME)
#CSVSEPARATOR=';'
#CSVDECIMAL=','
#%%#####################################
#               FUNCTIONS

#Export panda dataframe to .csv
def pd_to_csv(folder,file,pdata,nround=3):
    file = os.path.join(folder,file)
    if isinstance(pdata,pd.DataFrame):
        return pdata.round(nround).to_csv(file,sep=CSVSEPARATOR,decimal=CSVDECIMAL,index='False')
#Transfor dictionary in panda dataframe (for specific dic structure)
def dict_to_pd(dic):
    if dic != {}:
        #creates dataframe from nested dic, moves left column to index, and first index to column
        pdata=pd.concat({k: pd.DataFrame(v).T for k, v in dic.items()}, axis=0)
        if not pdata.empty:
            pdata=pdata.stack().unstack(0)
            return pdata #WARNING: indented so that empty dataframes are ignored, might result in problems if all result exclude that component
        else:
            return pd.DataFrame(columns=pdata.index.levels[0])

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
    def diff_panda(pdat,pdatref):
        # if not pdatref.empty:
        pdatf=pdat.subtract(pdatref)
        pdatf=pdatf.fillna(pdat)
        return pdatf
        # else:
        #     return pdat
    
    #generate elist if not passed
    if elist==0:
        elist=[]
        for scen in scenarios:
            for key in vardic[scen].keys():
                if key != 'AlCULAREAXXXX': #solve troubles
                    elist.append(key)        
        elist=set(elist)

    for elem in elist:
        print(elem)
        #scenarios that contain that element
        if keytype=='tuple':
            scenpresent=[scen for scen in scenarios if elem in vardic[scen].keys()] #scenarios that contain that element
        elif keytype=='nested': #ADD could be deleted as vardic[scen][elem] is usually not empty even if no data but is the years
            scenpresent=[scen for scen in scenarios if vardic[scen][elem]!={}] #scenarios that contain that element
            
        if refscen==0:
            #assemble different scenarios
            if keytype=='tuple':
                frames=[pd.Series(vardic[scen][elem]) for scen in scenpresent] 
            elif keytype=='nested':
                frames=[dict_to_pd(vardic[scen][elem]) for scen in scenpresent]
        else:
            #assemble different scenarios and subtract ref scenario
            if keytype=='tuple':
                pdat=pd.Series(vardic[scen][elem])
            elif keytype=='nested':
                pdat=dict_to_pd(vardic[scen][elem])
            #if refscen[scen] in vardic.keys():
            if keytype=='tuple':
                pdatref=pd.Series(vardic[refscen[scen]][elem])
            elif keytype=='nested':
                pdatref=dict_to_pd(vardic[refscen[scen]][elem])
            # else:
            #     pdatref=pd.DataFrame()
            frames=[diff_panda(pdat,pdatref) for scen in scenpresent]
        
        #discard scenarios that do not have the same amount of index levels (e.g. different model options - avoid concatenation problems) - WARNING- they will not appear in the final result
        nlevels=[pdat.index.nlevels for pdat in frames] #number of levels in each dataframe
        nlevelmax=max(set(nlevels), key = nlevels.count) #most occuring number of levels
        framestokeep=[frames[k] for k in range(len(frames)) if nlevels[k]==nlevelmax] #only keep frames with most common number of levels
        scentokeep=[scenpresent[k] for k in range(len(frames)) if nlevels[k]==nlevelmax] #only keep scenarios which frame has the most common number of levels
        #concatane to single dataframe
        if refscen==0:
            pdata=pd.concat(framestokeep,keys=[scen for scen in scentokeep])
        else:
            pdata=pd.concat(framestokeep,keys=[scen+'_'+refscen[scen] 
                                               #if refscen[scen]==refscen[scen] else scen 
                                               for scen in scentokeep])
        if keytype=='tuple':
            pdata=pdata.to_frame() 
        #elif keytype=='nested':
         #   pdata=pd.concat(frames,keys=[scen for scen in scenarios if elem in vardic[scen].keys()])
        
        if not pdata.empty: #don t export e
            #rename columns and indexes
            if indexname != 0:
                indexname[elem].insert(0,'scenario')
                pdata.index.names=indexname[elem]
                if renamecol==1:
                    pdata.columns=[elem]
        elif pdata.empty and keytype=='nested': #when crop market or power market off
            if indexname != 0:
                indexname[elem].insert(0,'scenario')
                for k in range(len(indexname[elem])):
                    pdata.insert(k,indexname[elem][k],0)
        #export to csv
        filename=str(elem)+'.csv'
        pd_to_csv(outpath,filename,pdata=pdata)
    
    #Create empty datasets for not exported Decision variables
    if indexname != 0:
        for elem in indexname.keys():
            if elem not in elist:
                indexname[elem].insert(0,'scenario')
                pdata=pd.DataFrame(columns=indexname[elem]+[elem])
                filename=str(elem)+'.csv'
                pdata.to_csv(os.path.join(outpath,filename),sep=CSVSEPARATOR,decimal=CSVDECIMAL,index=False)
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
DVpath=os.path.join(result_path,'IFPRI_IDX')

#Varaibles index names  #ADD: Variables can have different indexs (eg CULAREA)
VarIndex={
 'Hydropower_prod_GWh':['nyear', 'ntime', 'ncountry','nhpp'],
 'Hydropower_val_Md':['nyear', 'ntime', 'ncountry','nhpp'],
 'Power_price_dpkWh':['nyear', 'ntime', 'ncountry'],
 'Power_newprod_GWh':['nyear', 'ntime', 'ncountry','nptech'],
 'Crop_area_kha':['nyear', 'ncountry', 'ncrop', 'ntype'],
 'Crop_prod_kt':['nyear', 'ncountry', 'ncrop', 'ntype'],
 'Crop_val_Md':['nyear', 'ncountry', 'ncrop', 'ntype'],
 'Crop_price_dpt':['nyear', 'ncountry', 'ncrop'],
 'Cul_Dem_mm':['nyear', 'ncountry', 'nculture'],
 'Cul_Rain_mm':['nyear', 'ncountry', 'nculture'],
 'Runoff_Mm3':['nyear', 'ntime'],
 'ET_mm':['nyear', 'ntime'],
 'P_mm':['nyear', 'ntime'],
 'Cul_Cons_Mm3':['nyear', 'ntime', 'ncountry', 'ncrop']
 }

#load all decision variables (DV) from model runs
scen_to_load=set([s for s in refscen.keys()]+[s for s in refscen.values() if s==s]) 
ScenarioDV={scen:pickle.load(open(os.path.join(result_path,scen+'_ifpri_IDX.txt'),"rb")) for scen in scen_to_load}
#Assemble and Export to csv
REF = 0 if DIFFMODE == 0 else refscen
aggregate_scenarios_to_csv(scenarios,ScenarioDV,DVpath,keytype='tuple',indexname=VarIndex,elist=0,refscen=REF,renamecol=1)

        
