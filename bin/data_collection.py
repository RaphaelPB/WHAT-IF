# -*- coding: utf-8 -*-
"""
    Copyright 2018 Technical University of Denmark
    Copyright 2018 COWI A/S

    This file is part of WHAT-IF.

    WHAT-IF is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 (GPLv3) as published by the Free Software Foundation,

    WHAT-IF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with WHAT-IF. If not, see http://www.gnu.org/licenses/.

"""

#Import libraries
import os
import pandas as pd
import pickle
import locale
#set decimal and separator (for handling of csv files)
langlocale = locale.getdefaultlocale()[0]
locale.setlocale(locale.LC_ALL, langlocale)
CSVDECIMAL = ','#locale.localeconv()['decimal_point']
CSVSEPARATOR = ','
if CSVDECIMAL == ',':
    CSVSEPARATOR = ';'
    
class Database():   
            
    def __init__(self,update=0,DataFile=''):    
        if update ==1:
            Data=pickle.load(open(DataFile,"rb"))
            self.val=Data.val
            self.info=Data.info
            self.memo=Data.memo
            self.scen=Data.scen
            self.csv={}
            self.update=1
        else:    
            self.val={}
            self.info={}
            self.memo={}
            self.scen={}
            self.csv={}
            self.update=0            
        
    def harvest_sheet(self,Path,SheetName,Header=5,Index=1,DataType=0,
                      Scenario=None,MultiIndexName=None,MatrixName=None,
                      ColIndexName=None,OnlyCols=None):
        #Harvests one excel sheet, starting at row "Header"
        if Scenario != None and Scenario == Scenario:
            Index=Index+1
        Data = pd.read_excel(Path,sheet_name=SheetName, skiprows=Header, 
                             index_col=list(range(Index)), usecols=OnlyCols, engine='xlrd')
        if DataType==0: #Info data
            #Data.fillna('nan')
            for k in Data.axes[1]:
                self.info[k]=Data.to_dict()[k]
            for k in range(Index):
                self.info[Data.index.names[k]]=Data.index.values                        
        else: #Collect indexes
            self.collect_index_data(Data,Index,SheetName,MultiIndexName)                                        
        if DataType==1: #Classic column data
            self.collect_column_data(Data,Scenario)
        if DataType==2: #Matrix data
           self.collect_matrix_data(Data,Scenario,MatrixName,Index,ColIndexName,SheetName)
    
    def collect_index_data(self,Data,Index,SheetName,MultiIndexName):
        if Index==1:
            self._check_coherence(Data.index.name,Data.index.values,SheetName)
            self.val[Data.index.name]=Data.index.values
            self.memo[Data.index.name]=SheetName #remembers sheetname for check_coherence
        else:                
            for k in range(Index):                
                self._check_coherence(Data.index.names[k],Data.index.levels[k].values,SheetName)
                self.val[Data.index.names[k]]=Data.index.levels[k].values
                self.memo[Data.index.names[k]]=SheetName #remembers sheetname for check_coherence
            if MultiIndexName != None and MultiIndexName != 'nan':
                self.val[MultiIndexName]=Data.index.values
                    
    def collect_column_data(self,Data,Scenario,onlyscen=0):
        if Scenario != None and Scenario == Scenario: #if data has scenario, scenario is as separate index
            for k in Data.axes[1]:
                self.val[k]={scen:Data.loc[scen].to_dict()[k] 
                             for scen in Data.index.levels[0].values if onlyscen==0 or scen==onlyscen}
                self.scen[k]=Scenario
        else:            
            for k in Data.axes[1]:
                self.val[k]=Data.to_dict()[k]
                    
    def collect_matrix_data(self,Data,Scenario,MatrixName,Index,ColIndexName,SheetName,onlyscen=0):
    #Collect matrix data
        #drop "wrong" columns 
        todrop=[col for col in Data.columns if type(col) is str and 'Unnamed:' in col]
        Data=Data.drop(todrop,axis=1)
        if Scenario != None and Scenario == Scenario:
            self.val[MatrixName]={scen:Data.xs(scen,level=0).stack(dropna=False).to_dict() 
                                  for scen in Data.index.levels[0].values 
                                  if onlyscen==0 or scen==onlyscen}
            self.scen[MatrixName]=Scenario
        else:
            self.val[MatrixName]=Data.stack(dropna=False).to_dict()

        #Collect matrix column index
        if ColIndexName != None and ColIndexName != 'nan':
            self.val[ColIndexName]=[k for k in Data.axes[1] if type(k) is not str or 'Unnamed:' not in k]
            self.memo[ColIndexName]=SheetName
            
    def harvest_all(self,Paths):
        #Harvests all sheets in all files provided in Paths, files MUST have an "Info" sheet    
        for path in Paths:
            self.harvest_sheet(path,'Info')            
            for k in self.info['nsheet']:
                print(self.info['SheetName'][k])
                if self.info['OnlyCols'][k] != self.info['OnlyCols'][k]: #sign of no value provided
                    self.info['OnlyCols'][k]=None
                self.info['Index'][k]=int(self.info['Index'][k]) #make sur int options are int
                self.info['Header'][k]=int(self.info['Header'][k])
                self.info['DataType'][k]=int(self.info['DataType'][k])
                if '.csv' in self.info['SheetName'][k]: #Data is provided as csv file and not in excel
                    self.csv[self.info['SheetName'][k]]={ #sheetname is csv sheet, all other parameters work the same as in excel files
                            'Index':self.info['Index'][k],
                            'Update':self.info['Update'][k],
                            'DataType':self.info['DataType'][k],
                            'Scenario':self.info['Scenario'][k],
                            'MatrixName':self.info['MatrixName'][k],
                            'ColIndexName':self.info['ColIndexName'][k],
                            'MultiIndexName':self.info['MultiIndexName'][k]}
                    warning='LOAD CSV: Data will be loaded from '+self.info['SheetName'][k]
                    print(warning)
                else:
                    if self.update==0 or (self.update==1 and self.info['Update'][k]==1):
                        self.harvest_sheet(path,self.info['SheetName'][k], 
                                           Header=self.info['Header'][k], 
                                           Index=self.info['Index'][k],
                                           DataType=self.info['DataType'][k], 
                                           Scenario=self.info['Scenario'][k], 
                                           MultiIndexName=self.info['MultiIndexName'][k],
                                           MatrixName=self.info['MatrixName'][k], 
                                           ColIndexName=self.info['ColIndexName'][k], 
                                           OnlyCols=self.info['OnlyCols'][k])
                    if self.update==1 and self.info['Update'][k]==0:
                        warning= 'NOT UPDATED: Sheet "' + self.info['SheetName'][k] +'" was not updated'
                        print(warning)
                    
    def _check_coherence(self,index,newvalue,currentsheet):
        #Verifies if an index is already defined, if yes it verifies that it has the same value, indicating for coherence in the data structure
        if index in self.val and index==index and index!= None:
            if not set(self.val[index])==set(newvalue):
                if index in self.memo.keys():
                    warning= 'WARNING: Index "' + index + '" in "'+ currentsheet + '" is not coherent with "' + index +'" in "' + self.memo[index] +'"'
                else:
                    warning= 'WARNING: Index "' + index + '" in "'+ currentsheet + '" is not coherent with existing value'
                print(warning)
                
    def save(self,DataFile):
        pickle.dump(self,open(DataFile,"wb")) #python36
    
    def load_hydrology(self,datafolder,paths={},scenarios=[]):
        #this is a specific version of load_csv ... to be removed in the future ?
        if paths == {}:
            #default path names
            paths={'wRunOff':'wRunOff.csv', 'wRainFall':'wRainFall.csv', 'wET0':'wET0.csv'}
        for param in paths.keys(): #iterate through parameters
            if param not in self.val.keys():
                self.val[param]={} #initiate parameter if not existing
            #load file #WARNING sep and decimal might vary depending on how CSV file was created
            data = pd.read_csv(os.path.join(datafolder,paths[param]),
                               sep=CSVSEPARATOR,decimal=CSVDECIMAL,index_col=[0,1])
            #load scenarios
            for scen in data.index.levels[0]:
                if scenarios==[] or scen in scenarios: #only load selected scenarios
                    self.val[param][scen]={(t,c):data.loc[(scen,t),c] for t in data.loc[scen].index for c in data.loc[scen].columns}
            #verify if all requested scenarios where found
            for scen in scenarios:
                if scen not in self.val[param].keys():
                    print('WARNING: scenario '+scen+' of parameter '+param+' not found in '+paths[param])
    
    def load_csv(self,datafolder,nscenario):        
        for pfile in self.csv.keys(): #iterate through parameters
            opt=self.csv[pfile] #reading options of the csv file
            onlyscen=0
            iindex=opt['Index']
            if opt['Scenario'] != None and opt['Scenario']==opt['Scenario']:
                #onlyscen=self.val[opt['Scenario']][nscenario] #read only that scenario
                iindex+=1 #add scenario index
            data = pd.read_csv(os.path.join(datafolder,pfile),sep=CSVSEPARATOR,decimal=CSVDECIMAL,
                               index_col=list(range(iindex)))
            self.collect_index_data(data,opt['Index'],pfile,opt['MultiIndexName'])                                       
            if opt['DataType']==1: #Classic column data
                self.collect_column_data(data,opt['Scenario'],onlyscen=onlyscen)
            if opt['DataType']==2: #Matrix data
                self.collect_matrix_data(data,opt['Scenario'],opt['MatrixName'],iindex,opt['ColIndexName'],pfile,onlyscen=onlyscen)
        
    def read_param(self,ParamName,option=1,index=0,time='all',scenario='',directparam=0):
        #hard coded growing parameters
        GrowingParameters=['wUserDem',
                           'eEngyDem','eFuelCost','eCO2Val','eCAPEX', 'eOppCost', 'eOppCap',
                           'aCropDem','aCropVal','aTransCost','aFarmVal',
                           'aCulYield','aLandCap','aCulMax']
        if directparam!=0 and ParamName in directparam.keys(): #Assign parameter directly and not from parameters object 
            return directparam[ParamName]
        if ParamName not in self.val.keys() or option == 0: #parameter is not activated          
            if index==1:
                return [] #returns empty set
            return {} #returns empty parameter
        if ParamName in GrowingParameters:
            Param = self.grow_param(ParamName,scenario=scenario)
        else:
            if ParamName in self.scen.keys():
                datascen=self.val[self.scen[ParamName]][scenario]
                Param = self.val[ParamName][datascen]
            else:
                Param = self.val[ParamName]
        if time != 'all': #Only selected time steps 
            Param={k:Param[k] for k in Param.keys() if ((type(k)==tuple and k[0] in time) or (type(k)!=tuple and k in time))}
        return Param
    
    def grow_param(self,ParamName,scenario=''):
        #linear growth functiom
        def LinearGrowAtoB(ValuesA,ValuesB,yearA,yearB,yyears):
            if type(list(ValuesA.keys())[0]) is tuple: #keys need to be formulated differently if they are uni or multidimensional
                return {(y,)+keys:ValuesA[keys] + (ValuesB[keys]-ValuesA[keys])/(yearB-yearA)*(y-yearA) 
                        for y in yyears for keys in ValuesA.keys()}
            else:
                return {(y,keys): ValuesA[keys] + (ValuesB[keys]-ValuesA[keys])/(yearB-yearA)*(y-yearA) 
                        for y in yyears for keys in ValuesA.keys()}
        #Generate parameters
        optionON = self.val['Options']['Demand growth',self.val['sOptions'][scenario]]
        years=self.val['year']
        scen=self.val[self.scen[ParamName]][scenario] if ParamName in self.scen.keys() else ''
        if optionON==1 and scen in self.read_param('ngrowthscen'): #Linear interpolation from values A to B to C (optional)
            ValuesA=self.val[ParamName][self.read_param('StatusA')[scen]]
            ValuesB=self.val[ParamName][self.read_param('StatusB')[scen]]                                       
            yearA=int(self.read_param('yearA')[scen])
            yearB=int(self.read_param('yearB')[scen])     
            DicAB=LinearGrowAtoB(ValuesA,ValuesB,yearA,yearB,years)
            if self.read_param('StatusC')[scen] == self.read_param('StatusC')[scen]: #Value is not a NAN
                ValuesC=self.val[ParamName][self.read_param('StatusC')[scen]]
                yearC=int(self.read_param('yearC')[scen])
                yyears=[y for y in years if y>=yearB]
                DicBC=LinearGrowAtoB(ValuesB,ValuesC,yearB,yearC,yyears)
                DicAB.update(DicBC)                
        else: #Constant parameter through years 
            if ParamName in self.scen.keys():                    
                DicAB=LinearGrowAtoB(self.val[ParamName][scen],self.val[ParamName][scen],1,2,years) #1,2 are dummy values
            else:
                DicAB=LinearGrowAtoB(self.val[ParamName],self.val[ParamName],1,2,years)
        return DicAB
    
    def read_time(self,scenario):
        #Reads time parameters from options
        stime=self.val['sOptions'][scenario]
        tini=self.val['Options']['tini',stime]
        tfin=self.val['Options']['tfin',stime]
        t0,m0,y0=self.val['Options']['TimeMonthYear',stime].split('#')
        mo=self.val['monthorder']
        moi0=min(k for k in mo.keys()) #first index of ordered months (avoid i0 =1 or 0 problems)
        monthorder=[mo[k+moi0] for k in range(len(mo))] #months as ordered list
        nm0=[k for k in range(len(monthorder)) if monthorder[k]==m0][0] #initial month index
        ttime = [t for t in range(tini,tfin+1)] #considered time steps
        if 'wRunOff' not in self.scen.keys(): #no scenarios on runoff
            self.val['ntime']=[k[0] for k in self.val['wRunOff'].keys()]
        else: #scenarios on runoff
            scen=self.val[self.scen['wRunOff']][scenario] #select scenario
            self.val['ntime']=[k[0] for k in self.val['wRunOff'][scen].keys()]
        
        self.val['t_year'] = {t:int(y0)+(t-int(t0))//len(monthorder) for t in self.val['ntime']} #year of time steps
        self.val['t_month'] = {t:monthorder[(t-int(t0)+nm0)%len(monthorder)] for t in self.val['ntime']} #month of time steps
        t_prev = {t:t-1 for t in ttime} #previous time step       
        year = [self.val['t_year'][t] for t in range(tini,tfin+1,len(monthorder))]  
        self.val['year']=year
        if self.val['Options']['Initial time step',stime]== 0: #Close the loop for cyclic model
            t_prev[tini]=tfin
        return ttime,t_prev,year