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

class Database():   
            
    def __init__(self,update=0,DataFile=''):    
        if update ==1:
            Data=pickle.load(open(DataFile,"rb"))
            self.val=Data.val
            self.info=Data.info
            self.memo=Data.memo
            self.scen=Data.scen
            self.update=1
        else:    
            self.val={}
            self.info={}
            self.memo={}
            self.scen={}
            self.update=0            
        
    def harvest_sheet(self,Path,SheetName,Header=5,Index=1,DataType=0,
                      Scenario=None,MultiIndexName=None,MatrixName=None,
                      ColIndexName=None,OnlyCols=None):
        #Harvests one excel sheet, starting at row "Header"
        Data = pd.read_excel(Path,sheet_name=SheetName, skiprows=Header, 
                             index_col=list(range(Index)), usecols=OnlyCols, engine='xlrd')
        if DataType==0: #Info data
            #Data.fillna('nan')
            for k in Data.axes[1]:
                self.info[k]=Data.to_dict()[k]
            for k in range(Index):
                self.info[Data.index.names[k]]=Data.index.values                
        
        else: #Collect indexes
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
        
        if DataType==1: #Classic column data
            if Scenario != None and Scenario == Scenario: #if data has scenario, scenario is as separate index
                for k in Data.axes[1]: 
                    self.val[k]={scen:Data.loc[scen].to_dict()[k] for scen in Data.index.levels[0].values}
                    self.scen[k]=Scenario
            else:            
                for k in Data.axes[1]:
                    self.val[k]=Data.to_dict()[k]
        
        if DataType==2: #Matrix data
            #Collect matrix data
            if Scenario != None and Scenario == Scenario:
                self.val[MatrixName]={scen:{} for scen in Data.index.levels[0].values}
                self.scen[MatrixName]=Scenario
            else:
                self.val[MatrixName]={}            
            for i in Data.index.values:
                for j in Data.axes[1]:
                    if type(j) is not str or 'Unnamed:' not in j:
                        if Index==1: #2D matrix SPECIAL CASE NOT NECESSARY
                            self.val[MatrixName][i,j]=Data.to_dict()[j][i]
                        else: #>2D matrix
                            index=list(i)
                            index.append(j)
                            if Scenario != None and Scenario == Scenario: #if data has scenario, scenario is as separate index                            
                                self.val[MatrixName][index[0]][tuple(index[1:])]=Data.to_dict()[j][i]
                            else:    
                                self.val[MatrixName][tuple(index)]=Data.to_dict()[j][i]
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
                if self.info['OnlyCols'][k] != self.info['OnlyCols'][k]:
                    onlycols=None
                else:
                    onlycols=self.info['OnlyCols'][k]
                if self.update==0 or (self.update==1 and self.info['Update'][k]==1):
                    self.harvest_sheet(path,self.info['SheetName'][k], Header=int(self.info['Header'][k]), Index=int(self.info['Index'][k]),
                                       DataType=int(self.info['DataType'][k]), Scenario=self.info['Scenario'][k], MultiIndexName=self.info['MultiIndexName'][k],
                                       MatrixName=self.info['MatrixName'][k], ColIndexName=self.info['ColIndexName'][k], OnlyCols=onlycols)
                if self.update==1 and self.info['Update'][k]==0:
                    warning= 'NOT UPDATED: Sheet "' + self.info['SheetName'][k] +'" was not updated'
                    print(warning)
                    
    def _check_coherence(self,index,newvalue,currentsheet):
        #Verifies if an index is already defined, if yes it verifies that it has the same value, indicating for coherence in the data structure
        if index in self.val:
            if not set(self.val[index])==set(newvalue):
                warning= 'WARNING: Index "' + index + '" in "'+ currentsheet + '" is not coherent with "' + index +'" in "' + self.memo[index] +'"'
                print(warning)
                
    def save(self,DataFile):
        pickle.dump(self,open(DataFile,"wb")) #python36
    
    def load_hydrology(self,datafolder,paths={},scenarios=[]):
        if paths == {}:
            #default path names
            paths={'wRunOff':'wRunOff.csv', 'wRainFall':'wRainFall.csv', 'wET0':'wET0.csv'}
        for param in paths.keys(): #iterate through parameters
            if param not in self.val.keys():
                self.val[param]={} #initiate parameter if not existing
            #load file #WARNING sep and decimal might vary depending on how CSV file was created
            data = pd.read_csv(os.path.join(datafolder,paths[param]),sep=';',decimal=',',index_col=[0,1])
            #load scenarios
            for scen in data.index.levels[0]:
                if scenarios==[] or scen in scenarios: #only load selected scenarios
                    self.val[param][scen]={(t,c):data.loc[(scen,t),c] for t in data.loc[scen].index for c in data.loc[scen].columns}
            #verify if all requested scenarios where found
            for scen in scenarios:
                if scen not in self.val[param].keys():
                    print('WARNING: scenario '+scen+' of parameter '+param+' not found in '+paths[param])
    
        
        