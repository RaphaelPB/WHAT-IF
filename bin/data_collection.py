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
#the automatic set up of decimals seems to not work from servers
langlocale = locale.getdefaultlocale()[0]
locale.setlocale(locale.LC_ALL, langlocale)
CSVDECIMAL = locale.localeconv()['decimal_point']
CSVSEPARATOR = ','
if CSVDECIMAL == ',':
    CSVSEPARATOR = ';'
    
class Database():   
    #Object containing all model parameters, collected from excel and csv files        
    def __init__(self,update=0,DataFile=''):    
        if update ==1:
            #update from an existing data file (and update parameters as indicated in excel info files)
            Data=pickle.load(open(DataFile,"rb"))
            self.val=Data.val #store paraeters
            self.info=Data.info #store info data (how data is shaped)
            self.memo=Data.memo #to check consistency between parameters
            self.scen=Data.scen #store scenarios of data
            self.csv={} #data to collect fro csv files
            self.update=1 #update existing object (1) or create new (0)
        else:
            #collect all parameters from scratch (csv+xlsx files)
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
        #Parameters are collected from excel "info" sheet 
        i0=0 #start index to collect indices
        if Scenario != None and Scenario == Scenario:
            Index=Index+1 #add scenario as index
            i0=1#avoids collect scenario index  
        Data = pd.read_excel(Path,sheet_name=SheetName, skiprows=Header, 
                             index_col=list(range(Index)), usecols=OnlyCols, engine='openpyxl')
        #frame data depending on type (info,columns,or matrix)
        if DataType==0: #Info data
            self.info={}
            for k in Data.axes[1]:
                self.info[k]=Data.to_dict()[k]
            for k in range(Index):
                self.info[Data.index.names[k]]=Data.index.values             
        else: #Collect indexes
            self.collect_index_data(Data,Index,SheetName,MultiIndexName,i0=i0)                                        
        if DataType==1: #Classic column data
            self.collect_column_data(Data,Scenario)
        if DataType==2: #Matrix data
           self.collect_matrix_data(Data,Scenario,MatrixName,Index,ColIndexName,SheetName)
    
    def collect_index_data(self,Data,Index,SheetName,MultiIndexName,i0=0):
        #collect indexes
        if Index==1:
            self._check_coherence(Data.index.name,Data.index.values,SheetName)
            self.val[Data.index.name]=Data.index.values
            self.memo[Data.index.name]=SheetName #remembers sheetname for check_coherence
        else:
            for k in range(i0,Index):                
                self._check_coherence(Data.index.names[k],Data.index.levels[k].values,SheetName)
                self.val[Data.index.names[k]]=Data.index.levels[k].values
                self.memo[Data.index.names[k]]=SheetName #remembers sheetname for check_coherence
            if MultiIndexName != None and MultiIndexName == MultiIndexName:
                self.val[MultiIndexName]=Data.index.values
                    
    def collect_column_data(self,Data,Scenario,onlyscen=0):
        #collect column shaped data (first column(s) are indexes, next are parameters)
        if Scenario != None and Scenario == Scenario: #if data has scenario, scenario is as separate index
            for k in Data.axes[1]:
                tempdic={}
                scens=Data.index.levels[0].values
                for scen in scens:
                    if onlyscen==0 or scen==onlyscen:
                        if '*' in str(scen): 
                            #scenario is an update, name has the form "scenup*scenor"
                            scenup,scenor=scen.split('*')
                            #specific case where the scenario indicator is a number
                            if scenor.isdigit() and scenor not in scens: 
                                scenor=int(scenor)
                            #extract reference scenario
                            val=Data[k].xs(scenor,level=0).to_dict()
                            #update with modifications to orginal scenario
                            val.update(Data[k].xs(scen,level=0).to_dict())
                            tempdic[scenup]=val
                        else:
                            tempdic[str(scen)]=Data[k].xs(scen,level=0).to_dict()
                self.val[k]=tempdic
                self.scen[k]=Scenario
        else:            
            for k in Data.axes[1]:
                self.val[k]=Data.to_dict()[k]
                    
    def collect_matrix_data(self,Data,Scenario,MatrixName,Index,ColIndexName,
                            SheetName,onlyscen=0):
    #Collect matrix shaped data, rows and columns are indexes of a single parameter
        #drop "wrong" columns (happens if there is extra text in the excel sheet) 
        todrop=[col for col in Data.columns if type(col) is str and 'Unnamed:' in col]
        Data=Data.drop(todrop,axis=1)
        if Scenario != None and Scenario == Scenario:
            tempdic={}
            scens=Data.index.levels[0].values
            for scen in scens:
                if onlyscen==0 or scen==onlyscen:
                    if '*' in str(scen): #the scenario is updating another
                        scenup,scenor=scen.split('*') #get original and updated scenario
                        if scenor.isdigit() and scenor not in scens: #specific case where the scenario indicator is a number
                            scenor=int(scenor)
                        #extract reference scenario
                        val=Data.xs(scenor,level=0).stack(dropna=False).to_dict()
                        #update with modifications to orginal scenario
                        val.update(Data.xs(scen,level=0).stack(dropna=False).to_dict())
                        tempdic[scenup]=val
                    else:
                        tempdic[str(scen)]=Data.xs(scen,level=0).stack(dropna=False).to_dict()
            self.val[MatrixName]=tempdic
            self.scen[MatrixName]=Scenario
        else:
            self.val[MatrixName]=Data.stack(dropna=False).to_dict()

        #Collect matrix column index
        if ColIndexName != None and ColIndexName == ColIndexName:
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
                self.info['Header'][k]=int(self.info['Header'][k]) #how many rows to jump to data
                self.info['DataType'][k]=int(self.info['DataType'][k]) #matricx or columns
                if '.csv' in self.info['SheetName'][k]: #Data is provided as csv file and not in excel
                    self.csv[self.info['SheetName'][k]]={key:self.info[key][k] 
                                                         for key in self.info.keys() if key!='nsheet'}
                    if 'Onlyscen' not in self.info.keys():
                        self.csv[self.info['SheetName'][k]]['OnlyScen']=None
                            # 'Index':self.info['Index'][k],
                            # 'Update':self.info['Update'][k],
                            # 'DataType':self.info['DataType'][k],
                            # 'Scenario':self.info['Scenario'][k],
                            # 'OnlyScen':self.info['OnlyScen'][k] if 'OnlyScen' in self.info.keys() else None,
                            # 'MatrixName':self.info['MatrixName'][k],
                            # 'ColIndexName':self.info['ColIndexName'][k],
                            # 'MultiIndexName':self.info['MultiIndexName'][k]}
                    info='LOAD CSV: Data will be loaded from '+self.info['SheetName'][k]
                    print(info)
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
                        info= 'NOT UPDATED: Sheet "' + self.info['SheetName'][k] +'" was not updated'
                        print(info)
                    
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
        #Save parameters to a text file (can be directly called by next run)
        pickle.dump(self,open(DataFile,"wb")) #python36
    
    def load_csv(self,datafolder,nscenario):
        #function called when data is to be load from csv file (when parameter sheet in excel is .csv)
        #works the same as other parameter, can be columns or matrix        
        for pfile in self.csv.keys(): #iterate through parameters
            opt=self.csv[pfile] #reading options of the csv file
            onlyscen=0
            iindex=opt['Index']
            i0=0
            if opt['Scenario'] != None and opt['Scenario']==opt['Scenario']:
                if opt['OnlyScen'] == 1:
                    onlyscen=self.val[opt['Scenario']][nscenario] #read only that scenario
                iindex+=1 #add scenario index
                i0=1
            data = pd.read_csv(os.path.join(datafolder,pfile),sep=CSVSEPARATOR,decimal=CSVDECIMAL,
                               index_col=list(range(iindex)))
            self.collect_index_data(data,opt['Index'],pfile,opt['MultiIndexName'],i0=i0)                                       
            if opt['DataType']==1: #Classic column data
                self.collect_column_data(data,opt['Scenario'],onlyscen=onlyscen)
            if opt['DataType']==2: #Matrix data
                self.collect_matrix_data(data,opt['Scenario'],opt['MatrixName'],
                                         iindex,opt['ColIndexName'],pfile,onlyscen=onlyscen)
        
    def read_param(self,ParamName,option=1,index=0,time='all',scenario='',directparam=0):
        #reads parameter from database object, main feature is to return empty parameter
        #if parameter is not collected, and 
        #hard coded growing parameters
        GrowingParameters=['wUserDem',
                           'eEngyDem','eFuelCost','eCO2Val','eCAPEX', 'eOppCost', 'eOppCap',
                           'aCropDem','aCropVal','aTransCost','aFarmVal',
                           'aCulYield','aLandCap','aCulMax',
                           'iCAPEX']
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
                datascen=str(self.val[self.scen[ParamName]][scenario])
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
            ValuesA=self.val[ParamName][str(self.read_param('StatusA')[scen])]
            ValuesB=self.val[ParamName][str(self.read_param('StatusB')[scen])]                                 
            yearA=int(self.read_param('yearA')[scen])
            yearB=int(self.read_param('yearB')[scen])
            DicAB=LinearGrowAtoB(ValuesA,ValuesB,yearA,yearB,years)
            if self.read_param('StatusC')[scen] == self.read_param('StatusC')[scen]: #Value is not a NAN
                ValuesC=self.val[ParamName][str(self.read_param('StatusC')[scen])]
                yearC=int(self.read_param('yearC')[scen])
                yyears=[y for y in years if y>=yearB]
                DicBC=LinearGrowAtoB(ValuesB,ValuesC,yearB,yearC,yyears)
                DicAB.update(DicBC)              
        else: #Constant parameter through years 
            if ParamName in self.scen.keys():
                values=self.val[ParamName][str(scen)]                           
            else:
                values=self.val[ParamName]
            DicAB=LinearGrowAtoB(values,values,1,2,years) #1,2 are dummy values
        return DicAB
    
    def read_time(self,scenario):
        #Reads time parameters from options
        #ROOM for improvment here - a lot of the complexity comes from the MPC
        #framework that uses the "raw" parameters (e.g. parameters.val[] ...)
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