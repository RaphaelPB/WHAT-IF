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
from result_analysis import ResultAnalysis

#%%OPTIONS - MODIFY BY USER
#SCENFILE=os.path.join(dirname,'Data','Scenarios_to_compare.xlsx')

SHEET='main'
FOLDERNAME='WHATIF_main'
result_path = os.path.join(dirname,'Results',FOLDERNAME) 
SCENFILE='Scenarios_to_compare.xlsx'
#%%#####################################
#        COLLECT INFORMATION      
sceninfo=pd.read_excel(SCENFILE,sheet_name=SHEET, skiprows=1, index_col=[0], engine='xlrd')
refscen=sceninfo.to_dict()['refscen']
scenarios=[s for s in refscen.keys()]
#collect scenarios not in nscen but in refscen
tmpdic={}
for s in refscen.values():
    if s==s and s not in scenarios:
        scenarios.append(s)
        tmpdic[s]=s
refscen.update(tmpdic)

#%%#####################################
#        COMPARE SCENARIOS (to excel)

#load all results
parallelresults=[pickle.load(open(os.path.join(result_path,scen+'.txt'),"rb")) for scen in scenarios]
#Compare and Export results
result_analysis=ResultAnalysis()
result_analysis.export_scenario_analysis(scenarios,refscen,parallelresults,result_path)
