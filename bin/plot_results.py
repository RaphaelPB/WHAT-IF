# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 10:48:55 2018

@author: rapp
"""

#
from __future__  import division        
import matplotlib.pyplot as plt
import pickle
import os

#Result file
dirname = os.path.abspath(os.path.dirname(__file__))
ResultFile = os.path.join(dirname, 'SCENARIOS_results.txt')
ParamFile = os.path.join(dirname,'Parameters.txt') #parameters (python dictionnaries) saved as txt
results = pickle.load(open(ResultFile,"rb"))
parameters = pickle.load(open(ParamFile,"rb"))
ref_scen = parameters.val['sRefScen']

#HP
scenarios = ['HPcp0F','HPcp12F','HPF','HPcp50']
Legend = ['0','12','25','50']
#Reference=['Bcp0F','Bcp12F','BF','Bcp50']
Elem = 'EconomicBalance'
HPVal=[results[ss][Elem]['ToTal']-results[ref_scen[ss]][Elem]['ToTal'] for ss in scenarios]

plt.figure(1)
plt.subplot(211)
#fig, ax = plt.subplots()
plt.bar(ind, men_means, width, color='r')