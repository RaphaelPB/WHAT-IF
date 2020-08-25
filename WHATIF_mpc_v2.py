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
#Import python libraries: pandas, xlswriter, openpyxl, xlrd, pyomo, glpk 
import time
import os
import pickle
import sys
from pyomo.environ                import Suffix, SolverFactory, Var      #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
tmpbin = os.getenv('TMPBIN')
sys.path.append(os.path.join(dirname, 'bin'))
from data_collection              import Database                   #Stores all parameters from excel sheets in a python object
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Exports results to excel sheets
from SDP_nile_water_val import totalcostfunction, load_whatif_data, SDP_backwards, SDP_forwards
#For Model Predictive Control framework
#from joblib import Parallel, delayed
import multiprocessing as mp
import copy
import numpy as np
from pyomo.opt import TerminationCondition

#For block model
from pyomo.environ import ConcreteModel, Set, Objective, Constraint, Block

#%% DEFINE PARAMETERS
#OPTIONS
UPDATE = 0 #0 updates all parameters, 1 updates only selected parameters
PARALLEL_ensemble = 0 #Run ensemble forecasts in parallel
epll=20 #Maximum number of parallel runs for ensemble forecast (usefull if also using a parallel solver)
PARALLEL_scenario = 1 #Run scenarios in parallel (BOTH IS NOT POSSIBLE)
RESULTFOLDER = 'investX_ph24_mc3c50_nx2'
SCENITER=1 #0 normal, 1 generates scenarios based on keywords and keeps only key results, 2 same but keeps all results
#PH_all = 0 #Prediction horizon of MPC is until end of planning period
SDP_export = 0 #export water values for SDP
NRUNS=10 #number of years in the SDP backward runs
ICLASS = [(0,30),(30,60),(60,100)] #classes for SDP
FIXED_YEAR = 0 #MPC model always consider full years
NOFIXLASTY = 35
#%% knn_bootstrap KNN HYDROLOGY GENERATOR
def knn_bootstrap(parameters,TIME,k=20,nes=1,flen=23,ntime=0,climscen='base',weighted=0):
    #t=time.time()
    #parameters = Model parameters
    #TIME - ts = time step to start the forecast from
    #k = Number of nearest neighbors considered for resampling/weighted average
    #nes = Number of ensemble members desired
    #flen = Length of the desired forecast
    #ntime = time steps considered for knn_bootstrap analysis (0=use all time steps)
    #cs = climatic scenario of the hydrology data
    #weighted: 1=Instead of sampling uses all nearest neighbors, and saves also the probability/Kernel
    cT=0
    ts=TIME['t']
    if weighted==1:
        k=nes
    
    #Hydrology data
    Qo = parameters.val['wRunOff'][climscen]
    Po = parameters.val['wRainFall'][climscen]
    Eo = parameters.val['wET0'][climscen]
    #
    #climscen='nile'
    Q = parameters.val['wRunOff'][climscen]
    P = parameters.val['wRainFall'][climscen]
    E = parameters.val['wET0'][climscen]
    #Indexes
    def looptime(t):
        if t<=ntime[-1]:
            return t
        else:
            return cT+t-ntime[-1]    
    if ntime == 0:        
        #tfin = max(parameters.val['ntime'])
        ntime= np.array([t for t in parameters.val['ntime'] if t != ts and cT<t])
        #if PH_all == 1:
            #ntime = np.array([t for t in parameters.val['ntime'] if t != ts and cT<t<tfin])
    ncatch  = parameters.val['ncatch']    
    ncatchnile=['Upstream']
    #Create feature vector and observed pattern
    Di=np.array([sum(Qo[ts,c] for c in ncatchnile),
                 sum(Po[ts,c] for c in ncatchnile),
                 sum(Eo[ts,c] for c in ncatchnile)])
    d = len(Di)
    Dt = {t:[sum(Q[t,c] for c in ncatchnile),
             sum(P[t,c] for c in ncatchnile),
             sum(E[t,c] for c in ncatchnile)] 
         for t in ntime if (t-ts)%12==0}
    
    #Calculate weights and distance
    stdQ = np.std([sum(Q[t,c] for c in ncatchnile) for t in ntime if (t-ts)%12==0])
    #stdP = np.std([sum(P[t,c] for c in ncatch) for t in ntime if (t-ts)%12==0])
    #stdE = np.std([sum(E[t,c] for c in ncatch) for t in ntime if (t-ts)%12==0])
    
    w=[1/stdQ,0,0] #weights 
    mtime=np.array([t for t in ntime if (t-ts)%12==0])
    rt=np.array([sum(w[j]*(Di[j]-Dt[t][j])**2 for j in range(d))**0.5 for t in mtime]) #distances of the neighbors
    mtime_t=[mtime[k] for k in range(len(mtime)) if Dt[mtime[k]][0]==0 or rt[k]!=0] #Dt[mtime[k]][0]==0 specific case where total flow = 0 , otherwise remove as probably a cheat
    rt_t=[rt[k] for k in range(len(rt)) if Dt[mtime[k]][0]==0 or rt[k]!=0]
    mtime=np.array(mtime_t)
    rt=np.array(rt_t)
#    if len(mtime_t)<20:
#        print('time step: '+ str(ts))
#        print(mtime_t)
#        print(rt_t)
    #print('mtime lenght: '+str(len(mtime)) + ' mtime_t lenght: '+str(len(mtime_t)))
    #print('rt lenght: '+str(len(rt)) + ' mtime_t lenght: '+str(len(rt_t)))
    #Sort neighbours - and select k nearest
    time_sort = mtime[rt.argsort()] #sorted neighbors
    #rt_sort = np.sort(rt) #sorted weights
    
    #Resampling kernel 
    K = [1/l/sum(1/ll for ll in range(1,k+1)) for l in range(1,k+1)] #Sampled by ranking    
    #Percentiles method    
    #K = [1/rt_sort[l]/sum(1/rt_sort[ll] for ll in range(k)) for l in range(k)] #Sampled by the distance
    TQ = [sum(Q[looptime(tt),c] for tt in range(t,t+flen) for c in ncatchnile) for t in time_sort]
    #Ensemble forecastst
    forecast={es:{} for es in range(nes)}
    
    #Compute weighted average
    if nes==1: #Only average forecast
        forecast[0]['Q_ef'] = {(ts+t,c): sum(K[l]*Q[looptime(time_sort[l]+t),c] for l in range(k)) for t in range(1,flen+1) for c in ncatch}
        forecast[0]['P_ef'] = {(ts+t,c): sum(K[l]*P[looptime(time_sort[l]+t),c] for l in range(k)) for t in range(1,flen+1) for c in ncatch}
        forecast[0]['E_ef'] = {(ts+t,c): sum(K[l]*E[looptime(time_sort[l]+t),c] for l in range(k)) for t in range(1,flen+1) for c in ncatch}
        
    #Real Ensemble forecast
    else:
        NEWM=1
        fskill=5
        if weighted == '3class':
            ccuts=TIME['CCUTS']
            def cbounds(cl):
                lb=np.percentile(TQ,ccuts[cl-1]) if cl != 0 else -0.001
                up=np.percentile(TQ,ccuts[cl])
                return [lb,up]
            forecast={cl:{} for cl in range(len(ccuts))}
            rclass=[[l for l in range(k) if cbounds(cl)[0]<TQ[l]<=cbounds(cl)[1]] for cl in range(len(ccuts))]
            for cl in range(len(ccuts)):
                if NEWM != 1:
                    forecast[cl]['Q_ef'] = {(ts+t,c): sum(K[l]/sum(K[l] for l in rclass[cl])*Q[looptime(time_sort[l]+t),c] for l in rclass[cl]) for t in range(1,flen+1) for c in ncatch}
                    forecast[cl]['P_ef'] = {(ts+t,c): sum(K[l]/sum(K[l] for l in rclass[cl])*P[looptime(time_sort[l]+t),c] for l in rclass[cl]) for t in range(1,flen+1) for c in ncatch}
                    forecast[cl]['E_ef'] = {(ts+t,c): sum(K[l]/sum(K[l] for l in rclass[cl])*E[looptime(time_sort[l]+t),c] for l in rclass[cl]) for t in range(1,flen+1) for c in ncatch}
                else:
                    forecast[cl]['Q_ef'] = {(ts+t,c): sum(K[l]*Q[looptime(time_sort[l]+t),c] for l in range(k)) for t in range(1,fskill+1) for c in ncatch}
                    forecast[cl]['Q_ef'].update({(ts+t,c): sum(1/len(rclass[cl])*Q[looptime(time_sort[l]+t),c] for l in rclass[cl]) for t in range(fskill+1,flen+1) for c in ncatch})
                    forecast[cl]['P_ef'] = {(ts+t,c): sum(K[l]*P[looptime(time_sort[l]+t),c] for l in range(k)) for t in range(1,fskill+1) for c in ncatch}
                    forecast[cl]['P_ef'].update({(ts+t,c): sum(1/len(rclass[cl])*P[looptime(time_sort[l]+t),c] for l in rclass[cl]) for t in range(fskill+1,flen+1) for c in ncatch})
                    forecast[cl]['E_ef'] = {(ts+t,c): sum(K[l]*E[looptime(time_sort[l]+t),c] for l in range(k)) for t in range(1,fskill+1) for c in ncatch}
                    forecast[cl]['E_ef'].update({(ts+t,c): sum(1/len(rclass[cl])*E[looptime(time_sort[l]+t),c] for l in rclass[cl]) for t in range(fskill+1,flen+1) for c in ncatch})
                    
            if NEWM != 1:
                K=[sum(K[l] for l in rclass[cl]) for cl in range(len(ccuts))]
            else:
                K=[1/len(ccuts) for cl in range(len(ccuts))]
        else:
            for es in range(nes):
                #sample neighbor(=time step) with probability K
                if weighted == 0:
                    st = np.random.choice([time_sort[l] for l in range(k)],1,K)[0]
                else: #Use all neighbors (not a Sample)
                    st = time_sort[es]
                #historical time series of sampled neighbor
                forecast[es]['Q_ef']={(ts+t,c): Q[looptime(st+t),c] for t in range(1,flen+1) for c in ncatch} 
                forecast[es]['P_ef']={(ts+t,c): P[looptime(st+t),c] for t in range(1,flen+1) for c in ncatch}
                forecast[es]['E_ef']={(ts+t,c): E[looptime(st+t),c] for t in range(1,flen+1) for c in ncatch}
            
    #print('bootsrap time',time.time()-t)
    if weighted in [1,'3class']:
        forecast['Kernel']=K
    return forecast
   
#%% FUNCTIONS FOR MODEL PREDICTIVE CONTROL
    
def updatehydrology(TIME,parameters,forecast): #UPDATES HYDROLOGY 
    cT=0
    #creates a dictionary with hydrology parameters according to the forecast and MPC parameters
    #the output can then be passed as directparam=output to the hydroeconomic model
    #forecast is one of the ensemble members from the output of the knn_bootstrapping function
    pftime = TIME['t']+TIME['PFtime']-1 #perfect forecast time (absolute time) - pftime is to simulate a better forecast
    aftime = TIME['AFtime'] #average forecast time (relative) - past pftime+aftime, climatology is used instead of forecast
    ptime  = range(TIME['tini'],TIME['tfin']+1) #prediction horizon time steps to update
    ntime  = parameters.val['ntime'] #all time steps in the time series
    ncatch = parameters.val['ncatch'] #catchments
    if max(ptime) > max(ntime): #maximum time step of data timeseries
        return print('ERROR: trying to access time step further than available time series')
    climscen= parameters.val['sClimate'][TIME['ss']]
    #observed hydrology
    Qo = {(t,c):parameters.val['wRunOff'][climscen][t,c] for t in ptime for c in ncatch if t <= pftime}
    Po = {(t,c):parameters.val['wRainFall'][climscen][t,c] for t in ptime for c in ncatch if t <= pftime}
    Eo = {(t,c):parameters.val['wET0'][climscen][t,c] for t in ptime for c in ncatch if t <= pftime}
    #forecasted hydrology
    Qf = {(t,c):forecast['Q_ef'][t,c] for t in ptime for c in ncatch if pftime<t<=pftime+aftime}    
    Pf = {(t,c):forecast['P_ef'][t,c] for t in ptime for c in ncatch if pftime<t<=pftime+aftime}    
    Ef = {(t,c):forecast['E_ef'][t,c] for t in ptime for c in ncatch if pftime<t<=pftime+aftime}
    #average hydrology (climatology)
    Qa = {(t,c):sum(parameters.val['wRunOff'][climscen][tt,c] for tt in ntime if (tt-t)%12==0 and tt>cT)/(len(ntime)-cT)*12 for t in ptime for c in ncatch if t>pftime+aftime}    
    Pa = {(t,c):sum(parameters.val['wRainFall'][climscen][tt,c] for tt in ntime if (tt-t)%12==0 and tt>cT)/(len(ntime)-cT)*12 for t in ptime for c in ncatch if t>pftime+aftime}    
    Ea = {(t,c):sum(parameters.val['wET0'][climscen][tt,c] for tt in ntime if (tt-t)%12==0 and tt>cT)/(len(ntime)-cT)*12 for t in ptime for c in ncatch if t>pftime+aftime}
    #assamble
    paramup={}
    paramup['wRunOff'] = {**Qo,**Qf,**Qa}
    paramup['wRainFall'] = {**Po,**Pf,**Pa} #
    paramup['wET0'] = {**Eo,**Ef,**Ea} #
    return paramup
                
def save_vars(model,varlist,time=[],sw=0,se=0,block=False):
    #Saves variables from a model towards a dictionary
    ##model: model to save variables from
    ##varlist: list of variables to save
    ##time: only saves variables with index in these time steps
    ##sw, se: specifies reference (=original) scaling factor used in the model, 
    #as safe_solve() might have modified scaling factors in model
    ##block : True looks into blocks to find variables
    vardic={}
    for varname in varlist:
        if block==False:
            var=model.find_component(varname)
        else: #look up in arbitraly first ensemble model
            var=model.find_component('ensemble[0].'+varname)
        if var!=None:
            CF=1
            if varname in ['WwRSTORAGE'] and sw != 0:
                CF=model.sw.value/sw if block==False else model.ensemble[0].sw.value/sw
            if varname in ['EeGENCAP'] and se != 0:
                CF=model.se.value/se if block==False else model.ensemble[0].se.value/se
            if time == []:
                vardic[varname]={index:var[index].value*CF for index in var}
            else:
                vardic[varname]={index:var[index].value*CF for index in var if index[0] in time}
    return vardic

def average_vars(diclist,weights=0):
    #From a list of dictionaries containing the ensemble model variables
    #this function returns one dictionary with the average variables values
    #weights : weights the individual ensemble models (are normalized if not already the case)
    if len(diclist)==1: #only one model run, directly keep values
        return diclist[0]
    if weights==0: #creates average weights
        weights=[1/len(diclist) for s in range(len(diclist))] 
    if sum(weights) != 1: #normalizes the weights
        weights=[k/sum(weights) for k in weights]
    avvardic={}
    for varname in diclist[0].keys():
        avvardic[varname]={key:sum(diclist[s][varname][key]*weights[s] for s in range(len(diclist)))
                           for key in diclist[0][varname].keys()}
    return avvardic

def load_vars(model,vardic,time=[],fixed=0,cularea=[],gencap=[],onlyvar=[],reltime=0,block=None):
    #load_var loads (previously saved) variables into a model
    #model : model to load variables in
    #vardic : dictionary with keys the names of the decision variables
    #time : loads only the variables with the index within this time range
    #fixed : 1 fixes variables
    #cularea : [years] scpecial transformation from DV to parameter
    #gencap : [year] scpecial transformation towards main model
    #onlyvar : [varnames] only loads specified variables
    #reltime : n switches the time index by n (switch one month, n=1,switch a year n=12)
    #block: if number, the model to load is a block model, number specifies index of ensemble block
    varlist = [v for v in vardic.keys() if v in onlyvar] if onlyvar != [] else vardic.keys()
    for varname in varlist:
        if (cularea==[] or varname not in ['AlCULAREA']) and (gencap==[] or varname not in ['EeGENCAP']):
            if block == None:
                var=model.find_component(varname)
            else:
                var=model.find_component('ensemble['+str(block)+'].'+varname)
            if var != None:
                for index in vardic[varname].keys():
                    indexr=index
                    if reltime !=0: #Switch the loaded variable time by reltime (Warning cannot do months and years at the same time)
                        indexr=tuple([index[0]+reltime]+[index[k] for k in range(1,len(index))])
                    if indexr in var and (time==[] or indexr[0] in time): #ensures that the index exists in the model variable
                        var[indexr].value=vardic[varname][index] 
                        if fixed == 1: #fix variable
                            var[indexr].fixed=True
        elif cularea!=[] and varname in ['AlCULAREA']: #specific handling of cultivated area
            if block == None:
                var=model.CULAREA
                field_culture=model.field_culture
            else:
                var=model.ensemble[block].CULAREA
                field_culture=model.ensemble[block].field_culture
            for y,fz,cul in var:
                if y in cularea:
                    var[y,fz,cul].value=sum(vardic[varname][index] for index in vardic[varname].keys() if index[0]==y and index[1]==fz and field_culture[index[2],index[3]]==cul)
        elif gencap!=[] and varname in ['EeGENCAP']: #specific handling of generic capacity investments
            for y,pt,pm in vardic[varname].keys():
                if y in gencap:
                    PreviousInvestment=sum(model.EeGENCAP[ky,pt,pm].value for ky in model.nyear if ky < y and ky >= y-model.eLifeTime[pt,pm]) #IMPROVE only works if fixes current year
                    value=max(0,vardic[varname][y,pt,pm]-PreviousInvestment)
                    model.EeGENCAP[y,pt,pm].value=value
                    #model.EeGENCAP[y,pt,pm].fixed=True    
                    
def create_ensemble_block_model(TIME,parameters,ensforecast): 
    #Creates a single block model with all scenarios defined by ensemble forecasts
    #The first time step DV are linked and the objective function aggregated
    #TIME : time and other MPC options
    #parameters : parameters to create models
    #ensforecast : ensemble forecast generated by the knn_bootstrap function
    #ss : scenario name
    #obj_func : specifies if the ensemble block objective function is a weigthed average or a maxmin
    def create_block_ensemble(block,ens):
        paramup=updatehydrology(TIME,parameters,ensforecast[ens])
        block=HydroeconomicOptimization(parameters,scenario=TIME['ss'],directparam=paramup).model
        return block                   
    #Create ensemble model
    model=ConcreteModel()
    model.nens=Set(initialize=range(len(ensforecast.keys())-1),ordered=True)
    model.ensemble=Block(model.nens,rule=create_block_ensemble)
    #LINK ENSEMBLE MODELS VARIABLES TOGETHER:
    #define linking rule
    def link_block_variables(m,ens,nvar,nindex):
        #if ens == ens0: 
        #    return Constraint.Skip()                
        varname1='ensemble['+str(ens)+'].'+nvar
        varname2='ensemble['+str(max(0,ens-1))+'].'+nvar
        index=linkvar[nvar,nindex]
        return m.find_component(varname1)[index]==m.find_component(varname2)[index]
    #define multi index of variables+index that need to be linked to each other (=variable of this time step)
    ens0=min(model.nens) #arbitrary 'first' ensemble model
    linkvar={} #dictionary containing indexes of variables that need to be linked
    for var in model.ensemble[ens0].component_objects(Var):
        indexlist=[index for index in var if index[0]==TIME['t']] #WARNING if other index is same as time
        if var.local_name=='EeGENCAP':
            indexlist=[index for index in var if index[0]==TIME['year']] #Gencap investments need to be the same in each ensemble member
        for k in range(len(indexlist)):
            linkvar[(str(var.local_name),k)]=indexlist[k]    
    model.nlinkvar=Set(dimen=2,initialize=linkvar.keys())
    model.link=Constraint(model.nens,model.nlinkvar,rule=link_block_variables)
    return model
    
#AGGREGATE OBJECTIVE FUNCTIONS
def block_model_obj_function(model,ensforecast,obj_func='weighted'):
    model.ensemble[:].obj.deactivate() #deactivate existing objective functions
    if obj_func=='weighted':
        model.obj = Objective(expr=sum(model.ensemble[ens].obj.expr*ensforecast['Kernel'][ens] for ens in model.nens)) #weighted summed objective function
    elif obj_func=='maxmin':
        model.MAXMIN=Var()
        def max_min_obj(m,ens):
            return m.MAXMIN <= -m.ensemble[ens].obj.expr
        model.maxmin=Constraint(model.nens,rule=max_min_obj)
        model.obj=Objective(expr=-model.MAXMIN)

def generate_boundary_conditions(model,TIME,rstorage,gencap):
    #Fix initial reservoir storage to final storage from previous year
    #storage_value= 'hard' if model.Options['Reservoir Target'] != 'shadow' else 'shadow'
    DONOTFIXEND=False
    for res in model.nres:
        #Final storage target at end prediction horizon (sort of Rule-curve)
        if DONOTFIXEND or TIME['tfin']!=TIME['tfinal']:
            if TIME['storage_value']==0:
                target = model.Options['Reservoir Target']*model.wStorIni[res].value 
            else: #reservoir target is based on water shadowprice (no hard target)
                target = 0
            getattr(model,'wStorFin')[res] = target               
        #Initial storage = to last time step storage
        getattr(model,'wStorIni')[res] = rstorage[res]*model.sw.value               
    
    if TIME['storage_value'] != 0: #reservoir target is based on water shadowprice
        #empirical factor taking into account that shadowprices from Perfect Foresight runs are underestimated
        EF=float(model.Options['Reservoir Target'].split('w')[1]) if model.Options['Reservoir Target'] != 'shadow' else 1
        if DONOTFIXEND or TIME['tfin']!=TIME['tfinal']:
            model.obj.expr+=-sum(model.sw*model.WwRSTORAGE[t,r]*TIME['storage_value'][r,model.t_month[t]]*EF for t,r in model.WwRSTORAGE if t == TIME['tfin'])
    #Fix minimum (=existing) Generic cap investments
    #for pt in model.nptech: 
     #   for pm in model.npmarket:
      #      model.EeGENCAP[TIME['year'],pt,pm].setlb(gencap[pt,pm]) 
            
def safe_solve(model,TIME,gencap,solverstatus='unsolved',max_resolve=25,block=None):
    #Solves a model until a solution is found or the maximum trials is reached
    resolve=1
    def _rescale(m,spn,ov,rf,block=None):
        #m: model, spn:scaling parameter name, ov:original value, rf:rescaling factor
        if block==None:
            sp=m.find_component(spn)
            sp.value=ov*rf
            if spn == 'se': #scaling parameter is energy
                for pt in m.nptech: 
                    for pm in m.npmarket:
                        m.EeGENCAP[TIME['year'],pt,pm].setlb(gencap[pt,pm]/rf)
        else:
            for ens in m.ensemble:
                sp=m.find_component('ensemble['+str(ens)+'].'+spn)
                sp.value=ov*rf
                if spn == 'se': #scaling parameter is energy
                    for pt in m.ensemble[ens].nptech: 
                        for pm in m.ensemble[ens].npmarket:
                            m.ensemble[ens].EeGENCAP[TIME['year'],pt,pm].setlb(gencap[pt,pm]/rf)

    SE=model.se.value if block == None else model.ensemble[0].se.value
    BC=model.debugcost.value if block == None else model.ensemble[0].debugcost.value

    while solverstatus != TerminationCondition.optimal and resolve <= max_resolve:
        try:
            if resolve == 3:
                print('WARNING: re-scaling energy',TIME['ss'],TIME['t'])
                _rescale(model,'se',SE,0.1,block=block)
            elif resolve >=4:
                print('WARNING: emergency mode: random scaling parameters, resolve: '+str(resolve),TIME['ss'],TIME['t'])
                _rescale(model,'se',SE,np.random.choice([0.1,0.5,1,5,10]),block=block)
                _rescale(model,'debugcost',BC,np.random.choice([0.1,0.5,1,5,10]),block=block) 
            resolve+=1
            status=solver.solve(model)
            solverstatus=status.solver.termination_condition
#            if solverstatus == TerminationCondition.optimal:
#                print('Energy scaling',model.find_component('se').value,'Debugcost',model.find_component('debugcost').value)
        except:
            print('WARNING: model crashed, resolve: '+str(resolve),TIME['ss'],TIME['t'])
            resolve+=1
    if solverstatus != TerminationCondition.optimal:
        print('ERROR: Fail after '+str(max_resolve)+' resolve',TIME['ss'],TIME['t'],solverstatus)     


#%% Get reservoir water value
def storage_shadow_value(scenario,parameters,solver):
    save_mpc_param = parameters.val['Options']['MPC',parameters.val['sOptions'][scenario]] #save mpc parameter
    save_debug = parameters.val['Options']['Debug mode',parameters.val['sOptions'][scenario]]
    save_climate = parameters.val['sClimate'][scenario]
    parameters.val['Options']['MPC',parameters.val['sOptions'][scenario]] = 0 #activate perfect foresight
    parameters.val['Options']['Debug mode',parameters.val['sOptions'][scenario]] = 0 #Deactivate debugging (for clean shadowprices)
    #parameters.val['sClimate'][scenario]= 'nile'
    HOM = HydroeconomicOptimization(parameters,scenario=scenario) #perfect foresight model
    HOM.model.dual = Suffix(direction=Suffix.IMPORT) #load shadowprices 
    if solver.name == 'ipopt': #relax strict bounds use in the MPC framework 
        save_brf = solver.options['bound_relax_factor']
        solver.options['bound_relax_factor']=10**-8
    solver.solve(HOM.model) #solve model
    #extract shadowprices
    storage_value=  {(r,m):max(0,-sum(HOM.model.dual[HOM.model.water_waterbalance[t,c]]/len(HOM.model.nyear) 
                        for t in HOM.model.ntime for c in HOM.model.ncatch 
                        if HOM.model.res_catch[r]==c and HOM.model.t_month[t]==m))
                    for r in HOM.model.nres for m in HOM.model.nmonth}
    parameters.val['Options']['MPC',parameters.val['sOptions'][scenario]] = save_mpc_param #reassign original parameter value
    parameters.val['Options']['Debug mode',parameters.val['sOptions'][scenario]] = save_debug #reassign original parameter value
    parameters.val['sClimate'][scenario]= save_climate
    if solver.name == 'ipopt': #reassign original parameter value
        solver.options['bound_relax_factor']=save_brf
    return storage_value
#%% YEARLY ENSEMBLE FORECAST FUNCTION
def RunEnsembleModel_year(parameters,variables,forecast,TIME,solver,rstorage=0,gencap=0):
    #t0=time.time()
    #Update parameters
    parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]=TIME['cc'] #crop choice
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] #first time step of current year
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin'] #final step of model = prediction horizon
    paramup=updatehydrology(TIME,parameters,forecast) 
    #Construct model   
    model=HydroeconomicOptimization(parameters,scenario=TIME['ss'],directparam=paramup).model
    #Load variables, boundary conditions  
    if TIME['year'] != TIME['y0']:
        for yy in range(1,(TIME['tfin']-TIME['tini']+1)//12):
            load_vars(model,variables,onlyvar=['AlCULAREA','EeGENCAP','AcPROD',
                                               'AcTRANS','AcSUPPLY','AcEXTPROD'],reltime=yy)
            load_vars(model,variables,onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY',
                                               'WwOUTFLOW','EeGENPROD','EeHPPROD',
                                               'EeOPPROD','EeTRANS','EeSUPPLY'],reltime=12*yy)
        generate_boundary_conditions(model,TIME,rstorage,gencap)  

    #Solve mpc1 for cultivated area
    SW=model.sw.value
    SE=model.se.value
    safe_solve(model,TIME,gencap)   
    return save_vars(model,['AlCULAREA','EeGENCAP'],sw=SW,se=SE)

def BlockEnsembleModel_year(parameters,variables,ensforecast,TIME,solver,rstorage=0,gencap=0):
    #Update parameters
    parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]=TIME['cc'] #IMPROVE - put somewhere else than in TIME?
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] #first time step of current year
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin']
    #Generate ensemble block model
    model=create_ensemble_block_model(TIME,parameters,ensforecast)  
    #Load variables, boundary conditions  
    for ens in model.ensemble:
        if TIME['year'] != TIME['y0']:
            for yy in range(1,(TIME['tfin']-TIME['tini']+1)//12):
                load_vars(model,variables,onlyvar=['AlCULAREA','EeGENCAP','AcPROD',
                                                   'AcTRANS','AcSUPPLY','AcEXTPROD'],reltime=yy,block=ens)
                load_vars(model,variables,onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY',
                                                   'WwOUTFLOW','EeGENPROD','EeHPPROD',
                                                   'EeOPPROD','EeTRANS','EeSUPPLY'],reltime=12*yy,block=ens)
            generate_boundary_conditions(model.ensemble[ens],TIME,rstorage,gencap) 
    #Generate block model objective function
    block_model_obj_function(model,ensforecast,obj_func=TIME['obj_func'])
    #Solve model
    SE=model.ensemble[0].se.value
    safe_solve(model,TIME,gencap,block=True)
    return save_vars(model,['AlCULAREA','EeGENCAP'],se=SE,block=True)    

#%% MONTHLY ENSEMBLE FORECAST FUNCTION
def BlockEnsembleModel_month(parameters,variables,ensforecast,TIME,solver,rstorage=0,gencap=0):
    #Update parameters
    parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]='fixed'
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] 
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin']
    #Generate ensemble block model
    model=create_ensemble_block_model(TIME,parameters,ensforecast)
    #Load decision variables
    for ens in model.ensemble:
        if TIME['t']!=TIME['t0']:
        #    load_vars(model,variables,block=ens)
        #fix previous decisions
            load_vars(model,variables,fixed=1,time=range(TIME['tini'],TIME['t']),
                      onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY'],block=ens)            
        #load_vars(model,variables,cularea=model.ensemble[ens].nyear,onlyvar=['AlCULAREA'],block=ens)
                
        #Load reservoir levels and generic capacity
        if TIME['year'] != TIME['y0'] or (FIXED_YEAR == 0 and TIME['t']!=TIME['t0']):
            generate_boundary_conditions(model.ensemble[ens],TIME,rstorage,gencap)   
    #Generate block model objective function
    block_model_obj_function(model,ensforecast,obj_func=TIME['obj_func'])
    #Solve model
    SE=model.ensemble[0].se.value
    safe_solve(model,TIME,gencap,block=True)
    return save_vars(model,[str(v.local_name) for v in model.ensemble[0].component_objects(Var)],se=SE,block=True)            
    
def RunEnsembleModel_month(parameters,variables,forecast,TIME,solver,rstorage=0,gencap=0):
    #t0=time.time()
    #Update parameters
    parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]='fixed'
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] 
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin']
    paramup=updatehydrology(TIME,parameters,forecast)
    #Construct models    
    model=HydroeconomicOptimization(parameters,scenario=TIME['ss'],directparam=paramup).model
    #t1=time.time()
    #print('Creating base monthly model',t1-t0)
    #load all variables (to speed up sovler) + fix culture and gencap
    #if TIME['t']!=TIME['t0']:
    #    load_vars(model,variables)
    #load_vars(model,variables,cularea=model.nyear,onlyvar=['AlCULAREA'])
    #load_vars(model,variables,fixed=1,onlyvar=['EeGENCAP'])
    #t2=time.time()
    
    #print('Loading all variables',t2-t1)
    #fix previous decisions
    if TIME['t']!=TIME['tini']:
        load_vars(model,variables,fixed=1,time=range(TIME['tini'],TIME['t']),
                  onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY'])
    #Load reservoir levels and generic capacity
    if TIME['year'] != TIME['y0'] or (FIXED_YEAR == 0 and TIME['t']!=TIME['t0']):
        generate_boundary_conditions(model,TIME,rstorage,gencap) 
                          
    SW=model.sw.value
    SE=model.se.value
    safe_solve(model,TIME,gencap)
    return save_vars(model,[str(v) for v in model.component_objects(Var)],sw=SW,se=SE) # save all variables (to speed up solving time)
    #return save_vars(modelmpc2,['WwRSTORAGE','WwSUPPLY','AwSUPPLY'])#,time=[t for t in modelmpc2.ntime if t <= TIME['t']])
        
#%% Start MPC simulation
def SolveScenario(ss,parameters_in,solver):
# ---- Generate MPC ensemble models ----
    tt=time.time()
    #copy parameters
    parameters=copy.deepcopy(parameters_in)
    #Collect extra hydrology parameters
    if parameters.val['sClimate'][ss] not in parameters.val['wRunOff'].keys():
        parameters.load_hydrology(DataFolderPath,scenarios=[parameters.val['sClimate'][ss]])
    print('formulating scenario: ' + ss)
    #Option dictionnary
    TIME={}
    #Get water shadow price for reservoir hedging
    if 'shadow' in str(parameters.val['Options']['Reservoir Target',parameters.val['sOptions'][ss]]):
        TIME['storage_value']=storage_shadow_value(ss,parameters,solver)
        #print(ss,TIME['storage_value'])
    else:
        TIME['storage_value']=0
    #Generate main model
    #Define if model is MPC - if yes define MPC options
    MPC = parameters.val['Options']['MPC',parameters.val['sOptions'][ss]]
    if MPC != 0 and MPC != 'simulation' and MPC != 'SDPnile':
        TIME['ss']=ss #scenario - used in most MPC functions
        TIME['AFtime']=parameters.val['Options']['Average forecast',parameters.val['sOptions'][ss]]
        TIME['Phorizon']=int(MPC.split('#')[0]) #prediction horizon (=size of MPC model)
        TIME['PFtime']=int(MPC.split('#')[1]) #perfect forecast horizon (adds an horizon of pf to the forecast)
        nEM=int(MPC.split('#')[2]) #number of ensemble members
        ensoption = MPC.split('#')[3] if len(MPC.split('#'))>3 else 0
        weighted= 1 if ensoption in ['w','bm','bw'] else 0 #1: ensemble are not sampled randomly, but their output DV are weighted accordinf to kernel
        if ensoption in ['3c','3cm']:
            weighted= '3class' 
            TIME['CCUTS']=[int(k) for k in parameters.val['Options']['ClassCuts',parameters.val['sOptions'][ss]].split('#')]
        blocksolve= 1 if ensoption in ['bw','bm','3c','3cm'] else 0 #1: ensemble models are aggregated to a single block model
        TIME['obj_func'] = 'maxmin' if ensoption in ['bm','3cm'] else 'weighted' #ensemble model objective function for block model                             
        TIME['cc']=parameters.val['Options']['Crop choice',parameters.val['sOptions'][ss]] #crop choice option for main model
        parameters.val['Options']['Crop choice',parameters.val['sOptions'][ss]]='fixed' #crop choice is fixed for monthly mpc model
        climscen = parameters.val['sClimate'][ss] #climate scenario for forecasts
    #Define models
    HOM      = HydroeconomicOptimization(parameters,scenario=ss) #main model
    TIME['year'] = -1 #initialize value 
    TIME['y0'] = min([y for y in HOM.model.nyear]) #first simulated year
    TIME['t0'] = min([t for t in HOM.model.ntime]) #first simulated time step
    TIME['tfinal'] = HOM.model.Options['tfin'] #final time step of Planning Horizon
#%%SDP nile    
    if MPC == 'SDPnile':
        #options
        tini=HOM.model.Options['tini']
        tend=HOM.model.Options['tfin']
        cname=parameters.val['sClimate'][ss] if SCENITER == 0 or nSS == 0 else parameters.val['sClimate'][ss][:-1]
        ss_sdp={'sTransfer':parameters.val['sTransfer'][ss], 'sReservoir':parameters.val['sReservoir'][ss], 'sWaterDem':parameters.val['sWaterDem'][ss],
                   'sUserVal':parameters.val['sUserVal'][ss],'sHydropower':parameters.val['sHydropower'][ss],'sClimate':cname}
        #Load reservoir values or run sdp backwards
        outpath=0
        if SDP_export == 1:
            name=parameters.val['sUserVal'][ss]+'_'+cname+'_'+parameters.val['sWaterDem'][ss]+'_'+parameters.val['sTransfer'][ss]
            outpath=DataFolderPath + os.sep +'SDP'+ os.sep + name + '_sdpout.txt'
        if outpath!=0 and os.path.exists(outpath)==False:
            sdpout=SDP_backwards(min(parameters.val['ntime']),max(parameters.val['ntime']),ss_sdp,parameters,nruns=NRUNS,iclass=ICLASS,outpath=outpath)
        else:
            sdpout=pickle.load(open(outpath,"rb"))
        #Inflow time series for sdp forward
        runoff=parameters.val['wRunOff'][parameters.val['sClimate'][ss]]
        transopt = 1 if parameters.val['sTransfer'][ss]=='transfer' else 0
        It=[runoff[tini+t,'Upstream']+runoff[tini+t,'Transfer']*transopt for t in range(tend-tini+1)]
        sdpout['data']['It']=It
        sdpout['data']['time']=range(tini,tend+1)
        #Forward run
        Snow=HOM.model.wStorIni['Reservoir'].value
        sdpStorage=SDP_forwards(sdpout,Snow)
        for t in HOM.model.ntime:
            if t < tend-NOFIXLASTY:
                HOM.model.WwRSTORAGE[t,'Reservoir'].value=sdpStorage[t]
                HOM.model.WwRSTORAGE[t,'Reservoir'].fixed=True 
#%% SIMULATION
    if MPC == 'simulation':
        md=HOM.model
        Smin=3
        res='Reservoir'
        catch=['Upstream'] if parameters.val['sTransfer'][ss]!='transfer' else ['Upstream','Transfer']
        STORAGE=md.wStorIni[res].value
        for t in [t for t in md.ntime if t < max(md.ntime)-NOFIXLASTY]:
            INFLOW=sum(md.wRunOff[t,c].value for c in catch)
            DEMAND=sum(md.wUserDem[md.t_year[t],md.t_month[t],u] for u in md.nuser)
            if (t-7)%12 == 0 or t<TIME['t0']+7: 
                if STORAGE >= 60:
                    curt=0
                elif 55<=STORAGE<60:
                    curt=0.05
                elif 50<=STORAGE<55:
                    curt=0.10
                elif STORAGE<50: #45<=
                    curt=0.15
                HP = (2150 + (2850-2150)*(STORAGE-70)/(90-70))*md.m3pers_to_Mm3perMonth/1000
            if curt == 0:
                RELEASE = max(HP,DEMAND)
            else:
                RELEASE = DEMAND * (1-curt)
            stemp = STORAGE + INFLOW - RELEASE 
            if stemp > md.wStorCap[res]:
                RELEASE = RELEASE + (stemp-md.wStorCap[res])
                STORAGE = md.wStorCap[res]
            elif stemp < Smin and STORAGE+INFLOW>=Smin:
                RELEASE = STORAGE+INFLOW-Smin
                STORAGE = Smin
            elif STORAGE+INFLOW<Smin:
                RELEASE = 0
                STORAGE = STORAGE + INFLOW
            else:
                STORAGE = stemp
                
            md.WwOUTFLOW[t,catch[0]].value=RELEASE
            md.WwOUTFLOW[t,catch[0]].fixed=True
            md.WwRSTORAGE[t,res].value=STORAGE
            md.WwRSTORAGE[t,res].fixed=True           
#%%            
    # --- START OF MPC FRAMEWORK --- #
    if MPC != 0 and MPC != 'simulation' and MPC != 'SDPnile':       
    # ---- Iterate over time steps ----
        variables=[] #initiate variables
        for t in range(HOM.model.Options['tini'],HOM.model.Options['tfin']-NOFIXLASTY):
            #print(t)
            TIME['t']=t
            if FIXED_YEAR == 0:
                TIME['year']=HOM.model.t_year[t]
                TIME['tini']=t
                TIME['tfin']=min(t+TIME['Phorizon'],HOM.model.Options['tfin']) if TIME['Phorizon'] != 0 else HOM.model.Options['tfin']
                TIME['flen']=TIME['Phorizon'] #length of forecast
                rstorage = {res:HOM.model.WwRSTORAGE[t-1,res].value for res in HOM.model.nres} if TIME['t'] != TIME['t0'] else {}
                gencap   = {(pt,pm):sum(HOM.model.EeGENCAP[ky,pt,pm].value for ky in HOM.model.nyear if ky < TIME['year'] and ky >= TIME['year']-HOM.model.eLifeTime[pt,pm]) for pt in HOM.model.nptech for pm in HOM.model.npmarket}                    
            
            if FIXED_YEAR == 1 and HOM.model.t_year[t] != TIME['year']: #New year: Run crop choice model                        
                #New year
                TIME['year']=HOM.model.t_year[t]
                TIME['tini']=t
                TIME['tfin']=min(t+TIME['Phorizon'],HOM.model.Options['tfin']) if TIME['Phorizon'] != 0 else HOM.model.Options['tfin']
                TIME['flen']=TIME['Phorizon'] #length of forecast
#                print(ss,TIME['year'])                                                 
#                #Generate forecasts 
#                forecast = knn_bootstrap(parameters,TIME['t'],flen=TIME['flen'],climscen=climscen,nes=nEM,weighted=weighted)                
#                #Get previous storage levels and generic capacity investments    
                rstorage = {res:HOM.model.WwRSTORAGE[t-1,res].value for res in HOM.model.nres} if TIME['year'] != TIME['y0'] else {}
                gencap   = {(pt,pm):sum(HOM.model.EeGENCAP[ky,pt,pm].value for ky in HOM.model.nyear if ky < TIME['year'] and ky >= TIME['year']-HOM.model.eLifeTime[pt,pm]) for pt in HOM.model.nptech for pm in HOM.model.npmarket}
#                #Run ensemble members
#                if blocksolve ==1:
#                    ydecisions=[BlockEnsembleModel_year(parameters,variables,forecast,TIME,solver,rstorage=rstorage,gencap=gencap)]
#                elif PARALLEL_ensemble == 0 or min(epll,nEM)==1:
#                    ydecisions=[RunEnsembleModel_year(parameters,variables,forecast[em],TIME,solver,rstorage=rstorage,gencap=gencap) for em in range(nEM)]                                      
#                elif PARALLEL_ensemble == 1:
#                    pool = mp.Pool(min(epll,nEM))
#                    ydecisions = pool.starmap(RunEnsembleModel_year, [(parameters,variables,forecast[em],TIME,solver,rstorage,gencap) for em in range(nEM)])
#                    pool.close()
#                #Save and fix decision variables of cultivated area (year) for main model
#                #tt=time.time()
#                weights = 0 if weighted == 0 else forecast['Kernel'] #if activated are used to weight the DV from ensemble models
#                yavdecisions=average_vars(ydecisions,weights=weights) #calculate average of DV from ensemble models
#                load_vars(HOM.model,yavdecisions,cularea=[TIME['year']],onlyvar=['AlCULAREA']) #load crop choice in main model
#                #load_vars(HOM.model,yavdecisions,gencap=[TIME['year']],onlyvar=['EeGENCAP'])
#                variables=yavdecisions #save for mpc2 model
            
            #Run ensemble fixed crop choice monthly model           
            #Generate forecasts 
            if TIME['Phorizon'] == 0:
                TIME['flen'] = TIME['tfinal']-TIME['t']
            forecast = knn_bootstrap(parameters,TIME,flen=TIME['flen'],climscen=climscen,nes=nEM,weighted=weighted)
            #Run ensemble members
            if blocksolve == 1:
                mdecisions=[BlockEnsembleModel_month(parameters,variables,forecast,TIME,solver,rstorage=rstorage,gencap=gencap)]
            elif PARALLEL_ensemble == 0 or min(epll,nEM)==1:
                mdecisions=[RunEnsembleModel_month(parameters,variables,forecast[em],TIME,solver,rstorage=rstorage,gencap=gencap) for em in range(nEM)]
            elif PARALLEL_ensemble == 1:
                pool = mp.Pool(min(epll,nEM))
                mdecisions = pool.starmap(RunEnsembleModel_month, [(parameters,variables,forecast[em],TIME,solver,rstorage,gencap) for em in range(nEM)])
                pool.close()            
            #Fix decision variables of time step t in main model to average decision          
            weights = 0 if weighted == 0 else forecast['Kernel'] #if activated are used to weight the DV from ensemble models
            mavdecisions=average_vars(mdecisions,weights=weights) #calculate average of DV from ensemble models
            variables=mavdecisions #save DV for next runs
            load_vars(HOM.model,mavdecisions,time=[TIME['t']],fixed=1,onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY']) #fix water DV in main model
            #Save gencap investments #IMPROVE MOVE TO YEARLY MODEL
            if t == HOM.model.Options['tfin'] or HOM.model.t_year[t+1] != TIME['year']:
                load_vars(HOM.model,mavdecisions,time=[TIME['t']],gencap=[TIME['year']],onlyvar=['EeGENCAP'])
                
        #Set minimum end storage to 0 as it is not possible to garanty it will be reached with MPC
        #for res in HOM.model.nres:
            #getattr(HOM.model,'wStorFin')[res] = 0
    # --- END OF MPC FRAMEWORK --- #
    
    #Solve main model
    if solver.name == 'ipopt': #relax strict bounds use in the MPC framework 
        solver.options['bound_relax_factor']=10**-8
    HOM.model.dual  = Suffix(direction=Suffix.IMPORT) #shadowprices (REM: biased for MPC)
    solverstatus = 'unsolved'
    resolve=1
    while solverstatus != TerminationCondition.optimal and resolve <= 10:
        try:
            resolve+=1
            status=solver.solve(HOM.model)
            solverstatus=status.solver.termination_condition            
        except:
            pass
    
    #logger = logging.getLogger('pyomo.util.infeasible')
    #logger.setLevel(logging.INFO)
    #pyi.log_infeasible_constraints(HOM.model, tol=1E-5, logger=logger)
    #pyi.log_infeasible_bounds(HOM.model, tol=1E-5, logger=logger)
    #if status.solver.termination_condition != TerminationCondition.optimal:
        
    #print(ss)
    #print(status)

    #Read and export results  
    #exppath = ResultFolderPath + os.sep + ss + '.xlsx' 
    result_analysis = ResultAnalysis()
    #results=result_analysis.readresults(HOM.model,parameters,status,PRINT=0,scenario=ss) #   
    #result_analysis.export_to_excel(results,exppath,NEWSHEET=1) #export formated results to excel 
    #result_analysis.export_all_DV(HOM.model,ResultFolderPath,scenario=ss) #export all decision variables
    #result_analysis.export_mass_balances(results,ResultFolderPath,scenario=ss) #export mass balances (for powerBI)
    #result_analysis.export_index_mapping(HOM.model,os.path.join(ResultFolderPath,'Mapping')) #export index mapping/connections (for powerBI)
    
    #Save selected results for scenario analysis
    parameters.val['SCENITER']=SCENITER 
    output=result_analysis.selectedresults(HOM.model,parameters,scenario=ss)
    outputpath=ResultFolderPath + os.sep + ss + '.txt' 
    print(ss,'saving to text file')
    pickle.dump(output,open(outputpath,"wb"))
    print('scenario '+ ss +'solved in '+str(time.time()-tt)+' seconds')
    if SCENITER != 1:
        return output      

#%% Start parallel solving of problems
if __name__ == '__main__':
    #Paths
    DataFolderPath  = os.path.join(dirname, 'Data')
    ResultFolderPath= os.path.join(dirname, 'Results',RESULTFOLDER+time.strftime("%a%d_%m_%Y_%Hh%M"))    
    #Define path to data
    Main            = DataFolderPath + os.sep + 'MainFilenile.xlsx'
    Water           = DataFolderPath + os.sep + 'WaterModule.xlsx'
    Agriculture     = DataFolderPath + os.sep + 'AgricultureModule.xlsx'
    CropMarket      = DataFolderPath + os.sep + 'CropMarketModule.xlsx'
    Energy          = DataFolderPath + os.sep + 'EnergyModule.xlsx'
    Parameter       = DataFolderPath + os.sep + 'Parameters.txt' #parameters (python dictionnaries) saved as txt    
    #Collect parameters
    t=time.time()
    print('Harvesting parameters...')
    parameters      = Database(update=UPDATE,DataFile=Parameter)
    parameters.harvest_all([Main,Water,Energy]) #
    parameters.save(Parameter) #Save parameters
    #%%Generate scenarios
    if SCENITER in [1,2]:        
        pp=parameters.val
        idx0='nileiter'        
        prefix=pp['sOptions'][idx0].split('#')
        cini,cfin,cres,nSS=pp['sClimate'][idx0].split('#')
        cscenp=range(int(cini),int(cfin)+int(cres),int(cres))
        dini,dfin,dres=pp['sWaterDem'][idx0].split('#')
        dscenp=range(int(dini),int(dfin)+int(dres),int(dres))
        transfer=['','a'] if pp['sTransfer'][idx0]==1 else ['']
        iscen=[]
        nSS=int(nSS)
        ssyn=range(nSS) if nSS != 0 else ['']
        for p in prefix:
            for tr in transfer:
                for c in cscenp:
                    for d in dscenp:
                        for kss in ssyn:                
                            index=p+'c'+str(c)+str(kss)+'d'+str(d)+str(tr)
                            iscen.append(index)
                            pp['sClimate'][index]='c'+str(c)+str(kss)
                            pp['sWaterDem'][index]='d'+str(d)
                            if SCENITER == 1: 
                                pp['sRefScen'][index]=p+'c'+str(c)+str(kss)+'d'+str(max(d-10 if tr!='a' else d,int(dini)))
                            else:
                                pp['sRefScen'][index]='PF'+'c'+str(c)+str(kss)+'d'+str(d)+str(tr)
                            pp['sTransfer'][index]='base' if tr!='a' else 'transfer'
                            pp['sOptions'][index]=p
                            pp['sReservoir'][index]=pp['sReservoir'][idx0]
                            pp['sUserVal'][index]=pp['sUserVal'][idx0]
                            pp['sHydropower'][index]=pp['sHydropower'][idx0]
        pp['nscenario']=iscen
        #load synthetic hydro
        if nSS != 0:
            parameters.load_hydrology(DataFolderPath,scenarios=['nileiter'+str(kss) for kss in range(nSS)])            
        for cs in cscenp:
            cfactor={'Upstream':cs/100,'Downstream':1,'Transfer':1}
            for cparam in ['wRunOff','wET0','wRainFall']:
                parameters.val[cparam]['c'+str(cs)]={(t,c):max(0,parameters.val[cparam]['nileiter'][t,c])*cfactor[c] for t,c in parameters.val[cparam]['nileiter'].keys()}
                for kss in ssyn:
                    parameters.val[cparam]['c'+str(cs)+str(kss)]={(t,c):max(0,parameters.val[cparam]['nileiter'+str(kss)][t,c])*cfactor[c] for t,c in parameters.val[cparam]['nileiter'+str(kss)].keys()}
        for ds in dscenp:
            dfactor=ds/100
            parameters.val['wUserDem']['d'+str(ds)]={(u,m):parameters.val['wUserDem']['nile'][u,m]*dfactor for u,m in parameters.val['wUserDem']['nile'].keys()}
    print('*Parameters harvested in '+str(time.time()-t)+' seconds')
    #%%Create folder
    if not os.path.exists(ResultFolderPath):
        os.makedirs(ResultFolderPath)
    
    #Choose solver
    solver      = SolverFactory('ipopt',executable=os.path.join(tmpbin,'ipopt'))#'ipopt',executable="~/miniconda3/pkgs/ipopt_bin-3.7.1-10/bin/ipopt")#'cplex')#
    if solver.name == 'ipopt':
        solver.options['linear_solver']='ma97' #PARDISO
        #solver.options['mu_strategy']='adaptive'
        solver.options['bound_relax_factor']=10**-12
        
    #Run through scenarios
    scenarios = parameters.val['nscenario']
    if PARALLEL_scenario == 1:
        pool = mp.Pool(min(epll,len(scenarios)))    
        parallelresults = pool.starmap(SolveScenario, [(ss,parameters,solver) for ss in scenarios])
        pool.close()
    else:
        parallelresults=[SolveScenario(ss,parameters,solver) for ss in scenarios]
    if SCENITER != 1:
        result_analysis=ResultAnalysis()
        result_analysis.export_scenario_analysis(scenarios,parameters.val['sRefScen'],parallelresults,ResultFolderPath)

