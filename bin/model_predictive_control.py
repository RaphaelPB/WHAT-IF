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
from pyomo.environ                import Suffix, Var      #Pyomo library 
#Set local directory
dirname = os.path.abspath(os.path.dirname(__file__))
#Import own libraries from "bin" directory
sys.path.append(os.path.join(dirname, 'bin'))
from hydroeconomic_optimization   import HydroeconomicOptimization  #Generates the hydroeconomic optimization model
from result_analysis              import ResultAnalysis             #Exports results to excel sheets

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
PARALLEL_ensemble = 0 #Run ensemble forecasts in parallel
epll2=10 #Maximum number of parallel runs for ensemble forecast 
FIXPOWER=1
VARIANTMETHOD=0 #alternative considering a single weighted forecast for short term, and ensemble forecast for long-term with equal proabilities
fskill=6 #boundary between short term and long term, corresponds to forecasts skill
#%% knn_bootstrap KNN HYDROLOGY GENERATOR
def knn_bootstrap(parameters,TIME,k=20,nes=1,flen=23,ntime=0,weighted=0):
    #t=time.time()
    #parameters = Model parameters
    #ts = time step to start the forecast from
    #k = Number of nearest neighbors considered for resampling/weighted average
    #nes = Number of ensemble members desired
    #flen = Length of the desired forecast
    #ntime = time steps considered for knn_bootstrap analysis (0=use all time steps)
    #cs = climatic scenario of the hydrology data
    #weighted: 1=Instead of sampling uses all nearest neighbors, and saves also the probability/Kernel
    
    if weighted in [1,'ens_class']:
        k=nes
    #Hydrology data
    if TIME['climscen'] is not False: #select runoff scenario
        Q=parameters.val['wRunOff'][TIME['climscen']]
        P=parameters.val['wRainFall'][TIME['climscen']]
        E=parameters.val['wET0'][TIME['climscen']]     
    else: #no runoff scenario
        Q=parameters.val['wRunOff']
        P=parameters.val['wRainFall']
        E=parameters.val['wET0']

    #Indexes
    ts=TIME['t']    
    if ntime == 0:        
        tfin = max(parameters.val['ntime'])
        ntime= np.array([t for t in parameters.val['ntime'] if t != ts and t<tfin-flen and (TIME['stationary']==1 or t<ts)])
        # t = ts is not sampled as we use same time series for bootstrapping and planning - if non-stationary mode we use only past        
    ncatch  = parameters.val['ncatch']
            
    #Create feature vector and observed pattern
    Di=np.array([sum(Q[ts,c] for c in ncatch),
                 sum(P[ts,c] for c in ncatch),
                 sum(E[ts,c] for c in ncatch)])
    d = len(Di)
    Dt = {t:[sum(Q[t,c] for c in ncatch),
             sum(P[t,c] for c in ncatch),
             sum(E[t,c] for c in ncatch)] 
         for t in ntime if (t-ts)%12==0}
    
    #Calculate weights and distance
    stdQ = np.std([sum(Q[t,c] for c in ncatch) for t in ntime if (t-ts)%12==0])
    stdP = np.std([sum(P[t,c] for c in ncatch) for t in ntime if (t-ts)%12==0])
    stdE = np.std([sum(E[t,c] for c in ncatch) for t in ntime if (t-ts)%12==0])
    
    w=[1/stdQ,1/stdP,1/stdE] #weights 
    mtime=np.array([t for t in ntime if (t-ts)%12==0]) #time steps with the same month as ts
    k=min(k,len(mtime)) #reduce k if not enough occurences
    rt=np.array([sum(w[j]*(Di[j]-Dt[t][j])**2 for j in range(d))**0.5 for t in mtime]) #distances of the neighbors
    mtime_t=[mtime[k] for k in range(len(mtime)) if rt[k]!=0] 
    rt_t=[rt[k] for k in range(len(rt)) if rt[k]!=0]
    mtime=np.array(mtime_t)
    rt=np.array(rt_t)
    #Sort neighbours - and select k nearest
    time_sort = mtime[rt.argsort()] #sorted neighbors
    #rt_sort = np.sort(rt) #sorted weights
    
    #Resampling kernel 
    K = [1/l/sum(1/ll for ll in range(1,k+1)) for l in range(1,k+1)] #Sampled by ranking    
    #K = [1/rt_sort[l]/sum(1/rt_sort[ll] for ll in range(k)) for l in range(k)] #Sampled by the distance
    
    #Ensemble forecastst
    forecast={es:{} for es in range(nes)}
    
    #Compute weighted average - stored in ensemble member 0
    if nes==1: #Only average weighted forecast
        forecast[0]['Q_ef'] = {(ts+t,c): sum(K[l]*Q[time_sort[l]+t,c] for l in range(k)) for t in range(1,flen+1) for c in ncatch}
        forecast[0]['P_ef'] = {(ts+t,c): sum(K[l]*P[time_sort[l]+t,c] for l in range(k)) for t in range(1,flen+1) for c in ncatch}
        forecast[0]['E_ef'] = {(ts+t,c): sum(K[l]*E[time_sort[l]+t,c] for l in range(k)) for t in range(1,flen+1) for c in ncatch}
    
    #Real Ensemble forecast
    else:
        if weighted == 'ens_class': #ensemble forecast by grouping nearest neighbors in classes (by forecasted runoff)
            TQ = [sum(Q[tt,c] for tt in range(t,t+flen) for c in ncatch) for t in time_sort]
            ccuts=TIME['ens_class_cuts']
            def cbounds(cl):
                lb=np.percentile(TQ,ccuts[cl-1]) if cl != 0 else -0.001
                up=np.percentile(TQ,ccuts[cl])
                return [lb,up]
            forecast={cl:{} for cl in range(len(ccuts))}
            rclass=[[l for l in range(k) if cbounds(cl)[0]<TQ[l]<=cbounds(cl)[1]] for cl in range(len(ccuts))]
            for cl in range(len(ccuts)):
                if VARIANTMETHOD != 1:
                    forecast[cl]['Q_ef'] = {(ts+t,c): sum(K[l]/sum(K[l] for l in rclass[cl])*Q[time_sort[l]+t,c] for l in rclass[cl]) for t in range(1,flen+1) for c in ncatch}
                    forecast[cl]['P_ef'] = {(ts+t,c): sum(K[l]/sum(K[l] for l in rclass[cl])*P[time_sort[l]+t,c] for l in rclass[cl]) for t in range(1,flen+1) for c in ncatch}
                    forecast[cl]['E_ef'] = {(ts+t,c): sum(K[l]/sum(K[l] for l in rclass[cl])*E[time_sort[l]+t,c] for l in rclass[cl]) for t in range(1,flen+1) for c in ncatch}
                else:
                    forecast[cl]['Q_ef'] = {(ts+t,c): sum(K[l]*Q[time_sort[l]+t,c] for l in range(k)) for t in range(1,fskill+1) for c in ncatch}
                    forecast[cl]['Q_ef'].update({(ts+t,c): sum(1/len(rclass[cl])*Q[time_sort[l]+t,c] for l in rclass[cl]) for t in range(fskill+1,flen+1) for c in ncatch})
                    forecast[cl]['P_ef'] = {(ts+t,c): sum(K[l]*P[time_sort[l]+t,c] for l in range(k)) for t in range(1,fskill+1) for c in ncatch}
                    forecast[cl]['P_ef'].update({(ts+t,c): sum(1/len(rclass[cl])*P[time_sort[l]+t,c] for l in rclass[cl]) for t in range(fskill+1,flen+1) for c in ncatch})
                    forecast[cl]['E_ef'] = {(ts+t,c): sum(K[l]*E[time_sort[l]+t,c] for l in range(k)) for t in range(1,fskill+1) for c in ncatch}
                    forecast[cl]['E_ef'].update({(ts+t,c): sum(1/len(rclass[cl])*E[time_sort[l]+t,c] for l in rclass[cl]) for t in range(fskill+1,flen+1) for c in ncatch})
                    
            if VARIANTMETHOD != 1:
                K=[sum(K[l] for l in rclass[cl]) for cl in range(len(ccuts))]
            else:
                K=[1/len(ccuts) for cl in range(len(ccuts))]
                
        else: #ensemble forecast using directly nearest neighbors
            for es in range(nes):
                #sample neighbor(=time step) with probability K
                if weighted == 0:
                    st = np.random.choice([time_sort[l] for l in range(k)],1,K)[0]
                else: #Use all neighbors (not a Sample)
                    st = time_sort[es]
                #historical time series of sampled neighbor
                forecast[es]['Q_ef']={(ts+t,c): Q[st+t,c] for t in range(1,flen+1) for c in ncatch} 
                forecast[es]['P_ef']={(ts+t,c): P[st+t,c] for t in range(1,flen+1) for c in ncatch}
                forecast[es]['E_ef']={(ts+t,c): E[st+t,c] for t in range(1,flen+1) for c in ncatch}
    
    if weighted in [1,'ens_class']:
        forecast['Kernel']=K
    return forecast
   
#%% FUNCTIONS FOR MODEL PREDICTIVE CONTROL
    
def updatehydrology(TIME,parameters,forecast): #UPDATES HYDROLOGY 
    #creates a dictionary with hydrology parameters according to the forecast and MPC parameters
    #the output can then be passed as directparam=output to the hydroeconomic model
    #forecast is one of the ensemble members from the output of the knn_bootstrapping function
    pftime = TIME['t']+TIME['perfect_forecast']-1 #perfect forecast time (absolute time) - pftime is to simulate a better forecast
    aftime = TIME['average_forecast'] #average forecast time (relative) - past pftime+aftime, climatology is used instead of forecast
    ptime  = range(TIME['tini'],TIME['tfin']+1) #prediction horizon time steps to update        
    ctime  = parameters.val['ntime'] #time steps for which climatology is computed
    if TIME['stationary'] == 0:
        ctime = [t for t in ctime if TIME['t']-12*TIME['climate_horizon'] <= t < TIME['t']]
    cyears = sum(1 for t in ctime if (t-TIME['t'])%12==0)
    ncatch = parameters.val['ncatch'] #catchments
    if TIME['climscen'] is not False: #select runoff scenario
        Q=parameters.val['wRunOff'][TIME['climscen']]
        P=parameters.val['wRainFall'][TIME['climscen']]
        E=parameters.val['wET0'][TIME['climscen']]     
    else: #no runoff scenario
        Q=parameters.val['wRunOff']
        P=parameters.val['wRainFall']
        E=parameters.val['wET0']
    #observed hydrology
    Qo = {(t,c):Q[t,c] for t in ptime for c in ncatch if t <= pftime}
    Po = {(t,c):P[t,c] for t in ptime for c in ncatch if t <= pftime}
    Eo = {(t,c):E[t,c] for t in ptime for c in ncatch if t <= pftime}
    #forecasted hydrology
    Qf = {(t,c):forecast['Q_ef'][t,c] for t in ptime for c in ncatch if pftime<t<=pftime+aftime}    
    Pf = {(t,c):forecast['P_ef'][t,c] for t in ptime for c in ncatch if pftime<t<=pftime+aftime}    
    Ef = {(t,c):forecast['E_ef'][t,c] for t in ptime for c in ncatch if pftime<t<=pftime+aftime}
    #average hydrology (climatology)
    Qa = {(t,c):sum(Q[tt,c] for tt in ctime if (tt-t)%12==0)/cyears 
          for t in ptime for c in ncatch if t>pftime+aftime}    
    Pa = {(t,c):sum(P[tt,c] for tt in ctime if (tt-t)%12==0)/cyears 
          for t in ptime for c in ncatch if t>pftime+aftime}    
    Ea = {(t,c):sum(E[tt,c] for tt in ctime if (tt-t)%12==0)/cyears 
          for t in ptime for c in ncatch if t>pftime+aftime}

    #assemble
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
            if varname in ['WwRSTORAGE','WwOUTFLOW'] and sw != 0:
                CF=model.sw.value/sw if block==False else model.ensemble[0].sw.value/sw
            if varname in ['EeGENCAP','EeGENPROD','EeHPPROD','EeOPPROD','EeTRANS','EeSUPPLY'] and se != 0:
                CF=model.se.value/se if block==False else model.ensemble[0].se.value/se
            if time == []:
                vardic[varname]={index:var[index].value*CF for index in var}
            else:
                vardic[varname]={index:var[index].value*CF for index in var if time==[] or index[0] in time}
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
                    if FIXPOWER==1:
                        model.EeGENCAP[y,pt,pm].fixed=True    
                    
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
    model.nens=Set(initialize=range(len(ensforecast.keys())-1),ordered=True) #ensforecast contains all forecasts and 'Kernel' (hence len-1)
    model.ensemble=Block(model.nens,rule=create_block_ensemble)
    #LINK ENSEMBLE MODELS VARIABLES TOGETHER:
    #define linking rule
    def link_block_variables(m,ens,nvar,nindex):
        if ens == 0: 
            return Constraint.Skip               
        varname1='ensemble['+str(ens)+'].'+nvar
        varname2='ensemble[0].'+nvar #all models are linked to model 0
        index=linkvar[nvar,nindex]
        return m.find_component(varname1)[index]==m.find_component(varname2)[index]
    #define multi index of variables+index that need to be linked to each other (=variable of this time step)
    linkvar={} #dictionary containing indexes of variables that need to be linked
    for var in model.ensemble[0].component_objects(Var): #ensemble[0] arbitrary reference model
        indexlist=[index for index in var if index[0] in range(TIME['t'],TIME['t']+TIME['frequency'])] #WARNING if other index is same as time
        if var.local_name in ['EeGENCAP','AlCULAREA']: #
            indexlist=[index for index in var if index[0]==TIME['year']] #Yearly decisions of current year need to be the same in each ensemble member
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
    for res in model.nres:
        #Final storage target at end prediction horizon (hedging rule)
        if TIME['tfin']!=TIME['tfinal']:
            if TIME['storage_value']==0:
                target = model.Options['Reservoir Target']*model.wStorIni[res].value 
            else: #reservoir target is based on water shadowprice (no hard target)
                target = model.wStorFin[res].value*0
            getattr(model,'wStorFin')[res] = target               
        #Initial storage = to last time step storage
        getattr(model,'wStorIni')[res] = rstorage[res]*model.sw.value               
    
    if TIME['storage_value'] != 0: #reservoir target is based on water shadowprice
        #empirical factor taking into account that shadowprices from Perfect Foresight runs are underestimated
        EF=float(model.Options['Reservoir Target'].split('w')[1]) if model.Options['Reservoir Target'] != 'shadow' else 1
        if TIME['tfin']!=TIME['tfinal']:
            model.obj.expr+=-sum(model.sw*model.WwRSTORAGE[t,r]*TIME['storage_value'][r]*EF for t,r in model.WwRSTORAGE if t == TIME['tfin'])
    #Force existing Generic cap investments
    if FIXPOWER==1: #force the model with all the previous capacity investments as lower bound
        for pt in model.nptech: 
            for pm in model.npmarket:
                model.EeGENCAP[TIME['year'],pt,pm].setlb(gencap[pt,pm]) 
            
def safe_solve(model,TIME,gencap,solverstatus='unsolved',max_resolve=25,block=None):
    #Solves a model until a solution is found or the maximum trials is reached
    resolve=1
    def _rescale(m,spn,ov,rf,block=None):
        #m: model, spn:scaling parameter name, ov:original value, rf:rescaling factor
        #sometimes rescaling helps the model to find a solution - this is what _rescale does
        if block==None:
            sp=m.find_component(spn)
            sp.value=ov*rf
            if spn == 'se' and FIXPOWER==1: #scaling parameter is energy - hence lower bound to gen cap should be updated
                for pt in m.nptech: 
                    for pm in m.npmarket:
                        m.EeGENCAP[TIME['year'],pt,pm].setlb(gencap[pt,pm]/rf)
        else:
            for ens in m.ensemble:
                sp=m.find_component('ensemble['+str(ens)+'].'+spn)
                sp.value=ov*rf
                if spn == 'se' and FIXPOWER==1: #scaling parameter is energy                    
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
            solver=TIME['solver']
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
def storage_shadow_value(scenario,parameters,solver,TIME):
    #save options
    soption=parameters.val['sOptions'][scenario]
    save_debug = parameters.val['Options']['Debug mode',soption]
    save_cropchoice = parameters.val['Options']['Crop choice',soption]    
    parameters.val['Options']['Debug mode',soption] = 0 #Deactivate debugging (for clean shadowprices)
    parameters.val['Options']['Crop choice',soption] = TIME['crop_choice']   
    #Reframe time for non stationary-hydrology
    if TIME['stationary'] == 0:
        mint=min(t for t in parameters.val['ntime'] if (t-TIME['t'])%12==0)
        parameters.val['Options']['tini',soption]=max(mint,TIME['t']-TIME['climate_horizon']*12)
        parameters.val['Options']['tfin',soption]=TIME['t']-1    
    #Create model and solve
    HOM = HydroeconomicOptimization(parameters,scenario=scenario) #perfect foresight model
    HOM.model.dual = Suffix(direction=Suffix.IMPORT) #load shadowprices 
    if solver.name == 'ipopt': #relax strict bounds use in the MPC framework 
        save_brf = solver.options['bound_relax_factor']
        solver.options['bound_relax_factor']=10**-8
    solver.solve(HOM.model) #solve model    
    #Extract shadowprices average value of water at the month corresponding to end of prediction horizon
    storage_value=  {r:-sum(HOM.model.dual[HOM.model.water_waterbalance[t,c]]/len(HOM.model.nyear) 
                        for t in HOM.model.ntime for c in HOM.model.ncatch 
                        if HOM.model.res_catch[r]==c and (t-HOM.model.Options['tfin'])%12==0)
                    for r in HOM.model.nres}
    #Re-assign options
    parameters.val['Options']['Debug mode',soption] = save_debug #reassign original parameter value
    parameters.val['Options']['Crop choice',soption] = save_cropchoice #reassign original parameter value
    if solver.name == 'ipopt': #reassign original parameter value
        solver.options['bound_relax_factor']=save_brf
    return storage_value
#%% YEARLY ENSEMBLE FORECAST FUNCTION
def RunEnsembleModel_year(parameters,variables,forecast,TIME,solver,rstorage=0,gencap=0):
    #Update parameters
    parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]=TIME['crop_choice'] #crop choice
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] #initial time step of current year
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin'] #final time step of prediction horizon
    paramup=updatehydrology(TIME,parameters,forecast)
    #Construct model   
    model=HydroeconomicOptimization(parameters,scenario=TIME['ss'],directparam=paramup).model
    #Load variables, boundary conditions  
    if TIME['year'] != TIME['y0']:
        for yy in range(1,(TIME['tfin']-TIME['tini']+1)//12+1): #load "shifted" decision variables to speed up solving time
            load_vars(model,variables,onlyvar=['AlCULAREA','EeGENCAP','AcPROD',
                                               'AcTRANS','AcSUPPLY','AcEXTPROD'],reltime=yy)
            load_vars(model,variables,onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY',
                                               'WwOUTFLOW','EeGENPROD','EeHPPROD',
                                               'EeOPPROD','EeTRANS','EeSUPPLY'],reltime=12*yy)
        generate_boundary_conditions(model,TIME,rstorage,gencap)  
    #Solve yearly model and save all variables
    SE=model.se.value
    safe_solve(model,TIME,gencap)   
    #return save_vars(model,['AlCULAREA','EeGENCAP'],sw=SW,se=SE)
    return save_vars(model,[str(v.local_name) for v in model.component_objects(Var)],se=SE)

def BlockEnsembleModel_year(parameters,variables,ensforecast,TIME,solver,rstorage=0,gencap=0):
    #Update parameters
    parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]=TIME['crop_choice'] 
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] #first time step of current year
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin']
    #Generate ensemble block model
    model=create_ensemble_block_model(TIME,parameters,ensforecast)
    #Load variables, boundary conditions  
    for ens in model.ensemble:
        if TIME['year'] != TIME['y0']:
            for yy in range(1,(TIME['tfin']-TIME['tini']+1)//12+1): #load "shifted" decision variables to speed up solving time
                load_vars(model,variables,onlyvar=['AlCULAREA','EeGENCAP','AcPROD',
                                                   'AcTRANS','AcSUPPLY','AcEXTPROD'],reltime=yy,block=ens)
                load_vars(model,variables,onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY',
                                                   'WwOUTFLOW','EeGENPROD','EeHPPROD',
                                                   'EeOPPROD','EeTRANS','EeSUPPLY'],reltime=12*yy,block=ens)
            generate_boundary_conditions(model.ensemble[ens],TIME,rstorage,gencap) 
    #Generate block model objective function
    block_model_obj_function(model,ensforecast,obj_func=TIME['obj_func'])
    #Solve monthly model and save all variables
    SE=model.ensemble[0].se.value
    safe_solve(model,TIME,gencap,block=True)
    return save_vars(model,[str(v.local_name) for v in model.ensemble[0].component_objects(Var)],se=SE,block=True)
    #return save_vars(model,['AlCULAREA','EeGENCAP'],se=SE,block=True)    

#%% MONTHLY ENSEMBLE FORECAST FUNCTION
def BlockEnsembleModel_month(parameters,variables,ensforecast,TIME,solver,rstorage=0,gencap=0):
    #Update parameters
    if TIME['YWR'] != 'nonlinear':
        parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]='fixed'
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] 
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin']
    #Generate ensemble block model
    model=create_ensemble_block_model(TIME,parameters,ensforecast)
    #Load decision variables
    for ens in model.ensemble:
        load_vars(model,variables,block=ens) #preload variables from yearly model to speed up
        if TIME['t']!=TIME['tini']:            
        #fix previous decisions
            load_vars(model,variables,fixed=1,time=range(TIME['tini'],TIME['t']),
                      onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY'],block=ens)
        if TIME['YWR'] != 'nonlinear':            
            load_vars(model,variables,cularea=model.ensemble[ens].nyear,onlyvar=['AlCULAREA'],block=ens)
        else:
            load_vars(model,variables,time=[TIME['year']],onlyvar=['AlCULAREA'],block=ens,fixed=1)
        load_vars(model,variables,onlyvar=['EeGENCAP'],block=ens,fixed=1) 
                
        #Load reservoir levels and generic capacity
        if TIME['year'] != TIME['y0']:
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
    if TIME['YWR'] != 'nonlinear':
        parameters.val['Options']['Crop choice',parameters.val['sOptions'][TIME['ss']]]='fixed'
    parameters.val['Options']['tini',parameters.val['sOptions'][TIME['ss']]] = TIME['tini'] 
    parameters.val['Options']['tfin',parameters.val['sOptions'][TIME['ss']]] = TIME['tfin']
    paramup=updatehydrology(TIME,parameters,forecast)
    #Construct models    
    model=HydroeconomicOptimization(parameters,scenario=TIME['ss'],directparam=paramup).model
    #load all variables (to speed up sovler) + fix culture and gencap
    load_vars(model,variables)
    if TIME['t']!=TIME['tini']:
    #fix previous decisions
        load_vars(model,variables,fixed=1,time=range(TIME['tini'],TIME['t']),
                  onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY'])
    if TIME['YWR'] != 'nonlinear': 
        load_vars(model,variables,cularea=model.nyear,onlyvar=['AlCULAREA'])
    else:
        load_vars(model,variables,time=[TIME['year']],onlyvar=['AlCULAREA'],fixed=1)
    load_vars(model,variables,onlyvar=['EeGENCAP'],fixed=1)    

    #Load reservoir levels and generic capacity
    if TIME['year'] != TIME['y0']:
        generate_boundary_conditions(model,TIME,rstorage,gencap) 
                          
    SE=model.se.value
    safe_solve(model,TIME,gencap)
    return save_vars(model,[str(v) for v in model.component_objects(Var)],se=SE) # save all variables (to speed up solving time)
        
#%% Start MPC simulation
def SolveScenario(ss,parameters_in,solver,paths):
# ---- Generate MPC ensemble models ----
    tt=time.time()
    #copy parameters
    parameters=copy.deepcopy(parameters_in)
    #Collect extra parameters
    parameters.load_csv(paths['data'],ss)
    #if parameters.val['sClimate'][ss] not in parameters.val['wRunOff'].keys():
    #    parameters.load_hydrology(paths['data'],scenarios=[parameters.val['sClimate'][ss]])
    print('formulating scenario: ' + ss)    
           
    TIME={} #Option dictionnary

    #Define model
    #save yield water rersponse and crop choice parameter
    TIME['crop_choice']=parameters.val['Options']['Crop choice',parameters.val['sOptions'][ss]] #crop choice option for main model
    TIME['YWR']=parameters.val['Options']['Yield water response',parameters.val['sOptions'][ss]] #yield water response option for main model
    if TIME['YWR'] != 'nonlinear':
        parameters.val['Options']['Crop choice',parameters.val['sOptions'][ss]]='fixed'
    #Create main model
    HOM = HydroeconomicOptimization(parameters,scenario=ss) #main model
    
    #Define MPC options
    oopt=HOM.model.Options #shortcut for options
    MPC=oopt['MPC'] if 'MPC' in oopt.keys() else 0
    if MPC == MPC and MPC != 0: #
        #Model Predictive Control Parameters
        TIME['average_forecast']=oopt['Average forecast'] #Horizon at which climatology is used
        TIME['prediction_horizon']=int(MPC.split('#')[0]) #Prediction horizon (=size of MPC model)
        TIME['perfect_forecast']=int(MPC.split('#')[1]) #Perfect forecast horizon (real hydrology is used at this horizon)
        nEM=int(MPC.split('#')[2]) #Number of ensemble members (if =1 weighted forecast)
        ensoption = MPC.split('#')[3] if len(MPC.split('#'))>3 else 0 #Method to derive optimal control actions from forecast
        weighted= 1 if ensoption in ['w','bm','bw','ec'] else 0 #1: ensemble are not sampled randomly, but their output DV are weighted accordinf to kernel
        if ensoption == 'ec':
            weighted= 'ens_class'
            TIME['ens_class_cuts']=[int(k) for k in oopt['Ensemble Class Cuts'].split('#')]
        blocksolve= 1 if ensoption in ['bw','bm','ec'] else 0 #1: ensemble models are aggregated to a single block model
        TIME['obj_func'] = 'weighted' if ensoption in ['bw','ec'] else 'maxmin' #ensemble model objective function for block model                                     
        
        #Other relevant parameters
        if 'wRunOff' in parameters.scen.keys(): #scenarios on runoff exist
            TIME['climscen'] = parameters.val[parameters.scen['wRunOff']][ss] #climate scenario for forecast
        else:
            TIME['climscen'] = False
        TIME['ss']=ss #scenario 
        TIME['solver']=solver #solver
        
        #Define if framework considers stationarity of hydrology
        TIME['stationary']=1 #default is stationary hydrology
        stationary=oopt['Stationary'] if 'Stationary' in oopt.keys() else 1
        if stationary != 1:
            TIME['stationary']=0 #non-stationary hydrology
            TIME['climate_horizon']=int(stationary.split('#')[0]) #Number of past years used to compute average conditions
            TIME['update_shadow']=int(stationary.split('#')[1]) #Interval at which reservoir storage value is re-evaluated
        
        #Get water shadow price for reservoir hedging
        target_option=str(oopt['Reservoir Target']) if 'Reservoir Target' in oopt.keys() else 0
        TIME['storage_value']=0 #initialize 
        #Frequency of decision solving (needs to be lower or equal to perfect forecast)
        TIME['frequency']=oopt['Frequency'] if 'Frequency' in oopt.keys() else 1
        
        #Initialize time parameters
        TIME['year']=-1 #initialize value 
        TIME['y0']=min([y for y in HOM.model.nyear]) #first simulated year
        TIME['t0']=min([t for t in HOM.model.ntime]) #first simulated time step
        TIME['tfinal']=oopt['tfin'] #final time step of Planning Horizon
        
    # --- START OF MPC FRAMEWORK --- #
    if MPC != 0:       
    # ---- Iterate over time steps ----
        variables=[] #initiate variables
        for t in range(oopt['tini'],oopt['tfin']-oopt['Buffer'],TIME['frequency']):
            print(t)
            TIME['t']=t                    
            if HOM.model.t_year[t] != TIME['year']: #New year                      
                TIME['year']=HOM.model.t_year[t]
                TIME['tini']=t #start of prediction horizon
                TIME['tfin']=min(t+TIME['prediction_horizon']*12-1,oopt['tfin']) #end of prediction horizon (max end of planning period)
                TIME['flen']=TIME['prediction_horizon']*12-1 #length of forecast
                print(ss,TIME['year'])
                #Reservoir hedging rule - based on perfect foresight shadowprice
                if 'shadow' in target_option and (TIME['storage_value']==0 or (TIME['stationary']==0 and (TIME['year']-TIME['y0'])%TIME['update_shadow']==0)):
                    TIME['storage_value']=storage_shadow_value(ss,parameters,solver,TIME)
                    print(ss,TIME['storage_value'])       
                #Generate forecasts 
                forecast = knn_bootstrap(parameters,TIME,flen=TIME['flen'],nes=nEM,weighted=weighted)
                #Get previous storage levels and generic capacity investments    
                rstorage = {res:HOM.model.WwRSTORAGE[t-1,res].value for res in HOM.model.nres} if TIME['year'] != TIME['y0'] else {}
                gencap   = {(pt,pm):sum(HOM.model.EeGENCAP[ky,pt,pm].value for ky in HOM.model.nyear 
                                        if ky < TIME['year'] and ky >= TIME['year']-HOM.model.eLifeTime[pt,pm]) 
                            for pt in HOM.model.nptech for pm in HOM.model.npmarket}
                #Run ensemble members
                if blocksolve ==1:
                    ydecisions=[BlockEnsembleModel_year(parameters,variables,forecast,TIME,solver,rstorage=rstorage,gencap=gencap)]
                elif PARALLEL_ensemble == 0 or min(epll2,nEM)==1:
                    ydecisions=[RunEnsembleModel_year(parameters,variables,forecast[em],TIME,solver,rstorage=rstorage,gencap=gencap) for em in range(nEM)]                                      
                elif PARALLEL_ensemble == 1:
                    pool = mp.Pool(min(epll2,nEM))
                    ydecisions = pool.starmap(RunEnsembleModel_year, [(parameters,variables,forecast[em],TIME,solver,rstorage,gencap) for em in range(nEM)])
                    pool.close()
                #Save and fix decision variables for main model
                weights = 0 if weighted == 0 else forecast['Kernel'] #if activated are used to weight the DV from ensemble models (rmq won't apply to merged ensemble model)
                yavdecisions=average_vars(ydecisions,weights=weights) #calculate average of DV from ensemble models
                if TIME['YWR'] != 'nonlinear':
                    load_vars(HOM.model,yavdecisions,cularea=[TIME['year']],onlyvar=['AlCULAREA']) #load crop choice in main model
                else:
                    load_vars(HOM.model,yavdecisions,time=[TIME['year']],onlyvar=['AlCULAREA'],fixed=1) #load crop choice in main model
                load_vars(HOM.model,yavdecisions,gencap=[TIME['year']],onlyvar=['EeGENCAP']) #load generic capacity investments in main model
                variables=yavdecisions #save for monthly model
            
            #Run ensemble monthly model (with fixed yearly decision variables)          
            #Generate forecasts
            flen=TIME['flen']-(TIME['t']-TIME['tini']) #forecast lenght
            forecast = knn_bootstrap(parameters,TIME,flen=flen,nes=nEM,weighted=weighted)
            #Run ensemble members
            if blocksolve == 1:
                mdecisions=[BlockEnsembleModel_month(parameters,variables,forecast,TIME,solver,rstorage=rstorage,gencap=gencap)]
            elif PARALLEL_ensemble == 0 or min(epll2,nEM)==1:
                mdecisions=[RunEnsembleModel_month(parameters,variables,forecast[em],TIME,solver,rstorage=rstorage,gencap=gencap) for em in range(nEM)]
            elif PARALLEL_ensemble == 1:
                pool = mp.Pool(min(epll2,nEM))
                mdecisions = pool.starmap(RunEnsembleModel_month, [(parameters,variables,forecast[em],TIME,solver,rstorage,gencap) for em in range(nEM)])
                pool.close()            
            #Fix decision variables of time step t in main model to average decision          
            weights = 0 if weighted == 0 else forecast['Kernel'] #if activated are used to weight the DV from ensemble models (rmq won't apply to merged ensemble model)
            mavdecisions=average_vars(mdecisions,weights=weights) #calculate average of DV from ensemble models
            variables.update(mavdecisions) #update save variables with last monthly run
            load_vars(HOM.model,mavdecisions,time=range(TIME['t'],TIME['t']+TIME['frequency']),
                      fixed=1,onlyvar=['WwRSTORAGE','WwSUPPLY','AwSUPPLY']) #fix water DV in main model
    # --- END OF MPC FRAMEWORK --- #    
    #%%Solve main model
    if solver.name == 'ipopt': #relax strict bounds use in the MPC framework 
        solver.options['bound_relax_factor']=10**-8
    HOM.model.dual  = Suffix(direction=Suffix.IMPORT) #shadowprices (REM: biased for MPC)    
    
    #Safe solve for MPC
    solverstatus = 'unsolved'
    resolve=1
    while solverstatus != TerminationCondition.optimal and resolve <= 10:
        try:
            resolve+=1
            status=solver.solve(HOM.model)
            solverstatus=status.solver.termination_condition            
        except:
            pass
    
    #Resolve for investment
    if oopt['Investment module'] == 1: #Fix selected investment to switch to fully linear model and get dual values       
        for ip in HOM.model.ninvphase:
            for inv in HOM.model.ninvest:   
                HOM.model.IbINVEST[ip,inv].fixed = True
        print('...re-solving scenario ' + ss)
        solverstatus    = solver.solve(HOM.model)
        
    #Read and export results  
    result_analysis = ResultAnalysis()
    if paths['export'] in ['all','xlsx']:
        exppath = paths['output'] + os.sep + ss + '.xlsx' 
        results=result_analysis.readresults(HOM.model,parameters,status,PRINT=0,scenario=ss) #   
        result_analysis.export_to_excel(results,exppath,NEWSHEET=1) #export formated results to excel 
    if paths['export'] in ['all']:    
        result_analysis.export_all_DV(HOM.model,paths['output'],scenario=ss) #export all decision variables
        result_analysis.export_mass_balances(results,paths['output'],scenario=ss) #export mass balances (for powerBI)
        result_analysis.export_index_mapping(HOM.model,os.path.join(paths['output'],'Mapping')) #export index mapping/connections (for powerBI)
    
    #Save selected results for scenario analysis
    output=result_analysis.selectedresults(HOM.model,parameters,scenario=ss)
    outputpath=paths['output'] + os.sep + ss + '.txt' 
    print(ss,' saving to text file')
    pickle.dump(output,open(outputpath,"wb"))
    print('scenario '+ ss +' solved in '+str(time.time()-tt)+' seconds')
    return output      