# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 08:48:13 2020

@author: RAPY
"""
import os
tmpbin = os.getenv('TMPBIN')
import numpy as np
#import matplotlib.pyplot as pl
from pyomo.environ import ConcreteModel,Set,Var,NonNegativeReals,Reals,Param,Constraint,minimize,SolverFactory,Objective,Suffix,value
import pickle
#from data_collection import Database                   #Stores all parameters from excel sheets in a python object

#%%
#if __name__ == '__main__':
#    dirname = os.path.abspath(os.path.dirname(__file__))
#    #OPTIONS
#    outname='SDP_out_wI0_3cx_nx'
#    outpath=os.path.join(dirname,outname+'.txt')
#    stor_res=1 #resolution of storage cuts (BCM)
#    nruns=7
#    tini=277#+24
#    tend=588#-240
#    #iclass=[(0,20),(20,80),(80,100)]
#    iclass=[(0,30),(30,60),(60,100)] # as percentiles
#    #iclass=[(0,20),(20,40),(40,60),(60,80),(80,100)] # as percentiles
#    #iclass=[(0,15),(15,30),(30,50),(50,70),(70,85),(85,100)] # as percentiles
#    #scenarios
#    ss={'sTransfer':'base', 'sReservoir':'nile', 'sWaterDem':'nile',
#    'sUserVal':'nilex2','sHydropower':'nile','sClimate':'nile'}
#    #Setting initial storage
#    Snow = 50
#    #load WHAT-IF format data
#    #DataFile='C://Users//rapy//OneDrive - COWI//Scripts//SDP_nile//Parameters.txt'
#    DataFile=os.path.join(dirname,'Parameters.txt')
#    parameters=pickle.load(open(DataFile,"rb"))

#%% Optimization problem
def totalcostfunction(Sini,I,D,Seval,FCeval,SPeval,data):
    
    # This function calculates the total cost from any initial storage level (Sini), given any set of future cost cuts
    
    # Inputs:
    # Sini - initial storage 
    # CC - curtailment cost users in the current time step
    # I - inflow in the current time step
    # D - demand at users in the current time step
    # Seval - storage levels for which the total cost was evaluated in the period following the current period
    # FCeval - total cost values at evaluation points Seval from the period following the current period
    # SPeval - shadow prices at the evaluation points Seval from the period following the current period
    # data - various problem parameters (Storage cap, Hydropower turbines ...)
    
    # Outputs
    # totalcost - total cost for given Sini
    # SP - shadow price of reservoir storage constraint for given Sini
    # FC - future cost for given Sini
    # model - the pyomo model containing all optimized variables and shadow prices
    
    # Number of future cost points from previous evaluation
    nfc = range(0,len(Seval))
    # Convert previous eval results to dict
    Seval = dict(zip(nfc,Seval))
    FCeval = dict(zip(nfc,FCeval))
    SPeval = dict(zip(nfc,SPeval))
    
    # Create a pyomo model
    model = ConcreteModel()
    
    # Declare number of future cost constraints
    model.nfc = Set(initialize=nfc)
    model.nuser = Set(initialize=data['nuser'])

    # Declare decision variables
    model.A  = Var(model.nuser,within=NonNegativeReals) #allocation
    model.R   = Var(within=NonNegativeReals) #Release
    model.HR = Var(within=NonNegativeReals) #Hydropower release
    model.HP = Var(within=NonNegativeReals) #Hydropower production
    model.Send = Var(within=NonNegativeReals) #End storage
    model.FC = Var(within=Reals) #Future cost

    # Declare parameters
    #Reservoir storage
    model.Sini = Param(within=NonNegativeReals,initialize = Sini)
    model.Smax = Param(within=NonNegativeReals,initialize = data['Smax'])
    #Hydropower
    model.hmin = Param(within=NonNegativeReals,initialize = data['hmin']) #minimum head at MOL (m)
    model.hmax = Param(within=NonNegativeReals,initialize = data['hmax']) #maximum head at FSL (m)
    model.qmin = Param(within=NonNegativeReals,initialize = data['qmin']) #minimum turbine flow at MOL (m3/s)
    model.qmax = Param(within=NonNegativeReals,initialize = data['qmax']) #maximum turbine flow at FSL (m3/s)
    model.MOL  = Param(within=NonNegativeReals,initialize = data['MOL']) #Minimum Operational level (BCM)
    model.FSL  = Param(within=NonNegativeReals,initialize = data['FSL']) #Full supply level level (BCM)
    model.hpp  = Param(within=NonNegativeReals,initialize = data['hpp']) #Hydropower production factor (GWh/BCM)
    model.hpv  = Param(within=NonNegativeReals,initialize = data['hpv']) #Hydropower value (M$/GWh)
    #Inflow
    model.I    = Param(within=Reals,initialize = I) #Inflow
    #Demand
    model.D   = Param(model.nuser,within=NonNegativeReals,initialize = D)
    model.av = Param(model.nuser,within=NonNegativeReals,initialize = data['av'])
    #Values
    model.Seval = Param(model.nfc, within=NonNegativeReals,initialize = Seval)
    model.FCeval = Param(model.nfc, within=Reals,initialize = FCeval)
    model.SPeval = Param(model.nfc, within=Reals,initialize = SPeval)
    
    #Functions
    m3pers_to_Mm3perMonth = 1*3600*24*30.44/10**6 #1 m3/s * 3600 s/hour * 24 hours/day * 30.44 days/month / 10^6m3/Mm3
    def RelHeadVol(Vol,MOL,FSL,hmin,hmax):                
        return hmin/hmax + (1-hmin/hmax)*(Vol-MOL)/(FSL-MOL)
            
    def TurbFlowVol(Vol,MOL,FSL,qmin,qmax):
        return (qmin + (qmax-qmin)*(Vol-MOL)/(FSL-MOL))*m3pers_to_Mm3perMonth
    
    #Objective function
    def obj_rule(model):
        cost = -sum(model.av[u]*model.A[u] for u in model.nuser) - model.hpv*model.HP + model.FC   
        return cost
    model.obj = Objective(rule=obj_rule, sense = minimize)

    # Water availability constraint downstream
    def wa_down_c(model):
        return sum(model.A[u] for u in model.nuser) <= model.R
    model.wa_down = Constraint(rule=wa_down_c) 
    
    # Water demand constraint downstream
    def wd_down_c(model,user):
        return model.A[user] <= model.D[user]
    model.wd_down = Constraint(model.nuser, rule=wd_down_c)
    
    # Release constraint
    def rel_c(model):
        return model.R == model.Sini - model.Send + model.I
    model.rel = Constraint(rule=rel_c)
    
    # Max storage constraint
    def maxS_c(model):
        return model.Send <= model.Smax
    model.maxS = Constraint(rule=maxS_c)
    
    # Turbine Capacity constraint
    def maxTurbine(model):
        Vol = (model.Sini+model.Send)/2
        return model.HR <= TurbFlowVol(Vol,model.MOL,model.FSL,model.qmin,model.qmax)
    model.maxTurbine = Constraint(rule=maxTurbine)
    def maxTurbine2(model):
        return model.HR <= model.R
    model.maxTurbine2 = Constraint(rule=maxTurbine2)
    
    # Hydropower production constraint
    def hydroProd(model):
        Vol = (model.Sini+model.Send)/2
        return model.HP == model.HR*model.hpp*RelHeadVol(Vol,model.MOL,model.FSL,model.hmin,model.hmax)
    model.hydroProd = Constraint(rule=hydroProd)
    
    # Future cost constraints
    def fc_c(model, nfc):
        return model.FC >= model.SPeval[nfc]*(model.Send - model.Seval[nfc]) + model.FCeval[nfc]
    model.fc = Constraint(model.nfc, rule=fc_c)
    
    model.dual = Suffix(direction=Suffix.IMPORT)
    
    # Create a solver
    opt = SolverFactory('ipopt',executable=os.path.join(tmpbin,'ipopt'))
    #Solve
    results = opt.solve(model)     
    # Objective value and reservoir storage shadow price
    totalcost = value(model.obj)
    SP = model.dual[model.rel]
    #print(SP)
    FC = model.FC.value  
    results = model
    return totalcost,SP,FC,results

#%% Load Data Function
def load_whatif_data(tini,tend,iclass,ss,parameters):
    #Inflow classes
    mclass=[(kc[0]+kc[1])/2 for kc in iclass]
    nclass=len(iclass)
    t_m=parameters.val['t_month']
    #nmonth=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    nmonth=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    time=range(tini,tend+1) #time steps to be considered
    Inflow=parameters.val['wRunOff'][ss['sClimate']]
    Iclean={t:Inflow[t,'Upstream']+Inflow[t,'Transfer']*(ss['sTransfer']=='transfer') for t in time}
    Im_class=dict()
    for m in nmonth:
        #Im=np.array([max(0,Inflow[t,ca]) for t in time if t_m[t]==m])
        Im=np.array([Iclean[t] for t in time if t_m[t]==m])
        Im_class[m]={c:(np.percentile(Im,iclass[c][0]) if c != 0 else min(Im)-0.001, #min val of class - trick as 0% percentile is minimum value
                        np.percentile(Im,iclass[c][1]), #max val of class
                        np.percentile(Im,mclass[c])) #median val of class           
                for c in range(nclass)}
    #Transition probabilities, one per time step    
    T_count={m:{(c1,c2):sum(1 for t in time[:-1] if t_m[t]==m 
                    and Im_class[m][c1][0]<Iclean[t]<=Im_class[m][c1][1]
                    and Im_class[t_m[t+1]][c2][0]<Iclean[t+1]<=Im_class[t_m[t+1]][c2][1])
            for c1 in range(nclass) for c2 in range(nclass)} 
        for m in nmonth}
    T_prob={m:{(c1,c2):T_count[m][c1,c2]/max(1,sum(T_count[m][c1,ck] for ck in range(nclass))) for c1 in range(nclass) for c2 in range(nclass)} for m in nmonth}
    #Outputs
    TP={m:[[T_prob[m][c1,c2] for c2 in range(nclass)] for c1 in range(nclass)] for m in nmonth}    
    I1=[[Im_class[m][c][2] for m in nmonth] for c in range(nclass)] 
    It=[Iclean[t] for t in time]
    #Convert other data to lacal data style
    data=dict()
    data['Smax']=parameters.val['wStorCap'][ss['sReservoir']]['Reservoir'] #maximum storage level
    data['hmin']=parameters.val['wMinHead'][ss['sReservoir']]['Reservoir'] #minimum head at MOL
    data['hmax']=parameters.val['wMaxHead'][ss['sReservoir']]['Reservoir'] #maximum head at FSL
    data['qmin']=parameters.val['wMinTurb'][ss['sHydropower']]['hydropower'] #minimum flow at MOL
    data['qmax']=parameters.val['wMaxTurb'][ss['sHydropower']]['hydropower'] #maximum flow at FSL
    data['MOL']= parameters.val['wResMOL'][ss['sReservoir']]['Reservoir'] #Minimum Operational level
    data['FSL']= parameters.val['wResFSL'][ss['sReservoir']]['Reservoir'] #Maximum Operational level
    data['hpp']= parameters.val['eHppProd'][ss['sHydropower']]['hydropower']*parameters.val['eHppEff'][ss['sHydropower']]['hydropower'] #Hydropower production factor (Gwh/BCM)
    data['hpv']= parameters.val['eHppVal'][ss['sHydropower']]['hydropower'] #Hydropower value (M$/GWh)
    data['nuser']=parameters.val['nuser'] #users
    data['av']= parameters.val['wUserVal'][ss['sUserVal']] # water demand value [user] (M$/BCM)
    data['demand']= parameters.val['wUserDem'][ss['sWaterDem']] # water demand [month x user] (BCM)
    data['nmonth']=nmonth
    data['time']=time
    data['nclass']=nclass
    data['Im_class']=Im_class
    data['TP']=TP
    data['I1']=I1
    #data['It']=It
    return data,nclass,nmonth,TP,I1

#%% Backwards water value
def SDP_backwards(tini,tend,ss,parameters,iclass=[(0,30),(30,60),(60,100)],nruns=7,stor_res=1,outpath=0):
    #Load data
    data,nclass,nmonth,TP,I1=load_whatif_data(tini,tend,iclass,ss,parameters)
    nseval=round(data['Smax']/stor_res)
    # set this flag to True if you are evaluating the period before end of period and want to inititalize FC
    for run in range(nruns):
        #print('solving round'+str(run))
        if run==0:
            # Define dictionaries for totalcost and shadow prices
            tcall = dict()
            spall = dict()
            # Define S
            Sini = np.arange(stor_res, data['Smax']+1,stor_res, dtype=np.float) #Warning +1 stops working if sotr_res < 1
            Seval = Sini
            EFCeval = np.zeros((nseval,nclass), dtype=np.float)
            ESPeval = np.zeros((nseval,nclass), dtype=np.float)
            npriorruns = 0;        
        else:
            # use totalcost and SP from previous run for FC constraints
            EFCeval = intermediate
            ESPeval = intermediate2
            npriorruns = npriorruns + 1;
        # time step loop
        for m in range(len(nmonth)-1, -1, -1):
            D={u:data['demand'][nmonth[m],u] for u in data['nuser']} #this month's demand
            intermediate = np.zeros((nseval,nclass), dtype=np.float)
            intermediate2 = np.zeros((nseval,nclass), dtype=np.float)
            # inflow scenario loop
            for kclass in range(0,nclass):
                # Take weighted average of future cost and shadow prices
                EFCscen = np.zeros(nseval, dtype=np.float)
                ESPscen = np.zeros(nseval, dtype=np.float)
                for toclass in range (0,nclass):
                    EFCscen = EFCscen + TP[nmonth[m]][kclass][toclass]*EFCeval[:,toclass]
                    ESPscen = ESPscen + TP[nmonth[m]][kclass][toclass]*ESPeval[:,toclass]   
                # Initialize total cost
                totalcost = np.zeros(nseval)
                SP = np.zeros(nseval)
                FC = np.zeros(nseval)
                # Storage loop
                for ieval in range (0, Sini.size):
                    totalcost[ieval], SP[ieval], FC[ieval], results = totalcostfunction(Sini[ieval],I1[kclass][m],D,Seval.tolist(),EFCscen.tolist(),ESPscen.tolist(),data)
                    #print('solved time '+str(t)+' class '+str(k)+' to class '+str(toscen)+' res sclice '+str(i))
                intermediate[:,kclass] = totalcost
                intermediate2[:,kclass] = SP
            EFCeval = intermediate
            ESPeval = intermediate2
            tcall[m+npriorruns*len(nmonth)] = intermediate
            spall[m+npriorruns*len(nmonth)] = intermediate2
    
    #%%Save data
    #Prepare sim data
    spalla = dict()
    spalla4sim = dict()
    tcalla = dict()
    tcalla4sim = dict()
    for k in range (0,nclass):
        spalla[k] = np.zeros((nseval,npriorruns*len(nmonth)+len(nmonth)))
        tcalla[k] = np.zeros((nseval,npriorruns*len(nmonth)+len(nmonth)))
    for i in range(0,npriorruns*len(nmonth)+len(nmonth)):
        for k in range(0,nclass):
            spalla[k][:,i]=spall[i+(npriorruns-i//len(nmonth)*2)*len(nmonth)][:,k]
            tcalla[k][:,i]=tcall[i+(npriorruns-i//len(nmonth)*2)*len(nmonth)][:,k]
    for k in range(0,nclass):
        # Extracting steady-state total costs and water values for simulation
        helpvar = np.zeros((nseval,len(nmonth)))
        helpvar[:,0:len(nmonth)-1] = spalla[k][:,1:len(nmonth)]
        helpvar[:,len(nmonth)-1] = spalla[k][:,0]
        spalla4sim[k] = helpvar
        helpvar = np.zeros((nseval,len(nmonth)))
        helpvar[:,0:len(nmonth)-1] = tcalla[k][:,1:len(nmonth)]
        helpvar[:,len(nmonth)-1] = tcalla[k][:,0]
        tcalla4sim[k] = helpvar
    #send to output
    output={'tcalla4sim':tcalla4sim,'spalla4sim':spalla4sim, 'Seval':Seval, 'data':data,
            'nseval':nseval,'npriorruns':npriorruns,'Sini':Sini} #'tcall':tcall,'spall':spall,           
    #save
    if outpath != 0:
        pickle.dump(output,open(outpath,"wb"))
    return output
#%% SDP forward
def SDP_forwards(output,Snow):
    spsdp=output['spalla4sim'] 
    tcsdp=output['tcalla4sim'] 
    Seval=output['Seval'] 
    nseval=output['nseval']
    data=output['data'] 
    nclass=data['nclass'] 
    nmonth=data['nmonth']  
    TP=data['TP'] 
    Im_class=data['Im_class'] 
    inf=data['It']
    time=data['time']
    # Season and class
    infseason = [nmonth[t%12] for t in range(len(inf))]
    infclass = []
    for t in range(len(inf)):
        infclass.append(sum(c*(Im_class[nmonth[t%12]][c][0]<inf[t]<=Im_class[nmonth[t%12]][c][1]) for c in range(nclass)))    
    #%% time step loop
    totalcost = np.zeros(len(inf))
    SP = np.zeros(len(inf))
    FC = np.zeros(len(inf))
    realc = np.zeros(len(inf))
    Send = np.zeros(len(inf))
    SaveS = {}
    for t in range (0, len(inf)):
        # Demand
        D={u:data['demand'][nmonth[t%12],u] for u in data['nuser']} #this month's demand
        # Get inflow class
        k = infclass[t]
        # Get season
        seas = infseason[t]
        # Take weighted average of future cost and shadow prices
        EFCnow = np.zeros(nseval, dtype=np.float)
        ESPnow = np.zeros(nseval, dtype=np.float)
        for toscen in range (0,nclass):
            EFCnow = EFCnow + TP[seas][k][toscen]*tcsdp[toscen][:,t%12]
            ESPnow = ESPnow + TP[seas][k][toscen]*spsdp[toscen][:,t%12]
        totalcost[t], SP[t], FC[t], results= totalcostfunction(Snow,max(0,inf[t]),D,Seval.tolist(),EFCnow.tolist(),ESPnow.tolist(),data)
    
        #Update storage
        Snow = results.Send.value
        # Archive immediate cost and any other variable of interest
        realc[t] = totalcost[t] - FC[t]
        Send[t] = results.Send.value
        SaveS[time[t]] = results.Send.value
    return SaveS

#%% EXECUTE
#output=SDP_backwards(tini,tend,ss,parameters,iclass=iclass,nruns=nruns,outpath=outpath)    