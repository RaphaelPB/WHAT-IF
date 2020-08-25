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

from pyomo.environ import ConcreteModel, Set, Param, Var, Objective, Constraint, Expression
from pyomo.environ import NonNegativeReals, Binary, exp, value
#import math

class HydroeconomicOptimization():
    def __init__(self,parameters,scenario='WHATIF_main',directparam=0):
#KEYWORD LIST
        #Water Water, #Power Power, #Engy Energy, #Crop Crop, #Fuel Fuel
        #River River, #User User, #Res Reservoir, #Opp Other power plant, #Hpp Hydropower plant
        #Irg Irrigation, #Cul Cultivation, #Curt Curtailment, #Share Share/Percentage, 
        #Cap Capacity, #Eff Efficiency, #Cost Cost, #Val Value (WTP, marginal or shadow), #Dem Demand
        #loss Loss rate, #Evap Evaporation (water surface) or Evapotranspiration (land), #Prod Production
        #Trans for Water:Transfer, Crops:Transport, Electricity:Transmission 

#%%PRAMETER READING FUNCTION
        
        def read(name,time='all',option=1,index=0):
            return parameters.read_param(name,option=option,index=index,time=time,scenario=scenario,directparam=directparam)
        
#%%INTIALIZE MODEL
        md = ConcreteModel()
    #Options
        md.noption = Set(initialize=read('noption'))    #Configurations of the model
        Options = {opt:parameters.val['Options'][opt,parameters.val['sOptions'][scenario]] for opt in md.noption} #Select options of considered scenario
        md.Options = Param(md.noption, initialize=Options)
        #time framing
        ttime,t_prev,year=parameters.read_time(scenario)
    #Unit conversion
        md.mm_to_m3perha    = 1/1000*10000 # 0,001m * 10000m2/ha
        md.mm_to_m          = 1/1000# 0,001m * 10000m2/ha
        md.MW_to_GWhperMonth= 1*24*365/12/1000 # 1 MW * 24hours/day * 365/12days/month / 1000MW/GW 
        md.m3pers_to_Mm3perMonth = 1*3600*24*30.44/10**6 #1 m3/s * 3600 s/hour * 24 hours/day * 30.44 days/month / 10^6m3/Mm3
        md.kha_to_Mha       = 1/1000 # 1000 ha = 1/1000 Mha
        md.kt_to_Mt         = 1/1000 # 1000 t = 1/1000 Mt
    #Scaling parameters
        md.debugcost = Param(initialize=100000,mutable=True) #debug cost (variabes units $/m3 - 100$/t)   
        md.sw = Param(initialize=1000,mutable=True) # Scaling factor for WwRSTORAGE and WwOUTFLOW (1000=convert DV from Mm3 to km3)
        md.se = Param(initialize=1000,mutable=True) # Scaling factor for all energy decision variables
        
    # WETLAND MODULE UNDER DEVELOPMENT
        ANALYTICAL=1 if 'ANALYTICAL' not in md.Options.keys() else md.Options['ANALYTICAL']
        #if 1 uses continuous solution for wetland evaporation, any other value uses discrete evaporation
        OUTFLOWPOWER=1 if 'OUTFLOWPOWER' not in md.Options.keys() else md.Options['OUTFLOWPOWER']
        # Outflow(V) = wAlpha * V**OUTFLOWPOWER
        AREAPOWER=1 if 'AREAPOWER' not in md.Options.keys() else md.Options['AREAPOWER']
        # Area(V) = wkV* V**AREAPOWER + wResArea
        AREACONSTRAINT='minimum' if 'AREACONSTRAINT' not in md.Options.keys() else md.Options['AREACONSTRAINT']
        #%%###############
        #**Water Module**#
        ##################
        
    #Indices
        md.ncountry         = Set(initialize=read('ncountry'))                 #Countries [indice country] (id)
        md.ncatch           = Set(initialize=read('catch_ds').keys())          #Catchments [indice catchment] (id)
        md.ntime            = Set(initialize=ttime)                            #Time steps [indice time] (id)
        md.nyear            = Set(initialize=year)                             #Years [indice year] (id)
        md.nmonth           = Set(initialize=read('monthorder').values())      #Months [indice month] (id)         
        md.nres             = Set(initialize=read('wStorCap').keys())          #Reservoirs [indice reservoir] (id) 
        md.nuser            = Set(initialize=read('user_catch').keys())        #Water Users [indice user] (id)
        md.naquifer         = Set(initialize=read('wGwFlow',option=md.Options['Groundwater']).keys())  #Groundwater aquifers [indice aquifer] (id)           
    #Spatio-Temporal connections
        md.catch_ds         = Param(md.ncatch, initialize=read('catch_ds'))        #Downstream catchment of catchment [catchment] (id)
        md.catch_country    = Param(md.ncatch, initialize=read('catch_country'))   #Country of catchment [catchment] (id) #REM: DOES NOT ALWAYS MAKE SENSE
        md.res_catch        = Param(md.nres,  initialize=read('res_catch'))        #Reservoir's catchment [reservoir] (id)
        md.user_catch       = Param(md.nuser, initialize=read('user_catch'))       #User's catchment [user] (id)
        md.user_dscatch     = Param(md.nuser, initialize=read('user_dscatch'))     #User's downstream catchment [user] (id)
        md.user_country     = Param(md.nuser, initialize=read('user_country'))     #User's country [user] (id)
        md.t_month          = Param(md.ntime, initialize=read('t_month',time=ttime))#Month of time step [time] (id)
        md.t_year           = Param(md.ntime, initialize=read('t_year',time=ttime))#Year of time step [time] (id)
        md.t_prev           = Param(md.ntime, initialize=t_prev)                   #Previous time step [time] (id)       
            
    #Hydrology
        md.wRunOff          = Param(md.ntime, md.ncatch, mutable=True, initialize=read('wRunOff',time=ttime))   #Runoff [time x catchment] (Mm3/month)
        md.wRainFall        = Param(md.ntime, md.ncatch, mutable=True, initialize=read('wRainFall',time=ttime)) #Precipitation [time x catchment] (mm)
        md.wET0             = Param(md.ntime, md.ncatch, mutable=True, initialize=read('wET0',time=ttime))      #ET0 [time x catchment] (mm)
        md.wFlowLoss        = Param(md.ncatch, initialize=read('wFlowLoss'))                                   #Flow loss [catchment] (%)
    #Groundwater
        opt= md.Options['Groundwater']
        md.wGwRech          = Param(md.ntime, md.naquifer, initialize=read('wGwRech',option=opt,time=ttime))    #Groundwater recharge [time x catchment] (Mm3/month)
        md.wGwFlow          = Param(md.naquifer, initialize=read('wGwFlow',option=opt))        #Groundwater baseflow coeeficient [catchment] (-)
        md.wGwIni           = Param(md.naquifer, initialize=read('wGwIni',option=opt))         #Groundwater initial storage (if not cyclic) [catchment] (Mm3)
        md.wGwCost          = Param(md.naquifer, initialize=read('wGwCost',option=opt))        #Groundwater pumping cost [catchment] ($/m3)
    #Water Users
        md.wUserDem         = Param(md.nyear, md.nmonth, md.nuser, initialize=read('wUserDem',time=md.nyear))#User water demand [year x month x user] (Mm3/month)        
        md.wUserRturn       = Param(md.nuser, initialize=read('wUserRturn'))                   #Return rate of user water allocation [user] (%)
        md.wUserLoss        = Param(md.nuser, initialize=read('wUserLoss'))                    #Loss rate of user water allocation [user] (%)
        md.wUserVal         = Param(md.nuser, initialize=read('wUserVal'))                     #User marginal value for water [user] ($/m3)
        md.wSupCost         = Param(md.nuser, initialize=read('wSupCost'))                     #Water supply costs [user] ($/m3) 
    #Reservoirs
        md.wStorCap         = Param(md.nres, initialize=read('wStorCap'))                      #Storage Capacity [reservoir] (Mm3)
        if md.Options['Reservoirs'] == 'flood':
            md.wFloodRule   = Param(md.nres, md.nmonth, initialize=read('wFloodRule'))         #Flood rule curve [reservoir x month] (Mm3)
        if md.Options['Initial time step'] == 1:
            md.wStorIni     = Param(md.nres, mutable=True, initialize=read('wStorIni'))        #Initial storage [reservoir] (Mm3) - only used if initial time step is used
            md.wStorFin     = Param(md.nres, mutable=True, initialize=read('wStorFin'))        #Final storage [reservoir] (Mm3) - only used if initial time step is used
        md.wkV              = Param(md.nres, initialize=read('wkV'))                           #Area - volume coef, kV in "A = A0 + kV * V" [reservoir] (ha/m3)
        md.wResArea         = Param(md.nres, initialize=read('wResArea'))                      #Minimum area, A0 in "A = A0 + kV * V" [reservoir] (Mha)
    #Transfer schemes
        opt= md.Options['Transfers']
        md.ntransfer        = Set(initialize=read('wTransCap',option=opt).keys())              #Transfer schemes [indice transfer] (id)
        md.transfer_us      = Param(md.ntransfer, initialize=read('transfer_us',option=opt))   #Upstream catchment of transfer scheme [transfer] (id) 
        md.transfer_ds      = Param(md.ntransfer, initialize=read('transfer_ds',option=opt))   #Downstream catchment of transfer scheme [transfer] (id) 
        md.wTransCap        = Param(md.ntransfer, initialize=read('wTransCap',option=opt))     #Capacity of transfere scheme [transfer] (Mm3/month)
        md.wTransLoss       = Param(md.ntransfer, initialize=read('wTransLoss',option=opt))    #Loss rate of transfere scheme [transfer] (%)
    #Lakes and wetlands
        opt= md.Options['Lakes']
        md.nlake            = Set(initialize=[res for res in md.nres if read('res_type')[res]=='wetland'] if md.Options['Lakes']==1 else []) #Lakes are a sub-ensemble of reservoirs
        md.res_type         = Param(md.nres,  initialize=read('res_type',option=opt))      #Type of reservoir 'reservoir'=can be controlled [reservoir] (id)
        md.wAlpha           = Param(md.nres, initialize=read('wAlpha',option=opt))        #First order outflow coefficient Out = alpha * V [lake] (%)
        md.wkRO             = Param(md.nres, initialize=read('wkRO',option=opt))          #Share of runoff not diverted by wetland [lake] (%)
        md.wkUP             = Param(md.nres, initialize=read('wkUP',option=opt))          #Share of upstream flow not diverted by wetland [lake] (%)       
        
#%%**Environmental flow Module**#
        md.neflow           = Set(initialize=read('eflow_catch',option=md.Options['Eflows']).keys())               #Environmental constraint point [indice envconst] (id)
        md.eflow_catch      = Param(md.neflow, initialize=read('eflow_catch',option=md.Options['Eflows']))         #Catchment of environmental flow constraint [envconst] (id)
        md.eflow_hard       = Param(md.neflow, initialize=read('eflow_hard',option=md.Options['Eflows']))          #Eflow as hard or soft (limited to actual runoff) constraint [envconst] (binary)
        md.wEnvFlow         = Param(md.nmonth, md.neflow, initialize=read('wEnvFlow',option=md.Options['Eflows'])) #Env flow [month x envconst] (Mm3/month)
        
        #%%#####################
        #**Agriculture Module**#
        ########################  
        
    #Indices 
        opt=md.Options['Agriculture']
        md.ncrop            = Set(initialize=read('ncrop',option=opt,index=1))       #Crops [indice crop] (id)
        md.nypath           = Set(initialize=read('nypath',option=opt,index=1))      #Yield water response path [indice ypath] (id)
        md.nyphase          = Set(initialize=read('nyphase',option=opt,index=1))     #Yield water response phases [indice yphase] (id)
        md.nfzone           = Set(initialize=read('fzone_catch',option=opt).keys())  #Farm zones [indice fzone] (id)
        md.nifzone          = Set(initialize=[fz for fz in md.nfzone if read('aIrrigation')[read('fzone_type')[fz]]==1]) #Irrigated farming zones subindice of farming zones [indice ifzone] (id)
        md.nftype           = Set(initialize=read('nftype',option=opt,index=1))      #Farm types [indice ftype] (id)
        md.nculture         = Set(initialize=read('culture_crop',option=opt).keys()) #Culture -one crop at one time- [indice culture] (id)
        md.nfieldculture    = Set(dimen=2,initialize=read('nfieldculture',option=opt,index=1))#Field [field x culture] [indice field] (id)
    #Spatio-Temporal connections
        md.fzone_catch      = Param(md.nfzone, initialize=read('fzone_catch',option=opt))                #Farming zones's catchment [fzone] (id)
        md.fzone_cmarket    = Param(md.nfzone, initialize=read('fzone_cmarket',option=opt))              #Farming zones's crop market [fzone] (id)
        md.fzone_country    = Param(md.nfzone, initialize=read('fzone_country',option=opt))              #Farming zones's country [fzone] (id)
        md.fzone_dscatch    = Param(md.nfzone, initialize=read('fzone_dscatch',option=opt))              #Return flow recieving catchment of farming zone [fzone] (id)
        md.fzone_type       = Param(md.nfzone, initialize=read('fzone_type',option=opt))                 #Farming type of farming zone [fzone] (id)
        md.field_culture    = Param(md.nfieldculture, initialize=read('field_culture',option=opt))       #Culture of Field [field] (id)
        md.culture_crop     = Param(md.nculture, initialize=read('culture_crop',option=opt))             #Crop produced by culture [culture] (id) 
        md.yphase_month     = Param(md.nculture, md.nyphase, md.nmonth, initialize=read('yphase_month',option=opt)) # Phase month matrix, share of each phase in each month [culture x yphase x time]
    #Farming types
        md.aCulRturn        = Param(md.nftype, md.nculture, initialize=read('aCulRturn',option=opt))     #Water return rate from allocated [ftype x culture] (%) 
        md.aCulYield        = Param(md.nyear, md.nftype, md.nculture, initialize=read('aCulYield',option=opt,time=md.nyear))  #Maximum yield of culture [ftype x culture] (t/ha) 
        md.aCulCost         = Param(md.nftype, md.nculture, initialize=read('aCulCost',option=opt))      #Cultivation cost [ftype x culture] ($/ha)
        md.aIrrigation      = Param(md.nftype, initialize=read('aIrrigation',option=opt))                #Irrigation (1) or rainfed culture(0) [ftype] (binary) 
    #Farming Zones
        md.aLandCap         = Param(md.nyear, md.nfzone, initialize=read('aLandCap',option=opt,time=md.nyear)) #Land capacity for agriculture [year x fzone] (1000ha)
        md.aIrrgCost        = Param(md.nfzone, initialize=read('aIrrgCost',option=opt))            #Irrigation cost [fzone] ($/m3)
        md.aIrrgLoss        = Param(md.nfzone, initialize=read('aIrrgLoss',option=opt))            #Water loss rate from allocated [fzone] (%) 
    #Cultures
        if md.Options['KckY farmtype']==1: #Kc and kY are defined per farm type
            md.akY          = Param(md.nftype, md.nculture, md.nyphase, initialize=read('akY',option=opt))    #crop yield response phase-specific factor [ftype x culture x yphase](%) 
            md.aKc          = Param(md.nftype, md.nculture, md.nmonth, initialize=read('aKc',option=opt))     #crop water demand factor [ftypex culture x yphase](%) 
        else:
            md.akY          = Param(md.nculture, md.nyphase, initialize=read('akY',option=opt))    #crop yield response phase-specific factor [culture x yphase](%) 
            md.aKc          = Param(md.nculture, md.nmonth, initialize=read('aKc',option=opt))     #crop water demand factor [culture x yphase](%) 
        if md.Options['Yield water response'] != 'nonlinear':
            md.aYieldMat    = Param(md.nypath, md.nyphase, initialize=read('aYieldMat',option=opt))#yield response path matrix [ypath x yphase]
        if md.Options['Crop choice'] in ('max','once_max'):
            md.aCulMax      = Param(md.nyear, md.nfzone, md.nculture, initialize=read('aCulMax',option=opt, time=md.nyear)) #Maximum area per culture and farming zone [fzone x culture] (1000ha)
        
        #%%#####################
        #**Crop Market module**#
        ########################

    #Crop markets and demands
        opt=md.Options['Crop market']
        md.ncmarket         = Set(initialize=read('cmarket_country',option=opt).keys())                  #Crop markets [indice cmarket] (id)
        md.cmarket_country  = Param(md.ncmarket, initialize=read('cmarket_country',option=opt))          #Country of Crop Market [cmarket] (id)               
        md.nextcmarket      = Set(initialize=[cm for cm in md.ncmarket if read('cmarket_type')[cm]=='ExternalMarket']) #External crop markets with a crop supply curve (id)
        md.aCropDem         = Param(md.nyear, md.ncmarket, md.ncrop, initialize=read('aCropDem',option=opt,time=md.nyear)) #Crop demand [year x cmarket x crop] (1000t/year)
        md.aCropVal         = Param(md.nyear, md.ncmarket, md.ncrop, initialize=read('aCropVal',option=opt,time=md.nyear))       #Value of crop [cmarket x crop] ($/t)
        md.aMarkMarg        = Param(md.ncmarket, initialize=read('aMarkMarg',option=opt))                #Marketing margin per market [cmarket] (%of crop production costs)
    #Food security
        if md.Options['Minimum supply'] == 1:
            md.aMinDem      = Param(md.ncmarket, md.ncrop, initialize=read('aMinDem',option=opt))    #Crop minimum demand = food security [year x cmarket x crop] (1000t/year)
    #Crop demand own price elasticity
        if opt and md.Options['Crop demand elasticity'] == 'linearized':   #Parametrization of demand curve according to parameters
            #Demand steps - representing elasticity of crop demand
            aStepLen        = md.Options['aStepLen']
            aNumStep        = int(md.Options['aNumStep'])
            aDemEla         = read('aDemEla')
            md.ncdstep      = Set(initialize=range(2*aNumStep+1))         #Share of demand step for crop [indice cdstep] (id)
            aStepDem        = {(cm,cr,cds): max(0,(cds==0)*(1-aStepLen) + (cds!=0)*aStepLen/aNumStep) 
                                for cds in md.ncdstep for cr in md.ncrop for cm in md.ncmarket}
            aStepVal        = {(cm,cr,cds): max(0,1+(cds-aNumStep)*aStepLen/aNumStep/aDemEla[cm,cr])
                                            if aDemEla[cm,cr] != 0 else 1 
                                for cds in md.ncdstep for cr in md.ncrop for cm in md.ncmarket}
            md.aStepDem     = Param(md.ncmarket, md.ncrop, md.ncdstep, initialize=aStepDem)          #Demand share step percentage [cmarket x ncrop x cdstep] (%)
            md.aStepVal     = Param(md.ncmarket, md.ncrop, md.ncdstep, initialize=aStepVal)          #Demand share sclice marginal value coefficient [cmarket x ncrop x cdstep] (-)                                           
        elif opt and md.Options['Crop demand elasticity'] == 'nonlinear': # Value(D) = V0+1/aDemEla*(D-D0)/D0 where D is demand (AcSUPPLY in model)
            md.aDemEla      = Param(md.ncmarket, md.ncrop, initialize=read('aDemEla'))  #Crop demand elasticity to Own-price %Dchange/%PriceCHange
            md.ncdstep      = Set(initialize=[1])
        else:
            md.ncdstep      = Set(initialize=[1])       #Single demand point
            md.aStepDem     = Param(md.ncmarket, md.ncrop, md.ncdstep, 
                                    initialize={(cm,cr,1):1 for cr in md.ncrop for cm in md.ncmarket}) #No elasticity in demand
            md.aStepVal     = Param(md.ncmarket, md.ncrop, md.ncdstep, 
                                    initialize={(cm,cr,1):1 for cr in md.ncrop for cm in md.ncmarket}) #No elasticity in demand
    #Crop transport
        opt= opt and md.Options['Crop transport']
        md.nctrans          = Set(initialize=read('aTransLoss',option=opt).keys())                       #Crop transport roads [indice ctrans] (id)
        md.aTransIn         = Param(md.nctrans, initialize=read('aTransIn',option=opt))                  #Exporting crop market of crop transport road [ctrans] (id) 
        md.aTransOut        = Param(md.nctrans, initialize=read('aTransOut',option=opt))                 #Importing crop market of crop transport road [ctrans] (id)
        md.aTransLoss       = Param(md.nctrans, initialize=read('aTransLoss',option=opt))                #Crop transport loss matrix [ctrans] (%)
        md.aTransCost       = Param(md.nyear, md.nctrans, md.ncrop, initialize=read('aTransCost',option=opt,time=md.nyear))      #Crop transport cost matrix [year x ctrans x crop] ($/t)
        if md.Options['Crop market'] == 0:
            md.aFarmVal     = Param(md.nyear, md.ncountry, md.ncrop, initialize=read('aFarmVal',option=md.Options['Agriculture'],time=md.nyear))       #Value of crop at farm level [country x crop] ($/t)
        #%%###########################
        #**Energy production Module**#
        ##############################
        
        opt=md.Options['Energy production']
    #Hydropower
        md.nhpp             = Set(initialize=read('eHppCap',option=opt).keys())            #HydroPower Plants [indice hpp] (id)  
        md.hp_res           = Param(md.nhpp, initialize=read('hp_res',option=opt))         #Reservoir of HP turbines [hpp], value='ROR' if Run-Of-the-River HP
        md.hp_catch         = Param(md.nhpp, initialize=read('hp_catch',option=opt))       #Catchment of Run-Of-the-River HP turbines [hpp]
        md.hp_pmarket       = Param(md.nhpp, initialize=read('hp_pmarket',option=opt))     #Power pmarket of HP turbines [hpp]
        md.eHppProd         = Param(md.nhpp, initialize=read('eHppProd',option=opt))       #Production factor: energy produced per water flow of HP [hpp] (kWh/m3)
        md.eHppCap          = Param(md.nhpp, initialize=read('eHppCap',option=opt))        #Capacity of HP turbines [hpp] (MW)
        md.eHppEff          = Param(md.nhpp, initialize=read('eHppEff',option=opt))        #Efficiency of HP turbines [hpp] (%)
        md.eHppCost         = Param(md.nhpp, initialize=read('eHppCost',option=opt))       #Operational production cost of HP [hpp] ($/kWh)
    #Power plants and fuels
        opt= opt and md.Options['Power plants']
        md.nopp             = Set(initialize=read('op_pmarket',option=opt).keys())            #Other Power Plants (OPP) [indice opp] (id)
        md.op_pmarket       = Param(md.nopp, initialize=read('op_pmarket',option=opt))     #pmarket of OPP [opp] (id)
        md.eOppCost         = Param(md.nyear, md.nopp, initialize=read('eOppCost',option=opt,time=md.nyear)) #Operational production cost of OPP [opp] ($/kWh)
        md.eOppCap          = Param(md.nyear, md.nopp, initialize=read('eOppCap',option=opt,time=md.nyear))  #Production capacity of OPP [year x opp] (MW)        
        if md.Options['Power technologies'] == 1 or md.Options['Load capacity'] == 1:
            md.op_ptech     = Param(md.nopp, initialize=read('op_ptech',option=opt))       #Technology of OPP [opp] (id)
        if md.Options['Ramping'] == 1:
            md.eOppRamp     = Param(md.nopp, initialize=read('eOppRamp',option=opt))       #Ramping rate of OPP [opp] (%/load segment) 
        if md.Options['Fuels'] ==1:
            md.eOppEff      = Param(md.nopp, initialize=read('eOppEff',option=opt))        #Efficiency of OPP [opp] (kWh_net/kWh_fuel)            
            md.op_fuel      = Param(md.nopp, initialize=read('op_fuel',option=opt))        #Fuel of OPP [opp] (id)                 

#%% **Non linear hydro-power module**
        if md.Options['Hydropower'] == 'nonlinear' and md.Options['Energy production'] == 1:
            md.wResMOL      = Param(md.nres, initialize=read('wResMOL'))
            md.wResFSL      = Param(md.nres, initialize=read('wResFSL'))
            md.wMinTurb     = Param(md.nhpp, initialize=read('wMinTurb'))
            md.wMaxTurb     = Param(md.nhpp, initialize=read('wMaxTurb'))
            md.wMinHead     = Param(md.nres, initialize=read('wMinHead'))
            md.wMaxHead     = Param(md.nres, initialize=read('wMaxHead'))
            md.wFixHead     = Param(md.nres, initialize=read('wFixHead'))
            #Head-Volume relation (linear assumption)
            def _RelHeadVol(Vol,MOL,FSL,hmin,hmax):                
                return hmin/hmax + (1-hmin/hmax)*(Vol-MOL)/(FSL-MOL)
            #Turbine capacity - volume relation (linear assumption)
            def _TurbFlowVol(Vol,MOL,FSL,qmin,qmax):
                return qmin + (qmax-qmin)*(Vol-MOL)/(FSL-MOL)
      
        #%%#######################
        #**Energy market Module**#
        ##########################
        
        opt=md.Options['Energy production'] and md.Options['Energy market']
    #Power markets - Demands and value
        md.npmarket         = Set(initialize=read('pmarket_country',option=opt).keys())             #pmarkets [indice pmarket] (id)            
        md.pmarket_country  = Param(md.npmarket, initialize=read('pmarket_country',option=opt))     #Country of pmarket [pmarket] (id)            
        md.eSupLoss         = Param(md.npmarket, initialize=read('eSupLoss',option=opt))            #Local supply losses of pmarket [pmarket] (-)
        md.eEngyDem         = Param(md.nyear, md.nmonth, md.npmarket, initialize=read('eEngyDem',option=opt,time=md.nyear))  #Local power demand [month x pmarket] (GWh/year)
        md.eEngyVal         = Param(md.nmonth, md.npmarket, initialize=read('eEngyVal',option=opt)) #Curtailment cost (=Marginal value) of power [time x pmarket] ($/kWh)           
    #Transmission lines
        opt= opt and md.Options['Transmission']
        md.ntransline       = Set(initialize=read('eTransCap',option=opt).keys())              #Power transmission lines [indice tmsline] (id)
        md.eTransIn         = Param(md.ntransline, initialize=read('eTransIn',option=opt))     #"Exporting" power market [tmsline] (id)
        md.eTransOut        = Param(md.ntransline, initialize=read('eTransOut',option=opt))    #"Importing" power market [tmsline] (id)
        md.eTransLoss       = Param(md.ntransline, initialize=read('eTransLoss',option=opt))   #Transmition losses [tmsline] (%)
        md.eTransCap        = Param(md.ntransline, initialize=read('eTransCap',option=opt))    #Transmition capacity [tmsline] (MW)
        md.eTransCost       = Param(md.ntransline, initialize=read('eTransCost',option=opt))   #Transmition costs [tmsline] ($/kWh)        
    #Fuels and ressources
        opt= md.Options['Energy production'] and md.Options['Energy market'] and md.Options['Fuels']               
        md.nfuel            = Set(initialize=read('eFuelCO2',option=opt).keys())                   #Fuels [indice fuel] (id)
        md.eFuelCost        = Param(md.nyear, md.npmarket, md.nfuel, initialize=read('eFuelCost',option=opt,time=md.nyear)) #Cost of fuel [fuel] ($/kWh_fuel)
        md.eFuelCO2         = Param(md.nfuel, initialize=read('eFuelCO2',option=opt))              #CO2 emission factor of fuel [fuel] (t/kWh_fuel)
        md.eCO2Val          = Param(md.nyear, md.npmarket, initialize=read('eCO2Val',option=opt,time=md.nyear) 
                                                           if md.Options['Fuels'] == 1 else 
                                                           {(y,pm):0 for y in md.nyear for pm in md.npmarket})  #CO2 price [1] ($/t)
    #Generic capacity investment
        opt= md.Options['Energy production'] and md.Options['Energy market'] and md.Options['Power technologies']             
        md.nptech           = Set(initialize=set([key[0] for key in read('eVarOPEX',option=opt).keys()]))  #Power technologies [indice ptech] (id)
        md.eCAPEX           = Param(md.nyear, md.nptech, md.npmarket, initialize=read('eCAPEX',option=opt,time=md.nyear))#Cost of Generic power capacity [ptech x pmarket] (id) ($/MW)
        md.eVarOPEX         = Param(md.nptech, md.npmarket, initialize=read('eVarOPEX',option=opt))        #Variable Cost of Generic power production [ptech x pmarket] ($/kWh) 
        md.eFixOPEX         = Param(md.nptech, md.npmarket, initialize=read('eFixOPEX',option=opt))        #Fix Cost of Generic power production [ptech x pmarket] ($/[kWh/day] /year)
        md.eLifeTime        = Param(md.nptech, md.npmarket, initialize=read('eLifeTime',option=opt))       #Life time of power technology [ptech x pmarket] (years)
        md.eConstTime       = Param(md.nptech, md.npmarket, initialize=read('eConstTime',option=opt))      #Construction time of power technology [ptech x pmarket] (years)
        md.eMaxCap          = Param(md.nptech, md.npmarket, initialize=read('eMaxCap',option=opt))         #Maximum expandable capacity of power technology [ptech x pmarket] (kWh/day)
    #Fuels
        if md.Options['Fuels'] == 1: 
            md.eTechEff     = Param(md.nptech, md.npmarket, initialize=read('eTechEff',option=opt))     #Efficiency of power technology [ptech x pmarket] (%)
            md.ptech_fuel   = Param(md.nptech, md.npmarket, initialize=read('ptech_fuel',option=opt))   #Fuel of Technology [ptech x pmarket] (id)
        if md.Options['Ramping'] == 1:
            md.eTechRamp    = Param(md.nptech, md.npmarket, initialize=read('eTechRamp',option=opt))    #Ramping rate of Technology [ptech x pmarket] (%/load segment)           
    #Hydropower valuation if power market is off
        if md.Options['Energy market'] == 0:
            md.hp_country   = Param(md.nhpp, initialize=read('hp_country',option=md.Options['Energy production']))  #Country of HP if power market is off [hpp] (id)
            md.eHppVal      = Param(md.nhpp, initialize=read('eHppVal',option=md.Options['Energy production']))     #Value of HP production if power market is off [hpp] ($/kWh)
        else:
            md.hp_country   = Param(md.nhpp, initialize={hp:md.pmarket_country[md.hp_pmarket[hp]] for hp in md.nhpp})
    #Power load segments
        opt= md.Options['Energy production'] and md.Options['Energy market']
        md.npload           = Set(initialize=read('eLoadDem',option=opt).keys() if md.Options['Load'] == 1 or opt == 1 else [1])            #Slices of power demand per time: day week, day week end, night [indice pload] (id)
        md.eLoadDem         = Param(md.npload, initialize=read('eLoadDem',option=opt) if md.Options['Load'] == 1 or opt == 1 else {1:1})    #Share of demand per power load [pload] (-)
        md.eLoadTime        = Param(md.npload, initialize=read('eLoadTime',option=opt) if md.Options['Load'] == 1 or opt == 1 else {1:1})   #Share of time per power load [pload] (-)
        if md.Options['Load capacity'] == 1: #Assumes either "Other power plants" or "Power technologies" are ON
            if md.Options['Power technologies'] == 0:
                md.op_ptech = Param(md.nopp, initialize=read('op_ptech',option=opt))                      #Technology of Other power plant [pp] (id)
                md.nptech   = Set(initialize=set([read('op_ptech',option=opt)[op] for op in md.nopp]))    #Technologies are defined in Other power plants if Power technologies module is off [indice ptech] (id)
            md.eLoadCap     = Param(md.npload, md.nptech, md.nmonth, initialize=read('eLoadCap',option=opt))         #Load segment capacity factor of power technology [pload x ptech] (%)

        #%%#######################
            #**Activities**#
        ########################## 
        opt= md.Options['Activities'] 
        md.njactivity       = Set(initialize=read('jProdCap',option=opt).keys())
        md.j_country        = Param(md.njactivity,initialize=read('j_country',option=opt)) #country of activity
        md.j_catch          = Param(md.njactivity,initialize=read('j_catch',option=opt)) #catchment of activity (if relevant)
        md.j_pmarket        = Param(md.njactivity,initialize=read('j_pmarket',option=opt)) #power market of activity (if relevant)
        md.j_fzone          = Param(md.njactivity,initialize=read('j_fzone',option=opt)) #farming zone of activity (if relevant)
        md.j_cmarket        = Param(md.njactivity,initialize=read('j_cmarket',option=opt)) #crop market of activity (if relevant)
        md.j_cropin         = Param(md.njactivity,initialize=read('j_cropin',option=opt)) #crop as input to activity (if relevant)
        md.j_cropout        = Param(md.njactivity,initialize=read('j_cropout',option=opt)) #crop as output to activity (if relevant)
        md.jProdCap         = Param(md.njactivity,initialize=read('jProdCap',option=opt)) # Production capacity (units/month)
        md.jProdCost        = Param(md.njactivity,initialize=read('jProdCost',option=opt)) # Production costs (M$/month)
        #md.jProdVal         = Param(md.njactivity,initialize=read('jProdVal',option=opt)) # Production value (M$/month)
        md.jLandCons        = Param(md.njactivity,initialize=read('jLandCons',option=opt)) # Land consumption (1000 ha/month)
        md.jLandProd        = Param(md.njactivity,initialize=read('jLandProd',option=opt)) # Land production (1000 ha/month)
        md.jWatCons         = Param(md.njactivity,initialize=read('jWatCons',option=opt)) # Water consumptio (Mm3/month)
        md.jWatProd         = Param(md.njactivity,initialize=read('jWatProd',option=opt)) # Water production (Mm3/month)
        md.jPowCons         = Param(md.njactivity,initialize=read('jPowCons',option=opt)) # Power consumption (GWh/month)
        md.jPowProd         = Param(md.njactivity,initialize=read('jPowProd',option=opt)) # Power production (Gwh/month)
        md.jCropCons        = Param(md.njactivity,initialize=read('jCropCons',option=opt)) # Crop consumption (1000t/month)
        md.jCropProd        = Param(md.njactivity,initialize=read('jCropProd',option=opt)) # Crop production (1000t/month)
        
#%%DECLARE DECISION VARIABLES###       
    #Water Module#
        md.WwSUPPLY         = Var(md.ntime, md.nuser, within=NonNegativeReals)    #Water user allocation  [time x user] (Mm3)
        md.WwOUTFLOW        = Var(md.ntime, md.ncatch, within=NonNegativeReals)   #Catchment outflow [time x catchment] (km3)
        if ANALYTICAL == 1:
            md.WwRSTORAGE   = Var(md.ntime, md.nres, within=NonNegativeReals)     #Reservoir total storage [time x reservoir] (km3)            
        else:
            md.WwRSTORAGE   = Var(md.ntime, md.nres, bounds=(1e-20,None), initialize=1e-20)     #Reservoir total storage [time x reservoir] (km3)            
        md.WwGWSTORAGE      = Var(md.ntime, md.naquifer, within=NonNegativeReals)  #Groundwater reservoir total storage [time x aquifer] (Mm3)
        md.WwTRANSFER       = Var(md.ntime, md.ntransfer, within=NonNegativeReals) #Transfer scheme water transfer [time x transfer] (Mm3)                 
    
    #Agriculture Module#
        if md.Options['Crop choice'] == 'fixed': #Cultivated area is fixed (=parameter)
            md.CULAREA      = Param(md.nyear, md.nfzone, md.nculture, mutable=True, 
                                    initialize=read('aCulMax',option=md.Options['Agriculture'],time=md.nyear)) #Fixed cultivated area [year x fzone x culture] (1000ha)
            _xCULAREA       = lambda md,y,fz,cul : md.CULAREA[y,fz,cul]
            _xPHYAREA       = lambda md,y,fz : sum(md.CULAREA[y,fz,kcul] for kcul in md.nculture)
        elif md.Options['Yield water response'] == 'nonlinear':
            md.AlCULAREA    = Var(md.nyear, md.nfzone, md.nfieldculture, within=NonNegativeReals)   #Cultivated area  [year x fzone x field] (1000ha)
            _xCULAREA       = lambda md,y,fz,cul : sum(md.AlCULAREA[y,fz,kfd] for kfd in md.nfieldculture if md.field_culture[kfd]==cul)
            _xPHYAREA       = lambda md,y,fz : sum(md.AlCULAREA[y,fz,kfd] for kfd in md.nfieldculture if kfd[1]==1)
        else: #Cultivated area is linearized
            md.AlCULAREA    = Var(md.nyear, md.nfzone, md.nfieldculture, md.nypath, within=NonNegativeReals)    #Cultivated area and water demand satisfaction [year x fzone x field x ypath] (1000ha)
            _xCULAREA       = lambda md,y,fz,cul : sum(md.AlCULAREA[y,fz,kfd,kypt] for kfd in md.nfieldculture for kypt in md.nypath if md.field_culture[kfd]==cul)
            _xPHYAREA       = lambda md,y,fz : sum(md.AlCULAREA[y,fz,kfd,kypt] for kfd in md.nfieldculture for kypt in md.nypath if kfd[1]==1)
        md.xCULAREA         = Expression(md.nyear, md.nfzone, md.nculture, rule=_xCULAREA)          #Agriculture (harvested) area [year x fzone x culture] (kha)
        md.xPHYAREA         = Expression(md.nyear, md.nfzone, rule=_xPHYAREA)                       #Agriculture physical area [year x fzone] (kha)
        md.AwSUPPLY         = Var(md.ntime, md.nifzone, md.nculture, within=NonNegativeReals)       #Agricultural water allocation [time x ifzone x culture] (Mm3)
        md.AcPROD           = Var(md.nyear, md.nfzone, md.ncrop, within=NonNegativeReals)           #Crop Production [year x fzone x crop] (kt)       
        if md.Options['Groundwater'] == 1:
            md.AwGWSUPPLY   = Var(md.ntime, md.nifzone, md.nculture, within=NonNegativeReals)       #Agricultural ground water allocation [time x ifzone x culture] (Mm3)
            
    #Crop market module#
        md.AcTRANS          = Var(md.nyear, md.nctrans, md.ncrop, within=NonNegativeReals)              #Crop Transport [year x ctrans x crop] (kt)
        if md.Options['Crop demand elasticity'] != 'nonlinear':
            md.AcSUPPLY     = Var(md.nyear, md.ncmarket, md.ncrop, md.ncdstep, within=NonNegativeReals) #Crop Supply [year x cmarket x crop x cdstep] (kt) 
        else:
            md.AcSUPPLY     = Var(md.nyear, md.ncmarket, md.ncrop, within=NonNegativeReals)             #Crop Supply [year x cmarket x crop] (kt)     
        md.AcEXTPROD        = Var(md.nyear, md.nextcmarket, md.ncrop, within=NonNegativeReals)          #Crop external Production [year x cmarket x crop] (kt)
        
    #Energy production module#
        md.EwHPDISCHARGE    = Var(md.ntime, md.npload, md.nhpp, within=NonNegativeReals)    #Dicharge through hydropower turbine [time x pload x hpp] (Mm3)
        md.EeHPPROD         = Var(md.ntime, md.npload, md.nhpp, within=NonNegativeReals)    #Power production by hydropower turbine [time x pload x hpp] (GWh)
        md.EeOPPROD         = Var(md.ntime, md.npload, md.nopp, within=NonNegativeReals)    #Power production by power plant [time x pload x pp] (GWh)        
        
    #Energy market module#
        md.EeTRANS          = Var(md.ntime, md.npload, md.ntransline, within=NonNegativeReals)          #Power transmission between power markets [time x pload x tmsline] (GWh)
        md.EeSUPPLY         = Var(md.ntime, md.npload, md.npmarket, within=NonNegativeReals)            #Power supply to power market demand [time x pload x pmarket] (GWh)
        md.EeGENCAP         = Var(md.nyear, md.nptech, md.npmarket, within=NonNegativeReals)            #Generic Power technologies capacity [year x ptech x pmarket] (MW)
        md.EeGENPROD        = Var(md.ntime, md.npload, md.nptech, md.npmarket, within=NonNegativeReals) #Generic Power Tech Production [time x pload x ptech x pmarket] (GWh/load seg)            
    
    #Activities
        md.JjPROD           = Var(md.ntime,md.njactivity,within=NonNegativeReals) #Activities (units/month)
#%%DEBUG VARIABLES        
        if md.Options['Initial time step'] == 1:
            md.DUMYSTOR     = Var(md.nres, within=NonNegativeReals) #Dummy variable to avoid infeasible solutions regarding final storage objective
        if md.Options['Debug mode'] == 1:
            md.DUMYWATER    = Var(md.ntime, md.ncatch, within=NonNegativeReals) #DEBUG Variable adding water in the water balance (m3/month)
        if md.Options['Debug mode'] in [1,2]:
            md.DUMYCROP     = Var(md.nyear, md.nfzone, md.ncrop, within=NonNegativeReals)#, bounds=(0,1))   #Dummy Crop Production to compensate for negative prod [year x fzone x crop] (kt/year)       
            md.DUMYEFLOW    = Var(md.ntime, md.neflow, within=NonNegativeReals) #Tolerance to eflow objective if physically impossible (Mm3/month)

#%%INVESTMENT PLANNING MODULE
        #Discount rate (can also be used without investment planning module)
        if md.Options['Discount rate']==md.Options['Discount rate']: #option is defined
            DiscountFactor  = {y: (1-md.Options['Discount rate']/100)**(y-md.t_year[md.Options['tini']]) for y in md.nyear}
        else:
            DiscountFactor  = {y: 1 for y in md.nyear} #no discounting
        md.iDisFact         = Param(md.nyear, initialize=DiscountFactor)                     #Discounting factor [year] (-)
        
        if md.Options['Investment module'] in [1,'continuous']:                        
        #Declare indices
            md.ninvest      = Set(initialize=read('iCAPEX').keys())                          #Investments [indice inv] (id)
            md.ninvphase    = Set(initialize=read('invphase_t').keys())                      #Investment phases [indice invphase] (id)            
        #Declare parameters
            md.invphase_t   = Param(md.ninvphase, initialize=read('invphase_t'))              #time step of investment phase [invphase] (time)
            md.inv_after    = Param(md.ninvest, initialize=read('inv_after'))                 #Investment can occur only after another investment [inv] (id)
            md.inv_group    = Param(md.ninvest, initialize=read('inv_group'))                 #Investment group have to occur at the same time [inv] (id) #ADD: only works with 2 investments so far
            md.inv_replace  = Param(md.ninvest, initialize=read('inv_replace'))               #Investment replaces existing infrastructure [inv] (id) #ADD only active for hydropower now
            
            md.iMaxInv      = Param(md.ninvphase, initialize=read('iMaxInv'))                 #Maximum investment amount of investment phase [invphase] ($)
            md.iCAPEX       = Param(md.ninvest, initialize=read('iCAPEX'))                    #investment capital costs [inv] ($)
            md.iFixOPEX     = Param(md.ninvest, initialize=read('iFixOPEX'))                  #investment fix operational costs [inv] ($)
            md.iInvType     = Param(md.ninvest, initialize=read('iInvType'))                  #type of investment [invtype] (type)
            md.iInvCap      = Param(md.ninvest, initialize=read('iInvCap'))                   #capacity of investment [inv] (different units)
            md.iConstTime   = Param(md.ninvest, initialize=read('iConstTime'))                #implementation time of investment [inv] (month)
            md.iLifeTime    = Param(md.ninvest, initialize=read('iLifeTime'))                 #life time of investment after construction [inv] (month)
        
        #Declare Decision-Variables
            if md.Options['Investment module']==1:
                md.IbINVEST = Var(md.ninvphase, md.ninvest, within=Binary)                    #Investments [invphase x inv] (binary)        
            if md.Options['Investment module']=='continuous':
                md.IbINVEST = Var(md.ninvphase, md.ninvest, within=NonNegativeReals,bounds=(0,1.4))  #Investments [invphase x inv] (binary) 
                #ADD: bounds to continuous investment as binary
        ##Functions##
        def _InvestedCapacity(currentcap,iobject,itype,timestep,year=0):
            if md.Options['Investment module'] in [1,'continuous']: #Investment module is on
                if year==0:  #at monthly time step
                    InvestedCap=sum(md.iInvCap[inv]*md.IbINVEST[kip,inv] for kip in md.ninvphase for inv in md.ninvest 
                                    if  md.iInvType[inv]==itype and inv==iobject
                                    and md.invphase_t[kip]+md.iConstTime[inv] <= timestep 
                                    and md.invphase_t[kip]+md.iConstTime[inv]+md.iLifeTime[inv] >= timestep)                                       
                    ReplacedCap=sum(min(currentcap,md.iInvCap[kinv])*md.IbINVEST[kip,kinv] for kip in md.ninvphase for kinv in md.ninvest 
                                    if (md.iInvType[kinv]==itype and md.inv_replace[kinv]==iobject)
                                    and md.invphase_t[kip]+md.iConstTime[kinv] <= timestep
                                    and md.invphase_t[kip]+md.iConstTime[kinv]+md.iLifeTime[kinv] >= timestep)
                else:        #at yearly time step
                    InvestedCap=sum(md.iInvCap[inv]*md.IbINVEST[kip,inv] for kip in md.ninvphase for inv in md.ninvest 
                                    if md.iInvType[inv]==itype and inv==iobject
                                    and (md.invphase_t[kip]+md.iConstTime[inv] <= md.Options['tfin'] 
                                         and md.t_year[md.invphase_t[kip]+md.iConstTime[inv]] <= timestep) 
                                    and (md.invphase_t[kip]+md.iConstTime[inv]+md.iLifeTime[inv] >= md.Options['tfin'] 
                                         or md.t_year[md.invphase_t[kip]+md.iConstTime[inv]+md.iLifeTime[inv]] >= timestep))                                   
                    ReplacedCap=sum(min(currentcap,md.iInvCap[kinv])*md.IbINVEST[kip,kinv] for kip in md.ninvphase for kinv in md.ninvest 
                                    if  (md.iInvType[kinv]==itype and md.inv_replace[kinv]==iobject)
                                    and (md.invphase_t[kip]+md.iConstTime[kinv] <= md.Options['tfin'] 
                                         and md.t_year[md.invphase_t[kip]+md.iConstTime[kinv]] <= timestep)                                                                                                                     
                                    and (md.invphase_t[kip]+md.iConstTime[kinv]+md.iLifeTime[kinv] >= md.Options['tfin'] 
                                         or md.t_year[md.invphase_t[kip]+md.iConstTime[kinv]+md.iLifeTime[kinv]] >= timestep))
                return InvestedCap-ReplacedCap
            else:
                return 0                

#%%OBJECTIVE FUNCTION - ECONOMIC MODULE
        _ntime = lambda y : [kt for kt in md.ntime if md.t_year[kt]==y]
        #_nuser = lambda co : [ku for ku in md.nuser if md.user_country[ku]==co]
        #_nfzone = lambda co : [kfz for kfz in md.nfzone if md.fzone_country[kfz]==co]
        #_nifzone = lambda co : [kfz for kfz in md.nifzone if md.fzone_country[kfz]==co]
        
        #Intermediate expressions
        #Water users
        md.xUserBen=Expression(md.nyear,md.nuser, #User water supply benefits [M$]
                rule=lambda md,y,u:sum(md.WwSUPPLY[t,u] * md.wUserVal[u] for t in _ntime(y)))
        md.xUserSupCost=Expression(md.nyear,md.nuser, #User supply costs [M$]
                rule=lambda md,y,u:sum(md.WwSUPPLY[t,u]/(1-md.wUserLoss[u]) * md.wSupCost[u] for t in _ntime(y)))
        #Farming zones
        md.xFzCulCost=Expression(md.nyear,md.nfzone, #Cultivation costs [M$]
                rule=lambda md,y,fz:sum(md.kha_to_Mha*md.xCULAREA[y,fz,kcul]*md.aCulCost[md.fzone_type[fz],kcul] 
                                        for kcul in md.nculture))
        md.xFzIrrCost=Expression(md.nyear,md.nifzone, #Irrigation costs [M$]
                rule=lambda md,y,fz:sum(md.AwSUPPLY[t,fz,kcul]/(1-md.aIrrgLoss[fz])*md.aIrrgCost[fz] 
                                        for t in _ntime(y) for kcul in md.nculture))
        md.xFzPumpCost=Expression(md.nyear,md.nifzone, #Groundwater pumping costs [M$]
                rule=lambda md,y,fz:sum(md.AwGWSUPPLY[t,fz,cul]*md.wGwCost[aq] 
                                        for t in _ntime(y) for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]==md.aqui_catch[aq]))
        md.xFzBen=Expression(md.nyear,md.nfzone, #Farm level benefits [M$] (=0 if crop markets are activated)
                rule= lambda md,y,fz:(sum(md.kt_to_Mt*md.AcPROD[y,fz,kcr] * md.aFarmVal[y,md.fzone_country[fz],kcr] 
                                          for kcr in md.ncrop) 
                                     if md.Options['Crop market'] == 0 else 0))
        
        #Crop markets
        md.xCmMarkMarg =Expression(md.nyear,md.ncmarket, #Marketing margin costs [M$]
                rule=lambda md,y,cm:sum(md.aMarkMarg[cm]*md.xFzCulCost[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm)
                                    +sum(md.aMarkMarg[cm]*(md.xFzIrrCost[y,fz]+md.xFzPumpCost[y,fz]) for fz in md.nifzone if md.fzone_cmarket[fz]==cm))
        md.xCmTransCost=Expression(md.nyear,md.ncmarket, #Crop transport costs [M$]
                rule=lambda md,y,cm:sum(md.kt_to_Mt*md.AcTRANS[y,ct,kcr] * md.aTransCost[y,ct,kcr] 
                                        for ct in md.nctrans for kcr in md.ncrop if md.aTransOut[ct]==cm))
        md.xCmProdCost=Expression(md.nyear,md.ncmarket, #Crop external production costs [M$]
                rule=lambda md,y,cm:(sum(md.kt_to_Mt*md.AcEXTPROD[y,cm,kcr] * md.aCropVal[y,cm,kcr] * 1.0001 
                                         for kcr in md.ncrop) if cm in md.nextcmarket else 0))
        md.xCmBen=Expression(md.nyear,md.ncmarket, #Crop supply benefits [M$]
                rule=lambda md,y,cm:(sum(md.kt_to_Mt*md.AcSUPPLY[y,cm,kcr,kcds] * md.aCropVal[y,cm,kcr] * md.aStepVal[cm,kcr,kcds]
                                         for kcr in md.ncrop for kcds in md.ncdstep)
                                    if md.Options['Crop demand elasticity']!='nonlinear' else 
                                     sum(md.aCropVal[y,cm,kcr]
                                         *(md.kt_to_Mt*md.AcSUPPLY[y,cm,kcr]*(1-(1/md.aDemEla[cm,kcr] if md.aDemEla[cm,kcr] != 0 else 0))
                                           +(1/md.aDemEla[cm,kcr]*1/2*md.kt_to_Mt*md.AcSUPPLY[y,cm,kcr]**2/md.aCropDem[y,cm,kcr] #scaling factor is not **2 because also applies to aCropDem
                                             if md.aDemEla[cm,kcr] != 0 else 0) 
                                           if md.aCropDem[y,cm,kcr] !=0 else 0)
                                         for kcr in md.ncrop)))
        #Economic Balance
        def obj_rule(md):
        #Water module: Supply benefits, supply costs            
            UserBalance         = {y:sum(md.xUserSupCost[y,u]-md.xUserBen[y,u] for u in md.nuser) for y in md.nyear}            
        #Agricultural module: Cultivation costs of cultures (labour, machinery, fertilizers, infrastructure....), Irrigation costs of cultures (pumping ...)     
            AgricultureBalance  = {y:sum(md.xFzCulCost[y,fz] - md.xFzBen[y,fz] for fz in md.nfzone)
                                    +sum(md.xFzIrrCost[y,fz] + md.xFzPumpCost[y,fz] for fz in md.nifzone) for y in md.nyear}        
        #Crop market module: Crop allocated to crop markets, Cost and benefits of imported crops (external market), Transport costs          
            CropMarketBalance   = {y:sum(md.xCmMarkMarg[y,cm] + md.xCmTransCost[y,cm] + md.xCmProdCost[y,cm] 
                                         - md.xCmBen[y,cm] for cm in md.ncmarket) for y in md.nyear}
      
        #Energy production module: Production costs of hydropower and power plants (O&M, fuel ...)    
            HpOMCost            = {y:sum(md.eHppCost[hp]*md.se*md.EeHPPROD[t,pld,hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}            
            OppOMCost           = {y:sum(md.eOppCost[y,opp]*md.se*md.EeOPPROD[t,pld,opp] for t in md.ntime for pld in md.npload for opp in md.nopp if md.t_year[t]==y) for y in md.nyear}    
            OppFuelCost         = {y:sum(md.se*md.EeOPPROD[t,pld,opp]/md.eOppEff[opp]*md.eFuelCost[y,md.op_pmarket[opp],fu] for t in md.ntime for pld in md.npload for opp in md.nopp for fu in md.nfuel if md.t_year[t]==y and md.op_fuel[opp]==fu) for y in md.nyear}
            OppCO2Cost          = {y:sum(md.se*md.EeOPPROD[t,pld,opp]/md.eOppEff[opp]*md.eFuelCO2[fu]*md.eCO2Val[y,md.op_pmarket[opp]] for t in md.ntime for pld in md.npload for opp in md.nopp for fu in md.nfuel if md.t_year[t]==y and md.op_fuel[opp]==fu) for y in md.nyear}
            GenCapCost          = {y:sum(md.se*md.EeGENCAP[y,pt,pm]*md.eCAPEX[y,pt,pm]*min(1,(md.t_year[md.Options['tfin']]-y+1)/md.eLifeTime[pt,pm])+sum(md.se*md.EeGENCAP[ky,pt,pm] for ky in md.nyear if ky <= y and ky >= y-md.eLifeTime[pt,pm])*md.eFixOPEX[pt,pm] for pt in md.nptech for pm in md.npmarket) for y in md.nyear}
            GenProdCost         = {y:sum(md.se*md.EeGENPROD[t,pld,pt,pm]*md.eVarOPEX[pt,pm] for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear}
            GenFuelCost         = {y:sum(md.se*md.EeGENPROD[t,pld,pt,pm]/md.eTechEff[pt,pm]*md.eFuelCost[y,pm,fu]             for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.t_year[t]==y and md.ptech_fuel[pt,pm]==fu) for y in md.nyear}
            GenCO2Cost          = {y:sum(md.se*md.EeGENPROD[t,pld,pt,pm]/md.eTechEff[pt,pm]*md.eFuelCO2[fu]*md.eCO2Val[y,pm]  for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.t_year[t]==y and md.ptech_fuel[pt,pm]==fu) for y in md.nyear}
            EnergyProdBalance   = {y:HpOMCost[y] + OppOMCost[y] + OppFuelCost[y] + OppCO2Cost[y] + GenCapCost[y] + GenProdCost[y] + GenFuelCost[y] + GenCO2Cost[y] for y in md.nyear}           
                
        #Energy market module: Energy distributed, Energy transmission costs, Energy curtailment costs    
            EngyBenefit         = {y:sum(md.se*md.EeSUPPLY[t,pld,pm] * md.eEngyVal[md.t_month[t],pm] for t in md.ntime for pld in md.npload for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear}
            if md.Options['Energy market'] == 0:
                EngyBenefit     = {y:sum(md.se*md.EeHPPROD[t,pld,hp]*md.eHppVal[hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}
            EngyTransCost       = {y:sum(md.se*md.EeTRANS[t,pld,tl] * md.eTransCost[tl]  for t in md.ntime for pld in md.npload for tl in md.ntransline if md.t_year[t]==y) for y in md.nyear}
            EnergyMarketBalance = {y:EngyTransCost[y]-EngyBenefit[y] for y in md.nyear} #            
        #Activities
            ActivityBalance     = {y:sum(md.JjPROD[t,j]*md.jProdCost[j] for t in _ntime(y) for j in md.njactivity) for y in md.nyear}
            #ActBenefit          = {y:sum(md.JjPROD[t,j]*md.jProdVal[j] for j in md.njactivity for t in _ntime(y)) for y in md.nyear}
        #DEBUG MODE
            DebugCost=0
            if md.Options['Debug mode'] == 1:
                DebugCost += sum(md.DUMYWATER[t,c]*md.debugcost for t in md.ntime for c in md.ncatch)
            if md.Options['Debug mode'] in [1,2]:
                DebugCost += sum(md.kt_to_Mt*md.DUMYCROP[y,fz,cr]*md.debugcost*100 for y in md.nyear for fz in md.nfzone for cr in md.ncrop)
                DebugCost += sum(md.DUMYEFLOW[t,ef]*md.debugcost/2 for t in md.ntime for ef in md.neflow)
                #DebugCost += -sum(md.sw*md.WwRSTORAGE[t,res]*md.debugcost/10**10 for t in md.ntime for res in md.nres if t == md.Options['tfin'])
            if md.Options['Initial time step'] == 1:
                DebugCost += sum(md.DUMYSTOR[res]*md.debugcost/3 for res in md.nres)
            if md.Options['Crop choice'] == 'force': #ADD: assumes full supply path is path 1 NOT FUNCTIONING SO FAR
                DebugCost += sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt]*md.debugcost for y in md.nyear for fz in md.nfzone for fd in md.nfieldculture for ypt in md.nypath if ypt != 1)
        #Investment module               
            if md.Options['Investment module'] in [1,'continuous']:
                #Investment costs + discounting - #ADD in this case the remaining value at the end is proportional to the use of the infrastructure.
                CAPEX   = {y:sum(md.IbINVEST[ip,inv]*md.iCAPEX[inv] for ip in md.ninvphase for inv in md.ninvest if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear}                   
                CAPEXrem= sum(sum(md.IbINVEST[ip,inv]*md.iCAPEX[inv]
                                  *max(0,1-(md.Options['tfin']-md.invphase_t[ip]-md.iConstTime[inv])/md.iLifeTime[inv]) 
                                  for ip in md.ninvphase for inv in md.ninvest if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear)
                CAPEX[md.t_year[md.Options['tfin']]] += -CAPEXrem #remaining capital value at end of simulation period
                fixOPEX = {y:sum(_InvestedCapacity(0,inv,md.iInvType[inv],y,year=1)/md.iInvCap[inv]*md.iFixOPEX[inv] for ip in md.ninvphase for inv in md.ninvest if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear}
                InvestmentBalance={y:CAPEX[y]+fixOPEX[y] for y in md.nyear}
            else:
                InvestmentBalance={y:0 for y in md.nyear}
        #Final Balance
            #Twisted objective coefficient - prioritizing one sector over others
            ObjCoef=[1,1,1]    
            if 'Objective_coef' in md.Options.keys() and md.Options['Objective_coef']==md.Options['Objective_coef']:
                ObjCoef=[float(k) for k in md.Options['Objective_coef'].split('#')] #in the order: water, power, agriculture
            
            return DebugCost + sum(md.iDisFact[y]*(+ InvestmentBalance[y]
                                                   +  UserBalance[y] * ObjCoef[0]                                                   
                                                   + (EnergyProdBalance[y] + EnergyMarketBalance[y]) * ObjCoef[1]
                                                   + (AgricultureBalance[y] + CropMarketBalance[y]) * ObjCoef[2]
                                                   +  ActivityBalance[y])
                                       for y in md.nyear)
        md.obj = Objective(rule=obj_rule)
        
#%%CONSTRAINTS
#%%-----------------------------Investment planning module---------------------
        ##Constraint catalogue##
        def invest_once(md,ninv):
            return sum(md.IbINVEST[kip,ninv] for kip in md.ninvphase) <= 1
       
        def invest_max(md,nip):
            Budget=sum(md.iMaxInv[kip] for kip in md.ninvphase if md.invphase_t[kip] <= md.invphase_t[nip])
            Invest=sum(md.IbINVEST[kip,kinv]*md.iCAPEX[kinv] for kinv in md.ninvest for kip in md.ninvphase 
                       if md.invphase_t[kip] <= md.invphase_t[nip])
            return Invest<=Budget
              
        def invest_after(md,nip,ninv):
            if md.inv_after[ninv] not in md.ninvest:
                return Constraint.Skip
            else:
                return md.IbINVEST[nip,ninv] <= sum(md.IbINVEST[kip,md.inv_after[ninv]] for kip in md.ninvphase if md.invphase_t[kip] <= md.invphase_t[nip])
        
        def invest_group(md,nip,ninv):
            if md.inv_group[ninv] not in md:
                return Constraint.Skip
            else:
                return md.IbINVEST[nip,ninv] == sum(md.IbINVEST[nip,kinv] for kinv in md.ninvest if md.inv_group[kinv]==md.inv_group[ninv])
        
        ##Create constraints##
        if md.Options['Investment module'] in [1,'continuous']:
            md.invest_max   = Constraint(md.ninvphase,              rule=invest_max)
            md.invest_after = Constraint(md.ninvphase, md.ninvest,  rule=invest_after)
        if md.Options['Investment module'] == 1:    
            md.invest_once  = Constraint(md.ninvest,                rule=invest_once)
#%%----------------------------Water module----------------------------
#        def WaterBalance(elem,nt,nc):
#            if elem=='NetInflow': #Net incoming flow from upstream catchments (Mm3)
#                return (1-md.wFlowLoss[nc]) * sum(md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc)
                     
        ##Functions
        def _storage(t,r,mode='average'): #Average or relative storage from t-1 to t
            S_end=md.sw*md.WwRSTORAGE[t,r]
            S_ini=md.wStorIni[r] if md.Options['Initial time step'] == 1 and md.Options['tini'] == t else md.sw*md.WwRSTORAGE[md.t_prev[t],r]
            if mode=='average':
                return (S_end+S_ini)/2  
            elif mode=='delta': 
                return S_end-S_ini
            elif mode=='previous':
                return S_ini
        _nres = lambda c : [kres for kres in md.nres if md.res_catch[kres]==c] #Reservoirs of a catchment
        _nlake = lambda c : [kla for kla in md.nlake if md.res_catch[kla]==c] #Lakes of a catchment
        _nifzone = lambda c : [kfz for kfz in md.nifzone if md.fzone_catch[kfz]==c] #Irrigated farming zones of a catchment
        _naquifer = lambda c : [kaq for kaq in md.naquifer if md.aqui_catch[kaq]==c] #Aquifers of a catchment
        
        ##Expressions##
        def xResEvap(md,nt,nres): #Reservoir evaporation
            nc=md.res_catch[nres]
            ET0=max(0,value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc]))*md.mm_to_m
            if md.Options['Lakes'] == 1 and md.res_type[nres] == 'wetland':
                if ANALYTICAL ==1:
                    Inflow=(1-md.wFlowLoss[nc])*sum((1-md.wkUP[nres])*md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc)
                    a=ET0*md.wkV[nres]+md.wAlpha[nres] #REM: equation is not valid if a=0 !
                    b=Inflow+(1-md.wkRO[nres])*md.wRunOff[nt,nc]-md.wResArea[nres]*ET0
                    return ET0*(md.wkV[nres]/a*((_storage(nt,nres,mode='previous')-b/a)*(1-exp(-a))+b)+md.wResArea[nres])
                else:
                    return ET0*(md.wkV[nres]*_storage(nt,nres,mode='average')**AREAPOWER+md.wResArea[nres])                    
            else:
                return ET0*(md.wkV[nres]*_storage(nt,nres,mode='average')+md.wResArea[nres])
        md.xResEvap=Expression(md.ntime,md.nres,rule=xResEvap)
        ##Constraint catalogue##
        def water_waterbalance(md,nt,nc): #Water Balance per Watershed [catchment x time]
            ActProd             = sum(md.JjPROD[nt,kj]*md.jWatProd[kj] for kj in md.njactivity if md.j_catch[kj]==nc)
            ActCons             = sum(md.JjPROD[nt,kj]*md.jWatCons[kj] for kj in md.njactivity if md.j_catch[kj]==nc)
            NetInflow           = (1-md.wFlowLoss[nc]) * sum(md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) #Net incoming flow from upstream catchments (Mm3)
            RunOff              = md.wRunOff[nt,nc]   #Catchment Runoff (Mm3) 
            GwStorage           = sum(md.WwGWSTORAGE[md.t_prev[nt],kaq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[kaq] for kaq in _naquifer(nc)) #Previous ground water storage
            Pumping             = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in _nifzone(nc) for kcul in md.nculture) if md.Options['Groundwater']==1 else 0 #Water abstraction from groundwater reservoir
            BaseFlow            = sum(GwStorage*(1-exp(-md.wGwFlow[kaq])) + (md.wGwRech[nt,kaq]-Pumping)*(1-(1-exp(-md.wGwFlow[kaq]))/md.wGwFlow[kaq]) for kaq in _naquifer(nc)) #Baseflow from groundwater reservoir
            UserCons            = sum(md.WwSUPPLY[nt,ku] * (1/(1-md.wUserLoss[ku])-md.wUserRturn[ku]) for ku in md.nuser if md.user_catch[ku] == nc)      #Consumption of non agricultural users (Mm3)
            AgrCons             = sum(md.AwSUPPLY[nt,kfz,kcu] * (1/(1-md.aIrrgLoss[kfz])-md.aCulRturn[md.fzone_type[kfz],kcu]) for kfz in _nifzone(nc) for kcu in md.nculture)  #Agricultural module: Consumptiom to agriculture/crops (Mm3)                
            DeltaStorage        = sum(_storage(nt,kres,mode='delta') for kres in _nres(nc)) #Storage in reservoir (Mm3)
            Evaporation         = sum(md.xResEvap[nt,kres] for kres in _nres(nc)) #Evaporation from reservoir (Mm3)
            TransNetInflow      = sum(md.WwTRANSFER[nt,ktrans] * (1-md.wTransLoss[ktrans])   for ktrans in md.ntransfer if md.transfer_ds[ktrans] == nc)  #Incoming flow from upstream transfer scheme- if exists (Mm3)
            TransGrossOutflow   = sum(md.WwTRANSFER[nt,ktrans]                               for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc)  #Allocation to downstream transfer scheme- if exists (Mm3)                   
            CatchmentOutflow    = md.sw*md.WwOUTFLOW[nt,nc]  #Catchment outflow : Allocation to downstream catchment (Mm3)           
            DEBUG=md.DUMYWATER[nt,nc] if md.Options['Debug mode'] == 1 else 0
            return UserCons + AgrCons + TransGrossOutflow + DeltaStorage + Evaporation + CatchmentOutflow + ActCons == NetInflow + TransNetInflow + RunOff + BaseFlow + ActProd + DEBUG
                
        def water_waterbalance2(md,nt,nc): #Avoids allocation from downstream reservoir to upstream demand [catchment x time] 
            if  sum(1 for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc) + sum(1 for ku in md.nuser if md.user_catch[ku] == nc) + sum(1 for kfz in md.nfzone if md.fzone_catch[kfz]==nc) == 0:
                return Constraint.Skip #if there is no demand and no transfer project constraint would be empty and is therefore skipped            
            ActProd             = sum(md.JjPROD[nt,kj]*md.jWatProd[kj] for kj in md.njactivity if md.j_catch[kj]==nc)
            ActCons             = sum(md.JjPROD[nt,kj]*md.jWatCons[kj] for kj in md.njactivity if md.j_catch[kj]==nc)
            BaseFlow            = 0
            NetInflow           = (1-md.wFlowLoss[nc]) * sum(md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) #Net incoming flow from upstream catchments (Mm3)
            RunOff              = md.wRunOff[nt,nc] #Catchment Inflow (Mm3) 
            GwStorage           = sum(md.WwGWSTORAGE[md.t_prev[nt],kaq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[kaq] for kaq in _naquifer(nc)) #Previous ground water storage
            Pumping             = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in _nifzone(nc) for kcul in md.nculture) if md.Options['Groundwater']==1 else 0 #Water abstraction from groundwater reservoir
            BaseFlow            = sum(GwStorage*(1-exp(-md.wGwFlow[kaq])) + (md.wGwRech[nt,kaq]-Pumping)*(1-(1-exp(-md.wGwFlow[kaq]))/md.wGwFlow[kaq]) for kaq in _naquifer(nc)) #Baseflow from groundwater reservoir
            UserWithdrawal      = sum(md.WwSUPPLY[nt,ku] * 1/(1-md.wUserLoss[ku]) for ku in md.nuser if md.user_catch[ku] == nc) #Allocation to non agricultural users (Mm3)
            AgrWithdrawal       = sum(md.AwSUPPLY[nt,kfz,kcu]* 1/(1-md.aIrrgLoss[kfz]) for kfz in _nifzone(nc) for kcu in md.nculture ) #Agricultural module: Allocation to agriculture/crops (Mm3)    
            DeltaStorage        = sum(_storage(nt,kres,mode='delta') for kres in _nlake(nc)) #Storage in lakes and wetlands (Mm3)
            Evaporation         = sum(md.xResEvap[nt,kres] for kres in _nlake(nc)) #Evaporation from lakes and wetlands (Mm3)
            TransNetInflow      = sum(md.WwTRANSFER[nt,ktrans] * (1-md.wTransLoss[ktrans])   for ktrans in md.ntransfer if md.transfer_ds[ktrans] == nc)  #Incoming flow from upstream transfer scheme- if exists (Mm3)
            TransGrossOutflow   = sum(md.WwTRANSFER[nt,ktrans]                               for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc)  #Allocation to downstream transfer scheme- if exists (Mm3)                   
            DEBUG=md.DUMYWATER[nt,nc] if md.Options['Debug mode'] == 1 else 0
            return UserWithdrawal + AgrWithdrawal + TransGrossOutflow + DeltaStorage + Evaporation + ActCons <= NetInflow + TransNetInflow + RunOff + BaseFlow + ActProd +DEBUG
                
        def water_usermaxall(md,nt,nu): #Water allocation to non agricultural user limited by demand [time x user] 
            return md.WwSUPPLY[nt,nu] <= md.wUserDem[md.t_year[nt],md.t_month[nt],nu]
                        
        def water_transfercapacity(md,nt,ntrans): #Transfer scheme carrying capacity [time x transfer]          
            return md.WwTRANSFER[nt,ntrans] <= md.wTransCap[ntrans] + _InvestedCapacity(md.wTransCap[ntrans],ntrans,'transfer',nt)
        
        def water_lakebalance(md,nt,nla): #Storage controlled by first order reservoir or wetland model [time x reservoir]
            nc = md.res_catch[nla]
            Inflow      = (1-md.wFlowLoss[nc])*sum((1-md.wkUP[nla])*md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) 
            RunOff      = (1-md.wkRO[nla]) * md.wRunOff[nt,nc]            
            if ANALYTICAL==1: #Analytical version            
                ET0=max(0,value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc]))*md.mm_to_m
                a=ET0*md.wkV[nla]+md.wAlpha[nla]
                b=Inflow+RunOff-md.wResArea[nla]*ET0
                return md.sw*md.WwRSTORAGE[nt,nla] == _storage(nt,nla,mode='previous')*exp(-a) + b/a*(1-exp(-a))
            else:
                StorPrev    = _storage(nt,nla,mode='previous')
                AvStor      = _storage(nt,nla,mode='average')
                Evap        = md.xResEvap[nt,nla] 
                Outflow     = md.wAlpha[nla] * AvStor**OUTFLOWPOWER
                return md.sw*md.WwRSTORAGE[nt,nla] == StorPrev + Inflow + RunOff - Evap - Outflow #+ DEBUG   

        
        def water_groundwaterbalance(md,nt,naq):
            GwStorage       = md.WwGWSTORAGE[md.t_prev[nt],naq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[naq]
            Pumping         = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in _nifzone(md.aqui_catch[naq]) for kcul in md.nculture)
            if md.wGwFlow[naq] != 0:
                return md.WwGWSTORAGE[nt,naq] == GwStorage*exp(-md.wGwFlow[naq]) + (md.wGwRech[nt,naq]-Pumping)/md.wGwFlow[naq]*(1-exp(-md.wGwFlow[naq]))
            else:
                return md.WwGWSTORAGE[nt,naq] == GwStorage - Pumping
            
        def water_reservoircapacity(md,nt,nres): #Storage in reservoir limited by storage capacity [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir': #non man-operated reservoirs
                return Constraint.Skip
            return md.sw*md.WwRSTORAGE[nt,nres] <= md.wStorCap[nres] + _InvestedCapacity(md.wStorCap[nres],nres,'reservoir',nt)
        
        def water_floodrulecurve(md,nt,nres): #Storage in reservoir limited by flood rule curve [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir': #non man-operated reservoirs
                return Constraint.Skip
            return md.sw*md.WwRSTORAGE[nt,nres] <= md.wFloodRule[nres,md.t_month[nt]]
        
        def water_finalstorage(md,nres): #Storage in reservoir limited by storage capacity [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir':
                return Constraint.Skip
            Debug = md.DUMYSTOR[nres] if md.Options['Initial time step'] == 1 else 0 #Tolerance to constraint if physically impossible at very high cost
            return md.sw*md.WwRSTORAGE[md.Options['tfin'],nres] >= md.wStorFin[nres] - Debug
            
        ##Create constraints##
        md.water_waterbalance       = Constraint(md.ntime, md.ncatch, rule=water_waterbalance)
        md.water_waterbalance2      = Constraint(md.ntime, md.ncatch, rule=water_waterbalance2)   
        md.water_usermaxall         = Constraint(md.ntime, md.nuser,  rule=water_usermaxall)
        md.water_reservoircapacity  = Constraint(md.ntime, md.nres,   rule=water_reservoircapacity)
        md.water_groundwaterbalance = Constraint(md.ntime, md.naquifer, rule=water_groundwaterbalance)       
        md.water_transfercapacity   = Constraint(md.ntime, md.ntransfer, rule=water_transfercapacity)
        md.water_lakebalance        = Constraint(md.ntime, md.nlake,  rule=water_lakebalance)
        if md.Options['Reservoirs'] == 'flood':
            md.water_floodrulecurve = Constraint(md.ntime, md.nres,   rule=water_floodrulecurve)
        if md.Options['Initial time step'] == 1:
            md.water_finalstorage = Constraint(md.nres,  rule=water_finalstorage)
            
        #%%----------------------------Environmental module----------------------------
        #Compute natural outflow at catchment to use as minimum flow in upstream catchments or during softened years
        NatOutflow={t:{c:0 for c in md.ncatch} for t in md.ntime}
        for t in md.ntime:
            for c in md.ncatch:
                NatOutflow[t][c]=md.wRunOff[t,c].value
                upcatch=[c]
                while len(upcatch)>0:
                    UpInflow = sum(md.wRunOff[t,kc].value for kc in md.ncatch if md.catch_ds[kc] in upcatch)
                    NatOutflow[t][c] += UpInflow
                    upcatch=[kc for kc in md.ncatch if md.catch_ds[kc] in upcatch]        
        #Find most severe years to allow softening of constrain in these years
        if md.Options['Hard eflows'] != 1:         
            years       = [y for y in md.nyear]
            YearlyRunoff= [sum(md.wRunOff[t,c].value for t in md.ntime for c in md.ncatch if md.t_year[t]==y) for y in md.nyear]
            YearsRank   = sorted(range(len(YearlyRunoff)), key=lambda k: YearlyRunoff[k]) #Years ranked from minimal to maximal inflows
            nbSoftYears = int((1-md.Options['Hard eflows']) * len(md.nyear)) #Amount of years where the hard constraint can be softened in the optimization period
            SoftYears   = {years[ky]:ky in YearsRank[0:nbSoftYears] for ky in range(len(md.nyear))} #Years where the hard constraint is softened
        
        ##Constraint catalogue## 
        
        #minimum or average
        def env_minflow(md,nt,nef): #env minimum flow [envconstraint x time]
            #Constraint is on wetland area and not outlfow volume
            if md.eflow_hard[nef] == 'area': #
                if AREACONSTRAINT != 'minimum':
                    return Constraint.Skip
                nwet=[la for la in md.nlake if md.res_catch[la]==md.eflow_catch[nef] and md.res_type[la]=='wetland'][0] #find wetland of env constraint
                WetlandArea=(md.wkV[nwet]*_storage(nt,nwet,mode='average')**AREAPOWER+md.wResArea[nwet])
                Relax = md.DUMYEFLOW[nt,nef] if md.Options['Debug mode'] in [1,2] else 0 #Tolerance to constraint if physically impossible at very high cost
                return WetlandArea>=md.wEnvFlow[md.t_month[nt],nef] - Relax #Area > constraint - relax
            #Constraint is on outflow
            if md.eflow_hard[nef] == 0: #If no upstream reservoir is regulating the flow the env flow is limited to actual flow
                EflowConstraint = min(0.95*NatOutflow[nt][md.eflow_catch[nef]],md.wEnvFlow[md.t_month[nt],nef]) #Arbitrary 5% tolerance on consumption of available water
            elif md.Options['Hard eflows'] != 1 and SoftYears[md.t_year[nt]] == 1: #Soft year (=dry year) the constraint is relaxed to the natural outlfow
                EflowConstraint = min(0.95*NatOutflow[nt][md.eflow_catch[nef]],md.wEnvFlow[md.t_month[nt],nef]) #ADD 0.95 is arbitrary to allow some minimal losses on the way
            else:
                EflowConstraint = md.wEnvFlow[md.t_month[nt],nef]
            Relax = md.DUMYEFLOW[nt,nef] if md.Options['Debug mode'] in [1,2] else 0 #Tolerance to constraint if physically impossible at very high cost
            return md.sw*md.WwOUTFLOW[nt,md.eflow_catch[nef]] >= EflowConstraint - Relax  # river: env flow does not count transfer scheme water 
        
        def env_minarea_average(md,nef):
            if md.eflow_hard[nef] != 'area':
                return Constraint.Skip
            else:
                nwet=[la for la in md.nlake if md.res_catch[la]==md.eflow_catch[nef] and md.res_type[la]=='wetland'][0] #find wetland of env constraint
                TotArea=sum(md.wkV[nwet]*_storage(kt,nwet,mode='average')**AREAPOWER+md.wResArea[nwet] for kt in md.ntime)
                TotConst=sum(md.wEnvFlow[md.t_month[kt],nef] for kt in md.ntime)
                Relax = sum(md.DUMYEFLOW[kt,nef] for kt in md.ntime) if md.Options['Debug mode'] in [1,2] else 0
                return TotArea >= TotConst-Relax
        ##Create constraints##    
        md.env_minflow = Constraint(md.ntime, md.neflow, rule=env_minflow)
        if AREACONSTRAINT=='average':
            md.env_minarea = Constraint(md.neflow, rule=env_minarea_average)
        #%%----------------------------Agricultural module----------------------------
        ##Functions and Expressions
        _nfieldculture = lambda cul : [kfd for kfd in md.nfieldculture if md.field_culture[kfd]==cul]
        _nculture = lambda cr : [kcul for kcul in md.nculture if md.culture_crop[kcul]==cr]
        _ntime = lambda y : [kt for kt in md.ntime if md.t_year[kt]==y]
        
        #Crop coefficient (allows two resolutions in the parameter)
        md.xKc = Expression(md.nftype,md.nculture,md.nmonth,
                rule= lambda md,nft,ncul,nm: md.aKc[nft,ncul,nm] 
                                               if md.Options['KckY farmtype'] == 1
                                               else md.aKc[ncul,nm])
        #Yield water response coefficient (allows two resolutions in the parameter)
        md.xkY = Expression(md.nftype,md.nculture,md.nyphase,
                rule= lambda md,nft,ncul,nyps: md.akY[nft,ncul,nyps] 
                                               if md.Options['KckY farmtype'] == 1
                                               else md.akY[ncul,nyps])           
        #Culture water demand [mm/phase] 
        md.xCulDem = Expression(md.nyear,md.nfzone,md.nculture,md.nyphase,
                    rule=lambda md,y,fz,cul,yps:sum(md.xKc[md.fzone_type[fz],cul,md.t_month[kt]]()
                    *value(md.wET0[kt,md.fzone_catch[fz]])
                    *md.yphase_month[cul,yps,md.t_month[kt]] for kt in _ntime(y)))
        #Culture net rainfall [mm/phase] rainfall maximized by crop demand by month, so rainfall does not compensate from month to month       
        md.xCulRain = Expression(md.nyear,md.nfzone,md.nculture,md.nyphase,
                    rule=lambda md,y,fz,cul,yps:sum(min(value(md.wRainFall[kt,md.fzone_catch[fz]]),
                    md.xKc[md.fzone_type[fz],cul,md.t_month[kt]]()*value(md.wET0[kt,md.fzone_catch[fz]]))
                    *md.yphase_month[cul,yps,md.t_month[kt]] for kt in _ntime(y)))
        #Phase ratio [-] to determine in a month with two crop growth phases what share goes to what phase
        md.xPhaseRatio = Expression(md.ntime,md.nfzone,md.nculture,md.nyphase,
                    rule= lambda md,t,fz,cul,yps:(md.yphase_month[cul,yps,md.t_month[t]]
                    /sum(md.yphase_month[cul,kyps,md.t_month[t]] for kyps in md.nyphase)                                            
                    if sum(md.yphase_month[cul,kyps,md.t_month[t]] for kyps in md.nyphase) != 0 else 0))
        #Surface water allocation [Mm3/phase] includes return flow but not losses
        md.xCulSurf = Expression(md.nyear,md.nfzone,md.nculture,md.nyphase,
                    rule= lambda md,y,fz,cul,yps:(sum(md.AwSUPPLY[kt,fz,cul]*md.xPhaseRatio[kt,fz,cul,yps] for kt in _ntime(y)) 
                    if md.aIrrigation[md.fzone_type[fz]] == 1 else 0)) 
        #Groundwater allocation [Mm3/phase] includes return flow but not losses
        md.xCulPump = Expression(md.nyear,md.nfzone,md.nculture,md.nyphase,
                    rule= lambda md,y,fz,cul,yps:(sum(md.AwGWSUPPLY[kt,fz,cul]*md.xPhaseRatio[kt,fz,cul,yps] for kt in _ntime(y)) 
                    if md.Options['Groundwater']==1 else 0))
        #Crop production [1000t/year]
        def _xCulProd(md,y,fz,cul):
            ft=md.fzone_type[fz]
            #Water to crops [Mm3/phase]
            CulWater={yps:(md.xCulSurf[y,fz,cul,yps]+md.xCulPump[y,fz,cul,yps])*(1-md.aCulRturn[ft,cul]) #Surface + pumping - return flows Mm3
                           +md.kha_to_Mha*md.xCULAREA[y,fz,cul]*md.xCulRain[y,fz,cul,yps]*md.mm_to_m3perha for yps in md.nyphase} # Rainfall Mm3
            #Penalty area according to yield water response function [1000ha]
            PenaltyArea=sum(md.xkY[ft,cul,kyps]*(md.xCULAREA[y,fz,cul] 
                        -CulWater[kyps]/(md.xCulDem[y,fz,cul,kyps]*md.mm_to_m3perha*md.kha_to_Mha))
                        if md.xCulDem[y,fz,cul,kyps]() > md.xCulRain[y,fz,cul,kyps]() 
                        and md.xCulDem[y,fz,cul,kyps]() != 0 else 0 for kyps in md.nyphase)
            #Crop production [1000t/y] rmq: kha_to_Mha, and kt_to_Mt could be used as scaling factors for crop production and cultivated area           
            return md.aCulYield[y,md.fzone_type[fz],cul]*(md.xCULAREA[y,fz,cul] - PenaltyArea)*md.kha_to_Mha/md.kt_to_Mt #1000t/y
        md.xCulProd = Expression(md.nyear,md.nfzone,md.nculture, rule=_xCulProd)

        ##Constraint catalogue## 
        def agr_waterall(md,ny,nfz,ncul,nyps):  
            nft=md.fzone_type[nfz]
            #Rainfall accounts net rainfall i.e. rainfall within demand (so rainfall from another month cannot compensate)
            RainFall = md.xCulRain[ny,nfz,ncul,nyps]*md.mm_to_m3perha #Rainfall in m3/ha
            MaxDemand = md.xCulDem[ny,nfz,ncul,nyps]*md.mm_to_m3perha #Demand in m3/ha                                 
            if value(md.xCulDem[ny,nfz,ncul,nyps]) == 0: #no demand=skip constraint
                return Constraint.Skip
            Surface = md.xCulSurf[ny,nfz,ncul,nyps] #Mm3/phase
            Pumping = md.xCulPump[ny,nfz,ncul,nyps] #Mm3/phase
            PhaseDemand = MaxDemand*sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd,kypt]*md.aYieldMat[kypt,nyps] #Water allocated according to path matrix
                                        for kypt in md.nypath for kfd in _nfieldculture(ncul))                            
            return (Surface + Pumping)*(1-md.aCulRturn[nft,ncul]) + RainFall*md.kha_to_Mha*md.xCULAREA[ny,nfz,ncul] >= PhaseDemand #Mm3/phase
            
        def agr_maxland(md,ny,nfz): #Land limit constraint on cultivated area [year x fzone] [1000ha/year]
            ActProd = sum(md.JjPROD[kt,kj]*md.jLandProd[kj] for kt in _ntime(ny) for kj in md.njactivity if md.j_fzone[kj]==nfz) #land produced by activities
            ActCons = sum(md.JjPROD[kt,kj]*md.jLandCons[kj] for kt in _ntime(ny) for kj in md.njactivity if md.j_fzone[kj]==nfz) #land consumed by activities
            Invested = _InvestedCapacity(md.aLandCap[ny,nfz],nfz,'farming',ny,year=1) #land investments
            return md.xPHYAREA[ny,nfz] <= md.aLandCap[ny,nfz] + Invested + ActProd - ActCons            
                
        def agr_fieldarea(md,ny,nfz,nfd,nfdc): #Every culture in a same field has the same area [year x fzone x ftype x (field x culture)] 
            if nfdc==1:
                return Constraint.Skip
            else:
                if md.Options['Yield water response'] == 'nonlinear':
                    return md.AlCULAREA[ny,nfz,nfd,nfdc] == md.AlCULAREA[ny,nfz,nfd,1] #area of field[1] = 2,3...(if exist) have to match the one of field[1] = 1           
                else:
                    return sum(md.AlCULAREA[ny,nfz,nfd,nfdc,kypt] for kypt in md.nypath) == sum(md.AlCULAREA[ny,nfz,nfd,1,kypt] for kypt in md.nypath) #area of field[1] = 2,3...(if exist) have to match the one of field[1] = 1           
                
        def agr_maxcularea(md,ny,nfz,ncul): #Limits the cultivated area by an upper value (usually observed culture area) [year x fzone x culture]
            RF=1 #relative factor if area is not fully developed
            if md.Options['Investment module']== 'continuous' and nfz in md.ninvest:
                RF=_InvestedCapacity(md.aLandCap[ny,nfz],nfz,'farming',ny,year=1)/md.iInvCap[nfz]
                
            return md.xCULAREA[ny,nfz,ncul] <= md.aCulMax[ny,nfz,ncul]*RF 
        
        def agr_onecropchoice(md,ny,nfz,nfd,nfdc): #One culture area choice for the whole optimization period [year x fzone x ftype x (field x culture)] 
            if md.Options['Yield water response'] == 'nonlinear':
                return md.AlCULAREA[ny,nfz,nfd,nfdc] == md.AlCULAREA[md.t_year[md.Options['tini']],nfz,nfd,nfdc]
            else:
                return sum(md.AlCULAREA[ny,nfz,nfd,nfdc,kypt] for kypt in md.nypath) == sum(md.AlCULAREA[md.t_year[md.Options['tini']],nfz,nfd,nfdc,kypt] for kypt in md.nypath)
        
        def agr_cropprodlinearized(md,ny,nfz,ncr): #Crop production [year x fzone x crop]
            nft=md.fzone_type[nfz]
            #linearized yield response: crop production is linked to yield path choice        
            FzoneCropProd = sum(md.aCulYield[ny,nft,md.field_culture[kfd]] # t/ha
            * sum(md.AlCULAREA[ny,nfz,kfd,kypt]  # 1000 ha/year
            * (1-sum(md.xkY[nft,md.field_culture[kfd],kyps]*(1-md.aYieldMat[kypt,kyps]) for kyps in md.nyphase)) 
            for kypt in md.nypath) for kfd in md.nfieldculture if md.culture_crop[md.field_culture[kfd]]==ncr) 
            return md.AcPROD[ny,nfz,ncr] == FzoneCropProd #1000t/year
        
        def agr_maxirrigmonth(md,nt,nfz,ncul): #Irrigation is limited to max phase water demand to avoid compensation from one phase to another - Only for irrigated farming zones
            ny=md.t_year[nt]
            nm=md.t_month[nt]
            nft=md.fzone_type[nfz]
            nc=md.fzone_catch[nfz]                        
            MonthDem= sum(md.yphase_month[ncul,kyps,nm]*md.mm_to_m3perha 
                          *max(0,md.xKc[nft,ncul,nm]()*value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc])) 
                          for kyps in md.nyphase)
            Surface = md.AwSUPPLY[nt,nfz,ncul]  
            Pumping = md.AwGWSUPPLY[nt,nfz,ncul] if md.Options['Groundwater']==1 else 0
            return (Surface + Pumping) * (1-md.aCulRturn[nft,ncul]) <= MonthDem*md.kha_to_Mha*md.xCULAREA[ny,nfz,ncul] 
        
        def agr_cropprodexplicit(md,ny,nfz,ncr): #ADD THIS CAN LEAD TO NEGATIVE YIELDS ! 
            CropProd=sum(md.xCulProd[ny,nfz,kcul] for kcul in _nculture(ncr))
            Debug= md.DUMYCROP[ny,nfz,ncr] if md.Options['Debug mode'] in [1,2] else 0 #DUMMY crop if negative yield
            return md.AcPROD[ny,nfz,ncr] == CropProd + Debug  #1000t/y

        ##Create constraints##
        if md.Options['Crop choice'] != 'fixed' and md.Options['Yield water response'] != 'nonlinear':
            md.agr_waterall     = Constraint(md.nyear, md.nfzone, md.nculture, md.nyphase, rule=agr_waterall)
            md.agr_cropprod     = Constraint(md.nyear, md.nfzone, md.ncrop,         rule=agr_cropprodlinearized)
            md.agr_maxland      = Constraint(md.nyear, md.nfzone,                   rule=agr_maxland)
            md.agr_fieldarea    = Constraint(md.nyear, md.nfzone, md.nfieldculture, rule=agr_fieldarea)
        else:
            md.agr_cropprod     = Constraint(md.nyear, md.nfzone, md.ncrop,         rule=agr_cropprodexplicit)        
        md.agr_maxirrigmonth    = Constraint(md.ntime, md.nifzone, md.nculture,     rule=agr_maxirrigmonth)    
        if md.Options['Yield water response'] == 'nonlinear' and md.Options['Crop choice'] != 'fixed':
            md.agr_maxland      = Constraint(md.nyear, md.nfzone,                   rule=agr_maxland)
            md.agr_fieldarea    = Constraint(md.nyear, md.nfzone, md.nfieldculture, rule=agr_fieldarea)
        if md.Options['Crop choice'] in ['max','once_max']:
            md.agr_maxcularea   = Constraint(md.nyear, md.nfzone, md.nculture,      rule=agr_maxcularea)
        if md.Options['Crop choice'] in ['once','once_max']: 
            md.agr_onecropchoice= Constraint(md.nyear, md.nfzone, md.nfieldculture,   rule=agr_onecropchoice)
        
        #%%----------------------------Crop market module----------------------------
        ##Constraint catalogue##
        
        def agr_cropbalance(md,ny,ncm,ncr): #Crop mass balance (kt) [year x cmarket x crop]
            #Crop consummed and produced by activities
            ActProd     = sum(md.JjPROD[kt,kj]*md.jCropProd[kj] for kt in _ntime(ny) for kj in md.njactivity 
                                 if md.j_cmarket[kj]==ncm and md.j_cropout[kj]==ncr)
            ActCons     = sum(md.JjPROD[kt,kj]*md.jCropCons[kj] for kt in _ntime(ny) for kj in md.njactivity 
                                 if md.j_cmarket[kj]==ncm and md.j_cropin[kj]==ncr)
            #Crops produced in farming zones within the crop market
            Production  = sum(md.AcPROD[ny,kfz,ncr] for kfz in md.nfzone if md.fzone_cmarket[kfz]==ncm)
            #Imported and exported crops
            Import      = sum(md.AcTRANS[ny,kct,ncr]*(1-md.aTransLoss[kct]) for kct in md.nctrans if md.aTransOut[kct]==ncm)
            Export      = sum(md.AcTRANS[ny,kct,ncr] for kct in md.nctrans if md.aTransIn[kct]==ncm)
            #Crops consummed within the crop markets
            if md.Options['Crop demand elasticity']!='nonlinear':
                Supply  = sum(md.AcSUPPLY[ny,ncm,ncr,kcds] for kcds in md.ncdstep) 
            else:
                Supply  = md.AcSUPPLY[ny,ncm,ncr]          
            #Crops produced by external markets (no land and water constraints)
            if ncm in md.nextcmarket: 
                Production = md.AcEXTPROD[ny,ncm,ncr] #crop produced in international/external market (kt)
            return Supply <= Production + Import - Export + ActProd - ActCons
            # '<=' instead of '==' to allow crop waste in case of over production
               
        def agr_cropdemand(md,ny,ncm,ncr,ncds): #Crop market level crop step demand [year x cmarket x crop x cdstep]
            CropStepDemand  = md.aCropDem[ny,ncm,ncr] * md.aStepDem[ncm,ncr,ncds] 
            return  md.AcSUPPLY[ny,ncm,ncr,ncds] <= CropStepDemand 
        
        def agr_cropmindemand(md,ny,ncm,ncr): #Crop market level crop minimum demand ensuring food security [year x cmarket x crop]
            if md.Options['Crop demand elasticity']!='nonlinear':
                CropSupply  = sum(md.AcSUPPLY[ny,ncm,ncr,kcds] for kcds in md.ncdstep)
            else:
                CropSupply  = md.AcSUPPLY[ny,ncm,ncr]
            return  CropSupply >= md.aMinDem[ncm,ncr] 
                
        ##Create constraints##
        md.agr_cropbalance      = Constraint(md.nyear, md.ncmarket, md.ncrop, rule=agr_cropbalance)
        if md.Options['Crop demand elasticity']!='nonlinear':
            md.agr_cropdemand   = Constraint(md.nyear, md.ncmarket, md.ncrop, md.ncdstep, rule=agr_cropdemand)
        if md.Options['Minimum supply'] == 1:
            md.agr_cropmindemand= Constraint(md.nyear, md.ncmarket, md.ncrop, rule=agr_cropmindemand)
        
        #%%----------------------------Energy production module----------------------------
        #Constraint catalogue
        def engy_turbinecapacity(md,nt,npld,nhp):            
            if md.hp_res[nhp] == 'ROR':
                return md.EwHPDISCHARGE[nt,npld,nhp] <= md.wMaxTurb[nhp] * md.m3pers_to_Mm3perMonth * md.eLoadTime[npld]
            elif md.Options['Hydropower'] == 'nonlinear':
                nres = md.hp_res[nhp]
                if md.wFixHead[nres]!=md.wFixHead[nres]: #parameter is not defined
                    AvStorage = _storage(nt,nres,mode='average')
                    MaxTurbineFlow = _TurbFlowVol(AvStorage,md.wResMOL[nres],md.wResFSL[nres],md.wMinTurb[nhp],md.wMaxTurb[nhp])
                else:
                    MaxTurbineFlow = md.wFixHead[nres]*md.wMaxTurb[nhp]
                return md.EwHPDISCHARGE[nt,npld,nhp] <= MaxTurbineFlow * md.m3pers_to_Mm3perMonth * md.eLoadTime[npld]
            return Constraint.skip()
        
        def engy_hpdischarge(md,nt,nres): #Hydropower discharge [time x reservoir]
            if nres not in [md.hp_res[khp] for khp in md.nhpp]: #No hydropower in reservoir
                return Constraint.Skip
            else:
                HppDischarge = sum(md.EwHPDISCHARGE[nt,kpld,khp] for kpld in md.npload for khp in md.nhpp if md.hp_res[khp] == nres)
            return HppDischarge <= md.sw*md.WwOUTFLOW[nt,md.res_catch[nres]]
                
        def engy_hpdischarge2(md,nt,npld,nhp): #Hydropower discharge for ROR hydropowers [hydropower]
            if md.hp_res[nhp]=='ROR': #Run-Of-the-River Hydropower
                return md.EwHPDISCHARGE[nt,npld,nhp] <= md.sw*md.WwOUTFLOW[nt,md.hp_catch[nhp]] * md.eLoadTime[npld]
            else:
                return Constraint.Skip
                
        def engy_hpprod(md,nt,npld,nhp): #Hydropower production limited by downstream allocation [time x pload x hpp]            
            RelHead=1 #By default full head is available
            if md.Options['Hydropower'] == 'nonlinear' and md.hp_res[nhp] != 'ROR':
                nres = md.hp_res[nhp]
                if md.wFixHead[nres]!=md.wFixHead[nres]: #parameter is not defined
                    AvStorage = _storage(nt,nres,mode='average')
                    RelHead = _RelHeadVol(AvStorage,md.wResMOL[nres],md.wResFSL[nres],md.wMinHead[nres],md.wMaxHead[nres])
                else:
                    RelHead = md.wFixHead[nres]
            return md.se*md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * md.eHppProd[nhp] * md.EwHPDISCHARGE[nt,npld,nhp] * RelHead
                
        def engy_hpcapacity(md,nt,npld,nhp): #Hydropower production per load segment limited by hydropower capacity [time x pload x hpp]
            CapacityFactor = md.MW_to_GWhperMonth * md.eLoadTime[npld]
            return md.se*md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * (md.eHppCap[nhp] + _InvestedCapacity(md.eHppCap[nhp],nhp,'hydro',nt)) * CapacityFactor                 
        
        def engy_oppcapacity(md,nt,npld,nop): #Powerplant production per load segment limited by capacity [time x pload x opp]
            ny=md.t_year[nt]
            nm=md.t_month[nt]
            CapacityFactor = md.MW_to_GWhperMonth * md.eLoadTime[npld] #No particular restriction if not specified
            if md.Options['Load capacity'] == 1:
                CapacityFactor = CapacityFactor * md.eLoadCap[npld,md.op_ptech[nop],nm]                
            return md.se*md.EeOPPROD[nt,npld,nop] <= (md.eOppCap[ny,nop] + _InvestedCapacity(md.eOppCap[ny,nop],nop,'power',nt)) * CapacityFactor
        
        def engy_gencapacity(md,nt,npld,npt,npm): #Capacity constraint of Generic capacity [time x pload x ptech x pmarket]
            CapacityFactor = md.MW_to_GWhperMonth * md.eLoadTime[npld] #No particular restriction if not specified
            if md.Options['Load capacity'] == 1:
                CapacityFactor = CapacityFactor * md.eLoadCap[npld,npt,md.t_month[nt]]
            years= [ky for ky in md.nyear if ky <= md.t_year[nt]-md.eConstTime[npt,npm] and ky >= md.t_year[nt]-md.eLifeTime[npt,npm]]
            return md.EeGENPROD[nt,npld,npt,npm] <= sum(md.EeGENCAP[ky,npt,npm] for ky in years) * CapacityFactor 

        def engy_maxgeninvest(md,npt,npm): #Maximum increase in Generic capacity
            if md.eMaxCap[npt,npm] == 'unlimited': #hard coded key word for unlimited capacity investment
                return Constraint.Skip
            else:
                return sum(md.se*md.EeGENCAP[ky,npt,npm] for ky in md.nyear) <= md.eMaxCap[npt,npm]
        
        def engy_genramping(md,nt,nplda,npldb,npt,npm): #Ramping constraints for Generic capacity
            if nplda == npldb:
                return Constraint.Skip
            else:
                ProdRateSegA    = md.EeGENPROD[nt,nplda,npt,npm] / md.eLoadTime[nplda] #[kWh/day]
                ProdRateSegB    = md.EeGENPROD[nt,npldb,npt,npm] / md.eLoadTime[npldb] #[kWh/day]
                RampingCap      = sum(md.se*md.EeGENCAP[ky,npt,npm] for ky in md.nyear if ky <= md.t_year[nt] and ky >= md.t_year[nt]
                                      -md.eLifeTime[npt,npm]) * md.eTechRamp[npt,npm] #[kWh/day]
                return  ProdRateSegA <= ProdRateSegB + RampingCap   #Power plant's Load rate from load segment to other is not allow to vary more than ramp rate
        
        def engy_oppramping(md,nt,nplda,npldb,nop): #Ramping constraints for other power plants
            ny=md.t_year[nt]
            if nplda == npldb:
                return Constraint.Skip
            else:
                ProdRateSegA    = md.se*md.EeOPPROD[nt,nplda,nop] / md.eLoadTime[nplda] #[kWh/day]
                ProdRateSegB    = md.se*md.EeOPPROD[nt,npldb,nop] / md.eLoadTime[npldb] #[kWh/day]
                RampingCap      = md.eOppCap[ny,nop] * md.eOppRamp[nop] #[kWh/day]
                return  ProdRateSegA <= ProdRateSegB + RampingCap    #Power plant's Load rate from load segment to other is not allow to vary more than ramp rate
        
        #Create Constraints
        md.engy_hpdischarge  = Constraint(md.ntime, md.nres,            rule=engy_hpdischarge)
        md.engy_hpdischarge2 = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_hpdischarge2)
        md.engy_hpprod       = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_hpprod)
        md.engy_hpcapacity   = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_hpcapacity)        
        if md.Options['Hydropower'] == 'nonlinear':
            md.engy_turbinecapacity  = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_turbinecapacity)
        md.engy_oppcapacity     = Constraint(md.ntime, md.npload, md.nopp, rule=engy_oppcapacity)                         
        md.engy_gencapacity     = Constraint(md.ntime, md.npload, md.nptech, md.npmarket, rule=engy_gencapacity)
        md.engy_maxgeninvest    = Constraint(md.nptech, md.npmarket, rule=engy_maxgeninvest)
        if md.Options['Ramping'] == 1:
            md.engy_oppramping  = Constraint(md.ntime, md.npload, md.npload, md.nopp, rule=engy_oppramping)  
            md.engy_genramping  = Constraint(md.ntime, md.npload, md.npload, md.nptech, md.npmarket, rule=engy_genramping)
            
            
        #%%----------------------------Power market module----------------------------
        #Constraint catalogue    
        def engy_balance(md,nt,npld,npm): #Power supply per time slice [pmarket x time x pload]
            ActProd     = 1/md.se*sum(md.JjPROD[nt,kj]*md.jPowProd[kj]*md.eLoadTime[npld] for kj in md.njactivity if md.j_pmarket[kj]==npm)    #Energy production by activities
            ActCons     = 1/md.se*sum(md.JjPROD[nt,kj]*md.jPowCons[kj]*md.eLoadTime[npld] for kj in md.njactivity if md.j_pmarket[kj]==npm)    #Energy consumption by activities
            HpProd      = sum(md.EeHPPROD[nt,npld,khp]                          for khp in md.nhpp if md.hp_pmarket[khp]==npm)      #Hydropower energy prod
            OppProd     = sum(md.EeOPPROD[nt,npld,kop]                          for kop in md.nopp if md.op_pmarket[kop]==npm)      #Power plant energy prod
            GenProd     = sum(md.EeGENPROD[nt,npld,kpt,npm]                     for kpt in md.nptech)                               #Generic Power plant energy prod
            NetImport   = sum(md.EeTRANS[nt,npld,ktl] * (1-md.eTransLoss[ktl])  for ktl in md.ntransline if md.eTransOut[ktl]==npm) 
            GrossExport = sum(md.EeTRANS[nt,npld,ktl]                           for ktl in md.ntransline if md.eTransIn[ktl]==npm)
            GrossSupply = md.EeSUPPLY[nt,npld,npm] / (1-md.eSupLoss[npm])                                           
            return GrossSupply + GrossExport + ActCons == HpProd + OppProd + GenProd + NetImport + ActProd #
                
        def engy_transmissioncap(md,nt,npld,ntl): #Transmission limited to line capacity [time x pload x tmsline] 
            MaxCap=(md.eTransCap[ntl] + _InvestedCapacity(md.eTransCap[ntl],ntl,'transmission',nt)) * md.MW_to_GWhperMonth * md.eLoadTime[npld]
            return md.se*md.EeTRANS[nt,npld,ntl] <= MaxCap
                
        def engy_demand(md,nt,npld,npm): #Power demand satisfaction
            #ny0         = md.t_year[md.Options['tini']]
            Supply      = md.se*md.EeSUPPLY[nt,npld,npm] 
            LoadSegment = md.eEngyDem[md.t_year[nt],md.t_month[nt],npm] * md.eLoadDem[npld] # * md.eDemYear[npm]**(md.t_year[nt]-ny0) Demand of power market for power load segment          
            return Supply <= LoadSegment
                
        #Create constraints
        md.engy_balance         = Constraint(md.ntime, md.npload, md.npmarket, rule=engy_balance)
        md.engy_transmissioncap = Constraint(md.ntime, md.npload, md.ntransline, rule=engy_transmissioncap)
        md.engy_demand          = Constraint(md.ntime, md.npload, md.npmarket, rule=engy_demand)
        
        #%%----------------------------Activities----------------------------
        #Constraint catalogue
        #activity is constrained by capacity
        def activity_capacity(md,nt,nj):
            return md.JjPROD[nt,nj]<=md.jProdCap[nj]
        md.activity_capacity    = Constraint(md.ntime, md.njactivity, rule=activity_capacity)
#%%SAVE MODEL
        self.model=md


