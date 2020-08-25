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

from pyomo.environ import ConcreteModel, Set, Param, Var, Objective, Constraint
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
        
#%%PRAMETER READING FUNCTIONS
        #parameters that have a yearly value (to model growth) - This list is hardcoded as impacts variable naming
        GrowingParameters=['wUserDem','eEngyDem','eFuelCost','eCO2Val','eCAPEX','aCropDem','aCulYield','aLandCap','aCulMax'] 
        def ReadParam(ParamName,option=1,index=0,time='all'):
            if directparam!=0 and ParamName in directparam.keys(): #Assign parameter directly and not from parameters object 
                return directparam[ParamName]
            if ParamName not in parameters.val.keys() or option == 0: #parameter is not activated
                if index==1:
                    return [] #returns empty set
                return {} #returns empty parameter
            if ParamName in GrowingParameters:
                Param = GrowParam(ParamName)
            else:
                if ParamName in parameters.scen.keys():                
                    Param = parameters.val[ParamName][parameters.val[parameters.scen[ParamName]][scenario]]
                else:
                    Param = parameters.val[ParamName]
            if time != 'all': #Only selected time steps 
                Param={k:Param[k] for k in Param.keys() if ((type(k)==tuple and k[0] in time) or (type(k)!=tuple and k in time))}
            return Param
        
        def LinearGrowAtoB(ValuesA,ValuesB,nyAB,Keys,yyears):
            if type([key for key in Keys][0]) is tuple: #keys need to be formulated differently if they are uni or multidimensional
                return {(y,)+keys:ValuesA[keys] + (ValuesB[keys]-ValuesA[keys])/nyAB*(y-yyears[0]) for y in yyears for keys in Keys}
            else:
                return {(y,keys): ValuesA[keys] + (ValuesB[keys]-ValuesA[keys])/nyAB*(y-yyears[0]) for y in yyears for keys in Keys}
        
        def GrowParam(ParamName):
            years=sorted(parameters.val['nyear'])
            scen=parameters.val[parameters.scen[ParamName]][scenario] if ParamName in parameters.scen.keys() else ''
            if md.Options['Demand growth']==1 and scen in ReadParam('ngrowthscen'): #Linear interpolation from values A to B to C (optional)
                ValuesA=parameters.val[ParamName][ReadParam('StatusA')[scen]]
                ValuesB=parameters.val[ParamName][ReadParam('StatusB')[scen]]                                       
                nyAB=int(ReadParam('nyearsAtoB')[scen])                               
                DicAB=LinearGrowAtoB(ValuesA,ValuesB,nyAB,ValuesA.keys(),years[0:nyAB])
                if ReadParam('StatusC')[scen] == ReadParam('StatusC')[scen]: #Value is not a NAN
                    ValuesC=parameters.val[ParamName][ReadParam('StatusC')[scen]]
                    nyBC=int(ReadParam('nyearsBtoC')[scen])
                    DicBC=LinearGrowAtoB(ValuesB,ValuesC,nyBC,ValuesA.keys(),years[nyAB:nyAB+nyBC])
                    DicAB.update(DicBC)                
            else: #Constant parameter through years 
                if ParamName in parameters.scen.keys():                    
                    DicAB=LinearGrowAtoB(parameters.val[ParamName][scen],parameters.val[ParamName][scen],1,parameters.val[ParamName][scen].keys(),years)
                else:
                    DicAB=LinearGrowAtoB(parameters.val[ParamName],parameters.val[ParamName],1,parameters.val[ParamName].keys(),years)
            return DicAB
        
#%%INTIALIZE MODEL
        md = ConcreteModel()
    #Options
        md.nconfig          = Set(initialize=ReadParam('nconfig'))             #Configurations of the model
        Options= {opt:parameters.val['Options'][opt,parameters.val['sOptions'][scenario]] for opt in md.nconfig} #Select options of considered scenario
        md.Options          = Param(md.nconfig, initialize=Options)
        #time framing
        time = [t for t in range(md.Options['tini'],md.Options['tfin']+1)]
        t_year = ReadParam('t_year',time=time)
        year = [t_year[t] for t in time]       
        if md.Options['Initial time step'] == 0: #Close the loop for cyclic model
            parameters.val['t_prev'][md.Options['tini']]=md.Options['tfin'] #ADD: assumes 't_prev' parameter is not scenario dependent
    #Unit conversion
        md.mm_to_m3perha       = 1/1000*10000 # 0,001m * 10000m2/ha
        md.MW_to_GWhperMonth   = 1*24*365/12/1000 # 1 MW * 24hours/day * 365/12days/month / 1000MW/GW 
        md.m3pers_to_Mm3perMonth = 1*3600*24*30.44/10**6 #1 m3/s * 3600 s/hour * 24 hours/day * 30.44 days/month / 10^6m3/Mm3
        md.kha_to_Mha       = 1/1000 # 1000 ha = 1/1000 Mha
        md.kt_to_Mt         = 1/1000 # 1000 t = 1/1000 Mt
        #md.m3_to_Mm3        = 10**-6 # 1 m3 = 10**-6 Mm3        
        #md.t_to_Mt          = 10**-6 # 1 t = 10**-6 Mt
    #Scaling parameters
        md.debugcost = Param(initialize=1000000,mutable=True) #debug cost (variabes units $/m3 - 100$/t)   
        md.sw = Param(initialize=1,mutable=True) # Scaling factor for WwRSTORAGE and WwOUTFLOW (1000=convert DV from Mm3 to km3)
        md.se = Param(initialize=1,mutable=True) # Scaling factor for all energy decision variables
        #md.sa = 10**-0
        #md.sc = 10**-0
#%%**Water Module**#       
    #Indices
        md.ncountry         = Set(initialize=ReadParam('ncountry'))                 #Countries [indice country] (id)
        md.ncatch           = Set(initialize=ReadParam('catch_ds').keys())          #Catchments [indice catchment] (id)
        md.ntime            = Set(initialize=time)                                  #Time steps [indice time] (id)
        md.nyear            = Set(initialize=set(year))                             #Years [indice year] (id)
        md.nmonth           = Set(initialize=ReadParam('nmonth'))                   #Months [indice month] (id)         
        md.nres             = Set(initialize=ReadParam('wStorCap').keys())          #Reservoirs [indice reservoir] (id) 
        md.nuser            = Set(initialize=ReadParam('user_catch').keys())        #Water Users [indice user] (id)
        md.naquifer         = Set(initialize=ReadParam('wGwFlow',option=md.Options['Groundwater']).keys())  #Groundwater aquifers [indice aquifer] (id)           
    #Spatio-Temporal connections
        md.catch_ds         = Param(md.ncatch, initialize=ReadParam('catch_ds'))                   #Downstream catchment of catchment [catchment] (id)
        md.catch_country    = Param(md.ncatch, initialize=ReadParam('catch_country'))              #Country of catchment [catchment] (id) #REM: DOES NOT ALWAYS MAKE SENSE
        md.res_catch        = Param(md.nres,  initialize=ReadParam('res_catch'))                   #Reservoir's catchment [reservoir] (id)
        md.user_catch       = Param(md.nuser, initialize=ReadParam('user_catch'))                  #User's catchment [user] (id)
        md.user_dscatch     = Param(md.nuser, initialize=ReadParam('user_dscatch'))                #User's downstream catchment [user] (id)
        md.user_country     = Param(md.nuser, initialize=ReadParam('user_country'))                #User's country [user] (id)
        md.t_month          = Param(md.ntime, initialize=ReadParam('t_month',time=time))           #Month of time step [time] (id)
        md.t_year           = Param(md.ntime, initialize=ReadParam('t_year',time=time))            #Year of time step [time] (id)
        md.t_prev           = Param(md.ntime, initialize=ReadParam('t_prev',time=time))            #Previous time step [time] (id)       
            
    #Hydrology
        md.wRunOff          = Param(md.ntime, md.ncatch, mutable=True, initialize=ReadParam('wRunOff',time=time))   #Runoff [time x catchment] (Mm3/month)
        md.wRainFall        = Param(md.ntime, md.ncatch, mutable=True, initialize=ReadParam('wRainFall',time=time)) #Precipitation [time x catchment] (mm)
        md.wET0             = Param(md.ntime, md.ncatch, mutable=True, initialize=ReadParam('wET0',time=time))      #ET0 [time x catchment] (mm)
        md.wFlowLoss        = Param(md.ncatch, initialize=ReadParam('wFlowLoss'))                                   #Flow loss [catchment] (%)
    #Groundwater
        opt= md.Options['Groundwater']
        md.wGwRech          = Param(md.ntime, md.naquifer, initialize=ReadParam('wGwRech',option=opt,time=time))    #Groundwater recharge [time x catchment] (Mm3/month)
        md.wGwFlow          = Param(md.naquifer, initialize=ReadParam('wGwFlow',option=opt))        #Groundwater baseflow coeeficient [catchment] (-)
        md.wGwIni           = Param(md.naquifer, initialize=ReadParam('wGwIni',option=opt))         #Groundwater initial storage (if not cyclic) [catchment] (Mm3)
        md.wGwCost          = Param(md.naquifer, initialize=ReadParam('wGwCost',option=opt))        #Groundwater pumping cost [catchment] ($/m3)
    #Water Users
        md.wUserDem         = Param(md.nyear, md.nmonth, md.nuser, initialize=ReadParam('wUserDem',time=md.nyear))#User water demand [year x month x user] (Mm3/month)        
        md.wUserRturn       = Param(md.nuser, initialize=ReadParam('wUserRturn'))                   #Return rate of user water allocation [user] (%)
        md.wUserLoss        = Param(md.nuser, initialize=ReadParam('wUserLoss'))                    #Loss rate of user water allocation [user] (%)
        md.wUserVal         = Param(md.nuser, initialize=ReadParam('wUserVal'))                     #User marginal value for water [user] ($/m3)
        md.wSupCost         = Param(md.nuser, initialize=ReadParam('wSupCost'))                     #Water supply costs [user] ($/m3) 
    #Reservoirs
        md.wStorCap         = Param(md.nres, initialize=ReadParam('wStorCap'))                      #Storage Capacity [reservoir] (Mm3)
        if md.Options['Reservoirs'] == 'flood':
            md.wFloodRule   = Param(md.nres, md.nmonth, initialize=ReadParam('wFloodRule'))         #Flood rule curve [reservoir x month] (Mm3)
        if md.Options['Initial time step'] == 1:
            md.wStorIni     = Param(md.nres, mutable=True, initialize=ReadParam('wStorIni'))        #Initial storage [reservoir] (Mm3) - only used if initial time step is used
            md.wStorFin     = Param(md.nres, mutable=True, initialize=ReadParam('wStorFin'))        #Final storage [reservoir] (Mm3) - only used if initial time step is used
        md.wkV              = Param(md.nres, initialize=ReadParam('wkV'))                           #Area - volume coef, kV in "A = A0 + kV * V" [reservoir] (ha/m3)
        md.wResArea         = Param(md.nres, initialize=ReadParam('wResArea'))                      #Minimum area, A0 in "A = A0 + kV * V" [reservoir] (Mha)
    #Transfer schemes
        opt= md.Options['Transfers']
        md.ntransfer        = Set(initialize=ReadParam('wTransCap',option=opt).keys())              #Transfer schemes [indice transfer] (id)
        md.transfer_us      = Param(md.ntransfer, initialize=ReadParam('transfer_us',option=opt))   #Upstream catchment of transfer scheme [transfer] (id) 
        md.transfer_ds      = Param(md.ntransfer, initialize=ReadParam('transfer_ds',option=opt))   #Downstream catchment of transfer scheme [transfer] (id) 
        md.wTransCap        = Param(md.ntransfer, initialize=ReadParam('wTransCap',option=opt))     #Capacity of transfere scheme [transfer] (Mm3/month)
        md.wTransLoss       = Param(md.ntransfer, initialize=ReadParam('wTransLoss',option=opt))    #Loss rate of transfere scheme [transfer] (%)
    #Lakes and wetlands
        opt= md.Options['Lakes']
        md.nlake            = Set(initialize=[res for res in md.nres if ReadParam('res_type')[res]=='lake'] if md.Options['Lakes']==1 else []) #Lakes are a sub-ensemble of reservoirs
        md.res_type         = Param(md.nres,  initialize=ReadParam('res_type',option=opt))      #Type of reservoir 'reservoir'=can be controlled [reservoir] (id)
        md.wAlpha           = Param(md.nlake, initialize=ReadParam('wAlpha',option=opt))        #First order outflow coefficient Out = alpha * V [lake] (%)
        md.wkRO             = Param(md.nlake, initialize=ReadParam('wkRO',option=opt))          #Share of runoff not diverted by wetland [lake] (%)
        md.wkUP             = Param(md.nlake, initialize=ReadParam('wkUP',option=opt))          #Share of upstream flow not diverted by wetland [lake] (%)       
        
#%%**Environmental flow Module**#
        md.neflow           = Set(initialize=ReadParam('eflow_catch',option=md.Options['Eflows']).keys())               #Environmental constraint point [indice envconst] (id)
        md.eflow_catch      = Param(md.neflow, initialize=ReadParam('eflow_catch',option=md.Options['Eflows']))         #Catchment of environmental flow constraint [envconst] (id)
        md.eflow_hard       = Param(md.neflow, initialize=ReadParam('eflow_hard',option=md.Options['Eflows']))          #Eflow as hard or soft (limited to actual runoff) constraint [envconst] (binary)
        md.wEnvFlow         = Param(md.nmonth, md.neflow, initialize=ReadParam('wEnvFlow',option=md.Options['Eflows'])) #Env flow [month x envconst] (Mm3/month)
        
#%%**Agriculture Module**#        
    #Indices 
        opt=md.Options['Agriculture']
        md.ncrop            = Set(initialize=ReadParam('ncrop',option=opt,index=1))       #Crops [indice crop] (id)
        md.nypath           = Set(initialize=ReadParam('nypath',option=opt,index=1))      #Yield water response path [indice ypath] (id)
        md.nyphase          = Set(initialize=ReadParam('nyphase',option=opt,index=1))     #Yield water response phases [indice yphase] (id)
        md.nfzone           = Set(initialize=ReadParam('fzone_catch',option=opt).keys())  #Farm zones [indice fzone] (id)
        md.nifzone          = Set(initialize=[fz for fz in md.nfzone if ReadParam('aIrrigation')[ReadParam('fzone_type')[fz]]==1]) #Irrigated farming zones subindice of farming zones [indice ifzone] (id)
        md.nftype           = Set(initialize=ReadParam('nftype',option=opt,index=1))      #Farm types [indice ftype] (id)
        md.nculture         = Set(initialize=ReadParam('culture_crop',option=opt).keys()) #Culture -one crop at one time- [indice culture] (id)
        md.nfieldculture    = Set(dimen=2,initialize=ReadParam('nfieldculture',option=opt,index=1))#Field [field x culture] [indice field] (id)
    #Spatio-Temporal connections
        md.fzone_catch      = Param(md.nfzone, initialize=ReadParam('fzone_catch',option=opt))                #Farming zones's catchment [fzone] (id)
        md.fzone_cmarket    = Param(md.nfzone, initialize=ReadParam('fzone_cmarket',option=opt))              #Farming zones's crop market [fzone] (id)
        md.fzone_country    = Param(md.nfzone, initialize=ReadParam('fzone_country',option=opt))              #Farming zones's country [fzone] (id)
        md.fzone_dscatch    = Param(md.nfzone, initialize=ReadParam('fzone_dscatch',option=opt))              #Return flow recieving catchment of farming zone [fzone] (id)
        md.fzone_type       = Param(md.nfzone, initialize=ReadParam('fzone_type',option=opt))                 #Farming type of farming zone [fzone] (id)
        md.field_culture    = Param(md.nfieldculture, initialize=ReadParam('field_culture',option=opt))       #Culture of Field [field] (id)
        md.culture_crop     = Param(md.nculture, initialize=ReadParam('culture_crop',option=opt))             #Crop produced by culture [culture] (id) 
        md.yphase_month     = Param(md.nculture, md.nyphase, md.nmonth, initialize=ReadParam('yphase_month',option=opt)) # Phase month matrix, share of each phase in each month [culture x yphase x time]
    #Farming types
        md.aCulRturn        = Param(md.nftype, md.nculture, initialize=ReadParam('aCulRturn',option=opt))     #Water return rate from allocated [ftype x culture] (%) 
        md.aCulYield        = Param(md.nyear, md.nftype, md.nculture, initialize=ReadParam('aCulYield',option=opt,time=md.nyear))  #Maximum yield of culture [ftype x culture] (t/ha) 
        md.aCulCost         = Param(md.nftype, md.nculture, initialize=ReadParam('aCulCost',option=opt))      #Cultivation cost [ftype x culture] ($/ha)
        md.aIrrigation      = Param(md.nftype, initialize=ReadParam('aIrrigation',option=opt))                #Irrigation (1) or rainfed culture(0) [ftype] (binary) 
    #Farming Zones
        md.aLandCap         = Param(md.nyear, md.nfzone, initialize=ReadParam('aLandCap',option=opt,time=md.nyear)) #Land capacity for agriculture [year x fzone] (1000ha)
        md.aIrrgCost        = Param(md.nfzone, initialize=ReadParam('aIrrgCost',option=opt))            #Irrigation cost [fzone] ($/m3)
        md.aIrrgLoss        = Param(md.nfzone, initialize=ReadParam('aIrrgLoss',option=opt))            #Water loss rate from allocated [fzone] (%) 
    #Cultures
        md.akY              = Param(md.nculture, md.nyphase, initialize=ReadParam('akY',option=opt))    #crop yield response phase-specific factor [culture x yphase](%) 
        md.aKc              = Param(md.nculture, md.nmonth, initialize=ReadParam('aKc',option=opt))     #crop water demand factor [culture x yphase](%) 
        if md.Options['Yield water response'] != 'nonlinear':
            md.aYieldMat    = Param(md.nypath, md.nyphase, initialize=ReadParam('aYieldMat',option=opt))#yield response path matrix [ypath x yphase]
        if md.Options['Crop choice'] in ('max','once_max'):
            md.aCulMax      = Param(md.nyear, md.nfzone, md.nculture, initialize=ReadParam('aCulMax',option=opt, time=md.nyear)) #Maximum area per culture and farming zone [fzone x culture] (1000ha)
        
#%%**Crop Market module**#
    #Crop markets and demands
        opt=md.Options['Crop market']
        md.ncmarket         = Set(initialize=ReadParam('cmarket_country',option=opt).keys())                  #Crop markets [indice cmarket] (id)
        md.cmarket_country  = Param(md.ncmarket, initialize=ReadParam('cmarket_country',option=opt))          #Country of Crop Market [cmarket] (id)               
        md.nextcmarket      = Set(initialize=[cm for cm in md.ncmarket if ReadParam('cmarket_type')[cm]=='ExternalMarket']) #External crop markets with a crop supply curve (id)
        md.aCropDem         = Param(md.nyear, md.ncmarket, md.ncrop, initialize=ReadParam('aCropDem',option=opt,time=md.nyear)) #Crop demand [year x cmarket x crop] (1000t/year)
        md.aCropVal         = Param(md.ncmarket, md.ncrop, initialize=ReadParam('aCropVal',option=opt))       #Value of crop [cmarket x crop] ($/t)
    #Food security
        if md.Options['Minimum supply'] == 1:
            md.aMinDem      = Param(md.ncmarket, md.ncrop, initialize=ReadParam('aMinDem',option=opt))    #Crop minimum demand = food security [year x cmarket x crop] (1000t/year)
    #Demand steps - representing elasticity of crop demand
        if md.Options['Crop demand elasticity'] == 1:   #Parametrization of demand curve according to parameters
            aStepLen        = int(md.Options['aStepLen'])
            aNumStep        = int(md.Options['aNumStep'])
            aDemEla         = ReadParam('aDemEla')
            md.ncdstep      = Set(initialize=range(2*aNumStep+1))         #Share of demand step for crop [indice cdstep] (id)
            aStepDem        = {(cm,cr,cds): max(0,(cds==0)*(1-aStepLen) + (cds!=0)*aStepLen/aNumStep) for cds in md.ncdstep for cr in md.ncrop for cm in md.ncmarket}   
            aStepVal        = {(cm,cr,cds): max(0,1+(cds-aNumStep)*aStepLen/aNumStep/aDemEla[cm,cr]) for cds in md.ncdstep for cr in md.ncrop for cm in md.ncmarket}
            md.aStepDem     = Param(md.ncmarket, md.ncrop, md.ncdstep, initialize=aStepDem)          #Demand share step percentage [cmarket x ncrop x cdstep] (%)
            md.aStepVal     = Param(md.ncmarket, md.ncrop, md.ncdstep, initialize=aStepVal)          #Demand share sclice marginal value coefficient [cmarket x ncrop x cdstep] (-)                                           
        else:
            md.ncdstep      = Set(initialize=[1])       #Single demand point
            md.aStepDem     = Param(md.ncmarket, md.ncrop, md.ncdstep, initialize={(cm,cr,1):1 for cr in md.ncrop for cm in md.ncmarket}) #No elasticity in demand
            md.aStepVal     = Param(md.ncmarket, md.ncrop, md.ncdstep, initialize={(cm,cr,1):1 for cr in md.ncrop for cm in md.ncmarket}) #No elasticity in demand
    #Crop transport
        md.nctrans          = Set(initialize=ReadParam('aTransLoss',option=opt).keys())                       #Crop transport roads [indice ctrans] (id)
        md.aTransIn         = Param(md.nctrans, initialize=ReadParam('aTransIn',option=opt))                  #Exporting crop market of crop transport road [ctrans] (id) 
        md.aTransOut        = Param(md.nctrans, initialize=ReadParam('aTransOut',option=opt))                 #Importing crop market of crop transport road [ctrans] (id)
        md.aTransLoss       = Param(md.nctrans, initialize=ReadParam('aTransLoss',option=opt))                #Crop transport loss matrix [ctrans] (%)
        md.aTransCost       = Param(md.nctrans, md.ncrop, initialize=ReadParam('aTransCost',option=opt))      #Crop transport cost matrix [ctrans x crop] ($/t)
        if md.Options['Crop market'] == 0:
            md.aFarmVal     = Param(md.ncountry, md.ncrop, initialize=ReadParam('aFarmVal',option=md.Options['Agriculture']))       #Value of crop at farm level [country x crop] ($/t)
#%%**Energy production Module**#
        opt=md.Options['Energy production']
    #Hydropower
        md.nhpp             = Set(initialize=ReadParam('eHppCap',option=opt).keys())            #HydroPower Plants [indice hpp] (id)  
        md.hp_res           = Param(md.nhpp, initialize=ReadParam('hp_res',option=opt))         #Reservoir of HP turbines [hpp], value='ROR' if Run-Of-the-River HP
        md.hp_catch         = Param(md.nhpp, initialize=ReadParam('hp_catch',option=opt))       #Catchment of Run-Of-the-River HP turbines [hpp]
        md.hp_pmarket       = Param(md.nhpp, initialize=ReadParam('hp_pmarket',option=opt))     #Power pmarket of HP turbines [hpp]
        md.eHppProd         = Param(md.nhpp, initialize=ReadParam('eHppProd',option=opt))       #Production factor: energy produced per water flow of HP [hpp] (kWh/m3)
        md.eHppCap          = Param(md.nhpp, initialize=ReadParam('eHppCap',option=opt))        #Capacity of HP turbines [hpp] (MW)
        md.eHppEff          = Param(md.nhpp, initialize=ReadParam('eHppEff',option=opt))        #Efficiency of HP turbines [hpp] (%)
        md.eHppCost         = Param(md.nhpp, initialize=ReadParam('eHppCost',option=opt))       #Operational production cost of HP [hpp] ($/kWh)
    #Power plants and fuels
        opt= opt and md.Options['Power plants']
        md.nopp             = Set(initialize=ReadParam('eOppCap',option=opt).keys())            #Other Power Plants (OPP) [indice opp] (id)
        md.op_pmarket       = Param(md.nopp, initialize=ReadParam('op_pmarket',option=opt))     #pmarket of OPP [opp] (id)
        md.eOppCost         = Param(md.nopp, initialize=ReadParam('eOppCost',option=opt))       #Operational production cost of OPP [opp] ($/kWh)
        md.eOppCap          = Param(md.nopp, initialize=ReadParam('eOppCap',option=opt))        #Production capacity of OPP [opp] (MW)        
        if md.Options['Power technologies'] == 1 or md.Options['Load capacity'] == 1:
            md.op_ptech     = Param(md.nopp, initialize=ReadParam('op_ptech',option=opt))       #Technology of OPP [opp] (id)
        if md.Options['Ramping'] == 1:
            md.eOppRamp     = Param(md.nopp, initialize=ReadParam('eOppRamp',option=opt))       #Ramping rate of OPP [opp] (%/load segment) 
        if md.Options['Fuels'] ==1:
            md.eOppEff      = Param(md.nopp, initialize=ReadParam('eOppEff',option=opt))        #Efficiency of OPP [opp] (kWh_net/kWh_fuel)            
            md.op_fuel      = Param(md.nopp, initialize=ReadParam('op_fuel',option=opt))        #Fuel of OPP [opp] (id)                 
   
#%%**Energy market Module**#
        opt=md.Options['Energy production'] and md.Options['Energy market']
    #Power markets - Demands and value
        md.npmarket         = Set(initialize=ReadParam('pmarket_country',option=opt).keys())             #pmarkets [indice pmarket] (id)            
        md.pmarket_country  = Param(md.npmarket, initialize=ReadParam('pmarket_country',option=opt))     #Country of pmarket [pmarket] (id)            
        md.eSupLoss         = Param(md.npmarket, initialize=ReadParam('eSupLoss',option=opt))            #Local supply losses of pmarket [pmarket] (-)
        md.eEngyDem         = Param(md.nyear, md.nmonth, md.npmarket, initialize=ReadParam('eEngyDem',option=opt,time=md.nyear))  #Local power demand [month x pmarket] (GWh/year)
        md.eEngyVal         = Param(md.nmonth, md.npmarket, initialize=ReadParam('eEngyVal',option=opt)) #Curtailment cost (=Marginal value) of power [time x pmarket] ($/kWh)           
    #Transmission lines
        opt= opt and md.Options['Transmission']
        md.ntransline       = Set(initialize=ReadParam('eTransCap',option=opt).keys())              #Power transmission lines [indice tmsline] (id)
        md.eTransIn         = Param(md.ntransline, initialize=ReadParam('eTransIn',option=opt))     #"Exportin" power market [tmsline] (id)
        md.eTransOut        = Param(md.ntransline, initialize=ReadParam('eTransOut',option=opt))    #"Importing" power market [tmsline] (id)
        md.eTransLoss       = Param(md.ntransline, initialize=ReadParam('eTransLoss',option=opt))   #Transmition losses [tmsline] (%)
        md.eTransCap        = Param(md.ntransline, initialize=ReadParam('eTransCap',option=opt))    #Transmition capacity [tmsline] (MW)
        md.eTransCost       = Param(md.ntransline, initialize=ReadParam('eTransCost',option=opt))   #Transmition costs [tmsline] ($/kWh)        
    #Fuels and ressources
        opt= md.Options['Energy production'] and md.Options['Energy market'] and md.Options['Fuels']               
        md.nfuel            = Set(initialize=ReadParam('eFuelCO2',option=opt).keys())                   #Fuels [indice fuel] (id)
        md.eFuelCost        = Param(md.nyear, md.npmarket, md.nfuel, initialize=ReadParam('eFuelCost',option=opt,time=md.nyear)) #Cost of fuel [fuel] ($/kWh_fuel)
        md.eFuelCO2         = Param(md.nfuel, initialize=ReadParam('eFuelCO2',option=opt))              #CO2 emission factor of fuel [fuel] (t/kWh_fuel)
        md.eCO2Val          = Param(md.nyear, md.npmarket, initialize=ReadParam('eCO2Val',option=opt,time=md.nyear) if md.Options['Fuels'] == 1 else {(y,cm):0 for y in md.nyear for cm in md.ncmarket})  #CO2 price [1] ($/t)
    #Generic capacity investment
        opt= md.Options['Energy production'] and md.Options['Energy market'] and md.Options['Power technologies']             
        md.nptech           = Set(initialize=set([key[0] for key in ReadParam('eVarOPEX',option=opt).keys()]))  #Power technologies [indice ptech] (id)
        md.eCAPEX           = Param(md.nyear, md.nptech, md.npmarket, initialize=ReadParam('eCAPEX',option=opt,time=md.nyear))#Cost of Generic power capacity [ptech x pmarket] (id) ($/MW)
        md.eVarOPEX         = Param(md.nptech, md.npmarket, initialize=ReadParam('eVarOPEX',option=opt))        #Variable Cost of Generic power production [ptech x pmarket] ($/kWh) 
        md.eFixOPEX         = Param(md.nptech, md.npmarket, initialize=ReadParam('eFixOPEX',option=opt))        #Fix Cost of Generic power production [ptech x pmarket] ($/[kWh/day] /year)
        md.eLifeTime        = Param(md.nptech, md.npmarket, initialize=ReadParam('eLifeTime',option=opt))       #Life time of power technology [ptech x pmarket] (years)
        md.eMaxCap          = Param(md.nptech, md.npmarket, initialize=ReadParam('eMaxCap',option=opt))         #Maximum expandable capacity of power technology [ptech x pmarket] (kWh/day)
    #Fuels
        if md.Options['Fuels'] == 1: 
            md.eTechEff     = Param(md.nptech, md.npmarket, initialize=ReadParam('eTechEff',option=opt))     #Efficiency of power technology [ptech x pmarket] (%)
            md.ptech_fuel   = Param(md.nptech, md.npmarket, initialize=ReadParam('ptech_fuel',option=opt))   #Fuel of Technology [ptech x pmarket] (id)
        if md.Options['Ramping'] == 1:
            md.eTechRamp    = Param(md.nptech, md.npmarket, initialize=ReadParam('eTechRamp',option=opt))    #Ramping rate of Technology [ptech x pmarket] (%/load segment)           
    #Hydropower valuation if power market is off
        if md.Options['Energy market'] == 0:
            md.hp_country   = Param(md.nhpp, initialize=ReadParam('hp_country',option=md.Options['Energy production']))  #Country of HP if power market is off [hpp] (id)
            md.eHppVal      = Param(md.nhpp, initialize=ReadParam('eHppVal',option=md.Options['Energy production']))     #Value of HP production if power market is off [hpp] ($/kWh)
        else:
            md.hp_country   = Param(md.nhpp, initialize={hp:md.pmarket_country[md.hp_pmarket[hp]] for hp in md.nhpp})
    #Power load segments
        opt= md.Options['Energy production'] and md.Options['Energy market']
        md.npload           = Set(initialize=ReadParam('eLoadDem',option=opt).keys() if md.Options['Load'] == 1 or opt == 1 else [1])            #Slices of power demand per time: day week, day week end, night [indice pload] (id)
        md.eLoadDem         = Param(md.npload, initialize=ReadParam('eLoadDem',option=opt) if md.Options['Load'] == 1 or opt == 1 else {1:1})    #Share of demand per power load [pload] (-)
        md.eLoadTime        = Param(md.npload, initialize=ReadParam('eLoadTime',option=opt) if md.Options['Load'] == 1 or opt == 1 else {1:1})   #Share of time per power load [pload] (-)
        if md.Options['Load capacity'] == 1: #Assumes either "Other power plants" or "Power technologies" are ON
            if md.Options['Power technologies'] == 0:
                md.op_ptech = Param(md.nopp, initialize=ReadParam('op_ptech',option=opt))                      #Technology of Other power plant [pp] (id)
                md.nptech   = Set(initialize=set([ReadParam('op_ptech',option=opt)[op] for op in md.nopp]))    #Technologies are defined in Other power plants if Power technologies module is off [indice ptech] (id)
            md.eLoadCap     = Param(md.npload, md.nptech, initialize=ReadParam('eLoadCap',option=opt))         #Load segment capacity factor of power technology [pload x ptech] (%)

#%% **Non linear hydro-power module**
        if md.Options['Hydropower'] == 'nonlinear' and md.Options['Energy production'] == 1:
            md.wResMOL      = Param(md.nres, initialize=ReadParam('wResMOL'))
            md.wResFSL      = Param(md.nres, initialize=ReadParam('wResFSL'))
            md.wMinTurb     = Param(md.nhpp, initialize=ReadParam('wMinTurb'))
            md.wMaxTurb     = Param(md.nhpp, initialize=ReadParam('wMaxTurb'))
            md.wMinHead     = Param(md.nres, initialize=ReadParam('wMinHead'))
            md.wMaxHead     = Param(md.nres, initialize=ReadParam('wMaxHead'))
            
            def RelHeadVol(Vol,MOL,FSL,hmin,hmax):                
                return hmin/hmax + (1-hmin/hmax)*(Vol-MOL)/(FSL-MOL)
            
            def TurbFlowVol(Vol,MOL,FSL,qmin,qmax):
                return qmin + (qmax-qmin)*(Vol-MOL)/(FSL-MOL)
        
#%%DECLARE DECISION VARIABLES###       
    #Water Module#
        md.WwSUPPLY         = Var(md.ntime, md.nuser, within=NonNegativeReals)    #Water user allocation  [time x user] (Mm3)
        md.WwOUTFLOW        = Var(md.ntime, md.ncatch, within=NonNegativeReals)   #Catchment outflow [time x catchment] (km3)
        md.WwRSTORAGE       = Var(md.ntime, md.nres, within=NonNegativeReals)     #Reservoir total storage [time x reservoir] (km3)            
        md.WwGWSTORAGE      = Var(md.ntime, md.naquifer, within=NonNegativeReals)  #Groundwater reservoir total storage [time x aquifer] (Mm3)
        md.WwTRANSFER       = Var(md.ntime, md.ntransfer, within=NonNegativeReals) #Transfer scheme water transfer [time x transfer] (Mm3)                 
    
    #Agriculture Module#
        if md.Options['Crop choice'] == 'fixed': #Cultivated area is fixed
            md.CULAREA      = Param(md.nyear, md.nfzone, md.nculture, mutable=True, initialize=ReadParam('aCulMax',option=md.Options['Agriculture'],time=md.nyear)) #Fixed cultivated area [year x fzone x culture] (-)
            #md.CULAREA      = Param(md.nyear, md.nfzone, md.nculture, initialize=ReadParam('aCulMax',option=md.Options['Agriculture'],time=md.nyear))
        elif md.Options['Yield water response'] == 'nonlinear':
            md.AlCULAREA    = Var(md.nyear, md.nfzone, md.nfieldculture, within=NonNegativeReals)    #Cultivated area  [year x fzone x field] (-)
        else: #Cultivated area is a decision variable
            md.AlCULAREA    = Var(md.nyear, md.nfzone, md.nfieldculture, md.nypath, within=NonNegativeReals)    #Cultivated area and water demand satisfaction [year x fzone x field x ypath] (-)
        md.AwSUPPLY         = Var(md.ntime, md.nifzone, md.nculture, within=NonNegativeReals)                   #Agricultural water allocation [time x ifzone x culture] (Mm3)
        md.AcPROD           = Var(md.nyear, md.nfzone, md.ncrop, within=NonNegativeReals)                       #Crop Production [year x fzone x crop] (Mt)       
        if md.Options['Groundwater'] == 1:
            md.AwGWSUPPLY   = Var(md.ntime, md.nifzone, md.nculture, within=NonNegativeReals)                   #Agricultural ground water allocation [time x ifzone x culture] (Mm3)
            
    #Crop market module#
        md.AcTRANS          = Var(md.nyear, md.nctrans, md.ncrop, within=NonNegativeReals)              #Crop Transport [year x ctrans x crop] (Mt)
        md.AcSUPPLY         = Var(md.nyear, md.ncmarket, md.ncrop, md.ncdstep, within=NonNegativeReals) #Crop Supply [year x cmarket x crop x cdstep] (Mt) 
        md.AcEXTPROD        = Var(md.nyear, md.nextcmarket, md.ncrop, within=NonNegativeReals)          #Crop external Production [year x cmarket x crop] (Mt)
        
    #Energy production module#
        md.EwHPDISCHARGE    = Var(md.ntime, md.npload, md.nhpp, within=NonNegativeReals)    #Dicharge through hydropower turbine [time x pload x hpp] (Mm3)
        md.EeHPPROD         = Var(md.ntime, md.npload, md.nhpp, within=NonNegativeReals)    #Power production by hydropower turbine [time x pload x hpp] (GWh)
        md.EeOPPROD         = Var(md.ntime, md.npload, md.nopp, within=NonNegativeReals)    #Power production by power plant [time x pload x pp] (GWh)        
        md.EeFIRMHP         = Var(md.nyear, md.nhpp, within=NonNegativeReals)
        md.eHppFirm         = Param(md.nhpp, initialize=ReadParam('eHppFirm',option=md.Options['Energy production']))
    #Energy market module#
        md.EeTRANS          = Var(md.ntime, md.npload, md.ntransline, within=NonNegativeReals)          #Power transmission between power markets [time x pload x tmsline] (GWh)
        md.EeSUPPLY         = Var(md.ntime, md.npload, md.npmarket, within=NonNegativeReals)            #Power supply to power market demand [time x pload x pmarket] (GWh)
        md.EeGENCAP         = Var(md.nyear, md.nptech, md.npmarket, within=NonNegativeReals)            #Generic Power technologies capacity [year x ptech x pmarket] (MW)
        md.EeGENPROD        = Var(md.ntime, md.npload, md.nptech, md.npmarket, within=NonNegativeReals) #Generic Power Tech Production [time x pload x ptech x pmarket] (GWh/load seg)            

#%%DEBUG VARIABLES
        
        if md.Options['Initial time step'] == 1:
            md.DUMYSTOR     = Var(md.nres, within=NonNegativeReals) #Dummy variable to avoid infeasible solutions regarding final storage objective
        if md.Options['Debug mode'] == 1:
            md.DUMYWATER    = Var(md.ntime, md.ncatch, within=NonNegativeReals) #DEBUG Variable adding water in the water balance (m3/month)
        if md.Options['Debug mode'] in [1,2]:
            md.DUMYCROP     = Var(md.nyear, md.nfzone, md.ncrop, within=NonNegativeReals)#, bounds=(0,1))   #Dummy Crop Production to compensate for negative prod [year x fzone x crop] (kt/year)       
            md.DUMYEFLOW    = Var(md.ntime, md.neflow, within=NonNegativeReals) #Tolerance to eflow objective if physically impossible (Mm3/month)

#%%INVESTMENT PLANNING MODULE        
        if md.Options['Investment module'] == 1:                        
        #Declare indices
            md.ninvest      = Set(initialize=ReadParam('iCAPEX').keys())                          #Investments [indice inv] (id)
            md.ninvphase    = Set(initialize=ReadParam('invphase_t').keys())                      #Investment phases [indice invphase] (id)
            
        #Declare parameters
            md.iDisRate     = Param(md.nyear, initialize=ReadParam('iDisRate',time=md.nyear))      #Discounting factor [year] (-)
            md.invphase_t   = Param(md.ninvphase, initialize=ReadParam('invphase_t'))              #time step of investment phase [invphase] (time)
            md.inv_after    = Param(md.ninvest, initialize=ReadParam('inv_after'))                 #Investment can occur only after another investment [inv] (id)
            md.inv_group    = Param(md.ninvest, initialize=ReadParam('inv_group'))                 #Investment group have to occur at the same time [inv] (id) #ADD: only works with 2 investments so far
            md.inv_replace  = Param(md.ninvest, initialize=ReadParam('inv_replace'))               #Investment replaces existing infrastructure [inv] (id) #ADD only active for hydropower now
            
            md.iMaxInv      = Param(md.ninvphase, initialize=ReadParam('iMaxInv'))                 #Maximum investment amount of investment phase [invphase] ($)
            md.iCAPEX       = Param(md.ninvest, initialize=ReadParam('iCAPEX'))                    #investment capital costs [inv] ($)
            md.iFixOPEX     = Param(md.ninvest, initialize=ReadParam('iFixOPEX'))                  #investment fix operational costs [inv] ($)
            md.iInvType     = Param(md.ninvest, initialize=ReadParam('iInvType'))                  #type of investment [invtype] (type)
            md.iInvCap      = Param(md.ninvest, initialize=ReadParam('iInvCap'))                   #capacity of investment [inv] (different units)
            md.iConstTime   = Param(md.ninvest, initialize=ReadParam('iConstTime'))                #implementation time of investment [inv] (month)
            md.iLifeTime    = Param(md.ninvest, initialize=ReadParam('iLifeTime'))                 #life time of investment after construction [inv] (month)
        
        #Declare Decision-Variables    
            md.IbINVEST     = Var(md.ninvphase, md.ninvest, within=Binary)                         #Investments [invphase x inv] (binary)        
        
        ##general functions
        def InvestedCapacity(currentcap,iobject,itype,timestep,year=0):
            if md.Options['Investment module'] == 1: #Investment module is on
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
                                    and (md.invphase_t[kip]+md.iConstTime[inv] <= md.Options['tfin'] and md.t_year[md.invphase_t[kip]+md.iConstTime[inv]] <= timestep) 
                                    and (md.invphase_t[kip]+md.iConstTime[inv]+md.iLifeTime[inv] >= md.Options['tfin'] or md.t_year[md.invphase_t[kip]+md.iConstTime[inv]+md.iLifeTime[inv]] >= timestep))                                   
                    ReplacedCap=sum(min(currentcap,md.iInvCap[kinv])*md.IbINVEST[kip,kinv] for kip in md.ninvphase for kinv in md.ninvest 
                                    if  (md.iInvType[kinv]==itype and md.inv_replace[kinv]==iobject)
                                    and (md.invphase_t[kip]+md.iConstTime[kinv] <= md.Options['tfin'] and md.t_year[md.invphase_t[kip]+md.iConstTime[kinv]] <= timestep)                                                                                                                     
                                    and (md.invphase_t[kip]+md.iConstTime[kinv]+md.iLifeTime[kinv] >= md.Options['tfin'] or md.t_year[md.invphase_t[kip]+md.iConstTime[kinv]+md.iLifeTime[kinv]] >= timestep))
                return InvestedCap-ReplacedCap
            else:
                return 0                
#%%OBJECTIVE FUNCTION - ECONOMIC MODULE
        def obj_rule(md):
        #Water module            
            UserBenefit         = {y:sum(md.WwSUPPLY[t,u] * md.wUserVal[u] for u in md.nuser for t in md.ntime if md.t_year[t]==y) for y in md.nyear} #benefit from (non agricultural) water users
            UserSupCosts        = {y:sum(md.WwSUPPLY[t,u]/(1-md.wUserLoss[u]) * md.wSupCost[u] for u in md.nuser for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
            UserBalance         = {y:UserSupCosts[y]-UserBenefit[y] for y in md.nyear}
            
        #Agricultural module: Cultivation costs of cultures (labour, machinery, fertilizers, infrastructure....), Irrigation costs of cultures (pumping ...)     
            if md.Options['Crop choice'] == 'fixed': #ADD REM: for fixed area can be removed from OBJ function
                AgrCulCost      = {y:sum(md.kha_to_Mha*md.CULAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for fz in md.nfzone for cul in md.nculture for ypt in md.nypath) for y in md.nyear}
            elif md.Options['Yield water response'] == 'nonlinear':
                AgrCulCost      = {y:sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd] * md.aCulCost[md.fzone_type[fz],md.field_culture[fd]] for fz in md.nfzone for fd in md.nfieldculture) for y in md.nyear}
            else:
                AgrCulCost      = {y:sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt] * md.aCulCost[md.fzone_type[fz],md.field_culture[fd]] for fz in md.nfzone for fd in md.nfieldculture for ypt in md.nypath) for y in md.nyear}
            AgrIrrgCost         = {y:sum(md.AwSUPPLY[t,fz,cul]/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.t_year[t]==y) for y in md.nyear}
            GwPumpingCost       = {y:sum(md.AwGWSUPPLY[t,fz,cul] * md.wGwCost[aq] for t in md.ntime for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.t_year[t]==y and md.fzone_catch[fz]==md.aqui_catch[aq]) for y in md.nyear}
            AgricultureBalance  = {y:AgrCulCost[y] + AgrIrrgCost[y] + GwPumpingCost[y] for y in md.nyear}
        
        #Crop market module: Crop allocated to crop markets, Cost and benefits of imported crops (external market), Transport costs    
            SF=1
            if md.Options['Minimum supply'] == 'sufficiency': #ADD more elegant solution ?
                SF=10         
            AgrBenefit          = {y:sum(md.kt_to_Mt*md.AcSUPPLY[y,cm,cr,cds] * md.aCropVal[cm,cr] * md.aStepVal[cm,cr,cds] * SF for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep) for y in md.nyear}    
            if md.Options['Crop market'] == 0:
                AgrBenefit      = {y:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr] * md.aFarmVal[md.fzone_country[fz],cr] for fz in md.nfzone for cr in md.ncrop) for y in md.nyear}
            AgrTransCost        = {y:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr] * md.aTransCost[ct,cr] * SF for ct in md.nctrans for cr in md.ncrop) for y in md.nyear}
            ExtProdCost         = {y:sum(md.kt_to_Mt*md.AcEXTPROD[y,cm,cr] * md.aCropVal[cm,cr] * SF * 1.0001 for cm in md.nextcmarket for cr in md.ncrop) for y in md.nyear} #The assumption is that external market produce crops at their market price + 1% #ADD: crop supply curve
            CropMarketBalance   = {y:AgrTransCost[y] + ExtProdCost[y] - AgrBenefit[y] for y in md.nyear}
      
        #Energy production module: Production costs of hydropower and power plants (O&M, fuel ...)    
            HpOMCost            = {y:sum(md.eHppCost[hp]*md.se*md.EeHPPROD[t,pld,hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}            
            OppOMCost           = {y:sum(md.eOppCost[opp]*md.se*md.EeOPPROD[t,pld,opp] for t in md.ntime for pld in md.npload for opp in md.nopp if md.t_year[t]==y) for y in md.nyear}    
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
                Firmbenefit     = {y:sum(md.se*md.EeFIRMHP[y,hp]*md.eHppFirm[hp] for t in md.ntime for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}
                Scndbenefit     = {y:sum(md.se*(sum(md.EeHPPROD[t,pld,hp] for pld in md.npload)-md.EeFIRMHP[y,hp])*md.eHppVal[hp] for t in md.ntime for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}
                EngyBenefit     = {y:Firmbenefit[y]+Scndbenefit[y] for y in md.nyear}
            EngyTransCost       = {y:sum(md.se*md.EeTRANS[t,pld,tl] * md.eTransCost[tl]  for t in md.ntime for pld in md.npload for tl in md.ntransline if md.t_year[t]==y) for y in md.nyear}
            EnergyMarketBalance = {y:EngyTransCost[y]-EngyBenefit[y] for y in md.nyear} #            
        
        #DEBUG MODE
            #DEBUGCOST = 1000 #arbitrary high cost for additional water $/m3 or deviance from final storage objective
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
            if md.Options['Investment module'] == 1:
                #Investment costs + discounting - #ADD in this case the investments costs are assumed to be proportional to the use of the infrastructure.
                CAPEX   = {y:sum(md.IbINVEST[ip,inv]*md.iCAPEX[inv] for ip in md.ninvphase for inv in md.ninvest if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear}                   
                CAPEXrem= sum(sum(md.IbINVEST[ip,inv]*md.iCAPEX[inv]*max(0,1-(md.Options['tfin']-md.invphase_t[ip]-md.iConstTime[inv])/md.iLifeTime[inv]) for ip in md.ninvphase for inv in md.ninvest if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear)
                CAPEX[md.t_year[md.Options['tfin']]] += -CAPEXrem #remaining capital value at end of simulation period
                fixOPEX = {y:sum(InvestedCapacity(0,inv,md.iInvType[inv],y,year=1)/md.iInvCap[inv]*md.iFixOPEX[inv] for ip in md.ninvphase for inv in md.ninvest if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear}
                return DebugCost + sum(md.iDisRate[y]*(CAPEX[y] + fixOPEX[y] + UserBalance[y] + AgricultureBalance[y] + CropMarketBalance[y] + EnergyProdBalance[y] + EnergyMarketBalance[y]) for y in md.nyear)
            else: #No discounting               
                return DebugCost + sum(UserBalance[y] + AgricultureBalance[y] + CropMarketBalance[y] + EnergyProdBalance[y] + EnergyMarketBalance[y] for y in md.nyear)  
        md.obj = Objective(rule=obj_rule)
        
#%%CONSTRAINTS
#%%-----------------------------Investment planning module---------------------

        ##Constraint catalogue##
        def invest_once(md,ninv):
            return sum(md.IbINVEST[kip,ninv] for kip in md.ninvphase) <= 1
       
        def invest_max(md,nip):
            return sum(md.IbINVEST[kip,kinv]*md.iCAPEX[kinv] for kinv in md.ninvest for kip in md.ninvphase if md.invphase_t[kip] <= md.invphase_t[nip]) <= sum(md.iMaxInv[kip] for kip in md.ninvphase if md.invphase_t[kip] <= md.invphase_t[nip])
              
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
        if md.Options['Investment module'] == 1:
            md.invest_once  = Constraint(md.ninvest,                rule=invest_once)
            md.invest_max   = Constraint(md.ninvphase,              rule=invest_max)
            md.invest_after = Constraint(md.ninvphase, md.ninvest,  rule=invest_after)
#%%----------------------------Water module----------------------------
        ##Constraint catalogue##
        def water_waterbalance(md,nt,nc): #Water Balance per Watershed [catchment x time]
            nres=[kres for kres in md.nres if md.res_catch[kres]==nc]
            nifzone=[kfz for kfz in md.nifzone if md.fzone_catch[kfz]==nc]
            naquifer=[kaq for kaq in md.naquifer if md.aqui_catch[kaq]==nc]
            NetInflow           = (1-md.wFlowLoss[nc]) * sum(md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) #Net incoming flow from upstream catchments (Mm3)
            RunOff              = md.wRunOff[nt,nc]   #Catchment Runoff (Mm3) 
            GwStorage           = sum(md.WwGWSTORAGE[md.t_prev[nt],kaq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[kaq] for kaq in naquifer) #Previous ground water storage
            Pumping             = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in md.nifzone for kcul in md.nculture) if md.Options['Groundwater']==1 else 0 #Water abstraction from groundwater reservoir
            BaseFlow            = sum(GwStorage*(1-exp(-md.wGwFlow[kaq])) + (md.wGwRech[nt,kaq]-Pumping)*(1-(1-exp(-md.wGwFlow[kaq]))/md.wGwFlow[kaq]) for kaq in naquifer) #Baseflow from groundwater reservoir
            UserCons            = sum(md.WwSUPPLY[nt,ku] * (1/(1-md.wUserLoss[ku])-md.wUserRturn[ku]) for ku in md.nuser if md.user_catch[ku] == nc)      #Consumption of non agricultural users (Mm3)
            AgrCons             = sum(md.AwSUPPLY[nt,kfz,kcu] * (1/(1-md.aIrrgLoss[kfz])-md.aCulRturn[md.fzone_type[kfz],kcu]) for kfz in nifzone for kcu in md.nculture)  #Agricultural module: Consumptiom to agriculture/crops (Mm3)    
            if md.Options['Initial time step'] == 1 and md.Options['tini'] == nt: #Initial storage for initial time step
                DeltaStorage    = sum(md.sw*md.WwRSTORAGE[nt,kres] - md.wStorIni[kres] for kres in nres) #Storage in reservoir (Mm3)
                Evaporation     = max(0,value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc])) * md.mm_to_m3perha * sum(md.wkV[kres]*(md.sw*md.WwRSTORAGE[nt,kres]+md.wStorIni[kres])/2 + md.wResArea[kres] for kres in nres)
            else: #Cyclic run or not first time step 
                DeltaStorage    = sum(md.sw*md.WwRSTORAGE[nt,kres] - md.sw*md.WwRSTORAGE[md.t_prev[nt],kres] for kres in nres) #Storage in reservoir (Mm3)
                Evaporation     = max(0,value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc])) * md.mm_to_m3perha * sum(md.wkV[kres]*(md.sw*md.WwRSTORAGE[nt,kres]+md.sw*md.WwRSTORAGE[md.t_prev[nt],kres])/2 + md.wResArea[kres] for kres in nres)                
            TransNetInflow      = sum(md.WwTRANSFER[nt,ktrans] * (1-md.wTransLoss[ktrans])   for ktrans in md.ntransfer if md.transfer_ds[ktrans] == nc)  #Incoming flow from upstream transfer scheme- if exists (Mm3)
            TransGrossOutflow   = sum(md.WwTRANSFER[nt,ktrans]                               for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc)  #Allocation to downstream transfer scheme- if exists (Mm3)                   
            CatchmentOutflow    = md.sw*md.WwOUTFLOW[nt,nc]  #Catchment outflow : Allocation to downstream catchment (Mm3)           
            DEBUG=md.DUMYWATER[nt,nc] if md.Options['Debug mode'] == 1 else 0
            return  UserCons + AgrCons + TransGrossOutflow + DeltaStorage + Evaporation + CatchmentOutflow == NetInflow + TransNetInflow + RunOff + BaseFlow + DEBUG
                
        def water_waterbalance2(md,nt,nc): #Avoids allocation from downstream reservoir to upstream demand [catchment x time] 
            if  sum(1 for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc) + sum(1 for ku in md.nuser if md.user_catch[ku] == nc) + sum(1 for kfz in md.nfzone if md.fzone_catch[kfz]==nc) == 0:
                return Constraint.Skip #if there is no demand and no transfer project constraint would be empty and is therefore skipped            
            nifzone=[kfz for kfz in md.nifzone if md.fzone_catch[kfz]==nc]
            naquifer=[kaq for kaq in md.naquifer if md.aqui_catch[kaq]==nc]
            BaseFlow            = 0
            NetInflow           = (1-md.wFlowLoss[nc]) * sum(md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) #Net incoming flow from upstream catchments (Mm3)
            RunOff              = md.wRunOff[nt,nc] #Catchment Inflow (Mm3) 
            GwStorage           = sum(md.WwGWSTORAGE[md.t_prev[nt],kaq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[kaq] for kaq in naquifer) #Previous ground water storage
            Pumping             = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in md.nifzone for kcul in md.nculture) if md.Options['Groundwater']==1 else 0 #Water abstraction from groundwater reservoir
            BaseFlow            = sum(GwStorage*(1-exp(-md.wGwFlow[kaq])) + (md.wGwRech[nt,kaq]-Pumping)*(1-(1-exp(-md.wGwFlow[kaq]))/md.wGwFlow[kaq]) for kaq in naquifer) #Baseflow from groundwater reservoir
            UserWithdrawal      = sum(md.WwSUPPLY[nt,ku] * 1/(1-md.wUserLoss[ku]) for ku in md.nuser if md.user_catch[ku] == nc) #Allocation to non agricultural users (Mm3)
            AgrWithdrawal       = sum(md.AwSUPPLY[nt,kfz,kcu]* 1/(1-md.aIrrgLoss[kfz]) for kfz in nifzone for kcu in md.nculture ) #Agricultural module: Allocation to agriculture/crops (Mm3)    
            TransNetInflow      = sum(md.WwTRANSFER[nt,ktrans] * (1-md.wTransLoss[ktrans])   for ktrans in md.ntransfer if md.transfer_ds[ktrans] == nc)  #Incoming flow from upstream transfer scheme- if exists (Mm3)
            TransGrossOutflow   = sum(md.WwTRANSFER[nt,ktrans]                               for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc)  #Allocation to downstream transfer scheme- if exists (Mm3)                   
            DEBUG=md.DUMYWATER[nt,nc] if md.Options['Debug mode'] == 1 else 0
            return UserWithdrawal + AgrWithdrawal + TransGrossOutflow <= NetInflow + TransNetInflow + RunOff + BaseFlow + DEBUG
                
        def water_usermaxall(md,nt,nu): #Water allocation to non agricultural user limited by demand [time x user] 
            return md.WwSUPPLY[nt,nu] <= md.wUserDem[md.t_year[nt],md.t_month[nt],nu]
                        
        def water_transfercapacity(md,nt,ntrans): #Transfer scheme carrying capacity [time x transfer]          
            return md.WwTRANSFER[nt,ntrans] <= md.wTransCap[ntrans] + InvestedCapacity(md.wTransCap[ntrans],ntrans,'transfer',nt)
        
        def water_lakebalance(md,nt,nla): #Storage controlled by first order reservoir or wetland model [time x reservoir]
            nc = md.res_catch[nla]
            AvStor      = (md.sw*md.WwRSTORAGE[nt,nla] + md.sw*md.WwRSTORAGE[md.t_prev[nt],nla])/2 if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else (md.sw*md.WwRSTORAGE[nt,nla] + md.wStorIni[nla])/2
            StorPrev    = md.sw*md.WwRSTORAGE[md.t_prev[nt],nla] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wStorIni[nla]
            Inflow      = sum((1-md.wkUP[nla])*md.sw*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) * (1-md.wFlowLoss[nc])
            RunOff      = (1-md.wkRO[nla]) * md.wRunOff[nt,nc]
            Evap        = (md.wET0[nt,nc]-md.wRainFall[nt,nc]) * md.mm_to_m3perha * (md.wkV[nla]*AvStor + md.wResArea[nla]) #max evaporation is 2/kV to avoid impossible solution (which corresponds to the maximum evaporation if no inflow comes in)
            Outflow     = md.wAlpha[nla] * AvStor
            #Analytical version
            #a=(md.wET0[nt,nc]-md.wRainFall[nt,nc])*md.wkV[nla]+md.wAlpha[nla]
            #b=Inflow+RunOff-md.wResArea[nla]*(md.wET0[nt,nc]-md.wRainFall[nt,nc])
            return md.sw*md.WwRSTORAGE[nt,nla] == StorPrev + Inflow + RunOff - Evap - Outflow #+ DEBUG   
            #return md.sw*md.WwRSTORAGE[nt,nla] == StorPrev*exp(-a) + b/a*(1-exp(-a))
        
        def water_groundwaterbalance(md,nt,naq):
            GwStorage       = md.WwGWSTORAGE[md.t_prev[nt],naq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[naq]
            Pumping         = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in md.nifzone for kcul in md.nculture if md.fzone_catch[kfz]==md.aqui_catch[naq])
            if md.wGwFlow[naq] != 0:
                return md.WwGWSTORAGE[nt,naq] == GwStorage*exp(-md.wGwFlow[naq]) + (md.wGwRech[nt,naq]-Pumping)/md.wGwFlow[naq]*(1-exp(-md.wGwFlow[naq]))
            else:
                return md.WwGWSTORAGE[nt,naq] == GwStorage - Pumping
            
        def water_reservoircapacity(md,nt,nres): #Storage in reservoir limited by storage capacity [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir': #non man-operated reservoirs
                return Constraint.Skip
            return md.sw*md.WwRSTORAGE[nt,nres] <= md.wStorCap[nres] + InvestedCapacity(md.wStorCap[nres],nres,'reservoir',nt)
        
        def water_floodrulecurve(md,nt,nres): #Storage in reservoir limited by flood rule curve [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir': #non man-operated reservoirs
                return Constraint.Skip
            return md.sw*md.WwRSTORAGE[nt,nres] <= md.wFloodRule[nres,md.t_month[nt]]
        
        def water_finalstorage(md,nres): #Storage in reservoir limited by storage capacity [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir':
                return Constraint.Skip
            Debug = md.DUMYSTOR[nres] if md.Options['Initial time step'] == 1 else 0 #Tolerance to constraint if physically impossible at very high cost
            #if md.Options['Reservoirs'] == 'flood':
            #    if md.wStorFin[nres] > md.wFloodRule[nres,md.t_month[md.Options['tfin']]]:
            #        print('WARNING: unfeasible constraint as flood rule curve is lower than final storage target')
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

        if md.Options['Hard eflows'] != 1: 
        #Find most severe years to allow softening of constrain in these years
            years       = [y for y in md.nyear]
            YearlyRunoff= [sum(md.wRunOff[t,c].value for t in md.ntime for c in md.ncatch if md.t_year[t]==y) for y in md.nyear]
            YearsRank   = sorted(range(len(YearlyRunoff)), key=lambda k: YearlyRunoff[k]) #Years ranked from minimal to maximal inflows
            nbSoftYears = int((1-md.Options['Hard eflows']) * len(md.nyear)) #Amount of years where the hard constraint can be softened in the optimization period
            SoftYears   = {years[ky]:ky in YearsRank[0:nbSoftYears] for ky in range(len(md.nyear))} #Years where the hard constraint is softened
            
        def env_minflow(md,nt,nef): #env minimum flow [catchment x time]
            if md.eflow_hard[nef] == 0: #If no upstream reservoir is regulating the flow the env flow is limited to actual flow
                EflowConstraint = min(0.95*NatOutflow[nt][md.eflow_catch[nef]],md.wEnvFlow[md.t_month[nt],nef]) #Arbitrary 5% tolerance on consumption of available water
            elif md.Options['Hard eflows'] != 1 and SoftYears[md.t_year[nt]] == 1:
                EflowConstraint = min(0.95*NatOutflow[nt][md.eflow_catch[nef]],md.wEnvFlow[md.t_month[nt],nef]) #ADD 0.95 is arbitrary to allow some minimal losses on the way
            else:
                EflowConstraint = md.wEnvFlow[md.t_month[nt],nef]
            Debug = md.DUMYEFLOW[nt,nef] if md.Options['Debug mode'] in [1,2] else 0 #Tolerance to constraint if physically impossible at very high cost
            return md.sw*md.WwOUTFLOW[nt,md.eflow_catch[nef]] >= EflowConstraint - Debug  # river: env flow does not count transfer scheme water 
            
        md.env_minflow = Constraint(md.ntime, md.neflow, rule=env_minflow)
            
        #%%----------------------------Agricultural module----------------------------
        ##Constraint catalogue##
        
        def agr_waterall2(md,nt,nfz,ncul): #Water allocation limit depending on crop area [time x fzone x crop]
            nc=md.fzone_catch[nfz]
            nm=md.t_month[nt]
            ny=md.t_year[nt]
            nft=md.fzone_type[nfz]
            CulArea         = sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd,kypt] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
            CropWaterDem    = sum(md.yphase_month[ncul,kyps,nm] * md.wET0[nt,nc] * md.aKc[ncul,nm] * md.mm_to_m3perha     #Monthly crop water demand of phase                           
                            * sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd,kypt]*md.aYieldMat[kypt,kyps] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
                            for kyps in md.nyphase) 
            RainWater       = md.wRainFall[nt,nc] * md.mm_to_m3perha * CulArea * sum(md.yphase_month[ncul,kyps,nm] for kyps in md.nyphase)  
            if md.aIrrigation[nft] == 0: #Rainfed
                if sum(md.yphase_month[ncul,kyps,nm]*value(md.wET0[nt,nc])*md.aKc[ncul,nm] for kyps in md.nyphase) == 0: #No water demand of this culture at this month
                    return Constraint.Skip #As otherwise 0 >= 0 is not a permited constraint for the model
                else:
                    return  RainWater >= CropWaterDem # Water requierment has to be less than precipitation
            
            elif md.aIrrigation[nft] == 1: #Irrigated
                Pumping = sum(md.AwGWSUPPLY[nt,nfz,ncul] for aq in md.naquifer if md.fzone_catch[nfz] == md.aqui_catch[aq])                             
                return (md.AwSUPPLY[nt,nfz,ncul] + Pumping)*(1-md.aCulRturn[nft,ncul])  >= (CropWaterDem-RainWater)  #Water supply plus rainfall has to be over water dem 

        def agr_waterall(md,ny,nfz,ncul,nyps):  
            nft=md.fzone_type[nfz]
            nc=md.fzone_catch[nfz]
            #Rainfall accounts net rainfall i.e. rainfall within demand (so rainfall from another month cannot compensate)
            RainFall        = sum(min(value(md.wRainFall[kt,nc]),md.aKc[ncul,md.t_month[kt]]*value(md.wET0[kt,nc])) * md.mm_to_m3perha * md.yphase_month[ncul,nyps,md.t_month[kt]] for kt in md.ntime if md.t_year[kt]==ny)
            MaxDemand       = sum(md.aKc[ncul,md.t_month[kt]]*value(md.wET0[kt,nc])                                 * md.mm_to_m3perha * md.yphase_month[ncul,nyps,md.t_month[kt]] for kt in md.ntime if md.t_year[kt]==ny)
            if MaxDemand == 0:
                return Constraint.Skip
            #ADD: PhaseRatio assumes water supply in a month is distributed proportionally to demand between different growth phases of a same month (=> in a same month a phase can not be more curtailed than an other) 
            PhaseRatio      = {t: md.yphase_month[ncul,nyps,md.t_month[t]]/sum(md.yphase_month[ncul,kyps,md.t_month[t]] for kyps in md.nyphase)
                                  if sum(md.yphase_month[ncul,kyps,md.t_month[t]] for kyps in md.nyphase) != 0 else 0
                               for t in md.ntime if md.t_year[t]==ny}           
            SurfaceWater    = sum(md.AwSUPPLY[kt,nfz,ncul]   * PhaseRatio[kt] for kt in md.ntime if md.t_year[kt]==ny) if md.aIrrigation[nft] == 1 else 0  
            Pumping         = sum(md.AwGWSUPPLY[kt,nfz,ncul] * PhaseRatio[kt] for kt in md.ntime for aq in md.naquifer if md.fzone_catch[nfz] == md.aqui_catch[aq] and md.t_year[kt]==ny)
            CulArea         = sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd,kypt] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
            PhaseDemand     = MaxDemand * sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd,kypt]*md.aYieldMat[kypt,nyps] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
                            
            return (SurfaceWater + Pumping) * (1-md.aCulRturn[nft,ncul]) + RainFall*CulArea >= PhaseDemand
            
        def agr_maxland(md,ny,nfz): #Land limit constraint on cultivated area [year x fzone] 
            if md.Options['Yield water response'] == 'nonlinear':
                CultivatedLand  = sum(md.AlCULAREA[ny,nfz,kfd] for kfd in md.nfieldculture if kfd[1]==1) #field[1]=1 (first culture of field) is the reference for used area
            else:    
                CultivatedLand  = sum(md.AlCULAREA[ny,nfz,kfd,kypt] for kfd in md.nfieldculture for kypt in md.nypath if kfd[1]==1) #field[1]=1 (first culture of field) is the reference for used area
            
            return CultivatedLand <= md.aLandCap[ny,nfz] + InvestedCapacity(md.aLandCap[ny,nfz],nfz,'farming',ny,year=1)            
                
        def agr_fieldarea(md,ny,nfz,nfd,nfdc): #Every culture in a same field has the same area [year x fzone x ftype x (field x culture)] 
            if nfdc==1:
                return Constraint.Skip
            else:
                if md.Options['Yield water response'] == 'nonlinear':
                    return md.AlCULAREA[ny,nfz,nfd,nfdc] == md.AlCULAREA[ny,nfz,nfd,1] #area of field[1] = 2,3...(if exist) have to match the one of field[1] = 1           
                else:
                    return sum(md.AlCULAREA[ny,nfz,nfd,nfdc,kypt] for kypt in md.nypath) == sum(md.AlCULAREA[ny,nfz,nfd,1,kypt] for kypt in md.nypath) #area of field[1] = 2,3...(if exist) have to match the one of field[1] = 1           
                
        def agr_maxcularea(md,ny,nfz,ncul): #Limits the cultivated area by an upper value (usually observed culture area) [year x fzone x culture]
            if md.Options['Yield water response'] == 'nonlinear':
                return sum(md.AlCULAREA[ny,nfz,kfd] for kfd in md.nfieldculture if md.field_culture[kfd]==ncul) <= md.aCulMax[ny,nfz,ncul] 
            else:
                return sum(md.AlCULAREA[ny,nfz,kfd,kypt] for kfd in md.nfieldculture for kypt in md.nypath if md.field_culture[kfd]==ncul) <= md.aCulMax[ny,nfz,ncul] 
        
        def agr_onecropchoice(md,ny,nfz,nfd,nfdc): #One culture area choice for the whole optimization period [year x fzone x ftype x (field x culture)] 
            if md.Options['Yield water response'] == 'nonlinear':
                return md.AlCULAREA[ny,nfz,nfd,nfdc] == md.AlCULAREA[md.t_year[md.Options['tini']],nfz,nfd,nfdc]
            else:
                return sum(md.AlCULAREA[ny,nfz,nfd,nfdc,kypt] for kypt in md.nypath) == sum(md.AlCULAREA[md.t_year[md.Options['tini']],nfz,nfd,nfdc,kypt] for kypt in md.nypath)
        
        def agr_cropprodlinearized(md,ny,nfz,ncr): #Crop production [year x fzone x crop]
            nft=md.fzone_type[nfz]
            #linearized yield response: crop production is linked to yield path choice        
            FzoneCropProd       = sum(md.aCulYield[ny,nft,md.field_culture[kfd]] # t/ha
                                * sum(md.AlCULAREA[ny,nfz,kfd,kypt]  # 1000 ha/year
                                * (1-sum(md.akY[md.field_culture[kfd],kyps]*(1-md.aYieldMat[kypt,kyps]) for kyps in md.nyphase)) 
                                for kypt in md.nypath) for kfd in md.nfieldculture if md.culture_crop[md.field_culture[kfd]]==ncr) 
            return md.AcPROD[ny,nfz,ncr] == FzoneCropProd #1000t/year
        
        def agr_maxirrigmonth(md,nt,nfz,ncul): #Irrigation is limited to max phase water demand to avoid compensation from one phase to another - Only for irrigated farming zones
            ny=md.t_year[nt]
            nm=md.t_month[nt]
            nft=md.fzone_type[nfz]
            nc=md.fzone_catch[nfz]            
            
            if 1==2: #sum(md.yphase_month[ncul,kyps,nm] for kyps in md.nyphase) * max(0,md.aKc[ncul,nm]*value(md.wET0[nt,nc]) - value(md.wRainFall[nt,nc]))  == 0:
                return Constraint.Skip                
            else:                
                MonthDemand     = sum(md.yphase_month[ncul,kyps,nm] * max(0,md.aKc[ncul,nm]*value(md.wET0[nt,nc]) - value(md.wRainFall[nt,nc])) for kyps in md.nyphase)
                SurfaceWater    = md.AwSUPPLY[nt,nfz,ncul]  
                Pumping         = md.AwGWSUPPLY[nt,nfz,ncul] if md.Options['Groundwater']==1 else 0
                if md.Options['Crop choice'] == 'fixed':
                    CulArea     = md.kha_to_Mha*md.CULAREA[ny,nfz,ncul]
                elif md.Options['Yield water response'] == 'nonlinear':
                    CulArea     = sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd] for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
                else:
                    CulArea     = sum(md.kha_to_Mha*md.AlCULAREA[ny,nfz,kfd,kypt] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
                return (SurfaceWater + Pumping) * (1-md.aCulRturn[nft,ncul]) <= MonthDemand * md.mm_to_m3perha * CulArea  
        
        def agr_cropprodexplicit(md,ny,nfz,ncr): #ADD THIS CAN LEAD TO NEGATIVE YIELDS ! 
            nft=md.fzone_type[nfz]
            nc=md.fzone_catch[nfz]
            #Cultivated Area
            if md.Options['Yield water response'] == 'nonlinear' and md.Options['Crop choice'] != 'fixed':
                CulArea         = {cul:sum(md.AlCULAREA[ny,nfz,kfd] for kfd in md.nfieldculture if md.field_culture[kfd]==cul) for cul in md.nculture}
            else:
                CulArea         = {cul:md.CULAREA[ny,nfz,cul] for cul in md.nculture}
            #Rainfall accounts net rainfall i.e. rainfall within demand (so rainfall from another month cannot compensate)
            RainFall            = {(cul,yps): sum(min(value(md.wRainFall[kt,nc]),md.aKc[cul,md.t_month[kt]]*value(md.wET0[kt,nc])) * md.mm_to_m3perha * md.yphase_month[cul,yps,md.t_month[kt]] for kt in md.ntime if md.t_year[kt]==ny) for cul in md.nculture for yps in md.nyphase if md.culture_crop[cul]==ncr}
            Demand              = {(cul,yps): sum(md.aKc[cul,md.t_month[kt]]*value(md.wET0[kt,nc])                                 * md.mm_to_m3perha * md.yphase_month[cul,yps,md.t_month[kt]] for kt in md.ntime if md.t_year[kt]==ny) for cul in md.nculture for yps in md.nyphase if md.culture_crop[cul]==ncr}
            PhaseRatio          = {(cul,yps,t):      md.yphase_month[cul, yps,md.t_month[t]]
                                              /  sum(md.yphase_month[cul,kyps,md.t_month[t]] for kyps in md.nyphase)
                                              if sum(md.yphase_month[cul,kyps,md.t_month[t]] for kyps in md.nyphase) != 0 else 0
                                              for t in md.ntime for cul in md.nculture for yps in md.nyphase if md.t_year[t]==ny and md.culture_crop[cul]==ncr}
            #ADD: assumes water supply in a month is distributed proportionally to demand between different growth phases of a same month (=> in a same month a phase can not be more curtailed than an other) 
            SurfaceWater        = {(cul,yps): sum(md.AwSUPPLY[kt,nfz,cul]   * PhaseRatio[cul,yps,kt] for kt in md.ntime if md.t_year[kt]==ny and md.yphase_month[cul,yps,md.t_month[kt]] != 0) if md.aIrrigation[nft] == 1 else 0     for cul in md.nculture for yps in md.nyphase if md.culture_crop[cul]==ncr}             
            Pumping             = {(cul,yps): sum(md.AwGWSUPPLY[kt,nfz,cul] * PhaseRatio[cul,yps,kt] for kt in md.ntime if md.t_year[kt]==ny and md.yphase_month[cul,yps,md.t_month[kt]] != 0) if md.Options['Groundwater']==1 else 0 for cul in md.nculture for yps in md.nyphase if md.culture_crop[cul]==ncr}            
            Supply              = {(cul,yps): (SurfaceWater[cul,yps] + Pumping[cul,yps]) * (1-md.aCulRturn[nft,cul]) + RainFall[cul,yps]*md.kha_to_Mha*CulArea[cul] for cul in md.nculture for yps in md.nyphase if md.culture_crop[cul]==ncr}           
            FzoneCropProd       = sum(md.aCulYield[ny,nft,kcul] 
                                * (CulArea[kcul] - sum(md.akY[kcul,kyps]*(CulArea[kcul]-Supply[kcul,kyps]/Demand[kcul,kyps]/md.kha_to_Mha) 
                                                       if Demand[kcul,kyps] > RainFall[kcul,kyps] and Demand[kcul,kyps] != 0 else 0 for kyps in md.nyphase))                                  
                                for kcul in md.nculture if md.culture_crop[kcul]==ncr)                       
#            
            return md.AcPROD[ny,nfz,ncr] == FzoneCropProd + (md.DUMYCROP[ny,nfz,ncr] if md.Options['Debug mode'] in [1,2] else 0) #WARNING, #ADD: non linear dummy variable expression - will bug if used in linear formulation

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
            md.agr_onecropchoice   = Constraint(md.nyear, md.nfzone, md.nfieldculture,   rule=agr_onecropchoice)
        
        #%%----------------------------Crop market module----------------------------
        ##Constraint catalogue##
        def agr_cropbalance(md,ny,ncm,ncr): #Crop mass balance [year x cmarket x crop]
            CropLocalProd       = sum(md.AcPROD[ny,kfz,ncr] for kfz in md.nfzone if md.fzone_cmarket[kfz]==ncm)
            CropLocalImport     = sum(md.AcTRANS[ny,kct,ncr]*(1-md.aTransLoss[kct]) for kct in md.nctrans if md.aTransOut[kct]==ncm)
            CropExport          = sum(md.AcTRANS[ny,kct,ncr] for kct in md.nctrans if md.aTransIn[kct]==ncm) 
            CropSupply          = sum(md.AcSUPPLY[ny,ncm,ncr,kcds] for kcds in md.ncdstep)
            CropExtProd         = 0
            if ncm in md.nextcmarket: 
                CropExtProd     = md.AcEXTPROD[ny,ncm,ncr] #crop produced in international/external market (t)
            return CropSupply <= CropExtProd + CropLocalProd + CropLocalImport - CropExport # '<=' instead of '==' to allow crop waste in case of over production
               
        def agr_cropdemand(md,ny,ncm,ncr,ncds): #Crop market level crop step demand [year x cmarket x crop x cdstep]
            CropStepDemand      = md.aCropDem[ny,ncm,ncr] * md.aStepDem[ncm,ncr,ncds] 
            return  md.AcSUPPLY[ny,ncm,ncr,ncds] <= CropStepDemand 
        
        def agr_cropmindemand(md,ny,ncm,ncr): #Crop market level crop minimum demand ensuring food security [year x cmarket x crop]
            CropSupply  = sum(md.AcSUPPLY[ny,ncm,ncr,kcds] for kcds in md.ncdstep)
            return  CropSupply >= md.aMinDem[ncm,ncr] 
                
        ##Create constraints##
        md.agr_cropbalance      = Constraint(md.nyear, md.ncmarket, md.ncrop, rule=agr_cropbalance)
        md.agr_cropdemand       = Constraint(md.nyear, md.ncmarket, md.ncrop, md.ncdstep, rule=agr_cropdemand)
        if md.Options['Minimum supply'] == 1:
            md.agr_cropmindemand= Constraint(md.nyear, md.ncmarket, md.ncrop, rule=agr_cropmindemand)
        
        #%%----------------------------Energy production module----------------------------
        #Constraint catalogue
        def engy_firmhp(md,nt,nhp):
            return md.EeFIRMHP[md.t_year[nt],nhp]<=sum(md.EeHPPROD[nt,kpld,nhp] for kpld in md.npload)
        
        def engy_turbinecapacity(md,nt,npld,nhp):            
            if md.hp_res[nhp] == 'ROR':
                return md.EwHPDISCHARGE[nt,npld,nhp] <= md.wMaxTurb[nhp] * md.m3pers_to_Mm3perMonth * md.eLoadTime[npld]
            elif md.Options['Hydropower'] == 'nonlinear':
                nres = md.hp_res[nhp]
                AvStorage = (md.sw*md.WwRSTORAGE[nt,nres]+md.sw*md.WwRSTORAGE[md.t_prev[nt],nres])/2 if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else (md.sw*md.WwRSTORAGE[nt,nres]+md.wStorIni[nres])/2 
                MaxTurbineFlow = TurbFlowVol(AvStorage,md.wResMOL[nres],md.wResFSL[nres],md.wMinTurb[nhp],md.wMaxTurb[nhp])
                return md.EwHPDISCHARGE[nt,npld,nhp] <= MaxTurbineFlow * md.m3pers_to_Mm3perMonth * md.eLoadTime[npld]
        
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
            if md.Options['Hydropower'] == 'nonlinear' and md.hp_res[nhp] != 'ROR':
                nres = md.hp_res[nhp]
                AvStorage = (md.sw*md.WwRSTORAGE[nt,nres]+md.sw*md.WwRSTORAGE[md.t_prev[nt],nres])/2 if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else (md.sw*md.WwRSTORAGE[nt,nres]+md.wStorIni[nres])/2 
                RelHead = RelHeadVol(AvStorage,md.wResMOL[nres],md.wResFSL[nres],md.wMinHead[nres],md.wMaxHead[nres])                
                return md.se*md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * md.eHppProd[nhp] * md.EwHPDISCHARGE[nt,npld,nhp] * RelHead
            else:
                return md.se*md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * md.eHppProd[nhp] * md.EwHPDISCHARGE[nt,npld,nhp]
                
        def engy_hpcapacity(md,nt,npld,nhp): #Hydropower production per load segment limited by hydropower capacity [time x pload x hpp] 
            return md.se*md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * (md.eHppCap[nhp] + InvestedCapacity(md.eHppCap[nhp],nhp,'hydro',nt)) * md.MW_to_GWhperMonth * md.eLoadTime[npld]                 
        
        def engy_oppcapacity(md,nt,npld,nop): #Powerplant production per load segment limited by capacity [time x pload x opp]
            LoadCapacityFactor = md.eLoadCap[npld,md.op_ptech[nop]] if md.Options['Load capacity'] == 1 else 1                 
            return md.se*md.EeOPPROD[nt,npld,nop] <= (md.eOppCap[nop] + InvestedCapacity(md.eOppCap[nop],nop,'power',nt)) * md.MW_to_GWhperMonth * md.eLoadTime[npld] * LoadCapacityFactor
        
        def engy_gencapacity(md,nt,npld,npt,npm): #Capacity constraint of Generic capacity [time x pload x ptech x pmarket]
            LoadCapacityFactor = 1 #No particular restriction if not specified
            if md.Options['Load capacity'] == 1:
                LoadCapacityFactor = md.eLoadCap[npld,npt]
            return md.EeGENPROD[nt,npld,npt,npm] <= sum(md.EeGENCAP[ky,npt,npm] for ky in md.nyear if ky <= md.t_year[nt] and ky >= md.t_year[nt]-md.eLifeTime[npt,npm]) * md.MW_to_GWhperMonth * md.eLoadTime[npld] * LoadCapacityFactor 

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
                RampingCap      = sum(md.se*md.EeGENCAP[ky,npt,npm] for ky in md.nyear if ky <= md.t_year[nt] and ky >= md.t_year[nt]-md.eLifeTime[npt,npm]) * md.eTechRamp[npt,npm] #[kWh/day]
                return  ProdRateSegA <= ProdRateSegB + RampingCap   #Power plant's Load rate from load segment to other is not allow to vary more than ramp rate
        
        def engy_oppramping(md,nt,nplda,npldb,nop): #Ramping constraints for other power plants
            if nplda == npldb:
                return Constraint.Skip
            else:
                ProdRateSegA    = md.se*md.EeOPPROD[nt,nplda,nop] / md.eLoadTime[nplda] #[kWh/day]
                ProdRateSegB    = md.se*md.EeOPPROD[nt,npldb,nop] / md.eLoadTime[npldb] #[kWh/day]
                RampingCap      = md.eOppCap[nop] * md.eOppRamp[nop] #[kWh/day]
                return  ProdRateSegA <= ProdRateSegB + RampingCap    #Power plant's Load rate from load segment to other is not allow to vary more than ramp rate
        
        #Create Constraints
        md.engy_firmhp       = Constraint(md.ntime, md.nhpp, rule=engy_firmhp)
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
            OppProd=0
            GenProd=0
            HpProd      = sum(md.EeHPPROD[nt,npld,khp]                          for khp in md.nhpp if md.hp_pmarket[khp]==npm)      #Hydropower energy prod
            OppProd     = sum(md.EeOPPROD[nt,npld,kop]                          for kop in md.nopp if md.op_pmarket[kop]==npm)      #Power plant energy prod
            GenProd     = sum(md.EeGENPROD[nt,npld,kpt,npm]                     for kpt in md.nptech)                               #Generic Power plant energy prod
            NetImport   = sum(md.EeTRANS[nt,npld,ktl] * (1-md.eTransLoss[ktl])  for ktl in md.ntransline if md.eTransOut[ktl]==npm) 
            GrossExport = sum(md.EeTRANS[nt,npld,ktl]                           for ktl in md.ntransline if md.eTransIn[ktl]==npm)
            GrossSupply = md.EeSUPPLY[nt,npld,npm] / (1-md.eSupLoss[npm])                                           
            return GrossSupply + GrossExport == HpProd + OppProd + GenProd + NetImport
                
        def engy_transmissioncap(md,nt,npld,ntl): #Transmission limited to line capacity [time x pload x tmsline]   
            return md.se*md.EeTRANS[nt,npld,ntl] <= (md.eTransCap[ntl] + InvestedCapacity(md.eTransCap[ntl],ntl,'transmission',nt)) * md.MW_to_GWhperMonth * md.eLoadTime[npld]
                
        def engy_demand(md,nt,npld,npm): #Power demand satisfaction
            #ny0         = md.t_year[md.Options['tini']]
            Supply      = md.se*md.EeSUPPLY[nt,npld,npm] 
            LoadSegment = md.eEngyDem[md.t_year[nt],md.t_month[nt],npm] * md.eLoadDem[npld] # * md.eDemYear[npm]**(md.t_year[nt]-ny0) Demand of power market for power load segment          
            return Supply <= LoadSegment
                
        #Create constraints
        md.engy_balance         = Constraint(md.ntime, md.npload, md.npmarket, rule=engy_balance)
        md.engy_transmissioncap = Constraint(md.ntime, md.npload, md.ntransline, rule=engy_transmissioncap)
        md.engy_demand          = Constraint(md.ntime, md.npload, md.npmarket, rule=engy_demand)
        
#%%SAVE MODEL
        self.model=md


