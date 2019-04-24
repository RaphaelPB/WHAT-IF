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

from __future__ import division #As in python 2.7 division between integers gives an integer - avoid problems
from pyomo.environ import ConcreteModel, Set, Param, Var, Objective, Constraint
from pyomo.environ import NonNegativeReals, exp, value
#import math

class HydroeconomicOptimization():
    def __init__(self,parameters,scenario='WHATIF_main',MPC=0):
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

        def ReadParam(ParamName,option=1,index=0):
            if ParamName not in parameters.val.keys() or option == 0: #parameter is not activated
                if index==1:
                    return [] #returns empty set
                return {} #returns empty parameter
            if ParamName in GrowingParameters:
                return GrowParam(ParamName)
            else:
                if ParamName in parameters.scen.keys():                
                    return parameters.val[ParamName][parameters.val[parameters.scen[ParamName]][scenario]]
                else:
                    return parameters.val[ParamName]
        
        def LinearGrowAtoB(ValuesA,ValuesB,nyAB,Keys,yyears):
            if type(Keys[0]) is tuple: #keys need to be formulated differently if they are uni or multidimensional
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
        md.nconfig          = Set(initialize=ReadParam('Options').keys())             #Configurations of the model
        md.Options          = Param(md.nconfig, initialize=ReadParam('Options'))                    
    #Unit conversion
        md.mm_to_m3perha       = 1/1000*10000 # 0,001m * 10000m2/ha
        md.MW_to_GWhperMonth   = 1*24*30.44/1000 # 1 MW * 24hours/day * 30.44days/month / 1000MW/GW 
        md.m3pers_to_Mm3perMonth = 1*3600*24*30.44/10**6 #1 m3/s * 3600 s/hour * 24 hours/day * 30.44 days/month / 10^6m3/Mm3
#%%**Water Module**#       
    #Indices
        md.ncountry         = Set(initialize=ReadParam('ncountry'))                 #Countries [indice country] (id)
        md.ncatch           = Set(initialize=ReadParam('catch_ds').keys())          #Catchments [indice catchment] (id)
        md.ntime            = Set(initialize=ReadParam('t_month').keys())           #Time steps [indice time] (id)
        md.nyear            = Set(initialize=ReadParam('nyear'))                    #Years [indice year] (id)
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
        md.t_month          = Param(md.ntime, initialize=ReadParam('t_month'))                     #Month of time step [time] (id)
        md.t_year           = Param(md.ntime, initialize=ReadParam('t_year'))                      #Year of time step [time] (id)
        md.t_prev           = Param(md.ntime, initialize=ReadParam('t_prev'))                      #Previous time step [time] (id)        
    #Hydrology
        md.wRunOff          = Param(md.ntime, md.ncatch, mutable=True, initialize=ReadParam('wRunOff'))   #Runoff [time x catchment] (Mm3/month)
        md.wRainFall        = Param(md.ntime, md.ncatch, mutable=True, initialize=ReadParam('wRainFall')) #Precipitation [time x catchment] (mm)
        md.wET0             = Param(md.ntime, md.ncatch, mutable=True, initialize=ReadParam('wET0'))      #ET0 [time x catchment] (mm)
        md.wFlowLoss        = Param(md.ncatch, initialize=ReadParam('wFlowLoss'))                         #Flow loss [catchment] (%)
    #Groundwater
        opt= md.Options['Groundwater']
        md.wGwRech          = Param(md.ntime, md.naquifer, initialize=ReadParam('wGwRech',option=opt))     #Groundwater recharge [time x catchment] (Mm3/month)
        md.wGwFlow          = Param(md.naquifer, initialize=ReadParam('wGwFlow',option=opt))               #Groundwater baseflow coeeficient [catchment] (-)
        md.wGwIni           = Param(md.naquifer, initialize=ReadParam('wGwIni',option=opt))                #Groundwater initial storage (if not cyclic) [catchment] (Mm3)
        md.wGwCost          = Param(md.naquifer, initialize=ReadParam('wGwCost',option=opt))               #Groundwater pumping cost [catchment] ($/m3)
    #Water Users
        md.wUserDem         = Param(md.nyear, md.nmonth, md.nuser, initialize=ReadParam('wUserDem'))#User water demand [year x month x user] (Mm3/month)        
        md.wUserRturn       = Param(md.nuser, initialize=ReadParam('wUserRturn'))                   #Return rate of user water allocation [user] (%)
        md.wUserLoss        = Param(md.nuser, initialize=ReadParam('wUserLoss'))                    #Loss rate of user water allocation [user] (%)
        md.wUserVal         = Param(md.nuser, initialize=ReadParam('wUserVal'))                     #User marginal value for water [user] ($/m3)
        md.wSupCost         = Param(md.nuser, initialize=ReadParam('wSupCost'))                     #Water supply costs [user] ($/m3) 
    #Reservoirs
        md.wStorCap         = Param(md.nres, initialize=ReadParam('wStorCap'))                      #Storage Capacity [reservoir] (Mm3)
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
        md.aCulYield        = Param(md.nyear, md.nftype, md.nculture, initialize=ReadParam('aCulYield',option=opt))  #Maximum yield of culture [ftype x culture] (t/ha) 
        md.aCulCost         = Param(md.nftype, md.nculture, initialize=ReadParam('aCulCost',option=opt))      #Cultivation cost [ftype x culture] ($/ha)
        md.aIrrigation      = Param(md.nftype, initialize=ReadParam('aIrrigation',option=opt))                #Irrigation (1) or rainfed culture(0) [ftype] (binary) 
    #Farming Zones
        md.aLandCap         = Param(md.nyear, md.nfzone, initialize=ReadParam('aLandCap',option=opt))   #Land capacity for agriculture [year x fzone] (Mha)
        md.aIrrgCost        = Param(md.nfzone, initialize=ReadParam('aIrrgCost',option=opt))            #Irrigation cost [fzone] ($/m3)
        md.aIrrgLoss        = Param(md.nfzone, initialize=ReadParam('aIrrgLoss',option=opt))            #Water loss rate from allocated [fzone] (%) 
    #Cultures
        md.akY              = Param(md.nculture, md.nyphase, initialize=ReadParam('akY',option=opt))    #crop yield response phase-specific factor [culture x yphase](%) 
        md.aKc              = Param(md.nculture, md.nyphase, initialize=ReadParam('aKc',option=opt))    #crop water demand factor [culture x yphase](%) 
        md.aYieldMat        = Param(md.nypath, md.nyphase, initialize=ReadParam('aYieldMat',option=opt))#yield response path matrix [ypath x yphase]
        if md.Options['Culture max area'] in (1,'force'):
            md.aCulMax      = Param(md.nyear, md.nfzone, md.nculture, initialize=ReadParam('aCulMax',option=opt)) #Maximum area per culture and farming zone [fzone x culture] (Mha)
        
#%%**Crop Market module**#
    #Crop markets and demands
        opt=md.Options['Crop market']
        md.ncmarket         = Set(initialize=ReadParam('cmarket_country',option=opt).keys())                  #Crop markets [indice cmarket] (id)
        md.cmarket_country  = Param(md.ncmarket, initialize=ReadParam('cmarket_country',option=opt))          #Country of Crop Market [cmarket] (id)               
        md.nextcmarket      = Set(initialize=[cm for cm in md.ncmarket if ReadParam('cmarket_type')[cm]=='ExternalMarket']) #External crop markets with a crop supply curve (id)
        md.aCropDem         = Param(md.nyear, md.ncmarket, md.ncrop, initialize=ReadParam('aCropDem',option=opt)) #Crop demand [year x cmarket x crop] (Mt/year)
        md.aCropVal         = Param(md.ncmarket, md.ncrop, initialize=ReadParam('aCropVal',option=opt))       #Value of crop [cmarket x crop] ($/t)
    #Food security
        if md.Options['Minimum supply'] == 1:
            md.aMinDem      = Param(md.ncmarket, md.ncrop, initialize=ReadParam('aMinDem',option=opt))    #Crop minimum demand = food security [year x cmarket x crop] (Mt/year)
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
        md.eEngyDem         = Param(md.nyear, md.nmonth, md.npmarket, initialize=ReadParam('eEngyDem',option=opt))  #Local power demand [month x pmarket] (GWh/year)
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
        md.eFuelCost        = Param(md.nyear, md.npmarket, md.nfuel, initialize=ReadParam('eFuelCost',option=opt)) #Cost of fuel [fuel] ($/kWh_fuel)
        md.eFuelCO2         = Param(md.nfuel, initialize=ReadParam('eFuelCO2',option=opt))              #CO2 emission factor of fuel [fuel] (t/kWh_fuel)
        md.eCO2Val          = Param(md.nyear, md.npmarket, initialize=ReadParam('eCO2Val',option=opt) if md.Options['Fuels'] == 1 else {(y,cm):0 for y in md.nyear for cm in md.ncmarket})  #CO2 price [1] ($/t)
    #Generic capacity investment
        opt= md.Options['Energy production'] and md.Options['Energy market'] and md.Options['Power technologies']             
        md.nptech           = Set(initialize=set([key[0] for key in ReadParam('eVarOPEX',option=opt).keys()]))  #Power technologies [indice ptech] (id)
        md.eCAPEX           = Param(md.nyear, md.nptech, md.npmarket, initialize=ReadParam('eCAPEX',option=opt))#Cost of Generic power capacity [ptech x pmarket] (id) ($/MW)
        md.eVarOPEX         = Param(md.nptech, md.npmarket, initialize=ReadParam('eVarOPEX',option=opt))        #Variable Cost of Generic power production [ptech x pmarket] ($/kWh) 
        md.eFixOPEX         = Param(md.nptech, md.npmarket, initialize=ReadParam('eFixOPEX',option=opt))        #Fix Cost of Generic power production [ptech x pmarket] ($/[kWh/day] /year)
        md.eLifeTime        = Param(md.nptech, md.npmarket, initialize=ReadParam('eLifeTime',option=opt))       #Life time of power technology [ptech x pmarket] (years)
        md.eMaxCap          = Param(md.nptech, md.npmarket, initialize=ReadParam('eMaxCap',option=opt))         #Maximum expandable capacity of power technology [ptech x pmarket] (kWh/day)
    #Fuels
        if md.Options['Fuels'] == 1: 
            md.eTechEff     = Param(md.nptech, md.npmarket, initialize=ReadParam('eTechEff',option=opt))     #Efficiency of power technology [ptech x pmarket] (%)
            md.ptech_fuel   = Param(md.nptech, md.npmarket, initialize=ReadParam('ptech_fuel',option=opt))   #Fuel of Technology [ptech x pmarket] (id)
    #Hydropower valuation if power market is off
        if md.Options['Energy market'] == 0:
            md.hp_country   = Param(md.nhpp, initialize=ReadParam('hp_country',option=md.Options['Energy production']))  #Country of HP if power market is off [hpp] (id)
            md.eHppVal      = Param(md.nhpp, initialize=ReadParam('eHppVal',option=md.Options['Energy production']))     #Value of HP production if power market is off [hpp] ($/kWh)
    #Power load segments
        opt= md.Options['Energy production'] and md.Options['Energy market']
        md.npload           = Set(initialize=ReadParam('eLoadDem',option=opt).keys() if md.Options['Load'] == 1 or opt == 0 else [1])            #Slices of power demand per time: day week, day week end, night [indice pload] (id)
        md.eLoadDem         = Param(md.npload, initialize=ReadParam('eLoadDem',option=opt) if md.Options['Load'] == 1 or opt == 0 else {1:1})    #Percentage of demand per power load [pload] (%)
        md.eLoadTime        = Param(md.npload, initialize=ReadParam('eLoadTime',option=opt) if md.Options['Load'] == 1 or opt == 0 else {1:30})  #Length of power load [pload] (days)
        if md.Options['Load capacity'] == 1: #Assumes either "Other power plants" or "Power technologies" are ON
            if md.Options['Power technologies'] == 0:
                md.op_ptech = Param(md.nopp, initialize=ReadParam('op_ptech',option=opt))                      #Technology of Other power plant [pp] (id)
                md.nptech   = Set(initialize=set([ReadParam('op_ptech',option=opt)[op] for op in md.nopp]))    #Technologies are defined in Other power plants if Power technologies module is off [indice ptech] (id)
            md.eLoadCap     = Param(md.npload, md.nptech, initialize=ReadParam('eLoadCap',option=opt))         #Load segment capacity factor of power technology [pload x ptech] (%)

        
#%%DECLARE DECISION VARIABLES###       
    #Water Module#
        md.WwSUPPLY         = Var(md.ntime, md.nuser, within=NonNegativeReals)      #Water user allocation  [time x user] (Mm3)
        md.WwOUTFLOW        = Var(md.ntime, md.ncatch, within=NonNegativeReals)     #Catchment outflow [time x catchment] (Mm3)
        md.WwRSTORAGE       = Var(md.ntime, md.nres, within=NonNegativeReals)       #Reservoir total storage [time x reservoir] (Mm3)
        md.WwGWSTORAGE      = Var(md.ntime, md.naquifer, within=NonNegativeReals)   #Groundwater reservoir total storage [time x aquifer] (Mm3)
        md.WwTRANSFER       = Var(md.ntime, md.ntransfer, within=NonNegativeReals)  #Transfer scheme water transfer [time x transfer] (Mm3)                 
    
    #Agriculture Module#
        md.AlCULAREA        = Var(md.nyear, md.nfzone, md.nfieldculture, md.nypath, within=NonNegativeReals)    #Cultivated area and water demand satisfaction [year x fzone x field x ypath] (-)
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
        
    #Energy market module#
        md.EeTRANS          = Var(md.ntime, md.npload, md.ntransline, within=NonNegativeReals)          #Power transmission between power markets [time x pload x tmsline] (GWh)
        md.EeSUPPLY         = Var(md.ntime, md.npload, md.npmarket, within=NonNegativeReals)            #Power supply to power market demand [time x pload x pmarket] (GWh)
        md.EeGENCAP         = Var(md.nyear, md.nptech, md.npmarket, within=NonNegativeReals)            #Generic Power technologies capacity [year x ptech x pmarket] (MW)
        md.EeGENPROD        = Var(md.ntime, md.npload, md.nptech, md.npmarket, within=NonNegativeReals) #Generic Power Tech Production [time x pload x ptech x pmarket] (GWh/load seg)            

#%%DEBUG VARIABLES - IF ACTIVATED
        DEBUGCOST           = 1000 #arbitrary high cost for additional water $/m3 or deviance from final storage objective
        if md.Options['Initial time step'] == 1:
            md.DUMYSTOR     = Var(md.nres, within=NonNegativeReals) #Dummy variable to avoid infeasible solutions regarding final storage objective
        if md.Options['Debug mode'] == 1:
            md.DEBUG        = Var(md.ntime, md.ncatch, within=NonNegativeReals) #DEBUG Variable adding water in the water balance    
            md.DUMYCROP     = Var(md.nyear, md.nfzone, md.ncrop, within=NonNegativeReals)   #Dummy Crop Production to compensate for negative prod [year x fzone x crop] (Mt)       
            md.DUMYEFLOW    = Var(md.ntime, md.neflow, within=NonNegativeReals) #Tolerance to eflow objective if physically impossible

#%%OBJECTIVE FUNCTION - ECONOMIC MODULE
        def obj_rule(md):
        #Water module            
            UserBenefit         = {y:sum(md.WwSUPPLY[t,u] * md.wUserVal[u] for u in md.nuser for t in md.ntime if md.t_year[t]==y) for y in md.nyear} #benefit from (non agricultural) water users
            UserSupCosts        = {y:sum(md.WwSUPPLY[t,u]/(1-md.wUserLoss[u]) * md.wSupCost[u] for u in md.nuser for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
            UserBalance         = {y:UserSupCosts[y]-UserBenefit[y] for y in md.nyear}
            
        #Agricultural module: Cultivation costs of cultures (labour, machinery, fertilizers, infrastructure....), Irrigation costs of cultures (pumping ...)     
            AgrCulCost          = {y:sum(md.AlCULAREA[y,fz,fd,ypt] * md.aCulCost[md.fzone_type[fz],md.field_culture[fd]] for fz in md.nfzone for fd in md.nfieldculture for ypt in md.nypath) for y in md.nyear}
            AgrIrrgCost         = {y:sum(md.AwSUPPLY[t,fz,cul]/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.t_year[t]==y) for y in md.nyear}
            GwPumpingCost       = {y:sum(md.AwGWSUPPLY[t,fz,cul] * md.wGwCost[aq] for t in md.ntime for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.t_year[t]==y and md.fzone_catch[fz]==md.aqui_catch[aq]) for y in md.nyear}
            AgricultureBalance  = {y:AgrCulCost[y] + AgrIrrgCost[y] + GwPumpingCost[y] for y in md.nyear}
        
        #Crop market module: Crop allocated to crop markets, Cost and benefits of imported crops (external market), Transport costs          
            AgrBenefit          = {y:sum(md.AcSUPPLY[y,cm,cr,cds] * md.aCropVal[cm,cr] * md.aStepVal[cm,cr,cds] for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep) for y in md.nyear}    
            if md.Options['Crop market'] == 0:
                AgrBenefit      = {y:sum(md.AcPROD[y,fz,cr] * md.aFarmVal[md.fzone_country[fz],cr] for fz in md.nfzone for cr in md.ncrop) for y in md.nyear}
            AgrTransCost        = {y:sum(md.AcTRANS[y,ct,cr] * md.aTransCost[ct,cr]  for ct in md.nctrans for cr in md.ncrop) for y in md.nyear}
            ExtProdCost         = {y:sum(md.AcEXTPROD[y,cm,cr] * md.aCropVal[cm,cr] * 1.0001 for cm in md.nextcmarket for cr in md.ncrop) for y in md.nyear} #The assumption is that external market produce crops at their market price + 1% #ADD: crop supply curve
            CropMarketBalance   = {y:AgrTransCost[y] + ExtProdCost[y] - AgrBenefit[y] for y in md.nyear}
      
        #Energy production module: Production costs of hydropower and power plants (O&M, fuel ...)    
            HpOMCost            = {y:sum(md.eHppCost[hp]*md.EeHPPROD[t,pld,hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}            
            OppOMCost           = {y:sum(md.eOppCost[opp]*md.EeOPPROD[t,pld,opp] for t in md.ntime for pld in md.npload for opp in md.nopp if md.t_year[t]==y) for y in md.nyear}    
            OppFuelCost         = {y:sum(md.EeOPPROD[t,pld,opp]/md.eOppEff[opp]*md.eFuelCost[y,md.op_pmarket[opp],fu] for t in md.ntime for pld in md.npload for opp in md.nopp for fu in md.nfuel if md.t_year[t]==y and md.op_fuel[opp]==fu) for y in md.nyear}
            OppCO2Cost          = {y:sum(md.EeOPPROD[t,pld,opp]/md.eOppEff[opp]*md.eFuelCO2[fu]*md.eCO2Val[y,md.op_pmarket[opp]] for t in md.ntime for pld in md.npload for opp in md.nopp for fu in md.nfuel if md.t_year[t]==y and md.op_fuel[opp]==fu) for y in md.nyear}
            #ADD overestimates CAPEX if remaining time is higher than lifetime + obliges the model to reinvest in the same infrastructure at the end of the lifetime
            GenCapCost          = {y:sum(md.EeGENCAP[y,pt,pm]*md.eCAPEX[y,pt,pm]*(md.t_year[md.Options['tfin']]-y+1)/md.eLifeTime[pt,pm]+sum(md.EeGENCAP[ky,pt,pm] for ky in md.nyear if ky <= y)*md.eFixOPEX[pt,pm] for pt in md.nptech for pm in md.npmarket) for y in md.nyear}
            GenProdCost         = {y:sum(md.EeGENPROD[t,pld,pt,pm]*md.eVarOPEX[pt,pm] for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear}
            GenFuelCost         = {y:sum(md.EeGENPROD[t,pld,pt,pm]/md.eTechEff[pt,pm]*md.eFuelCost[y,pm,fu]             for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.t_year[t]==y and md.ptech_fuel[pt,pm]==fu) for y in md.nyear}
            GenCO2Cost          = {y:sum(md.EeGENPROD[t,pld,pt,pm]/md.eTechEff[pt,pm]*md.eFuelCO2[fu]*md.eCO2Val[y,pm]  for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.t_year[t]==y and md.ptech_fuel[pt,pm]==fu) for y in md.nyear}
            EnergyProdBalance   = {y:HpOMCost[y] + OppOMCost[y] + OppFuelCost[y] + OppCO2Cost[y] + GenCapCost[y] + GenProdCost[y] + GenFuelCost[y] + GenCO2Cost[y] for y in md.nyear}           
                
        #Energy market module: Energy distributed, Energy transmission costs, Energy curtailment costs    
            EngyBenefit         = {y:sum(md.EeSUPPLY[t,pld,pm] * md.eEngyVal[md.t_month[t],pm] for t in md.ntime for pld in md.npload for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear}
            if md.Options['Energy market'] == 0:
                EngyBenefit     = {y:-sum(md.EeHPPROD[t,pld,hp]*md.eHppVal[hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.t_year[t]==y) for y in md.nyear}
            EngyTransCost       = {y:sum(md.EeTRANS[t,pld,tl] * md.eTransCost[tl]  for t in md.ntime for pld in md.npload for tl in md.ntransline if md.t_year[t]==y) for y in md.nyear}
            EnergyMarketBalance = {y:EngyTransCost[y]-EngyBenefit[y] for y in md.nyear} #            
        
        #DEBUG MODE - Warning breaks shadowprices
            DebugCost=0    
            if md.Options['Debug mode'] == 1:
                DebugCost += sum(md.DEBUG[t,c]*DEBUGCOST for t in md.ntime for c in md.ncatch)
                DebugCost += sum(md.DUMYCROP[y,fz,cr]*DEBUGCOST*100 for y in md.nyear for fz in md.nfzone for cr in md.ncrop)
                DebugCost += sum(md.DUMYEFLOW[t,ef]*DEBUGCOST/2 for t in md.ntime for ef in md.neflow)
            if md.Options['Debug mode'] == 1 and md.Options['Initial time step'] == 1:
                DebugCost += sum(md.DUMYSTOR[res]*DEBUGCOST/2 for res in md.nres)
        
            return DebugCost + sum(UserBalance[y] + AgricultureBalance[y] + CropMarketBalance[y] + EnergyProdBalance[y] + EnergyMarketBalance[y] for y in md.nyear)  
        
        md.obj = Objective(rule=obj_rule)
        
#%%CONSTRAINTS

#%%----------------------------Water module----------------------------
        ##Constraint catalogue##
        def water_waterbalance(md,nt,nc): #Water Balance per Watershed [catchment x time]
            NetInflow           = (1-md.wFlowLoss[nc]) * sum(md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) #Net incoming flow from upstream catchments (m3)
            RunOff              = md.wRunOff[nt,nc]   #Catchment Runoff (m3) 
            StorPrev            = sum(md.WwGWSTORAGE[md.t_prev[nt],kaq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[kaq] for kaq in md.naquifer if md.aqui_catch[kaq]==nc) #Previous ground water storage
            Pumping             = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in md.nifzone for kcul in md.nculture for aq in md.naquifer if md.fzone_catch[kfz] == md.aqui_catch[aq] == nc) #Water abstraction from groundwater reservoir
            BaseFlow            = sum(StorPrev*(1-exp(-md.wGwFlow[kaq])) + (md.wGwRech[nt,kaq]-Pumping)*(1-(1-exp(-md.wGwFlow[kaq]))/md.wGwFlow[kaq]) for kaq in md.naquifer if md.aqui_catch[kaq]==nc) #Baseflow from groundwater reservoir
            UserCons            = sum(md.WwSUPPLY[nt,ku] * (1/(1-md.wUserLoss[ku])-md.wUserRturn[ku]) for ku in md.nuser if md.user_catch[ku] == nc)      #Consumption of non agricultural users (m3)
            AgrCons             = sum(md.AwSUPPLY[nt,kfz,kcu] * (1/(1-md.aIrrgLoss[kfz])-md.aCulRturn[md.fzone_type[kfz],kcu]) for kfz in md.nifzone for kcu in md.nculture if md.fzone_catch[kfz]==nc)  #Agricultural module: Consumptiom to agriculture/crops (m3)    
            if md.Options['Initial time step'] == 1 and md.Options['tini'] == nt: #Initial storage for initial time step
                DeltaStorage    = sum(md.WwRSTORAGE[nt,kres] - md.wStorIni[kres]                 for kres in md.nres    if md.res_catch[kres]==nc)       #Storage in reservoir (m3)
                Evaporation     = max(0,value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc])) * md.mm_to_m3perha * sum(md.wkV[kres]*(md.WwRSTORAGE[nt,kres]+md.wStorIni[kres])/2 + md.wResArea[kres] for kres in md.nres if md.res_catch[kres]==nc)
            else: #Cyclic run or not first time step 
                DeltaStorage    = sum(md.WwRSTORAGE[nt,kres] - md.WwRSTORAGE[md.t_prev[nt],kres] for kres in md.nres    if md.res_catch[kres]==nc)       #Storage in reservoir (m3)
                Evaporation     = max(0,value(md.wET0[nt,nc])-value(md.wRainFall[nt,nc])) * md.mm_to_m3perha * sum(md.wkV[kres]*(md.WwRSTORAGE[nt,kres]+md.WwRSTORAGE[md.t_prev[nt],kres])/2 + md.wResArea[kres] for kres in md.nres if md.res_catch[kres]==nc)                
            TransNetInflow      = sum(md.WwTRANSFER[nt,ktrans] * (1-md.wTransLoss[ktrans])   for ktrans in md.ntransfer if md.transfer_ds[ktrans] == nc)  #Incoming flow from upstream transfer scheme- if exists (m3)
            TransBrutOutflow    = sum(md.WwTRANSFER[nt,ktrans]                               for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc)  #Allocation to downstream transfer scheme- if exists (m3)                   
            CatchmentOutflow    = md.WwOUTFLOW[nt,nc]  #Catchment outflow : Allocation to downstream catchment (m3)           
            DEBUG = md.DEBUG[nt,nc] if md.Options['Debug mode'] == 1 else 0
            
            return  UserCons + AgrCons + TransBrutOutflow + DeltaStorage + Evaporation + CatchmentOutflow == NetInflow + TransNetInflow + RunOff + BaseFlow + DEBUG
                
        def water_waterbalance2(md,nt,nc): #Avoids allocation from downstream reservoir to upstream demand [catchment x time] 
            if  sum(1 for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc) + sum(1 for ku in md.nuser if md.user_catch[ku] == nc) + sum(1 for kfz in md.nfzone if md.fzone_catch[kfz]==nc) == 0:
                return Constraint.Skip #if there is no demand and no transfer project constraint would be empty and is therefore skipped
            BaseFlow            = 0
            NetInflow           = (1-md.wFlowLoss[nc]) * sum(md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) #Net incoming flow from upstream catchments (m3)
            RunOff              = md.wRunOff[nt,nc] #Catchment Inflow (m3) 
            GwStorPrev          = sum(md.WwGWSTORAGE[md.t_prev[nt],kaq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[kaq] for kaq in md.naquifer if md.aqui_catch[kaq]==nc) #Previous ground water storage
            Pumping             = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in md.nifzone for kcul in md.nculture for aq in md.naquifer if md.fzone_catch[kfz] == md.aqui_catch[aq] == nc) #Water abstraction from groundwater reservoir
            BaseFlow            = sum(GwStorPrev*(1-exp(-md.wGwFlow[kaq])) + (md.wGwRech[nt,kaq]-Pumping)*(1-(1-exp(-md.wGwFlow[kaq]))/md.wGwFlow[kaq]) for kaq in md.naquifer if md.aqui_catch[kaq]==nc) #Baseflow from groundwater reservoir
            UserBrutAll         = sum(md.WwSUPPLY[nt,ku] * 1/(1-md.wUserLoss[ku]) for ku in md.nuser if md.user_catch[ku] == nc) #Allocation to non agricultural users (m3)
            AgrBrutAll          = sum(md.AwSUPPLY[nt,kfz,kcu]* 1/(1-md.aIrrgLoss[kfz]) for kfz in md.nifzone for kcu in md.nculture if md.fzone_catch[kfz]==nc) #Agricultural module: Allocation to agriculture/crops (m3)    
            TransNetInflow      = sum(md.WwTRANSFER[nt,ktrans] * (1-md.wTransLoss[ktrans])   for ktrans in md.ntransfer if md.transfer_ds[ktrans] == nc)  #Incoming flow from upstream transfer scheme- if exists (m3)
            TransBrutOutflow    = sum(md.WwTRANSFER[nt,ktrans]                               for ktrans in md.ntransfer if md.transfer_us[ktrans] == nc)  #Allocation to downstream transfer scheme- if exists (m3)                   
            DEBUG=md.DEBUG[nt,nc] if md.Options['Debug mode'] == 1 else 0
            
            return UserBrutAll + AgrBrutAll + TransBrutOutflow <= NetInflow + TransNetInflow + RunOff + BaseFlow + DEBUG
                
        def water_usermaxall(md,nt,nu): #Water allocation to non agricultural user limited by demand [time x user] 
            return md.WwSUPPLY[nt,nu] <= md.wUserDem[md.t_year[nt],md.t_month[nt],nu]
                        
        def water_transfercapacity(md,nt,ntrans): #Transfer scheme carrying capacity [time x transfer]          
            return md.WwTRANSFER[nt,ntrans] <= md.wTransCap[ntrans]
        
        def water_lakebalance(md,nt,nla): #Storage controlled by first order reservoir or wetland model [time x reservoir]
            nc = md.res_catch[nla]
            AvStor      = (md.WwRSTORAGE[nt,nla] + md.WwRSTORAGE[md.t_prev[nt],nla])/2 if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else (md.WwRSTORAGE[nt,nla] + md.wStorIni[nla])/2
            StorPrev    = md.WwRSTORAGE[md.t_prev[nt],nla] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wStorIni[nla]
            Inflow      = sum((1-md.wkUP[nla])*md.WwOUTFLOW[nt,kc] for kc in md.ncatch if md.catch_ds[kc] == nc) * (1-md.wFlowLoss[nc])
            RunOff      = (1-md.wkRO[nla]) * md.wRunOff[nt,nc]
            Evap        = (md.wET0[nt,nc]-md.wRainFall[nt,nc]) * md.mm_to_m3perha * (md.wkV[nla]*AvStor + md.wResArea[nla]) #max evaporation is 2/kV to avoid impossible solution (which corresponds to the maximum evaporation if no inflow comes in)
            Outflow     = md.wAlpha[nla] * AvStor

            return md.WwRSTORAGE[nt,nla] == StorPrev + Inflow + RunOff - Evap - Outflow 

        
        def water_groundwaterbalance(md,nt,naq):
            StorPrev        = md.WwGWSTORAGE[md.t_prev[nt],naq] if not (md.Options['Initial time step'] == 1 and md.Options['tini'] == nt) else md.wGwIni[naq]
            Pumping         = sum(md.AwGWSUPPLY[nt,kfz,kcul] for kfz in md.nifzone for kcul in md.nculture if md.fzone_catch[kfz]==md.aqui_catch[naq])
            if md.wGwFlow[naq] != 0:
                return md.WwGWSTORAGE[nt,naq] == StorPrev*exp(-md.wGwFlow[naq]) + (md.wGwRech[nt,naq]-Pumping)/md.wGwFlow[naq]*(1-exp(-md.wGwFlow[naq]))
            else:
                return md.WwGWSTORAGE[nt,naq] == StorPrev - Pumping
            
        def water_reservoircapacity(md,nt,nres): #Storage in reservoir limited by storage capacity [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir':
                return Constraint.Skip
            return md.WwRSTORAGE[nt,nres] <= md.wStorCap[nres] 
        
        def water_finalstorage(md,nres): #Storage in reservoir limited by storage capacity [time x reservoir]
            if md.Options['Lakes'] == 1 and md.res_type[nres] != 'reservoir':
                return Constraint.Skip
            Debug = md.DUMYSTOR[nres] if md.Options['Debug mode'] == 1 else 0 #Tolerance to constraint if physically impossible at very high cost
            return md.WwRSTORAGE[md.Options['tfin'],nres] >= md.wStorFin[nres] - Debug
            
        ##Create constraints##
        md.water_waterbalance       = Constraint(md.ntime, md.ncatch, rule=water_waterbalance)
        md.water_waterbalance2      = Constraint(md.ntime, md.ncatch, rule=water_waterbalance2)   
        md.water_usermaxall         = Constraint(md.ntime, md.nuser,  rule=water_usermaxall)
        md.water_reservoircapacity  = Constraint(md.ntime, md.nres,   rule=water_reservoircapacity)
        md.water_groundwaterbalance = Constraint(md.ntime, md.naquifer, rule=water_groundwaterbalance)       
        md.water_transfercapacity   = Constraint(md.ntime, md.ntransfer, rule=water_transfercapacity)
        md.water_lakebalance        = Constraint(md.ntime, md.nlake,  rule=water_lakebalance)
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

        if ReadParam('Options')['Hard eflows'] != 1: 
        #Find most severe years to allow softening of constrain in these years
            years       = [y for y in md.nyear]
            YearlyRunoff= [sum(md.wRunOff[t,c].value for t in md.ntime for c in md.ncatch if md.t_year[t]==y) for y in md.nyear]
            YearsRank   = sorted(range(len(YearlyRunoff)), key=lambda k: YearlyRunoff[k]) #Years ranked from minimal to maximal inflows
            nbSoftYears = int((1-ReadParam('Options')['Hard eflows']) * len(md.nyear)) #Amount of years where the hard constraint can be softened in the optimization period
            SoftYears   = {years[ky]:ky in YearsRank[0:nbSoftYears] for ky in range(len(md.nyear))} #Years where the hard constraint is softened
            
        def env_minflow(md,nt,nef): #env minimum flow [catchment x time]
            if md.eflow_hard[nef] == 0: #If no upstream reservoir is regulating the flow the env flow is limited to actual flow
                EflowConstraint = min(0.95*NatOutflow[nt][md.eflow_catch[nef]],md.wEnvFlow[md.t_month[nt],nef]) #Arbitrary 5% tolerance on consumption of available water
            elif ReadParam('Options')['Hard eflows'] != 1 and SoftYears[md.t_year[nt]] == 1:
                EflowConstraint = min(0.95*NatOutflow[nt][md.eflow_catch[nef]],md.wEnvFlow[md.t_month[nt],nef]) #ADD 0.95 is arbitrary to allow some minimal losses on the way
            else:
                EflowConstraint = md.wEnvFlow[md.t_month[nt],nef]
            Debug = md.DUMYEFLOW[t,nef] if md.Options['Debug mode'] == 1 else 0 #Tolerance to constraint if physically impossible at very high cost
            return md.WwOUTFLOW[nt,md.eflow_catch[nef]] >= EflowConstraint - Debug  # river: env flow does not count transfer scheme water 
            
        if md.Options['Eflows'] == 1:
            md.env_minflow = Constraint(md.ntime, md.neflow, rule=env_minflow)
            
        #%%----------------------------Agricultural module----------------------------
        ##Constraint catalogue##
        
        def agr_waterall(md,nt,nfz,ncul): #Water allocation limit depending on crop area [time x fzone x crop]
            nc=md.fzone_catch[nfz]
            nm=md.t_month[nt]
            ny=md.t_year[nt]
            nft=md.fzone_type[nfz]
            CulArea         = sum(md.AlCULAREA[ny,nfz,kfd,kypt] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
            CropWaterDem    = sum(md.yphase_month[ncul,kyps,nm] * md.wET0[nt,nc] * md.aKc[ncul,kyps] * md.mm_to_m3perha     #Monthly crop water demand of phase                           
                            * sum(md.AlCULAREA[ny,nfz,kfd,kypt]*md.aYieldMat[kypt,kyps] for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==ncul)
                            for kyps in md.nyphase) 
            RainWater       = md.wRainFall[nt,nc] * md.mm_to_m3perha * CulArea * sum(md.yphase_month[ncul,kyps,nm] for kyps in md.nyphase)  
            if md.aIrrigation[nft] == 0: #Rainfed
                if sum(md.yphase_month[ncul,kyps,nm]*value(md.wET0[nt,nc])*md.aKc[ncul,kyps] for kyps in md.nyphase) == 0: #No water demand of this culture at this month
                    return Constraint.Skip #As otherwise 0 >= 0 is not a permited constraint for the model
                else:
                    return  RainWater >= CropWaterDem # Water requierment has to be less than precipitation
            
            elif md.aIrrigation[nft] == 1: #Irrigated
                Pumping = sum(md.AwGWSUPPLY[nt,nfz,ncul] for aq in md.naquifer if md.fzone_catch[nfz] == md.aqui_catch[aq])                             
                return (md.AwSUPPLY[nt,nfz,ncul] + Pumping)*(1-md.aCulRturn[nft,ncul])  >= (CropWaterDem-RainWater)  #Water supply plus rainfall has to be over water dem 
            
        def agr_maxland(md,ny,nfz): #Land limit constraint on cultivated area [year x fzone] 
            CultivatedLand  = sum(md.AlCULAREA[ny,nfz,kfd,kypt] for kfd in md.nfieldculture for kypt in md.nypath if kfd[1]==1) #field[1]=1 (first culture of field) is the reference for used area
            
            return CultivatedLand <= md.aLandCap[ny,nfz]          
                
        def agr_fieldarea(md,ny,nfz,nfd,nfdc): #Every culture in a same field has the same area [year x fzone x ftype x (field x culture)] 
            if nfdc==1:
                return Constraint.Skip
            else:
                return sum(md.AlCULAREA[ny,nfz,nfd,nfdc,kypt] for kypt in md.nypath) == sum(md.AlCULAREA[ny,nfz,nfd,1,kypt] for kypt in md.nypath) #area of field[1] = 2,3...(if exist) have to match the one of field[1] = 1           
        
        def agr_maxcularea(md,ny,nfz,ncul): #Limits the cultivated area by an upper value (usually observed culture area) [year x fzone x culture]
            return sum(md.AlCULAREA[ny,nfz,kfd,kypt] for kfd in md.nfieldculture for kypt in md.nypath if md.field_culture[kfd]==ncul) <= 1.05*md.aCulMax[ny,nfz,ncul] #ADD 5% tolerance is empirical
        
        def agr_fixedculture(md,ny,nfz,nfd,nfdc): #Forces one culture area choice for the whole optimization period [year x fzone x ftype x (field x culture)] 
            return sum(md.AlCULAREA[ny,nfz,nfd,nfdc,kypt] for kypt in md.nypath) == sum(md.AlCULAREA[md.t_year[md.Options['tini']],nfz,nfd,nfdc,kypt] for kypt in md.nypath)
        
        def agr_cropprod(md,ny,nfz,ncr): #Crop production [year x fzone x crop]
            nft=md.fzone_type[nfz]
            #linearized yield response: crop production is linked to yield path choice        
            FzoneCropProd       = sum(md.aCulYield[ny,nft,md.field_culture[kfd]] 
                                * sum(md.AlCULAREA[ny,nfz,kfd,kypt] 
                                * (1-sum(md.akY[md.field_culture[kfd],kyps]*(1-md.aYieldMat[kypt,kyps]) for kyps in md.nyphase)) 
                                for kypt in md.nypath) for kfd in md.nfieldculture if md.culture_crop[md.field_culture[kfd]]==ncr) 
            return md.AcPROD[ny,nfz,ncr] == FzoneCropProd
        

        ##Create constraints##
        md.agr_waterall     = Constraint(md.ntime, md.nfzone, md.nculture, rule=agr_waterall)
        md.agr_cropprod     = Constraint(md.nyear, md.nfzone, md.ncrop,         rule=agr_cropprod)
        md.agr_maxland      = Constraint(md.nyear, md.nfzone,                   rule=agr_maxland)
        md.agr_fieldarea    = Constraint(md.nyear, md.nfzone, md.nfieldculture, rule=agr_fieldarea)
        if md.Options['Culture max area'] == 1:
            md.agr_maxcularea   = Constraint(md.nyear, md.nfzone, md.nculture,      rule=agr_maxcularea)
        if md.Options['Fixed culture'] == 1 and md.Options['Culture max area'] != 'fixed' and MPC == 0:
            md.agr_fixedculture = Constraint(md.nyear, md.nfzone, md.nfieldculture, rule=agr_fixedculture)
        
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
            return CropSupply == CropExtProd + CropLocalProd + CropLocalImport - CropExport
               
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
        def engy_hpdischarge(md,nt,nres): #Hydropower discharge [time x reservoir]
            if nres not in [md.hp_res[khp] for khp in md.nhpp]: #No hydropower in reservoir
                return Constraint.Skip
            else:
                HppDischarge = sum(md.EwHPDISCHARGE[nt,kpld,khp] for kpld in md.npload for khp in md.nhpp if md.hp_res[khp] == nres)
            return HppDischarge <= md.WwOUTFLOW[nt,md.res_catch[nres]]
                
        def engy_hpdischarge2(md,nt,npld,nhp): #Hydropower discharge for ROR hydropowers [hydropower]
            if md.hp_res[nhp]=='ROR': #Run-Of-the-River Hydropower
                return md.EwHPDISCHARGE[nt,npld,nhp] <= md.WwOUTFLOW[nt,md.hp_catch[nhp]] * md.eLoadTime[npld]
            else:
                return Constraint.Skip
                
        def engy_hpprod(md,nt,npld,nhp): #Hydropower production limited by downstream allocation [time x pload x hpp]            
            return md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * md.eHppProd[nhp] * md.EwHPDISCHARGE[nt,npld,nhp]
                
        def engy_hpcapacity(md,nt,npld,nhp): #Hydropower production per load segment limited by hydropower capacity [time x pload x hpp] 
            return md.EeHPPROD[nt,npld,nhp] <= md.eHppEff[nhp] * md.eHppCap[nhp] * md.MW_to_GWhperMonth * md.eLoadTime[npld]                 
        
        def engy_oppcapacity(md,nt,npld,nop): #Powerplant production per load segment limited by capacity [time x pload x opp]
            LoadCapacityFactor = md.eLoadCap[npld,md.op_ptech[nop]] if md.Options['Load capacity'] == 1 else 1                 
            return md.EeOPPROD[nt,npld,nop] <= md.eOppCap[nop] * md.MW_to_GWhperMonth * md.eLoadTime[npld] * LoadCapacityFactor
        
        def engy_gencapacity(md,nt,npld,npt,npm): #Capacity constraint of Generic capacity [time x pload x ptech x pmarket]
            LoadCapacityFactor = 1 #No particular restriction if not specified
            if md.Options['Load capacity'] == 1:
                LoadCapacityFactor = md.eLoadCap[npld,npt]
            return md.EeGENPROD[nt,npld,npt,npm] <= sum(md.EeGENCAP[ky,npt,npm] for ky in md.nyear if ky <= md.t_year[nt]) * md.MW_to_GWhperMonth * md.eLoadTime[npld] * LoadCapacityFactor 

        def engy_maxgeninvest(md,npt,npm): #Maximum increase in Generic capacity
            if md.eMaxCap[npt,npm] == 'unlimited': #hard coded key word for unlimited capacity investment
                return Constraint.Skip
            else:
                return sum(md.EeGENCAP[y,npt,npm] for y in md.nyear) <= md.eMaxCap[npt,npm]            
        
        #Create Constraints
        md.engy_hpdischarge     = Constraint(md.ntime, md.nres,            rule=engy_hpdischarge)
        md.engy_hpdischarge2    = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_hpdischarge2)
        md.engy_hpprod          = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_hpprod)
        md.engy_hpcapacity      = Constraint(md.ntime, md.npload, md.nhpp, rule=engy_hpcapacity)    
        md.engy_oppcapacity     = Constraint(md.ntime, md.npload, md.nopp, rule=engy_oppcapacity)                         
        md.engy_gencapacity     = Constraint(md.ntime, md.npload, md.nptech, md.npmarket, rule=engy_gencapacity)
        md.engy_maxgeninvest    = Constraint(md.nptech, md.npmarket, rule=engy_maxgeninvest)

            
            
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
            return md.EeTRANS[nt,npld,ntl] <= md.eTransCap[ntl] * md.MW_to_GWhperMonth * md.eLoadTime[npld]
                
        def engy_demand(md,nt,npld,npm): #Power demand satisfaction
            Supply      = md.EeSUPPLY[nt,npld,npm] 
            LoadSegment = md.eEngyDem[md.t_year[nt],md.t_month[nt],npm] * md.eLoadDem[npld] 
            return Supply <= LoadSegment
                
        #Create constraints
        md.engy_balance         = Constraint(md.ntime, md.npload, md.npmarket, rule=engy_balance)
        md.engy_transmissioncap = Constraint(md.ntime, md.npload, md.ntransline, rule=engy_transmissioncap)
        md.engy_demand          = Constraint(md.ntime, md.npload, md.npmarket, rule=engy_demand)
        
#%%SAVE MODEL
        self.model=md


