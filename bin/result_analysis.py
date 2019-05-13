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

#RESULT ANALYSIS
from __future__       import division #As in python 2.7 division between integers gives an integer - avoid problems
import pandas as pd
import time
import pickle
from openpyxl         import load_workbook
from pyomo.environ    import Constraint

class ResultAnalysis():
    #%% Save all model parameters and solved decision variables
    def saveall(self,ExportFile,md,parameters,scenario='WHATIF_main'):
        self.model=md
        self.parameters=parameters
        self.scenario=scenario
        pickle.dump(self,open(ExportFile,"wb"))
    #%% Read only Selected result for scenario analysis (WHATIF_scenario.py)
    def selectedresults(self,md,parameters,scenario='WHATIF_main', MPC=0):
        #Length of index function
        def leni(index):
            return max(1,len(index))
        ##SHADOWPRICES##
        #Water shadowprice
        WaterValueShadow        = -sum(md.dual[md.water_waterbalance[t,c]] for t in md.ntime for c in md.ncatch)/(leni(md.ntime)*leni(md.ncatch))
        #Energy shadowprice
        if md.Options['Energy market'] == 1: 
            self.EnergyValue_co = {co:-sum(md.dual[md.engy_balance[t,pld,pm]]*md.EeSUPPLY[t,pld,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co)/sum(md.EeSUPPLY[t,pld,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co) for co in md.ncountry if sum(md.EeSUPPLY[t,pld,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co) != 0}
        else:
            self.EnergyValue_co = {co:-sum(md.eHppVal[hp] for hp in md.nhpp if md.hp_country[hp]==co)/sum(1 for hp in md.nhpp if md.hp_country[hp]==co) for co in md.ncountry if sum(1 for hp in md.nhpp if md.hp_country[hp]==co) != 0}
        #Crop shadowprice
        if md.Options['Crop market'] == 1:
            self.CropValue_co   = {co:-sum(md.dual[md.agr_cropbalance[y,cm,cr]]*sum(md.AcSUPPLY[y,cm,cr,cds].value for cds in md.ncdstep) for y in md.nyear for cr in md.ncrop for cm in md.ncmarket if md.cmarket_country[cm]==co)/sum(md.AcSUPPLY[y,cm,cr,cds].value for y in md.nyear for cr in md.ncrop for cm in md.ncmarket for cds in md.ncdstep if md.cmarket_country[cm]==co) for co in md.ncountry if sum(md.AcSUPPLY[y,cm,cr,cds].value for y in md.nyear for cr in md.ncrop for cm in md.ncmarket for cds in md.ncdstep if md.cmarket_country[cm]==co) != 0}      
        else:
            self.CropValue_co   = {co:-sum(md.aFarmVal[co,cr] for cr in md.ncrop)/leni(md.ncrop) for co in md.ncountry}
        ##WATER##
        if md.Options['Initial time step'] == 0: #Cyclic Run
            ReservoirEvaporation = sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[t,res].value+md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime for res in md.nres)/leni(md.nyear)
        else: #Initial time step
            ReservoirEvaporation = (sum(max(0,md.wET0[md.Options['tini'],md.res_catch[res]].value-1*md.wRainFall[md.Options['tini'],md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[md.Options['tini'],res].value + md.wStorIni[res].value)/2 + md.wResArea[res]) for res in md.nres )
                                 +  sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[t,res].value + md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime for res in md.nres if t != md.Options['tini'])
                                 /  leni(md.nyear))
        #Cultivated area
        GROSSAREA   = {(y,fz,cul):sum(md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if md.field_culture[fd]==cul) for cul in md.nculture for fz in md.nfzone for y in md.nyear}
        NETAREA     = {(y,fz):sum(md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if fd[1]==1) for fz in md.nfzone for y in md.nyear}
        
        ##ECONOMIC RESULTS##
        #User benefit
        UserBenefit_co          = {co:sum(md.WwSUPPLY[t,u].value * md.wUserVal[u] for t in md.ntime for u in md.nuser if md.user_country[u]==co)/leni(md.nyear) for co in md.ncountry}
        #Supply costs
        WaterSupplyCosts_co     = {co:sum(md.WwSUPPLY[t,u].value/(1-md.wUserLoss[u]) * md.wSupCost[u] for t in md.ntime for u in md.nuser if md.user_country[u]==co)/leni(md.nyear) for co in md.ncountry}
        #WaterBalance
        WaterBalance_co         = {co: UserBenefit_co[co] - WaterSupplyCosts_co[co] for co in md.ncountry}
        #Crop, cultivation and irrigation costs
        CropCultivationCost_co  = {co:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
        CropIrrigationCost_co   = {co:sum(md.AwSUPPLY[t,fz,cul].value/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
                                   
         
        #Crop transport costs
        CropTransportCost_co    = {co:sum(md.AcTRANS[y,ct,cr].value * md.aTransCost[ct,cr] for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co)/leni(md.nyear) for co in md.ncountry}
        #Crop import/export benefit and costs       
        CropExtProdCost_co      = {co:sum(md.AcEXTPROD[y,cm,cr].value * md.aCropVal[cm,cr] for y in md.nyear for cm in md.nextcmarket for cr in md.ncrop if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry}  #The assumption is that external market produce crops at their market price #ADD: crop supply curve
        self.CropExpBenefit_co  = {co:sum(md.AcTRANS[y,ct,cr].value * -md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransIn[ct]]==co)/leni(md.nyear) for co in md.ncountry}
        CropImpCost_co          = {co:sum(md.AcTRANS[y,ct,cr].value * -md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co)/leni(md.nyear) for co in md.ncountry}          
        #Crop benefits
        CropBenefit_co          = {co:sum(md.AcSUPPLY[y,cm,cr,cds].value * md.aStepVal[cm,cr,cds]*md.aCropVal[cm,cr] for y in md.nyear for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry}
        if md.Options['Crop market'] == 0:
            CropBenefit_co      = {co:sum(md.AcPROD[y,fz,cr].value * md.aFarmVal[md.fzone_country[fz],cr] for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
        #Crop production NET benefit    
        FarmBenefits_fz_y       = {fz:{y: sum(md.AcPROD[y,fz,cr].value * md.aCropVal[md.fzone_cmarket[fz],cr] for cr in md.ncrop)
                                            - sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for cul in md.nculture)
                                            - sum(md.AwSUPPLY[t,fz,cul].value/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz] for t in md.ntime for cul in md.nculture if md.t_year[t]==y and md.aIrrigation[md.fzone_type[fz]]==1)
                                            for y in md.nyear} 
                                       for fz in md.nfzone}
        #Consummer and Producer surplus
        self.AgrConsSurplus_co  = {co:sum(md.AcSUPPLY[y,cm,cr,cds].value * (md.aStepVal[cm,cr,cds]*md.aCropVal[cm,cr] - -md.dual[md.agr_cropbalance[y,cm,cr]])
                                      for y in md.nyear for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co)/leni(md.nyear) 
                                   for co in md.ncountry}
        self.AgrProdSurplus_co  = {co:sum(sum(md.AcPROD[y,fz,cr].value for fz in md.nfzone if md.fzone_cmarket[fz]==cm) * -md.dual[md.agr_cropbalance[y,cm,cr]] for y in md.nyear for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co)/leni(md.nyear) 
                                    - sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) 
                                    - sum(md.AwSUPPLY[t,fz,cul].value/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz]  for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) 
                                   for co in md.ncountry}    
    #Balances
        CropBalance_co          = {co:CropBenefit_co[co] - CropCultivationCost_co[co] - CropIrrigationCost_co[co] + self.CropExpBenefit_co[co] - CropExtProdCost_co[co] - CropImpCost_co[co] - CropTransportCost_co[co] for co in md.ncountry}
        CropBalance2_co         = {co:sum(FarmBenefits_fz_y[fz][y] for fz in md.nfzone for y in md.nyear if md.fzone_country[fz]==co) for co in md.ncountry}
        
    #Energy benefits        
        EnergyBenefit_co    = {co:sum(md.EeSUPPLY[t,pld,pm].value * md.eEngyVal[md.t_month[t],pm] for t in md.ntime for pm in md.npmarket for pld in md.npload if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}            
    #Energy import/export (rem: if 2 power markets are in the same country, the transmission will be counted as export and import - does not affect the balance)
        EnergyImpCost_co    = {co:sum(md.EeTRANS[t,pld,tl].value * -md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co)/leni(md.nyear) for co in md.ncountry}
        EnergyExpBenefit_co = {co:sum(md.EeTRANS[t,pld,tl].value * -md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransIn[tl]]==co)/leni(md.nyear) for co in md.ncountry}
        self.EnergyExports_co={co:sum(md.EeTRANS[t,pld,tl].value                                                    for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransIn[tl]]==co)/leni(md.nyear) for co in md.ncountry}
    #Energy transmission costs
        EnergyTransCost_co  = {co:sum(md.EeTRANS[t,pld,tl].value * md.eTransCost[tl] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co)/leni(md.nyear) for co in md.ncountry}
    #Operation costs of hydropower plants        
        HpOMCost_co         = {co:sum(md.EeHPPROD[t,pld,hp].value * md.eHppCost[hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.pmarket_country[md.hp_pmarket[hp]]==co)/leni(md.nyear) for co in md.ncountry}
    #Operation and fuel costs of other power plants                     
        OppOMCost_co        = {co:sum(md.EeOPPROD[t,pld,pp].value * md.eOppCost[pp] for t in md.ntime for pld in md.npload for pp in md.nopp if md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}                                                                  
        OppFuelCost_co      = {co:sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCost[md.t_year[t],md.op_pmarket[pp],fu] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}                 
        OppCO2Cost_co       = {co:sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCO2[fu]*md.eCO2Val[md.t_year[t],md.op_pmarket[pp]] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}
        self.CO2Emission_co = {co:sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCO2[fu]                                 for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}                             
    #Capacity, Operation and fuel costs of generic power technologies 
        GenCapCost_co       = {co:sum(md.EeGENCAP[y,pt,pm].value*md.eCAPEX[y,pt,pm]*(md.t_year[md.Options['tfin']]-y+1)/md.eLifeTime[pt,pm] + sum(md.EeGENCAP[ky,pt,pm].value for ky in md.nyear if ky <= y)*md.eFixOPEX[pt,pm] for y in md.nyear for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry} 
        GenProdCost_co      = {co:sum(md.EeGENPROD[t,pld,pt,pm].value*md.eVarOPEX[pt,pm] for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}                                                              
        GenFuelCost_co      = {co:sum(md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCost[md.t_year[t],pm,fu]            for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}  
        GenCO2Cost_co       = {co:sum(md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCO2[fu]*md.eCO2Val[md.t_year[t],pm] for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry} 
        self.CO2Emission_co = {co:self.CO2Emission_co[co]+sum(md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCO2[fu]     for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}     
    #CONSUMER AND PRODUCER SURPLUS 
        self.EgyConsSurplus_co  = {co:sum(md.EeSUPPLY[t,pld,pm].value * (md.eEngyVal[md.t_month[t],pm] - -md.dual[md.engy_balance[t,pld,pm]]/(1-md.eSupLoss[pm]))
                                      for t in md.ntime for pm in md.npmarket for pld in md.npload if md.pmarket_country[pm]==co)/leni(md.nyear) 
                                   for co in md.ncountry}
        self.EgyProdSurplus_co  = {co: sum(md.EeHPPROD[t,pld,hp].value * (-md.dual[md.engy_balance[t,pld,md.hp_pmarket[hp]]]) for t in md.ntime for pld in md.npload for hp in md.nhpp if md.pmarket_country[md.hp_pmarket[hp]]==co)/leni(md.nyear) 
                                    + sum(md.EeOPPROD[t,pld,pp].value * (-md.dual[md.engy_balance[t,pld,md.op_pmarket[pp]]]) for t in md.ntime for pld in md.npload for pp in md.nopp if md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) 
                                    + sum(md.EeGENPROD[t,pld,pt,pm].value * (-md.dual[md.engy_balance[t,pld,pm]]) for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear)
                                    - HpOMCost_co[co] 
                                    - OppOMCost_co[co] - OppFuelCost_co[co] - OppCO2Cost_co[co] 
                                    - GenCapCost_co[co] - GenProdCost_co[co] - GenFuelCost_co[co] - GenCO2Cost_co[co]
                                   for co in md.ncountry}
    #Balances
        EngyBalance_co      = {co:EnergyBenefit_co[co] + EnergyExpBenefit_co[co] - GenCapCost_co[co] - GenProdCost_co[co] - GenFuelCost_co[co] - GenCO2Cost_co[co]
                                - EnergyImpCost_co[co] - OppFuelCost_co[co] - OppCO2Cost_co[co] - HpOMCost_co[co] - OppOMCost_co[co] - EnergyTransCost_co[co] for co in md.ncountry}        
        
        EconomicBalance         = sum(UserBenefit_co[co] for co in md.ncountry) + sum(CropBalance_co[co] for co in md.ncountry) + sum(EngyBalance_co[co] for co in md.ncountry)
        self.EconomicBalance_co = {co: WaterBalance_co[co] + CropBalance_co[co] + EngyBalance_co[co] for co in md.ncountry}
        self.EconomicBalance    = {'Total':EconomicBalance,'Water':sum(WaterBalance_co[co] for co in md.ncountry),'Energy':sum(EngyBalance_co[co] for co in md.ncountry),'Agriculture':sum(CropBalance_co[co] for co in md.ncountry)}
        self.WaterEconomicBalance_co = UserBenefit_co
        self.EnergyEconomicBalance_co = EngyBalance_co
        self.AgricultureEconomicBalance_co = CropBalance_co
        self.AgricultureEconomicBalance2_co = CropBalance2_co

        ##OTHER KEY INDICATORS##
        #Downstream flow [Mm3]
        self.DownstreamFlow = sum(md.WwOUTFLOW[t,c].value for t in md.ntime for c in md.ncatch if md.catch_ds[c]=='outlet')/leni(md.nyear)
        #Irrigation net water consumption [Mm3]
        self.IrrigNetCons_co= {co:sum(md.AwSUPPLY[t,fz,cul].value * (1/(1-md.aIrrgLoss[fz]) - md.aCulRturn[md.fzone_type[fz],cul]) for t in md.ntime for cul in md.nculture for fz in md.nifzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
        #Net cultivated area [Mha]
        self.NetCulArea_co  = {co:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
        #Gross cultivated area [Mha]
        self.GrossCulArea_co= {co:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
        #Gross irrigated area [Mha]
        self.GrossIrrArea_co= {co:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co and md.aIrrigation[md.fzone_type[fz]]==1)/leni(md.nyear) for co in md.ncountry}
        #Net irrigated area [Mha]
        self.NetIrrArea_co= {co:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co and md.aIrrigation[md.fzone_type[fz]]==1)/leni(md.nyear) for co in md.ncountry}
        #Production [Mt]
        self.CropProd_co    = {co:sum(md.AcPROD[y,fz,cr].value for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}

        #Hydropower prod
        self.HydropowerProd_co  = {co:sum(md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload for hp in md.nhpp if md.pmarket_country[md.hp_pmarket[hp]]==co)/leni(md.nyear) for co in md.ncountry} 
        self.HydroFirmProd_co   = {co:sum(sorted([sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)/100*5))] for hp in md.nhpp if md.pmarket_country[md.hp_pmarket[hp]]==co) for co in md.ncountry}
        self.HydroFirmVal_co    = {co:sum(sorted([sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)/100*5))]*leni(md.nmonth)*0.06
                                              +(sum(md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload)/leni(md.nyear)-sorted([sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)/100*5))]*leni(md.nmonth))*0.02
                                              for hp in md.nhpp if md.pmarket_country[md.hp_pmarket[hp]]==co) for co in md.ncountry}
        #Energy exports
        self.EnergyExportBenefit_co = EnergyExpBenefit_co
        #Energy capacity investment #ADD WARNING THIS IS CASE SPECIFIC - Will not cause bug however !
        self.EnergyCapacityConst_co= {co:sum(md.EeGENCAP[y,pt,pm].value for y in md.nyear for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co and pt in ['CoalPP','GasPP']) for co in md.ncountry} 
        self.SolarCapacityConst_co = {co:sum(md.EeGENCAP[y,pt,pm].value for y in md.nyear for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co and pt =='SolarPV') for co in md.ncountry} 
        #Energy curtailment
        self.EnergyCurt_co      = {co:sum(md.eEngyDem[md.t_year[t],md.t_month[t],pm]-sum(md.EeSUPPLY[t,pld,pm].value for pld in md.npload) for t in md.ntime for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}
        self.OtherIndicators = {'Available runoff [Mm3]':        sum(md.wRunOff[t,c].value for t in md.ntime for c in md.ncatch)/leni(md.nyear),
                                'Downstream flow [Mm3]':         self.DownstreamFlow,
                                'Net irrigation cons [Mm3/year]':sum(self.IrrigNetCons_co[co] for co in md.ncountry),
                                'Gross Cultivated area [Mha]':  sum(self.GrossCulArea_co[co] for co in md.ncountry),
                                'Gross Irrigated area [Mha]':   sum(self.GrossIrrArea_co[co] for co in md.ncountry),
                                'Net Irrigated area [Mha]':     sum(self.NetIrrArea_co[co] for co in md.ncountry),
                                'Net Cultivated area [Mha]':    sum(self.NetCulArea_co[co] for co in md.ncountry),
                                'Crop exports benefits [M$/year]':sum(self.CropExpBenefit_co[co] for co in md.ncountry),
                                'Hydropower production [GWh]':  sum(self.HydropowerProd_co[co] for co in md.ncountry),
                                'Hydropower prod FIRM [GWh/month]':   sum(self.HydroFirmProd_co[co] for co in md.ncountry),
                                'Hydropower val FIRM [M$/year]':    sum(self.HydroFirmVal_co[co] for co in md.ncountry),
                                'Energy Curtailment [GWh]':     sum(self.EnergyCurt_co[co] for co in md.ncountry),
                                'Energy trade [GWh/year]':      sum(self.EnergyExports_co[co] for co in md.ncountry),
                                'Generic Power cap [MW]':       sum(self.EnergyCapacityConst_co[co] for co in md.ncountry),
                                'Generic Solar cap [MW]':       sum(self.SolarCapacityConst_co[co] for co in md.ncountry),
                                'CO2 emissions [Mt/year]':      sum(self.CO2Emission_co[co] for co in md.ncountry),
                                'Reservoir Evaporation [Mm3/year]':ReservoirEvaporation,
                                'Water shadowprice [$/m3]':     WaterValueShadow,
                                'Energy Shadowprice [$/kWh]':   sum(self.EnergyValue_co[co] for co in self.EnergyValue_co.keys() if co != 'SouthAfrica')/(leni(self.EnergyValue_co.keys())-1) if md.Options['Energy market'] == 1 else 0,
                                'Crop Shadowprice [$/t]':       sum(self.CropValue_co[co] for co in self.CropValue_co.keys() if co != 'World')/(leni(self.CropValue_co.keys())-1),
                                'Agr. cons. surplus [M$/year]': sum(self.AgrConsSurplus_co[co] for co in md.ncountry),
                                'Agr. prod. surplus [M$/year]': sum(self.AgrProdSurplus_co[co] for co in md.ncountry),
                                'Egy. cons. surplus [M$/year]': sum(self.EgyConsSurplus_co[co] for co in md.ncountry),
                                'Egy. prod. surplus [M$/year]': sum(self.EgyProdSurplus_co[co] for co in md.ncountry)} #ADD STUDY CASE SPECIFIC !!
        
#%% Read all results       
    def readresults(self,md,parameters,solverstatus,PRINT=0,VALIDATION=1,scenario='WHATIF_main',MPC=0):        
        
        def ReadParam(ParamName,option=1,index=0):
            if ParamName not in parameters.val.keys() or option == 0: #parameter is not activated
                if index==1:
                    return [] #returns empty set
                return {} #returns empty parameter
            if ParamName in parameters.scen.keys():                
                return parameters.val[ParamName][parameters.val[parameters.scen[ParamName]][scenario]]
            else:
                return parameters.val[ParamName]
#MAIN
    # Optimization Parameters
        self.MainTable={'Time and date of run': time.strftime("%c"),
                        'Number of years': len(md.nyear),
                        'Run time [s]': solverstatus.solver.Time,
                        'Scenario': scenario}                    

#%%SHADOWPRICES
        CONSTRAINTLIST = [const.name for const in md.component_objects(Constraint)]
        #Length of index function (return 1 for void index instead of 0 to avoid division by 0)
        def leni(index):
            return max(1,len(index))
        
        #Water shadow value
        self.WaterValueShadow_ca= {c: -sum(md.dual[md.water_waterbalance[t,c]]   for t in md.ntime)                     / leni(md.ntime)                 for c in md.ncatch}
        #Hydropower turbine capacity shadow value
        self.HpCapShadow_hp     = {hp:-sum(md.dual[md.engy_hpcapacity[t,pld,hp]] for t in md.ntime for pld in md.npload)/(leni(md.ntime)*leni(md.npload)) for hp in md.nhpp}
        #Reservoir storage capacity shadow value
        self.ResCapShadow_res   = {res:-sum(md.dual[md.water_reservoircapacity[t,res]] for t in md.ntime)               / leni(md.ntime)                 for res in md.nres if res not in md.nlake}    
        #Environmental flow
        self.EnvFlowShadow_m_ec = {m:{ec:sum(md.dual[md.env_minflow[t,ec]] for t in md.ntime if md.t_month[t]==m)/leni(md.nyear) for ec in md.neflow} for m in md.nmonth}        
        #Energy value
        EnergyValueShadow       = -sum(md.dual[md.engy_balance[t,pld,pm]]*md.eLoadDem[pld] for t in md.ntime for pld in md.npload for pm in md.npmarket)/(leni(md.ntime)*leni(md.npmarket))
        #Transmition capacity
        TransmitionCapShadow    = -sum(md.dual[md.engy_transmissioncap[t,pld,tl]] for t in md.ntime for pld in md.npload for tl in md.ntransline)/(leni(md.ntime)*leni(md.npload)*leni(md.ntransline))
        #Power plants capacity
        self.PpCapShadow_pp = {pp:-sum(md.dual[md.engy_oppcapacity[t,pld,pp]] for t in md.ntime for pld in md.npload)/(leni(md.ntime)*leni(md.npload)) for pp in md.nopp}                       
         #Land capacity shadow value
        self.AvLandShadow_fz    = {fz:-sum(md.dual[md.agr_maxland[y,fz]]         for y in md.nyear) /(leni(md.nyear)) if 'agr_maxland' in CONSTRAINTLIST else 0 for fz in md.nfzone}
        #Crop price is determined by cropbalance - of farm value if crop market option is off
        if md.Options['Crop market'] == 1:
            self.CropPrice_cr_y = {cr:{y: sum(-md.dual[md.agr_cropbalance[y,cm,cr]] for cm in md.ncmarket)/leni(md.ncmarket) for y in md.nyear} for cr in md.ncrop}
            self.CropPrice_cr_cm= {cr:{cm:sum(-md.dual[md.agr_cropbalance[y,cm,cr]] for y in md.nyear)/leni(md.nyear) for cm in md.ncmarket} for cr in md.ncrop}
        else: 
            self.CropPrice_cr_y = {cr:{y: sum(md.aFarmVal[co,cr] for co in md.ncountry)/leni(md.ncountry) for y in md.nyear}  for cr in md.ncrop}
            self.CropPrice_cr_cm= {cr:{co:sum(md.aFarmVal[co,cr] for y in md.nyear)    /leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}

        if PRINT==1:
            print('          <----------------SHADOWPRICES--------------->')
            #Water
            print('Water value [$/m3]')
            print(sum(self.WaterValueShadow_ca[ca] for ca in md.ncatch)/leni(md.ncatch))
            print('Storage capacity [$/m3]')
            print(sum(self.ResCapShadow_res[res] for res in md.nres)/leni(md.nres))
            print('Environmental flow [$/m3]')
            print(sum(self.EnvFlowShadow_m_ec[m][ec] for m in md.nmonth for ec in md.neflow)/leni(md.nmonth)/leni(md.neflow))
            #Agriculture
            print('Available Land [$/ha]')
            print(sum(self.AvLandShadow_fz[fz] for fz in md.nfzone)/leni(md.nfzone))
            print('Crop value [$/t]')
            print(sum(self.CropPrice_cr_y[cr][y] for cr in md.ncrop for y in md.nyear)/leni(md.ncrop)/leni(md.nyear))
            #Energy
            print('Energy value [$/kWh]')
            print(EnergyValueShadow)            
            print('Transmition lines capacity [$/kWh]')
            print(TransmitionCapShadow)
            print('Hydropower capacity [$/kWh]')
            print(sum(self.HpCapShadow_hp[hp] for hp in md.nhpp)/leni(md.nhpp))

#%%WATER MODULE
        #WATER BALANCES
        self.WaterBalance_t_c = {
                'RunOff [Mm3]':                 {t: {c: md.wRunOff[t,c].value for c in md.ncatch} for t in md.ntime},
                'Groundwater recharge [Mm3]':   {t: {c: sum(md.wGwRech[t,aq]              for aq in md.naquifer      if md.aqui_catch[aq]==c)        for c in md.ncatch} for t in md.ntime},
                'Upstream inflow [Mm3]':        {t: {c: sum(md.WwOUTFLOW[t,kc].value      for kc in md.ncatch        if md.catch_ds[kc] == c)        for c in md.ncatch} for t in md.ntime},
                'Transfer inflow [Mm3]':        {t: {c: sum(md.WwTRANSFER[t,ktrans].value for ktrans in md.ntransfer if md.transfer_ds[ktrans] == c) for c in md.ncatch} for t in md.ntime}, 
                'Transfer outflow [Mm3]':       {t: {c: sum(md.WwTRANSFER[t,ktrans].value for ktrans in md.ntransfer if md.transfer_us[ktrans] == c) for c in md.ncatch} for t in md.ntime},
                'Transfer loss [Mm3]':          {t: {c: sum(md.WwTRANSFER[t,ktrans].value * md.wTransLoss[ktrans] for ktrans in md.ntransfer if md.transfer_ds[ktrans] == c) for c in md.ncatch} for t in md.ntime}, 
                'User allocation [Mm3]':        {t: {c: sum(md.WwSUPPLY[t,u].value        for u in md.nuser          if md.user_catch[u]==c)         for c in md.ncatch} for t in md.ntime},
                'User return flow [Mm3]':       {t: {c: sum(md.WwSUPPLY[t,u].value * md.wUserRturn[u] for u in md.nuser if md.user_catch[u]==c)      for c in md.ncatch} for t in md.ntime},
                'User water loss [Mm3]':        {t: {c: sum(md.WwSUPPLY[t,u].value * md.wUserLoss[u]/(1-md.wUserLoss[u]) for u in md.nuser if md.user_catch[u]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture surface [Mm3]':    {t: {c: sum(md.AwSUPPLY[t,fz,cul].value   for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture groundwater [Mm3]':{t: {c: sum(md.AwGWSUPPLY[t,fz,cul].value for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]==md.aqui_catch[aq]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture return flow [Mm3]':{t: {c: sum(md.AwSUPPLY[t,fz,cul].value * md.aCulRturn[md.fzone_type[fz],cul]   for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture water loss [Mm3]' :{t: {c: sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgLoss[fz]/(1-md.aIrrgLoss[fz]) for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c) for c in md.ncatch} for t in md.ntime},
                'Outlet flow [Mm3]':            {t: {c: md.WwOUTFLOW[t,c].value for c in md.ncatch} for t in md.ntime},
                'River ET [Mm3]':               {t: {c: sum(md.WwOUTFLOW[t,kc].value*md.wFlowLoss[c] for kc in md.ncatch if md.catch_ds[kc] == c) for c in md.ncatch} for t in md.ntime},
                'Reservoir storage [Mm3]':      {t: {c: sum(md.WwRSTORAGE[t,res].value for res in md.nres if md.res_catch[res] == c) for c in md.ncatch} for t in md.ntime},
                'Groundwater storage [Mm3]':    {t: {c: sum(md.WwGWSTORAGE[t,aq].value for aq in md.naquifer if md.aqui_catch[aq] == c) for c in md.ncatch} for t in md.ntime}
                                 }
        if md.Options['Initial time step'] == 0: #Cyclic Run
            self.WaterBalance_t_c['Reservoir and lake ET [Mm3]']={t: {c: sum(max(0,md.wET0[t,c].value-1*md.wRainFall[t,c].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[t,res].value+md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for res in md.nres if md.res_catch[res]==c) for c in md.ncatch} for t in md.ntime}    
        else:
            self.WaterBalance_t_c['Reservoir and lake ET [Mm3]']={t: {c: sum(max(0,md.wET0[t,c].value-1*md.wRainFall[t,c].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[t,res].value+md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for res in md.nres if md.res_catch[res]==c) for c in md.ncatch} for t in md.ntime if t != md.Options['tini']}
            self.WaterBalance_t_c['Reservoir and lake ET [Mm3]'][md.Options['tini']]= {c: sum(max(0,md.wET0[md.Options['tini'],c].value-1*md.wRainFall[md.Options['tini'],c].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[md.Options['tini'],res].value + md.wStorIni[res].value)/2 + md.wResArea[res]) for res in md.nres if md.res_catch[res]==c) for c in md.ncatch}
                                 
        self.WaterBalance_t_c[' Balance = 0 [Mm3]']= {t: {c: + self.WaterBalance_t_c['RunOff [Mm3]'][t][c]
                                                             + self.WaterBalance_t_c['Groundwater recharge [Mm3]'][t][c]
                                                             + self.WaterBalance_t_c['User return flow [Mm3]'][t][c]
                                                             + self.WaterBalance_t_c['Agriculture return flow [Mm3]'][t][c]
                                                             + self.WaterBalance_t_c['Transfer inflow [Mm3]'][t][c]
                                                             + self.WaterBalance_t_c['Upstream inflow [Mm3]'][t][c]
                                                             + (self.WaterBalance_t_c['Reservoir storage [Mm3]'][md.t_prev[t]][c] if md.Options['Initial time step'] == 0 or t != md.Options['tini'] else sum(md.wStorIni[res].value for res in md.nres if md.res_catch[res]==c))
                                                             + (self.WaterBalance_t_c['Groundwater storage [Mm3]'][md.t_prev[t]][c] if md.Options['Initial time step'] == 0 or t != md.Options['tini'] else sum(md.wGwIni[aq].value for aq in md.naquifer if md.aqui_catch[aq]==c))
                                                             - self.WaterBalance_t_c['User allocation [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['User water loss [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Agriculture surface [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Agriculture water loss [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Agriculture groundwater [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Transfer outflow [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Transfer loss [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Reservoir storage [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Reservoir and lake ET [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Groundwater storage [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['River ET [Mm3]'][t][c]
                                                             - self.WaterBalance_t_c['Outlet flow [Mm3]'][t][c]
                                                             for c in md.ncatch} for t in md.ntime}
        #WATER BALANCES
        self.WaterBalance_ca = {keys:{ca:sum(self.WaterBalance_t_c[keys][t][ca] for t in md.ntime)/leni(md.nyear) for ca in md.ncatch} for keys in self.WaterBalance_t_c.keys()}
        self.WaterBalance_ca['Initial storage [Mm3]'] = {c:sum(md.wStorIni[res].value for res in md.nres if md.res_catch[res]==c) + sum(md.wGwIni[aq] for aq in md.naquifer if md.aqui_catch[aq]==c) for c in md.ncatch} if md.Options['Initial time step'] == 1 else {c:0 for c in md.ncatch}
        self.WaterBalance_ca['End storage [Mm3]']     = {c:sum(md.WwRSTORAGE[md.Options['tfin'],res].value for res in md.nres if md.res_catch[res]==c) + sum(md.WwGWSTORAGE[md.Options['tfin'],aq].value for aq in md.naquifer if md.aqui_catch[aq]==c) for c in md.ncatch} if md.Options['Initial time step'] == 1 else {c:0 for c in md.ncatch}
        self.WaterBalance_co = {keys:{co:sum(self.WaterBalance_ca[keys][c] for c in md.ncatch if md.catch_country[c] == co) for co in md.ncountry} for keys in self.WaterBalance_ca.keys()}
        self.WaterBalance    = {keys:sum(self.WaterBalance_ca[keys][c] for c in md.ncatch) for keys in self.WaterBalance_ca.keys()}                          
        self.WaterBalance_co['Upstream inflow [Mm3]'] = {co:sum(self.WaterBalance_ca['Outlet flow [Mm3]'][kc] for kc in md.ncatch if (md.catch_ds[kc] != 'outlet' and md.catch_country[md.catch_ds[kc]]==co) and md.catch_country[kc] != co) for co in md.ncountry}
        self.WaterBalance_co['Outlet flow [Mm3]']     = {co:sum(self.WaterBalance_ca['Outlet flow [Mm3]'][kc] for kc in md.ncatch if (md.catch_ds[kc] == 'outlet' or  md.catch_country[md.catch_ds[kc]]!=co) and md.catch_country[kc] == co) for co in md.ncountry} 
        self.WaterBalance['Outlet flow [Mm3]']     = self.WaterBalance['Outlet flow [Mm3]'] - self.WaterBalance['Upstream inflow [Mm3]']
        self.WaterBalance['Upstream inflow [Mm3]'] = 0
        
        #RESERVOIRS
        self.ReservoirTable  = {'AvStorage [Mm3]':               {res:sum(md.WwRSTORAGE[t,res].value/leni(md.ntime) for t in md.ntime) for res in md.nres},
                                'Capacity [Mm3]':                ReadParam('wStorCap'), 
                                'Capacity shadowprice [$/m3]':  self.ResCapShadow_res
                                }
        
        if md.Options['Initial time step'] != 1: #Cyclic Run
            self.ReservoirTable['Evaporation [Mm3/y]']  = {res: sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[t,res].value + md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime)/leni(md.nyear) for res in md.nres}
        else:
            self.ReservoirTable['Evaporation [Mm3/y]']  = {res: (max(0,md.wET0[md.Options['tini'],md.res_catch[res]].value-1*md.wRainFall[md.Options['tini'],md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[md.Options['tini'],res].value + md.wStorIni[res].value)/2 + md.wResArea[res])
                                                               +sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value)                                * md.mm_to_m3perha * (md.wkV[res]*(md.WwRSTORAGE[t,res].value + md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime if t != md.Options['tini']))
                                                               /leni(md.nyear) for res in md.nres}                                                                       
        self.ResStor_res_m   = {m:{res:sum(md.WwRSTORAGE[t,res].value for t in md.ntime if md.t_month[t]==m)/sum(1 for t in md.ntime if md.t_month[t]==m) for res in md.nres} for m in md.nmonth}
        self.ResStor_res_t   = {t:{res:md.WwRSTORAGE[t,res].value for res in md.nres} for t in md.ntime}
        
        #WATER USERS
        self.UserTable       = {'User consumption [Mm3/y]':      {u:sum(md.WwSUPPLY[t,u].value * (1-md.wUserRturn[u])                            for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'User loss [Mm3/y]':             {u:sum(md.WwSUPPLY[t,u].value * md.wUserLoss[u]/(1-md.wUserLoss[u])         for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'Net User abstraction [Mm3/y]':  {u:sum(md.WwSUPPLY[t,u].value * (1/(1-md.wUserLoss[u])-md.wUserRturn[u])   for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'User benefits [M$/y]':          {u:sum(md.WwSUPPLY[t,u].value * md.wUserVal[u]                                         for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'User curtailment [Mm3/y]':      {u:sum(-md.WwSUPPLY[t,u].value + md.wUserDem[md.t_year[t],md.t_month[t],u]                    for t in md.ntime)/leni(md.nyear) for u in md.nuser}
                                }

        #GROUNDWATER
        self.GwStor_c_m  = {m:{c:sum(md.WwGWSTORAGE[t,aq].value for t in md.ntime for aq in md.naquifer if md.t_month[t]==m and md.aqui_catch[aq]==c)/leni(md.nyear) for c in md.ncatch} for m in md.nmonth}
        self.GwStor_c_t  = {t:{c:sum(md.WwGWSTORAGE[t,aq].value for aq in md.naquifer if md.aqui_catch[aq]==c) for c in md.ncatch} for t in md.ntime}
        
        #%%PRINTS
        if PRINT ==1:
            print('          <----------------OPTIMAL DECISIONS--------------->')
            print('*************WATER*************')
            print('Inflow [Mm3]')
            print(self.WaterBalance['RunOff [Mm3]'])
            print('Water consumption per user [Mm3]')
            print('user:',sum(md.WwSUPPLY[t,u].value for t in md.ntime for u in md.nuser)/leni(md.nyear),'demand',sum(md.wUserDem[md.t_year[t],md.t_month[t],u] for t in md.ntime for u in md.nuser)/leni(md.nyear))
            print('Downstream flow [Mm3]')
            print(self.WaterBalance['Agriculture water loss [Mm3]'])
            print('Water losss [Mm3]')
            print('User:',self.WaterBalance['User water loss [Mm3]'])
            print('Agriculture:', self.WaterBalance['Outlet flow [Mm3]'])
            print('Return Flows [Mm3]')
            print(self.WaterBalance['User return flow [Mm3]']+self.WaterBalance['Agriculture return flow [Mm3]'])
            print('End Storage [Mm3]')
            print(self.WaterBalance['End storage [Mm3]'])
            print('Initial Storage [Mm3]')
            print(self.WaterBalance['Initial storage [Mm3]'])
            print('River and Reservoir ET')
            print(self.WaterBalance['River ET [Mm3]']+self.WaterBalance['Reservoir and lake ET [Mm3]'])
            print('Water Balance [Mm3]: InitialStorage + InFlow + ReturnFlow = UserAll + AgrAll + DownFlow + EndStorage + ET-Losses')
            print(self.WaterBalance[' Balance = 0 [Mm3]'])
            
#%% AGRICULTURE MODULE 
        #Available land, avergae cultivated land  
        GROSSAREA   = {(y,fz,cul):sum(md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if md.field_culture[fd]==cul) for cul in md.nculture for fz in md.nfzone for y in md.nyear}
        NETAREA     = {(y,fz):sum(md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if fd[1]==1) for fz in md.nfzone for y in md.nyear}      
        if md.Options['Crop market'] == 1:
            CropPrice   = {(y,fz,cr):-md.dual[md.agr_cropbalance[y,md.fzone_cmarket[fz],cr]] for cr in md.ncrop for fz in md.nfzone for y in md.nyear}
        else:
            CropPrice   = {(y,fz,cr): md.aFarmVal[md.fzone_country[fz],cr] for cr in md.ncrop for fz in md.nfzone for y in md.nyear}
        
        #CROP BALANCE
        self.CropBalance_y_fz={'Available land [Mha]':      {y: {fz:md.aLandCap[y,fz] for fz in md.nfzone} for y in md.nyear},
                               'Net Cultivated land [Mha]': {y: {fz:NETAREA[y,fz] for fz in md.nfzone} for y in md.nyear},
                               'Net Irrigated land [Mha]':  {y: {fz:NETAREA[y,fz] if md.aIrrigation[md.fzone_type[fz]]==1 else 0 for fz in md.nfzone} for y in md.nyear},
                               'Crop production [Mt]':      {y: {fz:sum(md.AcPROD[y,fz,cr].value  for cr in md.ncrop) for fz in md.nfzone} for y in md.nyear}}
        self.CropBalance_ca= {keys:{ca:sum(self.CropBalance_y_fz[keys][y][fz] for y in md.nyear for fz in md.nfzone if md.fzone_catch[fz]==ca)/leni(md.nyear) for ca in md.ncatch} for keys in self.CropBalance_y_fz.keys()}
        
        self.CropBalance_y_cm={'Available land [Mha]':      {y: {cm: sum(md.aLandCap[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Net Cultivated land [Mha]': {y: {cm: sum(NETAREA[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Net Irrigated land [Mha]':  {y: {cm: sum(NETAREA[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm and md.aIrrigation[md.fzone_type[fz]]==1) for cm in md.ncmarket} for y in md.nyear},
                               'Crop production [Mt]':      {y: {cm: sum(md.AcPROD[y,fz,cr].value for fz in md.nfzone for cr in md.ncrop if md.fzone_cmarket[fz]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Crop demand [Mt]':          {y: {cm: sum(md.aCropDem[y,cm,cr] for cr in md.ncrop) for cm in md.ncmarket} for y in md.nyear},
                               'Crop supply [Mt]':          {y: {cm: sum(md.AcSUPPLY[y,cm,cr,cds].value for cr in md.ncrop for cds in md.ncdstep) for cm in md.ncmarket} for y in md.nyear},
                               'Crop external prod [Mt]':   {y: {cm: sum(md.AcEXTPROD[y,cm,cr].value for cr in md.ncrop) if cm in md.nextcmarket else 0 for cm in md.ncmarket} for y in md.nyear},
                               'Crop net import [Mt]':      {y: {cm: sum(md.AcTRANS[y,ct,cr].value * (1-md.aTransLoss[ct]) for ct in md.nctrans for cr in md.ncrop if md.aTransOut[ct]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Crop gross export [Mt]':    {y: {cm: sum(md.AcTRANS[y,ct,cr].value                         for ct in md.nctrans for cr in md.ncrop if md.aTransIn[ct]==cm) for cm in md.ncmarket} for y in md.nyear}}
        self.CropBalance_y_cm[' Balance = 0 [Mt]']      = {y: {cm: self.CropBalance_y_cm['Crop supply [Mt]'][y][cm]
                                                                  +self.CropBalance_y_cm['Crop gross export [Mt]'][y][cm]
                                                                  -self.CropBalance_y_cm['Crop net import [Mt]'][y][cm]
                                                                  -self.CropBalance_y_cm['Crop production [Mt]'][y][cm]
                                                                  -self.CropBalance_y_cm['Crop external prod [Mt]'][y][cm] for cm in md.ncmarket} for y in md.nyear}
                   
        self.CropBalance_y = {keys:{y:sum(self.CropBalance_y_cm[keys][y][cm] for cm in md.ncmarket) for y in md.nyear} for keys in self.CropBalance_y_cm.keys()}
        self.CropBalance_co= {keys:{co:sum(self.CropBalance_y_cm[keys][y][cm] for y in md.nyear for cm in md.ncmarket if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry} for keys in self.CropBalance_y_cm.keys()}
        self.CropBalance   = {keys:sum(self.CropBalance_y_cm[keys][y][cm] for y in md.nyear for cm in md.ncmarket)/leni(md.nyear) for keys in self.CropBalance_y_cm.keys()}
        
        if md.Options['Crop market'] == 0: #BALANCE if crop market is off            
            self.CropBalance_y.update({keys:{y:sum(self.CropBalance_y_fz[keys][y][fz] for fz in md.nfzone) for y in md.nyear} for keys in self.CropBalance_y_fz.keys()})
            self.CropBalance_co.update({keys:{co:sum(self.CropBalance_y_fz[keys][y][fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry} for keys in self.CropBalance_y_fz.keys()})
            self.CropBalance.update({keys:sum(self.CropBalance_y_fz[keys][y][fz] for y in md.nyear for fz in md.nfzone)/leni(md.nyear) for keys in self.CropBalance_y_fz.keys()})
        
        #FARMING ZONES
        self.FarmZonesTable=  {'Crop Production [Mt/y]':        {fz:sum(md.AcPROD[y,fz,cr].value                                                    for y in md.nyear for cr in md.ncrop)    /leni(md.nyear) for fz in md.nfzone},                               
                               'Cultivation costs [M$/y]':      {fz:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul]                    for y in md.nyear for cul in md.nculture)/leni(md.nyear) for fz in md.nfzone},
                               'Irrigation costs [M$/y]':       {fz:sum(md.AwSUPPLY[t,fz,cul].value * 1/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz]     for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'GW pumping costs [M$/y]':       {fz:sum(md.AwGWSUPPLY[t,fz,cul].value * md.wGwCost[md.fzone_catch[fz]]              for t in md.ntime for cul in md.nculture)/leni(md.nyear) if md.Options['Groundwater']==1 else 0 for fz in md.nifzone},
                               'Cultivated land Net [Mha]':     {fz:sum(NETAREA[y,fz]                                                               for y in md.nyear)                       /leni(md.nyear) for fz in md.nfzone},
                               'Cultivated land Net [%]':       {fz:sum(NETAREA[y,fz] for y in md.nyear) / sum(md.aLandCap[y,fz] for y in md.nyear) for fz in md.nfzone if sum(md.aLandCap[y,fz] for y in md.nyear) != 0},
                               'Cultivated land Gross [%]':     {fz:sum(GROSSAREA[y,fz,cul] for y in md.nyear for cul in md.nculture) / sum(md.aLandCap[y,fz] for y in md.nyear) for fz in md.nfzone if sum(md.aLandCap[y,fz] for y in md.nyear) != 0},
                               'Irrig net abstraction [Mm3/y]': {fz:sum(md.AwSUPPLY[t,fz,cul].value * (1/(1-md.aIrrgLoss[fz]) - md.aCulRturn[md.fzone_type[fz],cul]) for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Irrig withdrawals [Mm3/y]':     {fz:sum(md.AwSUPPLY[t,fz,cul].value * 1/(1-md.aIrrgLoss[fz])                        for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Irrig losses [Mm3/y]':          {fz:sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgLoss[fz]/(1-md.aIrrgLoss[fz])         for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Irrig return flow [Mm3/y]':     {fz:sum(md.AwSUPPLY[t,fz,cul].value * md.aCulRturn[md.fzone_type[fz],cul]           for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Precipitation [mm/y]':          {fz:sum(md.wRainFall[t,md.fzone_catch[fz]].value                                    for t in md.ntime)                       /leni(md.nyear) for fz in md.nfzone},
                               'Land shadowprice [$/ha]':       self.AvLandShadow_fz,
                               'Production Value [M$/y]':       {fz:sum(md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for y in md.nyear for cr in md.ncrop)/leni(md.nyear) for fz in md.nfzone}}                       
        
        if md.Options['Culture max area'] == 'fixed' or MPC==1:
            self.FarmZonesTable['Crop dem. satisfaction [%]'] = {fz: 'Not available now' for fz in md.nfzone}
        else:
            self.FarmZonesTable['Crop dem. satisfaction [%]'] = {fz:sum(md.AlCULAREA[y,fz,fd,ypt].value*sum(md.aYieldMat[ypt,kyps] for kyps in md.nyphase)/leni(md.nyphase) for y in md.nyear for fd in md.nfieldculture for ypt in md.nypath)
                                                                   /sum(GROSSAREA[y,fz,cul] for y in md.nyear for cul in md.nculture) for fz in md.nfzone if sum(GROSSAREA[y,fz,cul] for y in md.nyear for cul in md.nculture) != 0}    
        #FARM BENEFITS 
        self.FarmBenefits_ca_ft = {c: {ft:sum(md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_type[fz]==ft and md.fzone_catch[fz]==c)/leni(md.nyear) 
                                        - sum(GROSSAREA[y,fz,cul] * md.aCulCost[ft,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_catch[fz]==c)/leni(md.nyear)
                                        - sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_type[fz]==ft and md.fzone_catch[fz]==c)/leni(md.nyear)
                                        for ft in md.nftype} 
                                   for c in md.ncatch}
        self.FarmBenefits_fz_y  = {fz:{y: sum(md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for cr in md.ncrop)
                                        - sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for cul in md.nculture)
                                        - sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz]  for t in md.ntime for cul in md.nculture if md.t_year[t]==y and md.aIrrigation[md.fzone_type[fz]]==1)
                                        for y in md.nyear} 
                                   for fz in md.nfzone}    
        
        #CROP TRANSPORT
        self.CropTransport_cm_cm = {cma:{cmb:sum(md.AcTRANS[y,ct,cr].value                      for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.aTransIn[ct]==cma and md.aTransOut[ct]==cmb)/leni(md.nyear) for cmb in md.ncmarket} for cma in md.ncmarket}
        self.CropTransport_ct_cr = {ct: {cr: sum(md.AcTRANS[y,ct,cr].value                      for y in md.nyear)/leni(md.nyear) for cr in md.ncrop} for ct in md.nctrans}
        self.CropTransport_ct    = {ct: sum(md.AcTRANS[y,ct,cr].value                           for y in md.nyear for cr in md.ncrop)/leni(md.nyear) for ct in md.nctrans}
        self.CropTransVal_ct     = {ct: sum(md.AcTRANS[y,ct,cr].value * -md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for y in md.nyear for cr in md.ncrop)/leni(md.nyear) for ct in md.nctrans}
        self.CropLoss_cm_cm      = {cma:{cmb:sum(md.AcTRANS[y,ct,cr].value * md.aTransLoss[ct]  for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.aTransIn[ct]==cma and md.aTransOut[ct]==cmb)/leni(md.nyear) for cmb in md.ncmarket} for cma in md.ncmarket}                                              
        self.CropImpShare_cr_co  = {cr: {co:(sum(md.AcTRANS[y,ct,cr].value * (1-md.aTransLoss[ct]) for y in md.nyear for ct in md.nctrans if md.cmarket_country[md.aTransOut[ct]]==co)
                                            -sum(md.AcTRANS[y,ct,cr].value                         for y in md.nyear for ct in md.nctrans if md.cmarket_country[md.aTransIn[ct]]==co))
                                            /sum(md.aCropDem[y,cm,cr]                           for y in md.nyear for cm in md.ncmarket if md.cmarket_country[cm]==co)
                                        for co in md.ncountry if sum(md.aCropDem[y,cm,cr] for y in md.nyear for cm in md.ncmarket if md.cmarket_country[cm]==co) != 0} 
                                    for cr in md.ncrop}
        
        #EXTERNAL MARKETS
        self.CropExtProd_cr_co  = {cr:{co:sum(md.AcEXTPROD[y,cm,cr].value for y in md.nyear for cm in md.nextcmarket if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}
        
        #CROP PRODUCTION
        self.CropProduction_cr_fz   = {cr:{fz:sum(md.AcPROD[y,fz,cr].value for y in md.nyear)/leni(md.nyear) for fz in md.nfzone} for cr in md.ncrop} 
        self.CropProduction_cr_co   = {cr:{co:sum(md.AcPROD[y,fz,cr].value for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}
        
        
        #YIELDS
        self.Yield_cul_ft       = {cul:{ft:sum(md.AlCULAREA[y,fz,kfd,kypt].value*(1-sum(md.akY[md.field_culture[kfd],kyps]*(1-md.aYieldMat[kypt,kyps]) for kyps in md.nyphase)) for y in md.nyear for fz in md.nfzone for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==cul and md.fzone_type[fz]==ft)
                                              /sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone if md.fzone_type[fz]==ft) 
                                               for ft in md.nftype if sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone if md.fzone_type[fz]==ft) != 0} 
                                               for cul in md.nculture}
        
        #LAND USE
        if md.Options['Culture max area'] == 'fixed' or MPC==1:
            self.LandUse_fd_ft      = {cul:{ft:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone if md.fzone_type[fz]==ft)/leni(md.nyear) for ft in md.nftype} for cul in md.nculture}
        else:
            self.LandUse_fd_ft      = {fd[0]:{ft:sum(md.AlCULAREA[y,fz,kfd,kypt].value for y in md.nyear for fz in md.nfzone for kfd in md.nfieldculture for kypt in md.nypath if md.fzone_type[fz]==ft and kfd[0]==fd[0])/leni(md.nyear) for ft in md.nftype} for fd in md.nfieldculture}
        self.LandUse_cul_fz         = {cul:{fz:sum(GROSSAREA[y,fz,cul] for y in md.nyear)/leni(md.nyear) for fz in md.nfzone} for cul in md.nculture}
        self.LandUse_cr_co          = {cr:{co:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co and md.culture_crop[cul]==cr)/leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}
        self.LandUse_cr_y           = {cr:{y: sum(GROSSAREA[y,fz,cul] for fz in md.nfzone for cul in md.nculture if md.culture_crop[cul]==cr) for y in md.nyear} for cr in md.ncrop}
        self.LandUse_co_ft          = {co:{ft:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co and md.fzone_type[fz]==ft)/leni(md.nyear) for ft in md.nftype} for co in md.ncountry}
        self.LandUse_co_irr         = {co:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co and md.aIrrigation[md.fzone_type[fz]]==1)/leni(md.nyear) for co in md.ncountry}
        self.LandUse_cr_ft          = {cr:{ft:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.culture_crop[cul]==cr and md.fzone_type[fz]==ft)/leni(md.nyear) for ft in md.nftype} for cr in md.ncrop}
        self.LandUse_ca   = {'Irrig Land Use [%]'  :{ca:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_catch[fz]==ca and md.aIrrigation[md.fzone_type[fz]] == 1)/leni(md.nyear)/ReadParam('CatchArea')[ca] for ca in md.ncatch},
                             'Irrig Used Cap [%]'  :{ca:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_catch[fz]==ca and md.aIrrigation[md.fzone_type[fz]] == 1)/sum(md.aLandCap[y,fz] for y in md.nyear for fz in md.nfzone if md.aIrrigation[md.fzone_type[fz]]==1 and md.fzone_catch[fz]==ca) for ca in md.ncatch if sum(md.aLandCap[y,fz] for y in md.nyear for fz in md.nfzone if md.aIrrigation[md.fzone_type[fz]]==1 and md.fzone_catch[fz]==ca)}}        
        
        #%%PRINTS
        if PRINT ==1:
            print('*************AGRICULTURE*************')
            print('Available Land [Mha]')
            print(self.CropBalance['Available land [Mha]'])
            print('Average cultivated area [Mha]')
            print(self.CropBalance['Net Cultivated land [Mha]'])
            print('Crop production [Mt]')
            print(self.CropBalance['Crop production [Mt]'])
            print('Net traded crop [Mt]')
            print(self.CropBalance['Crop net import [Mt]'])
            print('Crop transport losses [Mt]')
            print(self.CropBalance['Crop gross export [Mt]']-self.CropBalance['Crop net import [Mt]'])
            print('Crop external production [Mt]')
            print(sum(self.CropExtProd_cr_co[cr][co] for cr in md.ncrop for co in md.ncountry))
            print('Crop demand [Mt]')
            print(self.CropBalance['Crop demand [Mt]'])
            print('Crop balance [Mt]: CropProd = CropTrans - CropLoss')
            print(self.CropBalance[' Balance = 0 [Mt]'])
        

#%% ENERGY  
        #Hydropower plants
        #ADD: WARNING CASE SPECIFIC secondary power VALUE (HppVal/3)
        FIRMPOWER = 95/100 # Percentage of the time the power has to be generated to be firm
 
        self.HydropowerTable    = {'Power Prod [GWh/y]'     : {hp:sum(md.EeHPPROD[t,pld,hp].value  for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp} ,
                                   'Firm Power [GWh/month]' : {hp:sorted([sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)*(1-FIRMPOWER)))] for hp in md.nhpp},
                                   'Firm value [M$/year]'   : {hp:sorted([sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)*(1-FIRMPOWER)))]*leni(md.nmonth)*ReadParam('eHppVal')[hp]
                                                                   + (sum(md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload)/leni(md.nyear)-sorted([sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)*(1-FIRMPOWER)))]*leni(md.nmonth))*ReadParam('eHppVal')[hp]/3
                                                                   for hp in md.nhpp},
                                   'Capacity [MW]'          : ReadParam('eHppCap'), 
                                   'Discharge [Mm3/y]'      : {hp:sum(md.EwHPDISCHARGE[t,pld,hp].value  for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp},
                                   'AvStorage [Mm3/month]'  : {hp:sum(md.WwRSTORAGE[t,res].value for t in md.ntime for res in md.nres if md.hp_res[hp]==res)/leni(md.ntime) for hp in md.nhpp},
                                   'Operation costs [M$/y]' : {hp:sum(md.eHppCost[hp] * md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp},
                                   'Capacity shadowprice [$/GWh/month]':self.HpCapShadow_hp}
        
        self.HpProd_hp_m            = {m:{hp:sum(md.EeHPPROD[t,pld,hp].value        for t in md.ntime for pld in md.npload  if md.t_month[t]==m)                        /leni(md.nyear) for hp in md.nhpp}                           for m in md.nmonth}                      
        self.HpDis_hp_m             = {m:{hp:sum(md.EwHPDISCHARGE[t,pld,hp].value   for t in md.ntime for pld in md.npload  if md.t_month[t]==m)                        /leni(md.nyear) for hp in md.nhpp if md.hp_res[hp]!='ROR'}   for m in md.nmonth}
        self.HpStor_hp_m            = {m:{hp:sum(md.WwRSTORAGE[t,res].value         for t in md.ntime for res in md.nres    if md.hp_res[hp]==res and md.t_month[t]==m) /leni(md.nyear) for hp in md.nhpp}                           for m in md.nmonth}       
        self.HpProd_hp_y            = {y:{hp:sum(md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload if md.t_year[t]==y) for hp in md.nhpp} for y in md.nyear}        
        if md.Options['Energy market'] == 1:
            self.HydropowerTable['Energy value [M$/y]']     = {hp:sum(-md.EeHPPROD[t,pld,hp].value * md.dual[md.engy_balance[t,pld,md.hp_pmarket[hp]]] for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp}
            self.HpBenefits_hp_y    = {y:{hp:sum(md.EeHPPROD[t,pld,hp].value * (-md.dual[md.engy_balance[t,pld,md.hp_pmarket[hp]]] - md.eHppCost[hp]) for t in md.ntime for pld in md.npload if md.t_year[t]==y) for hp in md.nhpp} for y in md.nyear}
        else:
            self.HydropowerTable['Energy value [M$/y]'] = {hp:sum(-md.EeHPPROD[t,pld,hp].value * md.eHppVal[hp] for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp}
            self.HpBenefits_hp_y    = {y:{hp:sum(md.EeHPPROD[t,pld,hp].value * (md.eHppVal[hp]- md.eHppCost[hp]) for t in md.ntime for pld in md.npload if md.t_year[t]==y) for hp in md.nhpp} for y in md.nyear}
                        
        #Energy market prices
        self.EnergyVal_pm_m     = {pm:{m:sum(-md.dual[md.engy_balance[t,pld,pm]]*md.eLoadDem[pld] for t in md.ntime for pld in md.npload if md.t_month[t]==m)/leni(md.nyear) for m in md.nmonth} for pm in md.npmarket}
        
        #Energy Balance
        self.EnergyBalance_t_pm = {'Power demand [GWh]'             :{t:{pm:md.eEngyDem[md.t_year[t],md.t_month[t],pm] for pm in md.npmarket} for t in md.ntime},
                                   'Unserved demand [GWh]'          :{t:{pm:md.eEngyDem[md.t_year[t],md.t_month[t],pm]-sum(md.EeSUPPLY[t,pld,pm].value for pld in md.npload) for pm in md.npmarket} for t in md.ntime},
                                   'Hydropower production [GWh]'    :{t:{pm:sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload for hp in md.nhpp if md.hp_pmarket[hp]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Power transmission loss [GWh]'  :{t:{pm:sum(md.EeTRANS[t,pld,tl].value * md.eTransLoss[tl] for pld in md.npload for tl in md.ntransline if md.eTransIn[tl]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Power supply loss [GWh]'        :{t:{pm:sum(md.EeSUPPLY[t,pld,pm].value * md.eSupLoss[pm] / (1-md.eSupLoss[pm]) for pld in md.npload) for pm in md.npmarket} for t in md.ntime},
                                   'Net power export [GWh]'         :{t:{pm:sum(md.EeTRANS[t,pld,tl].value * (1-md.eTransLoss[tl]) for pld in md.npload for tl in md.ntransline if md.eTransIn[tl]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Net power import [GWh]'         :{t:{pm:sum(md.EeTRANS[t,pld,tl].value * (1-md.eTransLoss[tl]) for pld in md.npload for tl in md.ntransline if md.eTransOut[tl]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Energy price [$/kWh]'           :{t:{pm:sum(-md.dual[md.engy_balance[t,pld,pm]]*md.eLoadDem[pld] for pld in md.npload) for pm in md.npmarket} for t in md.ntime}}
        
        #Power Plants
        self.EnergyBalance_t_pm['Power plant prod [GWh]'] = {t:{pm:sum(md.EeOPPROD[t,pld,pp].value for pld in md.npload for pp in md.nopp if md.op_pmarket[pp]==pm) for pm in md.npmarket} for t in md.ntime}           
        self.PowerplantTable= {'Power Prod [GWh/y]'         : {pp:sum(md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for pp in md.nopp},
                               'Capacity [MW]'              : ReadParam('eOppCap'), 
                               'Fuel use [GWh/y]'           : {pp:sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp] for t in md.ntime for pld in md.npload)/leni(md.nyear) if md.Options['Fuels']==1 else 0 for pp in md.nopp},
                               'Operation costs [M$/y]'     : {pp:sum(md.eOppCost[pp] * md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for pp in md.nopp},
                               'Energy value [M$/y]'        : {pp:sum(-md.EeOPPROD[t,pld,pp].value * md.dual[md.engy_balance[t,pld,md.op_pmarket[pp]]] for t in md.ntime for pld in md.npload)/leni(md.nyear) for pp in md.nopp},
                               'Capacity shadowprice [$/kWh]':self.PpCapShadow_pp}
        self.PowerplantTable['Fuel cost [M$/y]']     = {pp:sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp] * md.eFuelCost[md.t_year[t],md.op_pmarket[pp],fu] for t in md.ntime for pld in md.npload for fu in md.nfuel if md.op_fuel[pp]==fu)/leni(md.nyear) for pp in md.nopp}
        self.PpProd_pp_m     = {m:{pp:sum(md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload  if md.t_month[t]==m)/leni(md.nyear) for pp in md.nopp} for m in md.nmonth}
        self.PpProd_pp_y     = {y:{pp:sum(md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload if md.t_year[t]==y) for pp in md.nopp} for y in md.nyear}
        self.PpBenefits_pp_y = {y:{pp:sum(md.EeOPPROD[t,pld,pp].value * (-md.dual[md.engy_balance[t,pld,md.op_pmarket[pp]]] - md.eOppCost[pp] - (md.Options['Fuels'] and md.eFuelCost[y,md.op_pmarket[pp],md.op_fuel[pp]]/md.eOppEff[pp])) for t in md.ntime for pld in md.npload if md.t_year[t]==y) for pp in md.nopp} for y in md.nyear}

        #Generic capacity
        self.EnergyBalance_t_pm['Gen power prod [GWh]'] = {t:{pm:sum(md.EeGENPROD[t,pld,pt,pm].value for pld in md.npload for pt in md.nptech) for pm in md.npmarket} for t in md.ntime}       
        self.EnergyBalance_t_pm['Gen power cap [MW]']   = {t:{pm:sum(md.EeGENCAP[md.t_year[t],pt,pm].value for pt in md.nptech)/leni(md.nmonth) for pm in md.npmarket} for t in md.ntime}  
        self.PowerGenCap_pt_pm  = {pt: {pm: sum(md.EeGENCAP[y,pt,pm].value for y in md.nyear) for pm in md.npmarket} for pt in md.nptech}
        self.PowerGenProd_pt_pm = {pt: {pm: sum(md.EeGENPROD[t,pld,pt,pm].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for pm in md.npmarket} for pt in md.nptech}
        self.PowerGenCap_pt_y   = {pt: {y:  sum(md.EeGENCAP[ky,pt,pm].value for ky in md.nyear for pm in md.npmarket if ky <= y) for y in md.nyear} for pt in md.nptech}
        self.PowerGenProd_pt_y  = {pt: {y:  sum(md.EeGENPROD[t,pld,pt,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear} for pt in md.nptech}

        #Transmission lines
        self.TmsLineTable  = {'Net transmission [GWh/y]'  : {tl:sum(md.EeTRANS[t,pld,tl].value * (1-md.eTransLoss[tl]) for t in md.ntime for pld in md.npload)/leni(md.nyear) for tl in md.ntransline},
                              'Power loss [GWh/y]'        : {tl:sum(md.EeTRANS[t,pld,tl].value * md.eTransLoss[tl] for t in md.ntime for pld in md.npload)/leni(md.nyear) for tl in md.ntransline},
                              'Cap. Shadow [$/GWh/y]'     : {tl:sum(-md.dual[md.engy_transmissioncap[t,pld,tl]]*md.eLoadTime[pld] for t in md.ntime for pld in md.npload)/leni(md.ntime) * leni(md.nmonth) for tl in md.ntransline}}             

        #Power markets and Energy balances
        self.EnergyBalance_t_pm['Balance = 0 [GWh]'] = {t:{pm:+ self.EnergyBalance_t_pm['Power demand [GWh]'][t][pm] 
                                                              - self.EnergyBalance_t_pm['Hydropower production [GWh]'][t][pm]
                                                              -(self.EnergyBalance_t_pm['Power plant prod [GWh]'][t][pm] if md.Options['Power plants'] == 1 else 0)
                                                              -(self.EnergyBalance_t_pm['Gen power prod [GWh]'][t][pm] if md.Options['Power technologies'] == 1 else 0)
                                                              - self.EnergyBalance_t_pm['Unserved demand [GWh]'][t][pm]
                                                              + self.EnergyBalance_t_pm['Power supply loss [GWh]'][t][pm]
                                                              + self.EnergyBalance_t_pm['Power transmission loss [GWh]'][t][pm]
                                                              + self.EnergyBalance_t_pm['Net power export [GWh]'][t][pm]
                                                              - self.EnergyBalance_t_pm['Net power import [GWh]'][t][pm]                                                                  
                                                              for pm in md.npmarket} for t in md.ntime}

        self.EnergyBalance      = {keys:sum(self.EnergyBalance_t_pm[keys][t][pm] for t in md.ntime for pm in md.npmarket)/leni(md.nyear) for keys in self.EnergyBalance_t_pm.keys()}
        self.EnergyBalance_co   = {keys:{co:sum(self.EnergyBalance_t_pm[keys][t][pm] for t in md.ntime for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry} for keys in self.EnergyBalance_t_pm.keys()}                       
        self.EnergyBalance_m    = {keys:{m:sum(self.EnergyBalance_t_pm[keys][t][pm] for t in md.ntime for pm in md.npmarket if md.t_month[t]==m)/leni(md.nyear) for m in md.nmonth} for keys in self.EnergyBalance_t_pm.keys()}    
        self.EnergyBalance_y    = {keys:{y:sum(self.EnergyBalance_t_pm[keys][t][pm] for t in md.ntime for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear} for keys in self.EnergyBalance_t_pm.keys()}
        #Correct intensive values (energy price)
        self.EnergyBalance['Energy price [$/kWh]']  = self.EnergyBalance['Energy price [$/kWh]']/leni(md.nmonth)/leni(md.npmarket)
        self.EnergyBalance_co['Energy price [$/kWh]']={co: self.EnergyBalance_co['Energy price [$/kWh]'][co]/leni(md.nmonth)/max(1,sum(1 for pm in md.npmarket if md.pmarket_country[pm]==co)) for co in md.ncountry}
        self.EnergyBalance_m['Energy price [$/kWh]']= {m : self.EnergyBalance_m['Energy price [$/kWh]'][m]/leni(md.npmarket) for m in md.nmonth}
        self.EnergyBalance_y['Energy price [$/kWh]']= {y : self.EnergyBalance_y['Energy price [$/kWh]'][y]/leni(md.npmarket)/leni(md.nmonth) for y in md.nyear}
        
        #Generic power tech
        for pt in md.nptech:
            self.EnergyBalance_y[pt+' production [GWh]']={y:self.PowerGenProd_pt_y[pt][y] for y in md.nyear}
        
        #Catchment scale energy balance
        self.EnergyBalance_ca={'Hydropower production [GWh]' : {c: sum(md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_catch[hp]==c)/leni(md.nyear) for c in md.ncatch}}

        #%%PRINTS
        if PRINT ==1:
            print('*************ENERGY*************')            
            print('Hydropower production [GWh]')
            print(self.EnergyBalance['Hydropower production [GWh]'])
            print('Power demand [GWh]')
            print(self.EnergyBalance['Power demand [GWh]'])
            print('Unserved power demand [GWh]')
            print(self.EnergyBalance['Unserved demand [GWh]'])
            print('Net power transmition [GWh]')
            print(self.EnergyBalance['Net power import [GWh]'])
            print('Power transmition losses [GWh]')
            print(self.EnergyBalance['Power transmission loss [GWh]'])
            print('Energy balance [GWh]: PowerDem = HpProd + PpProd - PowerTransLoss - UnsPowerDem')
            print(self.EnergyBalance['Balance = 0 [GWh]'])
            print('PowerPlant production [GWh]')
            print(self.EnergyBalance['Power plant prod [GWh]'])
            print('Generic capacity investment [kWh/month]')
            print(self.EnergyBalance['Gen power cap [MW]'])
            print('Generic capacity production [GWh]')
            print(self.EnergyBalance['Gen power prod [GWh]'])
                    
#%%ECONOMICS
        #Economic balance
        self.EconomicBalance_y_co = {'User benefits [M$]':       {y:{co:sum(md.wUserVal[u]*md.WwSUPPLY[t,u].value                       for t in md.ntime for u in md.nuser if md.user_country[u]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear},
                                'Water supply costs [M$]':       {y:{co:sum(md.WwSUPPLY[t,u].value/(1-md.wUserLoss[u]) * md.wSupCost[u] for t in md.ntime for u in md.nuser if md.user_country[u]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear},
                                'Crop cultivation costs [M$]':   {y:{co:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co) for co in md.ncountry} for y in md.nyear},
                                'Crop irrigation costs [M$]':    {y:{co:sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_country[fz]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear},
                                'Crop groundwater costs [M$]':   {y:{co:sum(md.AwGWSUPPLY[t,fz,cul].value * md.wGwCost[aq] for t in md.ntime for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]== md.aqui_catch[aq] and md.fzone_country[fz]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}}
        #Crop market
        self.EconomicBalance_y_co['Crop benefits [M$]']         = {y:{co:sum(md.AcSUPPLY[y,cm,cr,cds].value * md.aStepVal[cm,cr,cds]*md.aCropVal[cm,cr] for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co) for co in md.ncountry} for y in md.nyear}
        if md.Options['Crop market'] == 0:
            self.EconomicBalance_y_co['Crop benefits [M$]']     = {y:{co:sum(md.AcPROD[y,fz,cr].value * md.aFarmVal[md.fzone_country[fz],cr] for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Crop transport costs [M$]']  = {y:{co:sum(md.AcTRANS[y,ct,cr].value * md.aTransCost[ct,cr]                               for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Crop import costs [M$]']     = {y:{co:sum(md.AcTRANS[y,ct,cr].value * -md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Crop export benefits [M$]']  = {y:{co:sum(md.AcTRANS[y,ct,cr].value * -md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransIn[ct]]==co) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Crop ext prod costs [M$]']   = {y:{co:sum(md.AcEXTPROD[y,cm,cr].value * md.aCropVal[cm,cr] for cm in md.nextcmarket for cr in md.ncrop if md.cmarket_country[cm]==co) for co in md.ncountry} for y in md.nyear}
        
        #Energy markets       
        self.EconomicBalance_y_co['Energy benefits [M$]']       = {y:{co:sum(md.EeSUPPLY[t,pld,pm].value * md.eEngyVal[md.t_month[t],pm] for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        if md.Options['Energy market'] == 0:
            self.EconomicBalance_y_co['Energy benefits [M$]']   = {y:{co: sum(md.EeHPPROD[t,pld,hp].value*(md.eHppVal[hp]-md.eHppCost[hp]) for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_country[hp]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy import costs [M$]']   = {y:{co:sum(md.EeTRANS[t,pld,tl].value * -md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy export benefits [M$]']= {y:{co:sum(md.EeTRANS[t,pld,tl].value * -md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransIn[tl]]==co  and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy trans costs [M$]']    = {y:{co:sum(md.EeTRANS[t,pld,tl].value * md.eTransCost[tl] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy hydro O&M costs [M$]']= {y:{co:sum(md.EeHPPROD[t,pld,hp].value * md.eHppCost[hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.pmarket_country[md.hp_pmarket[hp]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy curtailment costs [M$]']={y:{co:sum(md.eEngyVal[md.t_month[t],pm]*(md.eEngyDem[md.t_year[t],md.t_month[t],pm]-sum(md.EeSUPPLY[t,pld,pm].value for pld in md.npload)) for t in md.ntime for pm in md.npmarket if md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        #Other power plants
        self.EconomicBalance_y_co['Energy power O&M costs [M$]']= {y:{co:sum(md.EeOPPROD[t,pld,pp].value*md.eOppCost[pp] for t in md.ntime for pld in md.npload for pp in md.nopp if md.pmarket_country[md.op_pmarket[pp]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy fuel costs [M$]']     = {y:{co:sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCost[md.t_year[t],md.op_pmarket[pp],fu] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        #Generic power technologies
        self.EconomicBalance_y_co['Energy Gen power O&M costs [M$]'] = {y:{co:sum(md.EeGENPROD[t,pld,pt,pm].value*md.eVarOPEX[pt,pm] for t in md.ntime for pld in md.npload for pm in md.npmarket for pt in md.nptech if md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy Gen cap costs [M$]']  = {y:{co:sum(md.EeGENCAP[y,pt,pm].value*md.eCAPEX[y,pt,pm]*(md.t_year[md.Options['tfin']]-y+1)/md.eLifeTime[pt,pm] + sum(md.EeGENCAP[ky,pt,pm].value for ky in md.nyear if ky <= y)*md.eFixOPEX[pt,pm] for pm in md.npmarket for pt in md.nptech if md.pmarket_country[pm]==co) for co in md.ncountry} for y in md.nyear}
        self.EconomicBalance_y_co['Energy Gen fuel costs [M$]'] = {y:{co:sum(md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCost[md.t_year[t],pm,fu] for t in md.ntime for pld in md.npload for pm in md.npmarket for pt in md.nptech for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}                   
        #CO2 emissions
        self.EconomicBalance_y_co['Energy CO2 costs [M$]']      = {y:{co: (sum(md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCO2[fu]*md.eCO2Val[md.t_year[t],pm] for t in md.ntime for pld in md.npload for pm in md.npmarket for pt in md.nptech for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co and md.t_year[t]==y))
                                                                         +(sum(md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCO2[fu]*md.eCO2Val[md.t_year[t],md.op_pmarket[pp]] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co and md.t_year[t]==y))
                                                                            for co in md.ncountry} for y in md.nyear}                
        #Balance
        self.EconomicBalance_y_co[' Balance [M$]'] = {y:{co:+ self.EconomicBalance_y_co['User benefits [M$]'][y][co]
                                                            - self.EconomicBalance_y_co['Water supply costs [M$]'][y][co]
                                                            + self.EconomicBalance_y_co['Crop benefits [M$]'][y][co]
                                                            - self.EconomicBalance_y_co['Crop cultivation costs [M$]'][y][co]
                                                            - self.EconomicBalance_y_co['Crop irrigation costs [M$]'][y][co]
                                                            - self.EconomicBalance_y_co['Crop groundwater costs [M$]'][y][co]       
                                                            - self.EconomicBalance_y_co['Crop transport costs [M$]'][y][co]       
                                                            - self.EconomicBalance_y_co['Crop import costs [M$]'][y][co]            
                                                            + self.EconomicBalance_y_co['Crop export benefits [M$]'][y][co]         
                                                            - self.EconomicBalance_y_co['Crop ext prod costs [M$]'][y][co]
                                                            + self.EconomicBalance_y_co['Energy benefits [M$]'][y][co]
                                                            - self.EconomicBalance_y_co['Energy import costs [M$]'][y][co]         
                                                            + self.EconomicBalance_y_co['Energy export benefits [M$]'][y][co]       
                                                            - self.EconomicBalance_y_co['Energy trans costs [M$]'][y][co]           
                                                            - self.EconomicBalance_y_co['Energy hydro O&M costs [M$]'][y][co]    
                                                            - self.EconomicBalance_y_co['Energy curtailment costs [M$]'][y][co]    
                                                            - self.EconomicBalance_y_co['Energy power O&M costs [M$]'][y][co]      
                                                            - self.EconomicBalance_y_co['Energy fuel costs [M$]'][y][co]            
                                                            - self.EconomicBalance_y_co['Energy Gen power O&M costs [M$]'][y][co]   
                                                            - self.EconomicBalance_y_co['Energy Gen cap costs [M$]'][y][co]       
                                                            - self.EconomicBalance_y_co['Energy Gen fuel costs [M$]'][y][co]      
                                                            - self.EconomicBalance_y_co['Energy CO2 costs [M$]'][y][co]               
                                                            for co in md.ncountry} for y in md.nyear}
        #Balance at other levels
        self.EconomicBalance      = {keys:     sum(self.EconomicBalance_y_co[keys][y][co] for y in md.nyear for co in md.ncountry)/leni(md.nyear) for keys in self.EconomicBalance_y_co.keys()}
        self.EconomicBalance_co   = {keys: {co:sum(self.EconomicBalance_y_co[keys][y][co] for y in md.nyear)/leni(md.nyear) for co in md.ncountry} for keys in self.EconomicBalance_y_co.keys()}
        self.EconomicBalance_y    = {keys: {y :sum(self.EconomicBalance_y_co[keys][y][co] for co in md.ncountry) for y in md.nyear} for keys in self.EconomicBalance_y_co.keys()}

        self.EconomicBalance_ca={'User benefits [M$]':           {c:sum(md.WwSUPPLY[t,u].value*md.wUserVal[u] for t in md.ntime for u in md.nuser if md.user_catch[u]==c)/leni(md.nyear) for c in md.ncatch},
                                 'Water supply costs [M$]':      {c:sum(md.WwSUPPLY[t,u].value/(1-md.wUserLoss[u]) * md.wSupCost[u] for t in md.ntime for u in md.nuser if md.user_catch[u]==c)/leni(md.nyear) for c in md.ncatch},                                  
                                 'Crop cultivation costs [M$]':  {c:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_catch[fz]==c)/leni(md.nyear)     for c in md.ncatch}, 
                                 'Crop irrigation costs [M$]':   {c:sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c)/leni(md.nyear)      for c in md.ncatch}}
        if md.Options['Crop market'] == 0:
            self.EconomicBalance_ca['Crop benefits [M$]']           = {c:sum(md.AcPROD[y,fz,cr].value * md.aFarmVal[md.fzone_country[fz],cr] for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_catch[fz]==c)/leni(md.nyear) for c in md.ncatch}

        self.EconomicBalance_ca['Crop groundwater costs [M$]']      = {c:sum(md.AwGWSUPPLY[t,fz,cul].value * md.wGwCost[md.fzone_catch[fz]]  for t in md.ntime for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]==md.aqui_catch[aq])/leni(md.nyear) for c in md.ncatch}
        if md.Options['Energy market'] == 0:
            self.EconomicBalance_ca['Energy benefits [M$]']         = {c:sum(md.EeHPPROD[t,pld,hp].value*(md.eHppVal[hp]-md.eHppCost[hp]) for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_catch[hp]==c)/leni(md.nyear) for c in md.ncatch} 
        #%%PRINTS    
        if PRINT==1:
            print('*************ECONOMIC*************')
            print('User benefit [M$]')
            print(self.EconomicBalance['User benefits [M$]'])
            print('Crop benefit [M$]')
            print(self.EconomicBalance['Crop benefits [M$]'])
            print('Cultivation cost [M$]')
            print(self.EconomicBalance['Crop cultivation costs [M$]'])
            print('Irrigation cost [M$]')
            print(self.EconomicBalance['Crop irrigation costs [M$]'])
            print('Crop import cost [M$]')
            print(self.EconomicBalance['Crop import costs [M$]'])
            print('Crop export benefits [M$]')
            print(self.EconomicBalance['Crop export benefits [M$]'])
            print('Crop transport cost [M$]')
            print(self.EconomicBalance['Crop transport costs [M$]'])
            print('Energy benefit [M$]')
            print(self.EconomicBalance['Energy benefits [M$]'])            
            print('Energy curtailment cost [M$]')
            print(self.EconomicBalance['Energy curtailment costs [M$]'])
            print('Hydropower operational cost [M$]')
            print(self.EconomicBalance['Energy hydro O&M costs [M$]'])
            print('Other Power plants operational cost [M$]')
            print(self.EconomicBalance['Energy power O&M costs [M$]'])
            print('Other Power plants fuel cost [M$]')
            print(self.EconomicBalance['Energy fuel costs [M$]'])
            print('Generic capacity construction costs [M$]')
            print(self.EconomicBalance['Energy Gen cap costs [M$]'])
            print('Generic capacity operation costs [M$]')
            print(self.EconomicBalance['Energy Gen power O&M costs [M$]'])
            print('Generic fuel cost [M$]')
            print(self.EconomicBalance['Energy Gen fuel costs [M$]'])
                
#%% VALIDATION - #ADD THIS SECTION IS SPECIFIC TO ZAMBEZI STUDY CASE !!!!! (SO FAR)
        if VALIDATION == 1:
            self.WaterAbs_fz    = {'Sim Consumption [Mm3/y]':{fz:self.FarmZonesTable['Irrig net abstraction [Mm3/y]'][fz] for fz in md.nifzone},
                                  'Obs Consumption [Mm3/y]':{fz:ReadParam('ObsWaterCons')[fz] for fz in md.nifzone}}
            self.WaterAbs_co    = {'Sim Consumption [Mm3/y]':{co:sum(self.FarmZonesTable['Irrig net abstraction [Mm3/y]'][fz] for fz in md.nifzone if md.fzone_country[fz]==co) for co in md.ncountry},
                                  'Obs Consumption [Mm3/y]':ReadParam('vAgrWaterCons')}
            self.Hydropower_hp  = {'Hp Prod sim [GWh/y]':self.HydropowerTable['Power Prod [GWh/y]'],
                                  'Hp Prod obs [GWh/y]':ReadParam('vObservedProd')}
            self.ObsCropProd_cr_co = {cr:{co:ReadParam('vCropProd')[co,cr]*ReadParam('vScaleFactor')[co] for co in md.ncountry} for cr in md.ncrop}
            

            self.CropValue_co  = {'Sim':{co:self.EconomicBalance_co['Crop benefits [M$]'][co]+self.EconomicBalance_co['Crop export benefits [M$]'][co]-self.EconomicBalance_co['Crop import costs [M$]'][co]-self.EconomicBalance_co['Crop transport costs [M$]'][co] for co in md.ncountry},
								  'Sim2':{co:sum(md.AcPROD[y,fz,cr].value * (md.aCropVal[md.fzone_cmarket[fz],cr] if md.Options['Crop market']==1 else md.aFarmVal[co,cr]) for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry},
                                  'Obs':{co:ReadParam('vProdValue')[co]*ReadParam('vScaleFactor')[co] for co in md.ncountry if co in ReadParam('vProdValue').keys()}}

            if md.Options['Energy market'] == 1:
                self.EnergyImExp_co= {'EngyExp sim [GWh/y]':self.EnergyBalance_co['Net power export [GWh]'],
                                      'EngyImp sim [GWh/y]':self.EnergyBalance_co['Net power import [GWh]'],
                                      'EngyExp obs [GWh/y]':ReadParam('PowerExp'),
                                      'EngyImp obs [GWh/y]':ReadParam('PowerImp')}
            else:
               self.EnergyImExp_co = {'No Data':{'No Data':0}} 
            
#%%EXPORT RESULT FUNCTIONS
    
    def exportsheet(self,writer,sheet,data,title,dataref=[],index=[],order={},total={},color={},cols=1,srowi=5,scoli=0,distance=4):    
        scol=scoli
        srow=srowi
        rowjump=0
        for k in range(len(data)):
            #Title of the data                                    
            ptitle = pd.DataFrame([title[k]]) 
            
            #To panda dataframes
            if k in index: #data is a single dictionary and needs orient='index'
                pdata    = pd.DataFrame.from_dict(data[k], orient='index')
                if dataref != []: #Substract reference data (dataref) 
                    pdataref = pd.DataFrame.from_dict(dataref[k], orient='index')
                    for obj in pdata.index:
                        if obj in pdataref.index: #only compare if object is also in reference results
                            pdata.loc[obj]=pdata.loc[obj]-pdataref.loc[obj]
            else:
                pdata    = pd.DataFrame.from_dict(data[k])
                if dataref != []: #Substract reference data (dataref) 
                    pdataref = pd.DataFrame.from_dict(dataref[k])
                    for obj in pdata.index:
                        if obj in pdataref.index: #only compare if object is also in reference results
                            pdata.loc[obj]=pdata.loc[obj]-pdataref.loc[obj]  #REM: supposes objects/projects are as rows                     
            
            #Reorder columns of dataframe
            if k in order.keys() and not pdata.empty: 
                pdata = pdata[order[k]]
            
            #Compute total 0=total of columns, 1=total of lines, 2=both            
            if k in total.keys() and not pdata.empty: 
                if total[k]==0:
                    pdata.loc['Total']=pdata.sum()
                elif total[k]==1:
                    pdata['Total']=pdata.sum(axis=1)
                elif total[k]==2:
                    pdata.loc['Total']=pdata.sum()
                    pdata['Total']=pdata.sum(axis=1)
                else:
                    print('WARNING: argument in "total" not valid')
            
            #Export to excel
            pdata.to_excel(writer, sheet_name=sheet, startrow=srow, startcol=scol)
            ptitle.to_excel(writer, sheet_name=sheet, index=False, header=False, startrow=srow, startcol=scol)
            
            #Move columns and lines for next data
            if (k+1)%cols == 0:
                scol = scoli
                srow= srow + max(rowjump,len(pdata)) + distance
                rowjump=0
            else:
                scol = scol + pdata.shape[1] + distance
                rowjump=max(rowjump,len(pdata))

    def exportresults(self,outpath,md,VALIDATION=1,NEWSHEET=0):
        months=['oct','nov','dec','jan','feb','mar','apr','may','jun','jul','aug','sep']        
        farmingzone=['Crop Production [Mt/y]','Production Value [M$/y]','Cultivation costs [M$/y]','Irrigation costs [M$/y]','Cultivated land Net [Mha]','Cultivated land Net [%]','Cultivated land Gross [%]','Irrig net abstraction [Mm3/y]','Irrig withdrawals [Mm3/y]','Irrig losses [Mm3/y]','Irrig return flow [Mm3/y]','Crop dem. satisfaction [%]','Land shadowprice [$/ha]']
        #sheets=['Investments','PowerNetwork','PowerPlants','Energy','Crops','FarmingZones','Reservoirs','Hydropowers','EnvFlow','Users','MassBalance','Validation']

        writer = pd.ExcelWriter(outpath, engine='openpyxl') #able to save formulas (graphs die however)
        if NEWSHEET==0:
            book = load_workbook(outpath)
            writer.book = book
            writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        
        #POWER NETWORK
        if md.Options['Energy market'] == 1:
            title=['Transmission lines']
            data=[self.TmsLineTable]
            self.exportsheet(writer,'PowerNetwork',data,title)
            
        #POWER PLANTS
        if md.Options['Power plants'] == 1:
            data=[self.PowerplantTable,self.PpProd_pp_m,self.PpProd_pp_y,self.PpBenefits_pp_y]
            title=['Power plant','Power production [kWh/month]','Power production [GWh/y]','Net Benefits* [M$/y]']
            self.exportsheet(writer,'PowerPlants',data,title,order={1:months})
            
        #ENERGY
        if md.Options['Energy market'] == 1:
            #ADD: if power tech not on some variables not defined
            data=[self.EnergyBalance,self.EnergyBalance_co,self.EnergyBalance_ca,self.EnergyBalance_m,self.EnergyBalance_y,self.EnergyVal_pm_m, self.PowerGenCap_pt_pm, self.PowerGenProd_pt_pm, self.PowerGenCap_pt_y, self.PowerGenProd_pt_y]
            title=['ENERGY BALANCE' ,'ENERGY BALANCE'     ,'ENERGY BALANCE'     ,'ENERGY BALANCE'    ,'ENERGY BALANCE'    ,'EneryShadowprice [$/kWh]', 'Generic capacity invest [MW]', 'Generic capacity prod [GWh/year]','Generic capacity invest [MW]', 'Generic capacity prod [GWh/year]']
            self.exportsheet(writer,'Energy',data,title,index=[0,1,2,3,4,5],order={3:months,5:months},total={3:1})
        
        #CROPS
        data=[self.CropBalance_y, self.CropPrice_cr_y,self.CropPrice_cr_cm,self.CropProduction_cr_co,self.LandUse_cr_co,self.CropImpShare_cr_co, self.CropExtProd_cr_co, self.LandUse_cr_y,self.LandUse_co_ft,self.LandUse_cr_ft,self.CropTransport_cm_cm, self.CropLoss_cm_cm, self.CropTransport_ct_cr]
        title=['Crop Balance','Crop Price [$/t]','Crop Price [$/t]','Crop Production [Mt]','Land Use [Mha]','Crop Import Share [%ofDemand]','Crop external production [Mt]','Land Use [Mha]','Land Use [Mha]','Land Use [Mha]','Crop Transport [Mt/y]', 'Crop Transport losses[Mt/y]', 'Crop transport [Mt/y]']
        self.exportsheet(writer,'Crops',data,title)
        #LAND USE
        #ADD PUT ALL LAND USE OUTPUTS HERE
        
        #FARMING ZONES
        data=[self.FarmZonesTable,self.LandUse_cul_fz,self.LandUse_fd_ft,self.Yield_cul_ft,self.CropProduction_cr_fz,self.FarmBenefits_ca_ft,self.FarmBenefits_fz_y]
        title=['Farming zone','Land Use [ha/y]','Land Use [ha/y]','Yield [%]','Crop production [Mt/y]','Net benefits [M$/y]','Net benefits [M$/y]']
        self.exportsheet(writer,'FarmingZones',data,title,order={0:farmingzone})        
        
        #RESERVOIRS
        data=[self.ReservoirTable,self.ResStor_res_m,self.ResStor_res_t]
        title=['Reservoir','Average storage [Mm3]','Reservoir storage [Mm3]']
        self.exportsheet(writer,'Reservoirs',data,title,order={1:months})
        
        #GROUNDWATER
        if md.Options['Groundwater']==1:
            self.exportsheet(writer,'Groundwater',[self.GwStor_c_m,self.GwStor_c_t],['Groundwater storage [Mm3]','Groundwater storage [Mm3]'],order={0:months})
        
        #HYDROPOWER
        data=[self.HydropowerTable,self.HpProd_hp_m,self.HpStor_hp_m,self.HpDis_hp_m,self.HpProd_hp_y,self.HpBenefits_hp_y]
        title=['Hydropower','Power production [GWh]','Storage* [Mm3]','Discharge [Mm3]','Power production [GWh/y]', 'Net Benefits* [M$/year]']
        self.exportsheet(writer,'Hydropower',data,title,order={1:months,2:months,3:months},total={1:2,3:2,4:2,5:2})
        
        #ENVIRONMENTAL FLOWS
        if md.Options['Eflows'] == 1:
            self.exportsheet(writer,'EnvFlow',[self.EnvFlowShadow_m_ec],['e-flow shadowprice [$/m3]'],order={0:months})
        
        #USERS
        self.exportsheet(writer,'Users',[self.UserTable],['Water Users'])
        
        #MASS BALANCES
        data=[self.EconomicBalance,self.EconomicBalance_co,self.EconomicBalance_ca,
              self.WaterBalance,self.WaterBalance_co,self.WaterBalance_ca,
              self.EnergyBalance,self.EnergyBalance_co,self.EnergyBalance_ca,
              self.CropBalance,self.CropBalance_co,self.CropBalance_ca]
        title=['ECONOMIC BALANCE','ECONOMIC BALANCE','ECONOMIC BALANCE',
               'WATER BALANCE','WATER BALANCE','WATER BALANCE',
               'ENERGY BALANCE','ENERGY BALANCE','ENERGY BALANCE',
               'CROP BALANCE','CROP BALANCE','CROP BALANCE']
        self.exportsheet(writer,'MassBalance',data,title,index=[0,1,2,3,4,5,6,7,8,9,10,11],cols=3)
        
        #MAIN
        DebugWater  = {'Debug is OFF':{0:'Debug is off'}}
        DebugCrop   = {'Debug is OFF':{0:'Debug is off'}}
        DebugEflow  = {'Debug is OFF':{0:'Debug is off'}}
        if md.Options['Debug mode'] == 1:
            DebugWater = {t:{c:md.DEBUG[t,c].value for c in md.ncatch} for t in md.ntime}
            DebugCrop  = {y:{fz:sum(md.DUMYCROP[y,fz,cr].value for cr in md.ncrop) for fz in md.nfzone} for y in md.nyear}
            DebugEflow = {t:{ef:md.DUMYEFLOW[t,ef].value for ef in md.neflow} for t in md.ntime}
        self.exportsheet(writer,'Main',[self.MainTable,md.Options,DebugWater,DebugCrop,DebugEflow],['Configurations','Options','DebugWater','DebugCrop','DebugEflow'],index=[0,1,2])
        
        #VALIDATION
        if VALIDATION == 1:
            data=[self.WaterAbs_fz, self.WaterAbs_co, self.CropProduction_cr_co, self.ObsCropProd_cr_co, self.CropValue_co, self.Hydropower_hp, self.EnergyImExp_co]
            title=['Water Abstraction', 'Water Abstraction', 'Sim Crop production [Mt/y]', 'Obs Crop production [Mt/y]', 'Value of crop Prod [M$/y]', 'Hydropower', 'Energy']
            self.exportsheet(writer,'Validation',data,title,total={0:0,1:0,2:0,3:0,4:0,5:0,6:0})
        
        #Save excel files
        writer.save()