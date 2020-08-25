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
import os
import pandas as pd
import time
import pickle
from openpyxl         import load_workbook
from pyomo.environ    import Constraint, Var

class ResultAnalysis():
    #%% Save all model parameters and solved decision variables
    def export_all_DV(self,md,expfolder,scenario='WHATIF_main'):
        #scaling factors (descale variables)
        def descaling(varname):
            if varname in ['WwRSTORAGE','WwOUTFLOW']:
                return md.sw.value
            elif varname in ['EeGENCAP','EeGENPROD', 'EeHPPROD', 'EeOPPROD', 'EeSUPPLY','EeTRANS']:
                return md.se.value
            else:
                return 1
        All_DV={str(var):{index:var[index].value*descaling(str(var)) for index in var} for var in md.component_objects(Var)} #Save all decision variables to dictionary
        if md.Options['Yield water response']=='linearized':
            All_DV['AlCULAREA']={[y,fz,fd,fc]:sum(md.AlCULAREA[y,fz,fd,fc,yp].value for yp in md.nypath) for y,fz,fd,fc, in md.AlCULAREA}
        #Constraints
        constrlist = [const.name for const in md.component_objects(Constraint)]
        All_DV['energy_shadow'] = {index:md.dual[md.engy_balance[index]] for index in md.engy_balance if 'engy_balance' in constrlist}
        All_DV['crop_shadow']   = {index:md.dual[md.agr_cropbalance[index]] for index in md.agr_cropbalance if 'agr_cropbalance' in constrlist}
        All_DV['water_shadow']  = {index:md.dual[md.water_waterbalance[index]] for index in md.water_waterbalance if 'water_waterbalance' in constrlist}
        if md.find_component('CULAREA') is not None:
            All_DV['CULAREA']       = {index:md.CULAREA[index].value for index in md.CULAREA} 
        #oo['parameters']=parameters
        #oo['scenario']=scenario
        expfile = expfolder + os.sep + scenario + '_DV.txt' 
        pickle.dump(All_DV,open(expfile,"wb"))
    #%% Read only Selected result for scenario analysis (WHATIF_scenario.py)
    def selectedresults(self,md,parameters,scenario='WHATIF_main'):
        #output dictionary
        oo={}
        #NILE CASE
        def TurbFlowVol(Vol,MOL,FSL,qmin,qmax):
            return qmin + (qmax-qmin)*(Vol-MOL)/(FSL-MOL)
        def RelHeadVol(Vol,MOL,FSL,hmin,hmax):                
            return hmin/hmax + (1-hmin/hmax)*(Vol-MOL)/(FSL-MOL)
        def _storage(t,r,mode='average'): #Average or relative storage from t-1 to t
            S_end=md.WwRSTORAGE[t,r].value
            S_ini=md.wStorIni[r].value if md.Options['Initial time step'] == 1 and md.Options['tini'] == t else md.WwRSTORAGE[md.t_prev[t],r].value
            return (S_end+S_ini)/2 if mode=='average' else S_end-S_ini
        oo['Runoff'] = {t: sum(md.wRunOff[t,c].value for c in md.ncatch) for t in md.ntime}
        oo['Reservoir storage'] = {t: sum(md.WwRSTORAGE[t,r].value for r in md.nres) for t in md.ntime}
        oo['Shadow value'] = {t: -sum(md.dual[md.water_waterbalance[t,c]] for c in md.ncatch) for t in md.ntime}
        oo['Turbine capacity'] = {t: sum(TurbFlowVol(_storage(t,md.hp_res[hp]),md.wResMOL[md.hp_res[hp]],md.wResFSL[md.hp_res[hp]],md.wMinTurb[hp],md.wMaxTurb[hp])*md.m3pers_to_Mm3perMonth for hp in md.nhpp) if md.Options['Hydropower']=='nonlinear' else 0 for t in md.ntime}
        oo['Head factor'] = {t: sum(md.eHppEff[hp]*md.eHppProd[hp]*RelHeadVol(_storage(t,md.hp_res[hp]),md.wResMOL[md.hp_res[hp]],md.wResFSL[md.hp_res[hp]],md.wMinHead[md.hp_res[hp]],md.wMaxHead[md.hp_res[hp]]) for hp in md.nhpp) if md.Options['Hydropower']=='nonlinear' else 0 for t in md.ntime}
        oo['Hydropower production'] = {t: sum(md.EeHPPROD[t,pld,hp].value for hp in md.nhpp for pld in md.npload) for t in md.ntime}
        oo['zHydropower firm'] = {t: sum(md.EeFIRMHP[md.t_year[t],hp].value for hp in md.nhpp) for t in md.ntime}
        oo['zHydropower sec'] = {t: oo['Hydropower production'][t]-oo['zHydropower firm'][t] for t in md.ntime}
        oo['Hydropower value'] = {t: sum(md.EeHPPROD[t,pld,hp].value*md.eHppVal[hp] for hp in md.nhpp for pld in md.npload)
                                    +sum(md.EeFIRMHP[md.t_year[t],hp].value*(md.eHppFirm[hp]-md.eHppVal[hp]) for hp in md.nhpp) for t in md.ntime}
        oo['Hydropower spill'] = {t: sum(md.WwOUTFLOW[t,md.res_catch[md.hp_res[hp]]].value-sum(md.EwHPDISCHARGE[t,pld,hp].value for pld in md.npload) for hp in md.nhpp) for t in md.ntime}
        oo['Hydropower spill val'] = {t: oo['Hydropower spill'][t]*oo['Head factor'][t]*sum(md.eHppVal[hp] for hp in md.nhpp)/len(md.nhpp) for t in md.ntime} #nile specific
        oo['Hydropower discharge'] = {t: sum(md.EwHPDISCHARGE[t,pld,hp].value for pld in md.npload for hp in md.nhpp) for t in md.ntime}
        oo['User curt.'] = {t: sum(md.wUserDem[md.t_year[t],md.t_month[t],u] for u in md.nuser)-sum(md.WwSUPPLY[t,u].value for u in md.nuser) for t in md.ntime}
        oo['User value'] = {t: sum(md.WwSUPPLY[t,u].value*md.wUserVal[u] for u in md.nuser) for t in md.ntime}
        oo['Benefits']={t: oo['User value'][t]+oo['Hydropower value'][t] for t in md.ntime}
        
        oo['yRunoff'] = {y: sum(oo['Runoff'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yReservoir storage'] = {y: sum(oo['Reservoir storage'][t] for t in md.ntime if md.t_year[t]==y and md.t_month[t]=='jul') for y in md.nyear}
        oo['yShadow value'] = {y: sum(oo['Shadow value'][t] for t in md.ntime if md.t_year[t]==y)/12 for y in md.nyear}
        oo['yTurbine capacity'] = {y: sum(oo['Turbine capacity'][t] for t in md.ntime if md.t_year[t]==y)/12 for y in md.nyear}
        oo['yHead factor'] = {y: sum(oo['Head factor'][t] for t in md.ntime if md.t_year[t]==y)/12 for y in md.nyear}
        oo['yHydropower production'] = {y: sum(oo['Hydropower production'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yzHydropower firm'] = {y: sum(oo['zHydropower firm'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yzHydropower sec'] = {y: sum(oo['zHydropower sec'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yHydropower value'] = {y: sum(oo['Hydropower value'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yHydropower spill'] = {y: sum(oo['Hydropower spill'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yHydropower spill val'] = {y: sum(oo['Hydropower spill val'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yHydropower discharge'] = {y: sum(oo['Hydropower discharge'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yUser curt.'] = {y: sum(oo['User curt.'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yUser value'] = {y: sum(oo['User value'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        oo['yBenefits']={y: sum(oo['Benefits'][t] for t in md.ntime if md.t_year[t]==y) for y in md.nyear}
        
        oo['SUM_NILE']={
        #'tRunoff': sum(oo['yRunoff'][y] for y in md.nyear)/len(md.nyear),
        'tReservoir storage': sum(oo['yReservoir storage'][y] for y in md.nyear)/len(md.nyear), 
        'tShadow value': sum(oo['yShadow value'][y] for y in md.nyear)/len(md.nyear),
        #'tTurbine capacity': sum(oo['yTurbine capacity'][y] for y in md.nyear)/len(md.nyear),
        'tHead factor': sum(oo['yHead factor'][y] for y in md.nyear)/len(md.nyear),
        'tHydropower production':  sum(oo['yHydropower production'][y] for y in md.nyear)/len(md.nyear), 
        #'tzHydropower firm':  sum(oo['yzHydropower firm'][y] for y in md.nyear)/len(md.nyear), 
        #'tzHydropower sec':  sum(oo['yzHydropower sec'][y] for y in md.nyear)/len(md.nyear), 
        'tHydropower value': sum(oo['yHydropower value'][y] for y in md.nyear)/len(md.nyear),
        'tHydropower spill': sum(oo['yHydropower spill'][y] for y in md.nyear)/len(md.nyear), 
        'tHydropower spill val': sum(oo['yHydropower spill val'][y] for y in md.nyear)/len(md.nyear), 
        #'tHydropower discharge':sum(oo['yHydropower discharge'][y] for y in md.nyear)/len(md.nyear),
        'tUser curt.': sum(oo['yUser curt.'][y] for y in md.nyear)/len(md.nyear), 
        'tUser value': sum(oo['yUser value'][y] for y in md.nyear)/len(md.nyear), 
        'tUser av value': sum(oo['yUser value'][y] for y in md.nyear)/sum(md.WwSUPPLY[t,u].value for t in md.ntime for u in md.nuser), 
        'tBenefits': sum(oo['yBenefits'][y] for y in md.nyear)/len(md.nyear),
        'FinalStorage':oo['Reservoir storage'][md.Options['tfin']],
        'tBenefitsm2':sum(oo['yBenefits'][y] for y in md.nyear if y<=max(md.nyear)-2)/(len(md.nyear)-2),
        }
#        #Remove indicators
        if parameters.val['SCENITER']==1:
            del oo['Runoff']
            del oo['Reservoir storage'] 
            del oo['Shadow value'] 
            del oo['Turbine capacity'] 
            del oo['Head factor'] 
            del oo['Hydropower production']
            del oo['zHydropower firm'] 
            del oo['zHydropower sec'] 
            del oo['Hydropower value'] 
            del oo['Hydropower spill'] 
            del oo['Hydropower spill val']
            del oo['Hydropower discharge'] 
            del oo['User curt.'] 
            del oo['User value']
            del oo['Benefits']
            #remove yearly values
            del oo['yRunoff']
            del oo['yReservoir storage'] 
            del oo['yShadow value'] 
            del oo['yTurbine capacity'] 
            del oo['yHead factor'] 
            del oo['yHydropower production']
            del oo['yzHydropower firm'] 
            del oo['yzHydropower sec'] 
            del oo['yHydropower value'] 
            del oo['yHydropower spill'] 
            del oo['yHydropower spill val']
            del oo['yHydropower discharge'] 
            del oo['yUser curt.'] 
            del oo['yUser value']
            del oo['yBenefits']
        if parameters.val['SCENITER']!=1:
            #Length of index function
            def leni(index):
                return max(1,len(index))
            ##SHADOWPRICES##
            #Water shadowprice
            WaterValueShadow        = -sum(md.dual[md.water_waterbalance[t,c]] for t in md.ntime for c in md.ncatch)/(leni(md.ntime)*leni(md.ncatch))
            #Energy shadowprice
            if md.Options['Energy market'] == 1: 
                oo['EnergyValue_co'] = {co:-sum(1/md.se.value*md.dual[md.engy_balance[t,pld,pm]]*md.se.value*md.EeSUPPLY[t,pld,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co)/sum(md.se.value*md.EeSUPPLY[t,pld,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co) for co in md.ncountry if sum(md.se.value*md.EeSUPPLY[t,pld,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co) != 0}
            else:
                oo['EnergyValue_co'] = {co:-sum(md.eHppVal[hp] for hp in md.nhpp if md.hp_country[hp]==co)/sum(1 for hp in md.nhpp if md.hp_country[hp]==co) for co in md.ncountry if sum(1 for hp in md.nhpp if md.hp_country[hp]==co) != 0}
            #Crop shadowprice
            if md.Options['Crop market'] == 1:
                oo['CropValue_co']   = {co:-sum(1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,cm,cr]]*sum(md.AcSUPPLY[y,cm,cr,cds].value for cds in md.ncdstep) for y in md.nyear for cr in md.ncrop for cm in md.ncmarket if md.cmarket_country[cm]==co)/sum(md.AcSUPPLY[y,cm,cr,cds].value for y in md.nyear for cr in md.ncrop for cm in md.ncmarket for cds in md.ncdstep if md.cmarket_country[cm]==co) for co in md.ncountry if sum(md.AcSUPPLY[y,cm,cr,cds].value for y in md.nyear for cr in md.ncrop for cm in md.ncmarket for cds in md.ncdstep if md.cmarket_country[cm]==co) != 0}      
            else:
                oo['CropValue_co']   = {co: sum(md.aFarmVal[co,cr] for cr in md.ncrop)/leni(md.ncrop) for co in md.ncountry}
            ##WATER##
            if md.Options['Initial time step'] == 0: #Cyclic Run
                ReservoirEvaporation = sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[t,res].value+md.sw.value*md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime for res in md.nres)/leni(md.nyear)
            else: #Initial time step
                ReservoirEvaporation = (sum(max(0,md.wET0[md.Options['tini'],md.res_catch[res]].value-1*md.wRainFall[md.Options['tini'],md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[md.Options['tini'],res].value + md.wStorIni[res].value)/2 + md.wResArea[res]) for res in md.nres )
                                     +  sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[t,res].value + md.sw.value*md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime for res in md.nres if t != md.Options['tini'])
                                     /  leni(md.nyear))
            #Cultivated area
            if md.Options['Crop choice'] == 'fixed':
                GROSSAREA   = {(y,fz,cul): md.kha_to_Mha*md.CULAREA[y,fz,cul].value for y in md.nyear for fz in md.nfzone for cul in md.nculture}
                NETAREA     = {(y,fz): sum(md.kha_to_Mha*md.CULAREA[y,fz,cul].value for cul in md.nculture) for fz in md.nfzone for y in md.nyear}
            elif md.Options['Yield water response'] == 'nonlinear':
                GROSSAREA   = {(y,fz,cul):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd].value for fd in md.nfieldculture if md.field_culture[fd]==cul) for cul in md.nculture for fz in md.nfzone for y in md.nyear}
                NETAREA     = {(y,fz):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd].value for fd in md.nfieldculture if fd[1]==1) for fz in md.nfzone for y in md.nyear}
            else:
                GROSSAREA   = {(y,fz,cul):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if md.field_culture[fd]==cul) for cul in md.nculture for fz in md.nfzone for y in md.nyear}
                NETAREA     = {(y,fz):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if fd[1]==1) for fz in md.nfzone for y in md.nyear}
            if md.Options['Crop market'] == 1:
                CropPrice   = {(y,fz,cr):-1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.fzone_cmarket[fz],cr]] for cr in md.ncrop for fz in md.nfzone for y in md.nyear}
            else:
                CropPrice   = {(y,fz,cr): md.aFarmVal[md.fzone_country[fz],cr] for cr in md.ncrop for fz in md.nfzone for y in md.nyear}
            ##INVESTMENT DECISIONS##
            if md.Options['Investment module'] == 1:
                oo['InvestmentPlan_inv_ip']  = {inv:{ip:md.IbINVEST[ip,inv].value for ip in md.ninvphase} for inv in md.ninvest}
                oo['InvestmentDecision_inv'] = {inv:sum(md.IbINVEST[ip,inv].value for ip in md.ninvphase) for inv in md.ninvest}
                oo['InvestmentPhase_inv']    = {inv:[ip for ip in md.ninvphase if round(md.IbINVEST[ip,inv].value)==1][0] if round(oo['InvestmentDecision_inv'][inv])==1 else 'no' for inv in md.ninvest}
            
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
            CropTransportCost_co    = {co:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * md.aTransCost[ct,cr] for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co)/leni(md.nyear) for co in md.ncountry}
            #Crop import/export benefit and costs       
            CropExtProdCost_co      = {co:sum(md.kt_to_Mt*md.AcEXTPROD[y,cm,cr].value * md.aCropVal[cm,cr] for y in md.nyear for cm in md.nextcmarket for cr in md.ncrop if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry}  #The assumption is that external market produce crops at their market price #ADD: crop supply curve
            oo['CropExpBenefit_co'] = {co:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransIn[ct]]==co)/leni(md.nyear) for co in md.ncountry}
            CropImpCost_co          = {co:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co)/leni(md.nyear) for co in md.ncountry}          
            #Crop benefits
            CropBenefit_co          = {co:sum(md.kt_to_Mt*md.AcSUPPLY[y,cm,cr,cds].value * md.aStepVal[cm,cr,cds]*md.aCropVal[cm,cr] for y in md.nyear for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry}
            if md.Options['Crop market'] == 0:
                CropBenefit_co      = {co:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * md.aFarmVal[md.fzone_country[fz],cr] for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
            #Crop production NET benefit    
            FarmBenefits_fz_y       = {fz:{y: sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for cr in md.ncrop)
                                                - sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for cul in md.nculture)
                                                - sum(md.AwSUPPLY[t,fz,cul].value/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz] for t in md.ntime for cul in md.nculture if md.t_year[t]==y and md.aIrrigation[md.fzone_type[fz]]==1)
                                                for y in md.nyear} 
                                           for fz in md.nfzone}
            #Consummer and Producer surplus
            oo['AgricultureConsSurplus_co']  = {co:sum(md.kt_to_Mt*md.AcSUPPLY[y,cm,cr,cds].value * (md.aStepVal[cm,cr,cds]*md.aCropVal[cm,cr] - -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,cm,cr]])
                                          for y in md.nyear for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co)/leni(md.nyear) 
                                       for co in md.ncountry}
            oo['AgricultureProdSurplus_co']  = {co:sum(sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value for fz in md.nfzone if md.fzone_cmarket[fz]==cm) * -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,cm,cr]] for y in md.nyear for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co)/leni(md.nyear) 
                                        - sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) 
                                        - sum(md.AwSUPPLY[t,fz,cul].value/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz]  for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) 
                                       for co in md.ncountry}    
        #Balances
            CropBalance_co          = {co:CropBenefit_co[co] - CropCultivationCost_co[co] - CropIrrigationCost_co[co] + oo['CropExpBenefit_co'][co] - CropExtProdCost_co[co] - CropImpCost_co[co] - CropTransportCost_co[co] for co in md.ncountry}
            CropBalance2_co         = {co:sum(FarmBenefits_fz_y[fz][y] for fz in md.nfzone for y in md.nyear if md.fzone_country[fz]==co) for co in md.ncountry}
            
        #Energy benefits        
            EnergyBenefit_co    = {co:sum(md.se.value*md.EeSUPPLY[t,pld,pm].value * md.eEngyVal[md.t_month[t],pm] for t in md.ntime for pm in md.npmarket for pld in md.npload if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}            
            if md.Options['Energy market'] == 0:
                Firmbenefit     = {co:sum(md.se.value*md.EeFIRMHP[md.t_year[t],hp].value*md.eHppFirm[hp] for t in md.ntime for hp in md.nhpp if md.hp_country[hp]==co)/leni(md.nyear) for co in md.ncountry}
                Scndbenefit     = {co:sum(md.se.value*(sum(md.EeHPPROD[t,pld,hp].value for pld in md.npload)-md.EeFIRMHP[md.t_year[t],hp].value)*md.eHppVal[hp] for t in md.ntime for hp in md.nhpp if md.hp_country[hp]==co)/leni(md.nyear) for co in md.ncountry}
                EnergyBenefit_co ={co:Firmbenefit[co]+Scndbenefit[co] for co in md.ncountry}
        #Energy import/export (rem: if 2 power markets are in the same country, the transmission will be counted as export and import - does not affect the balance)
            EnergyImpCost_co    = {co:sum(md.se.value*md.EeTRANS[t,pld,tl].value * -1/md.se.value*md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co)/leni(md.nyear) for co in md.ncountry}
            EnergyExpBenefit_co = {co:sum(md.se.value*md.EeTRANS[t,pld,tl].value * -1/md.se.value*md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransIn[tl]]==co)/leni(md.nyear) for co in md.ncountry}
            oo['EnergyExports_co']={co:sum(md.se.value*md.EeTRANS[t,pld,tl].value                                                    for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransIn[tl]]==co)/leni(md.nyear) for co in md.ncountry}
        #Energy transmission costs
            EnergyTransCost_co  = {co:sum(md.se.value*md.EeTRANS[t,pld,tl].value * md.eTransCost[tl] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co)/leni(md.nyear) for co in md.ncountry}
        #Operation costs of hydropower plants        
            HpOMCost_co         = {co:sum(md.se.value*md.EeHPPROD[t,pld,hp].value * md.eHppCost[hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_country[hp]==co)/leni(md.nyear) for co in md.ncountry}
        #Operation and fuel costs of other power plants                     
            OppOMCost_co        = {co:sum(md.se.value*md.EeOPPROD[t,pld,pp].value * md.eOppCost[pp] for t in md.ntime for pld in md.npload for pp in md.nopp if md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}                                                                  
            OppFuelCost_co      = {co:sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCost[md.t_year[t],md.op_pmarket[pp],fu] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}                 
            OppCO2Cost_co       = {co:sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCO2[fu]*md.eCO2Val[md.t_year[t],md.op_pmarket[pp]] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}
            oo['CO2Emission_co'] = {co:sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCO2[fu]                                 for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) for co in md.ncountry}                             
        #Capacity, Operation and fuel costs of generic power technologies 
            GenCapCost_co       = {co:sum(md.se.value*md.EeGENCAP[y,pt,pm].value*md.eCAPEX[y,pt,pm]*min(1,(md.t_year[md.Options['tfin']]-y+1)/md.eLifeTime[pt,pm]) + sum(md.se.value*md.EeGENCAP[ky,pt,pm].value for ky in md.nyear if ky <= y and ky >= y-md.eLifeTime[pt,pm])*md.eFixOPEX[pt,pm] for y in md.nyear for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry} 
            GenProdCost_co      = {co:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value*md.eVarOPEX[pt,pm] for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}                                                              
            GenFuelCost_co      = {co:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCost[md.t_year[t],pm,fu]            for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}  
            GenCO2Cost_co       = {co:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCO2[fu]*md.eCO2Val[md.t_year[t],pm] for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry} 
            oo['CO2Emission_co'] = {co:oo['CO2Emission_co'][co]+sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCO2[fu]     for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}     
        #CONSUMER AND PRODUCER SURPLUS 
            oo['EnergyConsSurplus_co']  = {co:sum(md.se.value*md.EeSUPPLY[t,pld,pm].value * (md.eEngyVal[md.t_month[t],pm] - -1/md.se.value*md.dual[md.engy_balance[t,pld,pm]]/(1-md.eSupLoss[pm]))
                                          for t in md.ntime for pm in md.npmarket for pld in md.npload if md.pmarket_country[pm]==co)/leni(md.nyear) 
                                       for co in md.ncountry}
            oo['EnergyProdSurplus_co']  = {co: sum(md.se.value*md.EeHPPROD[t,pld,hp].value * (-1/md.se.value*md.dual[md.engy_balance[t,pld,md.hp_pmarket[hp]]] if md.Options['Energy market'] == 1 else md.eHppVal[hp]) for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_country[hp]==co)/leni(md.nyear) 
                                        + sum(md.se.value*md.EeOPPROD[t,pld,pp].value * (-1/md.se.value*md.dual[md.engy_balance[t,pld,md.op_pmarket[pp]]]) for t in md.ntime for pld in md.npload for pp in md.nopp if md.pmarket_country[md.op_pmarket[pp]]==co)/leni(md.nyear) 
                                        + sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value * (-1/md.se.value*md.dual[md.engy_balance[t,pld,pm]]) for t in md.ntime for pld in md.npload for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear)
                                        - HpOMCost_co[co] 
                                        - OppOMCost_co[co] - OppFuelCost_co[co] - OppCO2Cost_co[co] 
                                        - GenCapCost_co[co] - GenProdCost_co[co] - GenFuelCost_co[co] - GenCO2Cost_co[co]
                                       for co in md.ncountry}
        #Balances
            EngyBalance_co      = {co:EnergyBenefit_co[co] + EnergyExpBenefit_co[co] - GenCapCost_co[co] - GenProdCost_co[co] - GenFuelCost_co[co] - GenCO2Cost_co[co]
                                    - EnergyImpCost_co[co] - OppFuelCost_co[co] - OppCO2Cost_co[co] - HpOMCost_co[co] - OppOMCost_co[co] - EnergyTransCost_co[co] for co in md.ncountry}        
            
            EconomicBalance         = sum(UserBenefit_co[co] for co in md.ncountry) + sum(CropBalance_co[co] for co in md.ncountry) + sum(EngyBalance_co[co] for co in md.ncountry)
            oo['EconomicBalance_co'] = {co: WaterBalance_co[co] + CropBalance_co[co] + EngyBalance_co[co] for co in md.ncountry}
            oo['EconomicBalance']    = {'Total':EconomicBalance,'Water':sum(WaterBalance_co[co] for co in md.ncountry),'Energy':sum(EngyBalance_co[co] for co in md.ncountry),'Agriculture':sum(CropBalance_co[co] for co in md.ncountry)}
            oo['WaterEconomicBalance_co'] = UserBenefit_co
            oo['EnergyEconomicBalance_co'] = EngyBalance_co
            oo['AgricultureEconomicBalance_co'] = CropBalance_co
            oo['AgricultureEconomicBalance2_co'] = CropBalance2_co
    
            ##OTHER KEY INDICATORS##
            #Downstream flow [Mm3]
            oo['DownstreamFlow'] = sum(md.sw.value*md.WwOUTFLOW[t,c].value for t in md.ntime for c in md.ncatch if md.catch_ds[c]=='outlet')/leni(md.nyear)
            #Irrigation net water consumption [Mm3]
            oo['IrrigNetCons_co']= {co:sum(md.AwSUPPLY[t,fz,cul].value * (1/(1-md.aIrrgLoss[fz]) - md.aCulRturn[md.fzone_type[fz],cul]) for t in md.ntime for cul in md.nculture for fz in md.nifzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
            #Net cultivated area [Mha]
            oo['NetCulArea_co']  = {co:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
            #Gross cultivated area [Mha]
            oo['GrossCultivatedArea_co']= {co:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
            #Gross irrigated area [Mha]
            oo['GrossIrrArea_co']= {co:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co and md.aIrrigation[md.fzone_type[fz]]==1)/leni(md.nyear) for co in md.ncountry}
            #Net irrigated area [Mha]
            oo['NetIrrArea_co']  = {co:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co and md.aIrrigation[md.fzone_type[fz]]==1)/leni(md.nyear) for co in md.ncountry}
            #Production [Mt]
            oo['CropProd_co']    = {co:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry}
    
            #Hydropower prod
            oo['HydropowerProduction_co']  = {co:sum(md.se.value*md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_country[hp]==co)/leni(md.nyear) for co in md.ncountry} 
            oo['HydroFirmProd_co']   = {co:sum(sorted([sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)/100*5))] for hp in md.nhpp if md.hp_country[hp]==co) for co in md.ncountry}
            oo['HydroFirmVal_co']    = {co:sum(sorted([sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)/100*5))]*leni(md.nmonth)*0.06
                                                  +(sum(md.se.value*md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload)/leni(md.nyear)-sorted([sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)/100*5))]*leni(md.nmonth))*0.02
                                                  for hp in md.nhpp if md.hp_country[hp]==co) for co in md.ncountry}
            #Energy exports
            oo['EnergyExportBenefit_co'] = EnergyExpBenefit_co
            #Energy capacity investment #ADD WARNING THIS IS CASE SPECIFIC - Will not cause bug however !
            oo['EnergyCapacityConst_co']= {co:sum(md.se.value*md.EeGENCAP[y,pt,pm].value for y in md.nyear for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co and pt in ['CoalPP','GasPP'])/leni(md.nyear) for co in md.ncountry} 
            oo['SolarCapacityConst_co'] = {co:sum(md.se.value*md.EeGENCAP[y,pt,pm].value for y in md.nyear for pt in md.nptech for pm in md.npmarket if md.pmarket_country[pm]==co and pt =='SolarPV')/leni(md.nyear) for co in md.ncountry} 
            #Energy curtailment
            oo['EnergyCurt_co']      = {co:sum(md.eEngyDem[md.t_year[t],md.t_month[t],pm]-sum(md.se.value*md.EeSUPPLY[t,pld,pm].value for pld in md.npload) for t in md.ntime for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry}
            oo['OtherIndicators'] = {'Available runoff [Mm3]':        sum(md.wRunOff[t,c].value for t in md.ntime for c in md.ncatch)/leni(md.nyear),
                                    'Downstream flow [Mm3]':         oo['DownstreamFlow'],
                                    'Net irrigation cons [Mm3/year]':sum(oo['IrrigNetCons_co'][co] for co in md.ncountry),
                                    'Gross Cultivated area [Mha]':  sum(oo['GrossCultivatedArea_co'][co] for co in md.ncountry),
                                    'Gross Irrigated area [Mha]':   sum(oo['GrossIrrArea_co'][co] for co in md.ncountry),
                                    'Net Irrigated area [Mha]':     sum(oo['NetIrrArea_co'][co] for co in md.ncountry),
                                    'Net Cultivated area [Mha]':    sum(oo['NetCulArea_co'][co] for co in md.ncountry),
                                    'Crop exports benefits [M$/year]':sum(oo['CropExpBenefit_co'][co] for co in md.ncountry),
                                    'Hydropower production [GWh]':  sum(oo['HydropowerProduction_co'][co] for co in md.ncountry),
                                    'Hydropower prod FIRM [GWh/month]':   sum(oo['HydroFirmProd_co'][co] for co in md.ncountry),
                                    'Hydropower val FIRM [M$/year]':    sum(oo['HydroFirmVal_co'][co] for co in md.ncountry),
                                    'Energy Curtailment [GWh]':     sum(oo['EnergyCurt_co'][co] for co in md.ncountry),
                                    'Energy trade [GWh/year]':      sum(oo['EnergyExports_co'][co] for co in md.ncountry),
                                    'Generic Power cap [MW]':       sum(oo['EnergyCapacityConst_co'][co] for co in md.ncountry),
                                    'Generic Solar cap [MW]':       sum(oo['SolarCapacityConst_co'][co] for co in md.ncountry), #ADD STUDY CASE SPECIFIC !!
                                    'CO2 emissions [Mt/year]':      sum(oo['CO2Emission_co'][co] for co in md.ncountry),
                                    'Reservoir Evaporation [Mm3/year]':ReservoirEvaporation,
                                    'Reservoir Final Storage [Mm3]':sum(md.WwRSTORAGE[md.Options['tfin'],res].value*md.sw.value for res in md.nres),
                                    'Eflow fail [Mm3/year]':        sum(md.DUMYEFLOW[t,ef].value for t in md.ntime for ef in md.neflow)/leni(md.nyear) if md.Options['Debug mode'] in [1,2] else 0,
                                    'Debug Water [Mm3/year]':       sum(md.DUMYWATER[t,c].value for t in md.ntime for c in md.ncatch)/leni(md.nyear) if md.Options['Debug mode'] in [1] else 0,
                                    'Water shadowprice [$/m3]':     WaterValueShadow,
                                    'Energy Shadowprice [$/kWh]':   sum(oo['EnergyValue_co'][co] for co in oo['EnergyValue_co'].keys() if co != 'SouthAfrica')/max(1,leni(oo['EnergyValue_co'].keys())-1) if md.Options['Energy market'] == 1 else 0,
                                    'Crop Shadowprice [$/t]':       sum(oo['CropValue_co'][co] for co in oo['CropValue_co'].keys() if co != 'World')/max(1,leni(oo['CropValue_co'].keys())-1),
                                    'Agr. cons. surplus [M$/year]': sum(oo['AgricultureConsSurplus_co'][co] for co in md.ncountry),
                                    'Agr. prod. surplus [M$/year]': sum(oo['AgricultureProdSurplus_co'][co] for co in md.ncountry),
                                    'Egy. cons. surplus [M$/year]': sum(oo['EnergyConsSurplus_co'][co] for co in md.ncountry),
                                    'Egy. prod. surplus [M$/year]': sum(oo['EnergyProdSurplus_co'][co] for co in md.ncountry),
                                    ' Water User Benefit [M$/year]':sum(UserBenefit_co[co] for co in md.ncountry),
                                    ' Water Supply Cost [M$/year]': sum(WaterSupplyCosts_co[co] for co in md.ncountry),
                                    ' Crop Supply Benefit [M$/year]':sum(CropBenefit_co[co] for co in md.ncountry),
                                    ' Crop Cul Cost [M$/year]':     sum(CropCultivationCost_co[co] for co in md.ncountry),
                                    ' Crop Irrig Cost [M$/year]':   sum(CropIrrigationCost_co[co] for co in md.ncountry),
                                    ' Crop Ext Prod Cost [M$/year]':sum(CropExtProdCost_co[co] for co in md.ncountry),
                                    ' Crop Trans Cost [M$/year]':   sum(CropTransportCost_co[co] for co in md.ncountry),
                                    ' Energy Supply Benefit [M$/year]': sum(EnergyBenefit_co[co] for co in md.ncountry),
                                    ' Energy Gen Cap Cost [M$/year]':   sum(GenCapCost_co[co] for co in md.ncountry),
                                    ' Energy Gen Prod Cost [M$/year]':  sum(GenProdCost_co[co] for co in md.ncountry),
                                    ' Energy Gen Fuel Cost [M$/year]':  sum(GenFuelCost_co[co] for co in md.ncountry),
                                    ' Energy Gen CO2 Cost [M$/year]':   sum(GenCO2Cost_co[co] for co in md.ncountry),
                                    ' Energy Opp Fuel Cost [M$/year]':  sum(OppFuelCost_co[co] for co in md.ncountry),
                                    ' Energy Opp OetM Cost [M$/year]':  sum(OppOMCost_co[co] for co in md.ncountry),                                
                                    ' Energy Opp CO2 Cost [M$/year]':   sum(OppCO2Cost_co[co] for co in md.ncountry),                                
                                    ' Energy Hpp OetM Cost [M$/year]':  sum(HpOMCost_co[co] for co in md.ncountry),
                                    ' Energy Trans Cost [M$/year]':     sum(EnergyTransCost_co[co] for co in md.ncountry)} 
            oo['Options']={key:md.Options[key] for key in md.Options}
        return oo 
#%% Read all results       
    def readresults(self,md,parameters,solverstatus,PRINT=0,VALIDATION=0,scenario='WHATIF_main'):        
        oo={} #initiate result dictionary
        oo['Options']={key:md.Options[key] for key in md.Options}
        #Read parameters function #warning: does not read growing parameters
        def ReadParam(ParamName):
            if ParamName in parameters.scen.keys():                
                return parameters.val[ParamName][parameters.val[parameters.scen[ParamName]][scenario]]
            else:
                return parameters.val[ParamName]
#MAIN
    # Optimization Parameters
        oo['MainTable']={'Time and date of run': time.strftime("%c"),
                        'Number of years': len(md.nyear),
                        'Run time [s]': solverstatus.solver.Time,
                        'Scenario': scenario}                    

#%%SHADOWPRICES
        CONSTRAINTLIST = [const.name for const in md.component_objects(Constraint)]
        #Length of index function (return 1 for void index instead of 0 to avoid division by 0)
        def leni(index):
            return max(1,len(index))
        
        #Water shadow value
        oo['WaterValueShadow_ca']= {c: -sum(md.dual[md.water_waterbalance[t,c]]   for t in md.ntime)                     / leni(md.ntime)                 for c in md.ncatch}
        #Hydropower turbine capacity shadow value
        oo['HpCapShadow_hp']     = {hp:-sum(md.dual[md.engy_hpcapacity[t,pld,hp]] for t in md.ntime for pld in md.npload)/(leni(md.ntime)*leni(md.npload)) for hp in md.nhpp}
        #Reservoir storage capacity shadow value
        oo['ResCapShadow_res']   = {res:-sum(md.dual[md.water_reservoircapacity[t,res]] for t in md.ntime)               / leni(md.ntime)                 for res in md.nres if res not in md.nlake}    
        #Environmental flow
        oo['EnvFlowShadow_m_ec'] = {m:{ec:sum(md.dual[md.env_minflow[t,ec]] for t in md.ntime if md.t_month[t]==m)/leni(md.nyear) for ec in md.neflow} for m in md.nmonth}        
        #Energy value
        EnergyValueShadow       = -sum(1/md.se.value*md.dual[md.engy_balance[t,pld,pm]]*md.eLoadDem[pld] for t in md.ntime for pld in md.npload for pm in md.npmarket)/(leni(md.ntime)*leni(md.npmarket))
        #Transmition capacity
        TransmitionCapShadow    = -sum(md.dual[md.engy_transmissioncap[t,pld,tl]] for t in md.ntime for pld in md.npload for tl in md.ntransline)/(leni(md.ntime)*leni(md.npload)*leni(md.ntransline))
        #Power plants capacity
        oo['PpCapShadow_pp'] = {pp:-sum(md.dual[md.engy_oppcapacity[t,pld,pp]] for t in md.ntime for pld in md.npload)/(leni(md.ntime)*leni(md.npload)) for pp in md.nopp}                       
         #Land capacity shadow value
        oo['AvLandShadow_fz']    = {fz:-sum(1/md.kha_to_Mha*md.dual[md.agr_maxland[y,fz]] for y in md.nyear) /(leni(md.nyear)) if 'agr_maxland' in CONSTRAINTLIST else 0 for fz in md.nfzone}
        #Crop price is determined by cropbalance - of farm value if crop market option is off
        if md.Options['Crop market'] == 1:
            oo['CropPrice_cr_y'] = {cr:{y: sum(-1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,cm,cr]] for cm in md.ncmarket)/leni(md.ncmarket) for y in md.nyear} for cr in md.ncrop}
            oo['CropPrice_cr_cm']= {cr:{cm:sum(-1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,cm,cr]] for y in md.nyear)/leni(md.nyear) for cm in md.ncmarket} for cr in md.ncrop}
        else: 
            oo['CropPrice_cr_y'] = {cr:{y: sum(md.aFarmVal[co,cr] for co in md.ncountry)/leni(md.ncountry) for y in md.nyear}  for cr in md.ncrop}
            oo['CropPrice_cr_cm']= {cr:{co:sum(md.aFarmVal[co,cr] for y in md.nyear)    /leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}

        if PRINT==1:
            print('          <----------------SHADOWPRICES--------------->')
            #Water
            print('Water value [$/m3]')
            print(sum(oo['WaterValueShadow_ca'][ca] for ca in md.ncatch)/leni(md.ncatch))
            print('Storage capacity [$/m3]')
            print(sum(oo['ResCapShadow_res'][res] for res in md.nres)/leni(md.nres))
            print('Environmental flow [$/m3]')
            print(sum(oo['EnvFlowShadow_m_ec'][m][ec] for m in md.nmonth for ec in md.neflow)/leni(md.nmonth)/leni(md.neflow))
            #Agriculture
            print('Available Land [$/ha]')
            print(sum(oo['AvLandShadow_fz'][fz] for fz in md.nfzone)/leni(md.nfzone))
            print('Crop value [$/t]')
            print(sum(oo['CropPrice_cr_y'][cr][y] for cr in md.ncrop for y in md.nyear)/leni(md.ncrop)/leni(md.nyear))
            #Energy
            print('Energy value [$/kWh]')
            print(EnergyValueShadow)            
            print('Transmition lines capacity [$/kWh]')
            print(TransmitionCapShadow)
            print('Hydropower capacity [$/kWh]')
            print(sum(oo['HpCapShadow_hp'][hp] for hp in md.nhpp)/leni(md.nhpp))

#%%WATER MODULE
        #WATER BALANCES
        oo['WaterBalance_t_c'] = {
                'RunOff [Mm3]':                 {t: {c: md.wRunOff[t,c].value for c in md.ncatch} for t in md.ntime},
                'Groundwater recharge [Mm3]':   {t: {c: sum(md.wGwRech[t,aq]              for aq in md.naquifer      if md.aqui_catch[aq]==c)        for c in md.ncatch} for t in md.ntime},
                'Upstream inflow [Mm3]':        {t: {c: sum(md.sw.value*md.WwOUTFLOW[t,kc].value      for kc in md.ncatch        if md.catch_ds[kc] == c)        for c in md.ncatch} for t in md.ntime},
                'Transfer inflow [Mm3]':        {t: {c: sum(md.WwTRANSFER[t,ktrans].value for ktrans in md.ntransfer if md.transfer_ds[ktrans] == c) for c in md.ncatch} for t in md.ntime}, 
                'Transfer outflow [Mm3]':       {t: {c: sum(md.WwTRANSFER[t,ktrans].value for ktrans in md.ntransfer if md.transfer_us[ktrans] == c) for c in md.ncatch} for t in md.ntime},
                'Transfer loss [Mm3]':          {t: {c: sum(md.WwTRANSFER[t,ktrans].value * md.wTransLoss[ktrans] for ktrans in md.ntransfer if md.transfer_ds[ktrans] == c) for c in md.ncatch} for t in md.ntime}, 
                'User allocation [Mm3]':        {t: {c: sum(md.WwSUPPLY[t,u].value        for u in md.nuser          if md.user_catch[u]==c)         for c in md.ncatch} for t in md.ntime},
                'User return flow [Mm3]':       {t: {c: sum(md.WwSUPPLY[t,u].value * md.wUserRturn[u] for u in md.nuser if md.user_catch[u]==c)      for c in md.ncatch} for t in md.ntime},
                'User water loss [Mm3]':        {t: {c: sum(md.WwSUPPLY[t,u].value * md.wUserLoss[u]/(1-md.wUserLoss[u]) for u in md.nuser if md.user_catch[u]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture irrigation [Mm3]': {t: {c: sum(md.AwSUPPLY[t,fz,cul].value   for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture groundwater [Mm3]':{t: {c: sum(md.AwGWSUPPLY[t,fz,cul].value for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]==md.aqui_catch[aq]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture return flow [Mm3]':{t: {c: sum(md.AwSUPPLY[t,fz,cul].value * md.aCulRturn[md.fzone_type[fz],cul]   for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c) for c in md.ncatch} for t in md.ntime},
                'Agriculture irrig. loss [Mm3]':{t: {c: sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgLoss[fz]/(1-md.aIrrgLoss[fz]) for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c) for c in md.ncatch} for t in md.ntime},
                'Outlet flow [Mm3]':            {t: {c: md.sw.value*md.WwOUTFLOW[t,c].value for c in md.ncatch} for t in md.ntime},
                'River ET [Mm3]':               {t: {c: sum(md.sw.value*md.WwOUTFLOW[t,kc].value*md.wFlowLoss[c] for kc in md.ncatch if md.catch_ds[kc] == c) for c in md.ncatch} for t in md.ntime},
                'Reservoir storage [Mm3]':      {t: {c: sum(md.sw.value*md.WwRSTORAGE[t,res].value for res in md.nres if md.res_catch[res] == c) for c in md.ncatch} for t in md.ntime},
                'Groundwater storage [Mm3]':    {t: {c: sum(md.WwGWSTORAGE[t,aq].value for aq in md.naquifer if md.aqui_catch[aq] == c) for c in md.ncatch} for t in md.ntime}
                                 }
        if md.Options['Initial time step'] == 0: #Cyclic Run
            oo['WaterBalance_t_c']['Reservoir evaporation [Mm3]']={t: {c: sum(max(0,md.wET0[t,c].value-1*md.wRainFall[t,c].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[t,res].value+md.sw.value*md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for res in md.nres if md.res_catch[res]==c) for c in md.ncatch} for t in md.ntime}    
        else:
            oo['WaterBalance_t_c']['Reservoir evaporation [Mm3]']={t: {c: sum(max(0,md.wET0[t,c].value-1*md.wRainFall[t,c].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[t,res].value+md.sw.value*md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for res in md.nres if md.res_catch[res]==c) for c in md.ncatch} for t in md.ntime if t != md.Options['tini']}
            oo['WaterBalance_t_c']['Reservoir evaporation [Mm3]'][md.Options['tini']]= {c: sum(max(0,md.wET0[md.Options['tini'],c].value-1*md.wRainFall[md.Options['tini'],c].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[md.Options['tini'],res].value + md.wStorIni[res].value)/2 + md.wResArea[res]) for res in md.nres if md.res_catch[res]==c) for c in md.ncatch}
                                 
        oo['WaterBalance_t_c'][' Balance = 0 [Mm3]']= {t: {c: + oo['WaterBalance_t_c']['RunOff [Mm3]'][t][c]
                                                             + oo['WaterBalance_t_c']['Groundwater recharge [Mm3]'][t][c]
                                                             + oo['WaterBalance_t_c']['User return flow [Mm3]'][t][c]
                                                             + oo['WaterBalance_t_c']['Agriculture return flow [Mm3]'][t][c]
                                                             + oo['WaterBalance_t_c']['Transfer inflow [Mm3]'][t][c]
                                                             + oo['WaterBalance_t_c']['Upstream inflow [Mm3]'][t][c]
                                                             + (oo['WaterBalance_t_c']['Reservoir storage [Mm3]'][md.t_prev[t]][c] if md.Options['Initial time step'] == 0 or t != md.Options['tini'] else sum(md.wStorIni[res].value for res in md.nres if md.res_catch[res]==c))
                                                             + (oo['WaterBalance_t_c']['Groundwater storage [Mm3]'][md.t_prev[t]][c] if md.Options['Initial time step'] == 0 or t != md.Options['tini'] else sum(md.wGwIni[aq].value for aq in md.naquifer if md.aqui_catch[aq]==c))
                                                             - oo['WaterBalance_t_c']['User allocation [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['User water loss [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Agriculture irrigation [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Agriculture irrig. loss [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Agriculture groundwater [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Transfer outflow [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Transfer loss [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Reservoir storage [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Reservoir evaporation [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Groundwater storage [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['River ET [Mm3]'][t][c]
                                                             - oo['WaterBalance_t_c']['Outlet flow [Mm3]'][t][c]
                                                             for c in md.ncatch} for t in md.ntime}
        #WATER BALANCES
        oo['WaterBalance_y'] = {keys:{y:sum(oo['WaterBalance_t_c'][keys][t][ca] for t in md.ntime for ca in md.ncatch if md.t_year[t]==y) for y in md.nyear} for keys in oo['WaterBalance_t_c'].keys()}
        oo['WaterBalance_ca'] = {keys:{ca:sum(oo['WaterBalance_t_c'][keys][t][ca] for t in md.ntime)/leni(md.nyear) for ca in md.ncatch} for keys in oo['WaterBalance_t_c'].keys()}
        oo['WaterBalance_ca']['Initial storage [Mm3]'] = {c:sum(md.wStorIni[res].value for res in md.nres if md.res_catch[res]==c) + sum(md.wGwIni[aq] for aq in md.naquifer if md.aqui_catch[aq]==c) for c in md.ncatch} if md.Options['Initial time step'] == 1 else {c:0 for c in md.ncatch}
        oo['WaterBalance_ca']['End storage [Mm3]']     = {c:sum(md.sw.value*md.WwRSTORAGE[md.Options['tfin'],res].value for res in md.nres if md.res_catch[res]==c) + sum(md.WwGWSTORAGE[md.Options['tfin'],aq].value for aq in md.naquifer if md.aqui_catch[aq]==c) for c in md.ncatch} if md.Options['Initial time step'] == 1 else {c:0 for c in md.ncatch}
        oo['WaterBalance_co'] = {keys:{co:sum(oo['WaterBalance_ca'][keys][c] for c in md.ncatch if md.catch_country[c] == co) for co in md.ncountry} for keys in oo['WaterBalance_ca'].keys()}
        oo['WaterBalance']    = {keys:sum(oo['WaterBalance_ca'][keys][c] for c in md.ncatch) for keys in oo['WaterBalance_ca'].keys()}                          
        oo['WaterBalance_co']['Upstream inflow [Mm3]'] = {co:sum(oo['WaterBalance_ca']['Outlet flow [Mm3]'][kc] for kc in md.ncatch if (md.catch_ds[kc] != 'outlet' and md.catch_country[md.catch_ds[kc]]==co) and md.catch_country[kc] != co) for co in md.ncountry}
        oo['WaterBalance_co']['Outlet flow [Mm3]']     = {co:sum(oo['WaterBalance_ca']['Outlet flow [Mm3]'][kc] for kc in md.ncatch if (md.catch_ds[kc] == 'outlet' or  md.catch_country[md.catch_ds[kc]]!=co) and md.catch_country[kc] == co) for co in md.ncountry} 
        oo['WaterBalance']['Outlet flow [Mm3]']     = oo['WaterBalance']['Outlet flow [Mm3]'] - oo['WaterBalance']['Upstream inflow [Mm3]']
        oo['WaterBalance']['Upstream inflow [Mm3]'] = 0
        
        #RESERVOIRS
        oo['ReservoirTable']  = {'AvStorage [Mm3]':               {res:sum(md.sw.value*md.WwRSTORAGE[t,res].value/leni(md.ntime) for t in md.ntime) for res in md.nres},
                                'Capacity [Mm3]':                ReadParam('wStorCap'), 
                                'Capacity shadowprice [$/m3]':  oo['ResCapShadow_res']
                                }
        
        if md.Options['Initial time step'] != 1: #Cyclic Run
            oo['ReservoirTable']['Evaporation [Mm3/y]']  = {res: sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[t,res].value + md.sw.value*md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime)/leni(md.nyear) for res in md.nres}
        else:
            oo['ReservoirTable']['Evaporation [Mm3/y]']  = {res: (max(0,md.wET0[md.Options['tini'],md.res_catch[res]].value-1*md.wRainFall[md.Options['tini'],md.res_catch[res]].value) * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[md.Options['tini'],res].value + md.wStorIni[res].value)/2 + md.wResArea[res])
                                                               +sum(max(0,md.wET0[t,md.res_catch[res]].value-1*md.wRainFall[t,md.res_catch[res]].value)                                * md.mm_to_m3perha * (md.wkV[res]*(md.sw.value*md.WwRSTORAGE[t,res].value + md.sw.value*md.WwRSTORAGE[md.t_prev[t],res].value)/2 + md.wResArea[res]) for t in md.ntime if t != md.Options['tini']))
                                                               /leni(md.nyear) for res in md.nres}                                                                       
        oo['ResStor_res_m']   = {m:{res:sum(md.sw.value*md.WwRSTORAGE[t,res].value for t in md.ntime if md.t_month[t]==m)/sum(1 for t in md.ntime if md.t_month[t]==m) for res in md.nres} for m in md.nmonth}
        oo['ResStor_res_t']   = {t:{res:md.sw.value*md.WwRSTORAGE[t,res].value for res in md.nres} for t in md.ntime}
        
        #WATER USERS
        oo['UserTable']       = {'User consumption [Mm3/y]':      {u:sum(md.WwSUPPLY[t,u].value * (1-md.wUserRturn[u])                            for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'User loss [Mm3/y]':             {u:sum(md.WwSUPPLY[t,u].value * md.wUserLoss[u]/(1-md.wUserLoss[u])         for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'Net User abstraction [Mm3/y]':  {u:sum(md.WwSUPPLY[t,u].value * (1/(1-md.wUserLoss[u])-md.wUserRturn[u])   for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'User benefits [M$/y]':          {u:sum(md.WwSUPPLY[t,u].value * md.wUserVal[u]                                         for t in md.ntime)/leni(md.nyear) for u in md.nuser},
                                'User curtailment [Mm3/y]':      {u:sum(-md.WwSUPPLY[t,u].value + md.wUserDem[md.t_year[t],md.t_month[t],u]                    for t in md.ntime)/leni(md.nyear) for u in md.nuser}
                                }

        #GROUNDWATER
        oo['GwStor_c_m']  = {m:{c:sum(md.WwGWSTORAGE[t,aq].value for t in md.ntime for aq in md.naquifer if md.t_month[t]==m and md.aqui_catch[aq]==c)/leni(md.nyear) for c in md.ncatch} for m in md.nmonth}
        oo['GwStor_c_t']  = {t:{c:sum(md.WwGWSTORAGE[t,aq].value for aq in md.naquifer if md.aqui_catch[aq]==c) for c in md.ncatch} for t in md.ntime}
        
        #%%PRINTS
        if PRINT ==1:
            print('          <----------------OPTIMAL DECISIONS--------------->')
            print('*************WATER*************')
            print('Inflow [Mm3]')
            print(oo['WaterBalance']['RunOff [Mm3]'])
            print('Water consumption per user [Mm3]')
            print('user:',sum(md.WwSUPPLY[t,u].value for t in md.ntime for u in md.nuser)/leni(md.nyear),'demand',sum(md.wUserDem[md.t_year[t],md.t_month[t],u] for t in md.ntime for u in md.nuser)/leni(md.nyear))
            print('Downstream flow [Mm3]')
            print(oo['WaterBalance']['Outlet flow [Mm3]'])
            print('Water losss [Mm3]')
            print('User:',oo['WaterBalance']['User water loss [Mm3]'])
            print('Agriculture:', oo['WaterBalance']['Agriculture irrig. loss [Mm3]'])
            print('Return Flows [Mm3]')
            print(oo['WaterBalance']['User return flow [Mm3]']+oo['WaterBalance']['Agriculture return flow [Mm3]'])
            print('End Storage [Mm3]')
            print(oo['WaterBalance']['End storage [Mm3]'])
            print('Initial Storage [Mm3]')
            print(oo['WaterBalance']['Initial storage [Mm3]'])
            print('River and Reservoir ET')
            print(oo['WaterBalance']['River ET [Mm3]']+oo['WaterBalance']['Reservoir evaporation [Mm3]'])
            print('Water Balance [Mm3]: InitialStorage + InFlow + ReturnFlow = UserAll + AgrAll + DownFlow + EndStorage + ET-Losses')
            print(oo['WaterBalance'][' Balance = 0 [Mm3]'])
            
#%% AGRICULTURE MODULE 
        #Available land, avergae cultivated land  
        if md.Options['Crop choice'] == 'fixed':
            GROSSAREA   = {(y,fz,cul): md.kha_to_Mha*md.CULAREA[y,fz,cul].value for y in md.nyear for fz in md.nfzone for cul in md.nculture}
            NETAREA     = {(y,fz): sum(md.kha_to_Mha*md.CULAREA[y,fz,cul].value for cul in md.nculture) for fz in md.nfzone for y in md.nyear}
        elif md.Options['Yield water response'] == 'nonlinear':
            GROSSAREA   = {(y,fz,cul):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd].value for fd in md.nfieldculture if md.field_culture[fd]==cul) for cul in md.nculture for fz in md.nfzone for y in md.nyear}
            NETAREA     = {(y,fz):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd].value for fd in md.nfieldculture if fd[1]==1) for fz in md.nfzone for y in md.nyear}
        else:
            GROSSAREA   = {(y,fz,cul):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if md.field_culture[fd]==cul) for cul in md.nculture for fz in md.nfzone for y in md.nyear}
            NETAREA     = {(y,fz):sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt].value for fd in md.nfieldculture for ypt in md.nypath if fd[1]==1) for fz in md.nfzone for y in md.nyear}      
        if md.Options['Crop market'] == 1:
            CropPrice   = {(y,fz,cr):-1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.fzone_cmarket[fz],cr]] for cr in md.ncrop for fz in md.nfzone for y in md.nyear}
        else:
            CropPrice   = {(y,fz,cr): md.aFarmVal[md.fzone_country[fz],cr] for cr in md.ncrop for fz in md.nfzone for y in md.nyear}
        
        #CROP BALANCE
        oo['CropBalance_y_fz']={'Available land [Mha]':      {y: {fz:md.kha_to_Mha*md.aLandCap[y,fz] for fz in md.nfzone} for y in md.nyear},
                               'Net Cultivated land [Mha]': {y: {fz:NETAREA[y,fz] for fz in md.nfzone} for y in md.nyear},
                               'Net Irrigated land [Mha]':  {y: {fz:NETAREA[y,fz] if md.aIrrigation[md.fzone_type[fz]]==1 else 0 for fz in md.nfzone} for y in md.nyear},
                               'Crop production [Mt]':      {y: {fz:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value  for cr in md.ncrop) for fz in md.nfzone} for y in md.nyear}}
        oo['CropBalance_ca']= {keys:{ca:sum(oo['CropBalance_y_fz'][keys][y][fz] for y in md.nyear for fz in md.nfzone if md.fzone_catch[fz]==ca)/leni(md.nyear) for ca in md.ncatch} for keys in oo['CropBalance_y_fz'].keys()}
        
        oo['CropBalance_y_cm']={'Available land [Mha]':      {y: {cm: sum(md.kha_to_Mha*md.aLandCap[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Net Cultivated land [Mha]': {y: {cm: sum(NETAREA[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Net Irrigated land [Mha]':  {y: {cm: sum(NETAREA[y,fz] for fz in md.nfzone if md.fzone_cmarket[fz]==cm and md.aIrrigation[md.fzone_type[fz]]==1) for cm in md.ncmarket} for y in md.nyear},
                               'Crop production [Mt]':      {y: {cm: sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value for fz in md.nfzone for cr in md.ncrop if md.fzone_cmarket[fz]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Crop demand [Mt]':          {y: {cm: sum(md.kt_to_Mt*md.aCropDem[y,cm,cr] for cr in md.ncrop) for cm in md.ncmarket} for y in md.nyear},
                               'Crop supply [Mt]':          {y: {cm: sum(md.kt_to_Mt*md.AcSUPPLY[y,cm,cr,cds].value for cr in md.ncrop for cds in md.ncdstep) for cm in md.ncmarket} for y in md.nyear},
                               'Crop external prod [Mt]':   {y: {cm: sum(md.kt_to_Mt*md.AcEXTPROD[y,cm,cr].value for cr in md.ncrop) if cm in md.nextcmarket else 0 for cm in md.ncmarket} for y in md.nyear},
                               'Crop net import [Mt]':      {y: {cm: sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * (1-md.aTransLoss[ct]) for ct in md.nctrans for cr in md.ncrop if md.aTransOut[ct]==cm) for cm in md.ncmarket} for y in md.nyear},
                               'Crop gross export [Mt]':    {y: {cm: sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value                         for ct in md.nctrans for cr in md.ncrop if md.aTransIn[ct]==cm) for cm in md.ncmarket} for y in md.nyear}}
        oo['CropBalance_y_cm'][' Balance = 0 [Mt]']      = {y: {cm: oo['CropBalance_y_cm']['Crop supply [Mt]'][y][cm]
                                                                  +oo['CropBalance_y_cm']['Crop gross export [Mt]'][y][cm]
                                                                  -oo['CropBalance_y_cm']['Crop net import [Mt]'][y][cm]
                                                                  -oo['CropBalance_y_cm']['Crop production [Mt]'][y][cm]
                                                                  -oo['CropBalance_y_cm']['Crop external prod [Mt]'][y][cm] for cm in md.ncmarket} for y in md.nyear}
                   
        oo['CropBalance_y'] = {keys:{y:sum(oo['CropBalance_y_cm'][keys][y][cm] for cm in md.ncmarket) for y in md.nyear} for keys in oo['CropBalance_y_cm'].keys()}
        oo['CropBalance_co']= {keys:{co:sum(oo['CropBalance_y_cm'][keys][y][cm] for y in md.nyear for cm in md.ncmarket if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry} for keys in oo['CropBalance_y_cm'].keys()}
        oo['CropBalance']   = {keys:sum(oo['CropBalance_y_cm'][keys][y][cm] for y in md.nyear for cm in md.ncmarket)/leni(md.nyear) for keys in oo['CropBalance_y_cm'].keys()}
        
        if md.Options['Crop market'] == 0: #BALANCE if crop market is off            
            oo['CropBalance_y'].update({keys:{y:sum(oo['CropBalance_y_fz'][keys][y][fz] for fz in md.nfzone) for y in md.nyear} for keys in oo['CropBalance_y_fz'].keys()})
            oo['CropBalance_co'].update({keys:{co:sum(oo['CropBalance_y_fz'][keys][y][fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry} for keys in oo['CropBalance_y_fz'].keys()})
            oo['CropBalance'].update({keys:sum(oo['CropBalance_y_fz'][keys][y][fz] for y in md.nyear for fz in md.nfzone)/leni(md.nyear) for keys in oo['CropBalance_y_fz'].keys()})
        
        #FARMING ZONES
        oo['FarmZonesTable']=  {'Crop Production [Mt/y]':       {fz:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value                                                    for y in md.nyear for cr in md.ncrop)    /leni(md.nyear) for fz in md.nfzone},                               
                               'Cultivation costs [M$/y]':      {fz:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul]                    for y in md.nyear for cul in md.nculture)/leni(md.nyear) for fz in md.nfzone},
                               'Irrigation costs [M$/y]':       {fz:sum(md.AwSUPPLY[t,fz,cul].value * 1/(1-md.aIrrgLoss[fz]) * md.aIrrgCost[fz]     for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'GW pumping costs [M$/y]':       {fz:sum(md.AwGWSUPPLY[t,fz,cul].value * md.wGwCost[md.fzone_catch[fz]]              for t in md.ntime for cul in md.nculture)/leni(md.nyear) if md.Options['Groundwater']==1 else 0 for fz in md.nifzone},
                               'Cultivated land Net [Mha]':     {fz:sum(NETAREA[y,fz]                                                               for y in md.nyear)                       /leni(md.nyear) for fz in md.nfzone},
                               'Cultivated land Net [%]':       {fz:sum(NETAREA[y,fz] for y in md.nyear) / sum(md.kha_to_Mha*md.aLandCap[y,fz] for y in md.nyear) for fz in md.nfzone if sum(md.aLandCap[y,fz] for y in md.nyear) != 0},
                               'Cultivated land Gross [%]':     {fz:sum(GROSSAREA[y,fz,cul] for y in md.nyear for cul in md.nculture) / sum(md.kha_to_Mha*md.aLandCap[y,fz] for y in md.nyear) for fz in md.nfzone if sum(md.aLandCap[y,fz] for y in md.nyear) != 0},
                               'Irrig net abstraction [Mm3/y]': {fz:sum(md.AwSUPPLY[t,fz,cul].value * (1/(1-md.aIrrgLoss[fz]) - md.aCulRturn[md.fzone_type[fz],cul]) for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Irrig withdrawals [Mm3/y]':     {fz:sum(md.AwSUPPLY[t,fz,cul].value * 1/(1-md.aIrrgLoss[fz])                        for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Irrig losses [Mm3/y]':          {fz:sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgLoss[fz]/(1-md.aIrrgLoss[fz])         for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Irrig return flow [Mm3/y]':     {fz:sum(md.AwSUPPLY[t,fz,cul].value * md.aCulRturn[md.fzone_type[fz],cul]           for t in md.ntime for cul in md.nculture)/leni(md.nyear) for fz in md.nifzone},
                               'Precipitation [mm/y]':          {fz:sum(md.wRainFall[t,md.fzone_catch[fz]].value                                    for t in md.ntime)                       /leni(md.nyear) for fz in md.nfzone},
                               'Land shadowprice [$/ha]':       oo['AvLandShadow_fz'],
                               'Production Value [M$/y]':       {fz:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for y in md.nyear for cr in md.ncrop)/leni(md.nyear) for fz in md.nfzone}}                       
        
        if md.Options['Crop choice'] == 'fixed' or md.Options['Yield water response'] == 'nonlinear':
            oo['FarmZonesTable']['Crop dem. satisfaction [%]'] = {fz: 'Not available now' for fz in md.nfzone}
        else:
            oo['FarmZonesTable']['Crop dem. satisfaction [%]'] = {fz:sum(md.kha_to_Mha*md.AlCULAREA[y,fz,fd,ypt].value*sum(md.aYieldMat[ypt,kyps] for kyps in md.nyphase)/leni(md.nyphase) for y in md.nyear for fd in md.nfieldculture for ypt in md.nypath)
                                                                   /sum(GROSSAREA[y,fz,cul] for y in md.nyear for cul in md.nculture) for fz in md.nfzone if sum(GROSSAREA[y,fz,cul] for y in md.nyear for cul in md.nculture) != 0}    
        #FARM BENEFITS 
        oo['FarmBenefits_ca_ft'] = {c: {ft:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_type[fz]==ft and md.fzone_catch[fz]==c)/leni(md.nyear) 
                                        - sum(GROSSAREA[y,fz,cul] * md.aCulCost[ft,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_catch[fz]==c)/leni(md.nyear)
                                        - sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_type[fz]==ft and md.fzone_catch[fz]==c)/leni(md.nyear)
                                        for ft in md.nftype} 
                                   for c in md.ncatch}
        oo['FarmBenefits_fz_y']  = {fz:{y: sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * CropPrice[y,fz,cr] for cr in md.ncrop)
                                        - sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for cul in md.nculture)
                                        - sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz]  for t in md.ntime for cul in md.nculture if md.t_year[t]==y and md.aIrrigation[md.fzone_type[fz]]==1)
                                        for y in md.nyear} 
                                   for fz in md.nfzone}    
        
        #CROP TRANSPORT
        oo['CropTransport_cm_cm'] = {cma:{cmb:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value                      for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.aTransIn[ct]==cma and md.aTransOut[ct]==cmb)/leni(md.nyear) for cmb in md.ncmarket} for cma in md.ncmarket}
        oo['CropTransport_ct_cr'] = {ct: {cr: sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value                      for y in md.nyear)/leni(md.nyear) for cr in md.ncrop} for ct in md.nctrans}
        oo['CropTransport_ct']    = {ct: sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value                           for y in md.nyear for cr in md.ncrop)/leni(md.nyear) for ct in md.nctrans}
        oo['CropTransVal_ct']     = {ct: sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for y in md.nyear for cr in md.ncrop)/leni(md.nyear) for ct in md.nctrans}
        oo['CropLoss_cm_cm']      = {cma:{cmb:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * md.aTransLoss[ct]  for y in md.nyear for ct in md.nctrans for cr in md.ncrop if md.aTransIn[ct]==cma and md.aTransOut[ct]==cmb)/leni(md.nyear) for cmb in md.ncmarket} for cma in md.ncmarket}                                              
        oo['CropImpShare_cr_co']  = {cr: {co:(sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * (1-md.aTransLoss[ct]) for y in md.nyear for ct in md.nctrans if md.cmarket_country[md.aTransOut[ct]]==co)
                                            -sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value                         for y in md.nyear for ct in md.nctrans if md.cmarket_country[md.aTransIn[ct]]==co))
                                            /sum(md.aCropDem[y,cm,cr]                           for y in md.nyear for cm in md.ncmarket if md.cmarket_country[cm]==co)
                                        for co in md.ncountry if sum(md.aCropDem[y,cm,cr] for y in md.nyear for cm in md.ncmarket if md.cmarket_country[cm]==co) != 0} 
                                    for cr in md.ncrop}
        
        #EXTERNAL MARKETS
        oo['CropExtProd_cr_co']  = {cr:{co:sum(md.kt_to_Mt*md.AcEXTPROD[y,cm,cr].value for y in md.nyear for cm in md.nextcmarket if md.cmarket_country[cm]==co)/leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}
        
        #CROP PRODUCTION
        oo['CropProduction_cr_fz']   = {cr:{fz:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value for y in md.nyear)/leni(md.nyear) for fz in md.nfzone} for cr in md.ncrop} 
        oo['CropProduction_cr_co']   = {cr:{co:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}
        
        
        #YIELDS
        if md.Options['Crop choice'] == 'fixed' or md.Options['Yield water response'] == 'nonlinear': #ADD: the two options are not really the equivalent - as option 1 does not differentiate rainfed and irrigated agriculture
            #ADD: this is not relative yield but absolute yield [t/ha] (in the excel file it is indicated as in [%])
            oo['Yield_cul_ft']       = {cr:{co:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co)/sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co and md.culture_crop[cul]==cr)
                                            if sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co and md.culture_crop[cul]==cr) != 0  else 0 for co in md.ncountry} for cr in md.ncrop}
        else:    
            oo['Yield_cul_ft']       = {cul:{ft:sum(md.kha_to_Mha*md.AlCULAREA[y,fz,kfd,kypt].value*(1-sum(md.akY[md.field_culture[kfd],kyps]*(1-md.aYieldMat[kypt,kyps]) for kyps in md.nyphase)) for y in md.nyear for fz in md.nfzone for kypt in md.nypath for kfd in md.nfieldculture if md.field_culture[kfd]==cul and md.fzone_type[fz]==ft)
                                              /sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone if md.fzone_type[fz]==ft) 
                                               for ft in md.nftype if sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone if md.fzone_type[fz]==ft) != 0} 
                                               for cul in md.nculture}
        
        #LAND USE
        if md.Options['Crop choice'] == 'fixed' or md.Options['Yield water response'] == 'nonlinear':
            oo['LandUse_fd_ft']      = {cul:{ft:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone if md.fzone_type[fz]==ft)/leni(md.nyear) for ft in md.nftype} for cul in md.nculture}
        else:
            oo['LandUse_fd_ft']      = {fd[0]:{ft:sum(md.kha_to_Mha*md.AlCULAREA[y,fz,kfd,kypt].value for y in md.nyear for fz in md.nfzone for kfd in md.nfieldculture for kypt in md.nypath if md.fzone_type[fz]==ft and kfd[0]==fd[0])/leni(md.nyear) for ft in md.nftype} for fd in md.nfieldculture}
        oo['LandUse_cul_fz']         = {cul:{fz:sum(GROSSAREA[y,fz,cul] for y in md.nyear)/leni(md.nyear) for fz in md.nfzone} for cul in md.nculture}
        oo['LandUse_cr_co']          = {cr:{co:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co and md.culture_crop[cul]==cr)/leni(md.nyear) for co in md.ncountry} for cr in md.ncrop}
        oo['LandUse_cr_y']           = {cr:{y: sum(GROSSAREA[y,fz,cul] for fz in md.nfzone for cul in md.nculture if md.culture_crop[cul]==cr) for y in md.nyear} for cr in md.ncrop}
        oo['LandUse_co_ft']          = {co:{ft:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co and md.fzone_type[fz]==ft)/leni(md.nyear) for ft in md.nftype} for co in md.ncountry}
        oo['LandUse_co_irr']         = {co:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_country[fz]==co and md.aIrrigation[md.fzone_type[fz]]==1)/leni(md.nyear) for co in md.ncountry}
        oo['LandUse_cr_ft']          = {cr:{ft:sum(GROSSAREA[y,fz,cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.culture_crop[cul]==cr and md.fzone_type[fz]==ft)/leni(md.nyear) for ft in md.nftype} for cr in md.ncrop}
        oo['LandUse_ca']   = {'Irrig Land Use [%]'  :{ca:sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_catch[fz]==ca and md.aIrrigation[md.fzone_type[fz]] == 1)/leni(md.nyear)/ReadParam('CatchArea')[ca] for ca in md.ncatch},
                             'Irrig Used Cap [%]'  :{ca: sum(NETAREA[y,fz] for y in md.nyear for fz in md.nfzone if md.fzone_catch[fz]==ca and md.aIrrigation[md.fzone_type[fz]] == 1)
                                                        /sum(md.kha_to_Mha*md.aLandCap[y,fz] for y in md.nyear for fz in md.nfzone if md.aIrrigation[md.fzone_type[fz]]==1 and md.fzone_catch[fz]==ca) for ca in md.ncatch 
                                                        if sum(md.aLandCap[y,fz] for y in md.nyear for fz in md.nfzone if md.aIrrigation[md.fzone_type[fz]]==1 and md.fzone_catch[fz]==ca) != 0}}        
        
        #%%PRINTS
        if PRINT ==1:
            print('*************AGRICULTURE*************')
            print('Available Land [Mha]')
            print(oo['CropBalance']['Available land [Mha]'])
            print('Average cultivated area [Mha]')
            print(oo['CropBalance']['Net Cultivated land [Mha]'])
            print('Crop production [Mt]')
            print(oo['CropBalance']['Crop production [Mt]'])
            print('Net traded crop [Mt]')
            print(oo['CropBalance']['Crop net import [Mt]'])
            print('Crop transport losses [Mt]')
            print(oo['CropBalance']['Crop gross export [Mt]']-oo['CropBalance']['Crop net import [Mt]'])
            print('Crop external production [Mt]')
            print(sum(oo['CropExtProd_cr_co'][cr][co] for cr in md.ncrop for co in md.ncountry))
            print('Crop demand [Mt]')
            print(oo['CropBalance']['Crop demand [Mt]'])
            print('Crop balance [Mt]: CropProd = CropTrans - CropLoss')
            print(oo['CropBalance'][' Balance = 0 [Mt]'])
        

#%% ENERGY  
        #Hydropower plants
        #ADD: WARNING CASE SPECIFIC secondary power VALUE (HppVal/3)
        FIRMPOWER = 95/100 # Percentage of the time the power has to be generated to be firm
 
        oo['HydropowerTable']    = {'Power Prod [GWh/y]'    : {hp:sum(md.se.value*md.EeHPPROD[t,pld,hp].value  for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp} ,
                                   'Firm Power [GWh/month]' : {hp:sorted([sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)*(1-FIRMPOWER)))] for hp in md.nhpp},
                                   'Firm value [M$/year]'   : {hp:sorted([sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)*(1-FIRMPOWER)))]*leni(md.nmonth)*ReadParam('eHppVal')[hp]
                                                                   + (sum(md.se.value*md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload)/leni(md.nyear)-sorted([sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload) for t in md.ntime])[int(round(leni(md.ntime)*(1-FIRMPOWER)))]*leni(md.nmonth))*ReadParam('eHppVal')[hp]/3
                                                                   for hp in md.nhpp},
                                   'Capacity [MW]'          : ReadParam('eHppCap'), 
                                   'Discharge [Mm3/y]'      : {hp:sum(md.EwHPDISCHARGE[t,pld,hp].value  for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp},
                                   'AvStorage [Mm3/month]'  : {hp:sum(md.sw.value*md.WwRSTORAGE[t,res].value for t in md.ntime for res in md.nres if md.hp_res[hp]==res)/leni(md.ntime) for hp in md.nhpp},
                                   'Operation costs [M$/y]' : {hp:sum(md.eHppCost[hp] * md.se.value*md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp},
                                   'Capacity shadowprice [$/GWh/month]':oo['HpCapShadow_hp']}
        
        oo['HpProd_hp_m']            = {m:{hp:sum(md.se.value*md.EeHPPROD[t,pld,hp].value        for t in md.ntime for pld in md.npload  if md.t_month[t]==m)                        /leni(md.nyear) for hp in md.nhpp}                           for m in md.nmonth}                      
        oo['HpDis_hp_m']             = {m:{hp:sum(md.EwHPDISCHARGE[t,pld,hp].value   for t in md.ntime for pld in md.npload  if md.t_month[t]==m)                        /leni(md.nyear) for hp in md.nhpp if md.hp_res[hp]!='ROR'}   for m in md.nmonth}
        oo['HpStor_hp_m']            = {m:{hp:sum(md.sw.value*md.WwRSTORAGE[t,res].value         for t in md.ntime for res in md.nres    if md.hp_res[hp]==res and md.t_month[t]==m) /leni(md.nyear) for hp in md.nhpp}                           for m in md.nmonth}       
        oo['HpProd_hp_y']            = {y:{hp:sum(md.se.value*md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload if md.t_year[t]==y) for hp in md.nhpp} for y in md.nyear}        
        if md.Options['Energy market'] == 1:
            oo['HydropowerTable']['Energy value [M$/y]']     = {hp:sum(-md.se.value*md.EeHPPROD[t,pld,hp].value * 1/md.se.value*md.dual[md.engy_balance[t,pld,md.hp_pmarket[hp]]] for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp}
            oo['HpBenefits_hp_y']    = {y:{hp:sum(md.se.value*md.EeHPPROD[t,pld,hp].value * (-1/md.se.value*md.dual[md.engy_balance[t,pld,md.hp_pmarket[hp]]] - md.eHppCost[hp]) for t in md.ntime for pld in md.npload if md.t_year[t]==y) for hp in md.nhpp} for y in md.nyear}
        else:
            oo['HydropowerTable']['Energy value [M$/y]'] = {hp:sum(-md.se.value*md.EeHPPROD[t,pld,hp].value * md.eHppVal[hp] for t in md.ntime for pld in md.npload)/leni(md.nyear) for hp in md.nhpp}
            oo['HpBenefits_hp_y']    = {y:{hp:sum(md.se.value*md.EeHPPROD[t,pld,hp].value * (md.eHppVal[hp]- md.eHppCost[hp]) for t in md.ntime for pld in md.npload if md.t_year[t]==y) for hp in md.nhpp} for y in md.nyear}
                        
        #Energy market prices
        oo['EnergyVal_pm_m']     = {pm:{m:sum(-1/md.se.value*md.dual[md.engy_balance[t,pld,pm]]*md.eLoadDem[pld] for t in md.ntime for pld in md.npload if md.t_month[t]==m)/leni(md.nyear) for m in md.nmonth} for pm in md.npmarket}
        
        #Energy Balance
        oo['EnergyBalance_t_pm'] = {'Power demand [GWh]'            :{t:{pm:md.eEngyDem[md.t_year[t],md.t_month[t],pm] for pm in md.npmarket} for t in md.ntime},
                                   'Unserved demand [GWh]'          :{t:{pm:md.eEngyDem[md.t_year[t],md.t_month[t],pm]-sum(md.se.value*md.EeSUPPLY[t,pld,pm].value for pld in md.npload) for pm in md.npmarket} for t in md.ntime},
                                   'Hydropower production [GWh]'    :{t:{pm:sum(md.se.value*md.EeHPPROD[t,pld,hp].value for pld in md.npload for hp in md.nhpp if md.hp_pmarket[hp]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Power transmission loss [GWh]'  :{t:{pm:sum(md.se.value*md.EeTRANS[t,pld,tl].value * md.eTransLoss[tl] for pld in md.npload for tl in md.ntransline if md.eTransIn[tl]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Power supply loss [GWh]'        :{t:{pm:sum(md.se.value*md.EeSUPPLY[t,pld,pm].value * md.eSupLoss[pm] / (1-md.eSupLoss[pm]) for pld in md.npload) for pm in md.npmarket} for t in md.ntime},
                                   'Net power export [GWh]'         :{t:{pm:sum(md.se.value*md.EeTRANS[t,pld,tl].value * (1-md.eTransLoss[tl]) for pld in md.npload for tl in md.ntransline if md.eTransIn[tl]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Net power import [GWh]'         :{t:{pm:sum(md.se.value*md.EeTRANS[t,pld,tl].value * (1-md.eTransLoss[tl]) for pld in md.npload for tl in md.ntransline if md.eTransOut[tl]==pm) for pm in md.npmarket} for t in md.ntime},
                                   'Energy price [$/kWh]'           :{t:{pm:sum(-1/md.se.value*md.dual[md.engy_balance[t,pld,pm]]*md.eLoadDem[pld] for pld in md.npload) for pm in md.npmarket} for t in md.ntime}}
        
        #Power Plants
        oo['EnergyBalance_t_pm']['Power plant prod [GWh]'] = {t:{pm:sum(md.se.value*md.EeOPPROD[t,pld,pp].value for pld in md.npload for pp in md.nopp if md.op_pmarket[pp]==pm) for pm in md.npmarket} for t in md.ntime}           
        oo['PowerplantTable']= {'Power Prod [GWh/y]'        : {pp:sum(md.se.value*md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for pp in md.nopp},
                               'Capacity [MW]'              : ReadParam('eOppCap'), 
                               'Fuel use [GWh/y]'           : {pp:sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp] for t in md.ntime for pld in md.npload)/leni(md.nyear) if md.Options['Fuels']==1 else 0 for pp in md.nopp},
                               'Operation costs [M$/y]'     : {pp:sum(md.eOppCost[pp] * md.se.value*md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for pp in md.nopp},
                               'Energy value [M$/y]'        : {pp:sum(-md.se.value*md.EeOPPROD[t,pld,pp].value * 1/md.se.value*md.dual[md.engy_balance[t,pld,md.op_pmarket[pp]]] for t in md.ntime for pld in md.npload)/leni(md.nyear) for pp in md.nopp},
                               'Capacity shadowprice [$/kWh]':oo['PpCapShadow_pp']}
        oo['PowerplantTable']['Fuel cost [M$/y]']     = {pp:sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp] * md.eFuelCost[md.t_year[t],md.op_pmarket[pp],fu] for t in md.ntime for pld in md.npload for fu in md.nfuel if md.op_fuel[pp]==fu)/leni(md.nyear) for pp in md.nopp}
        oo['PpProd_pp_m']     = {m:{pp:sum(md.se.value*md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload  if md.t_month[t]==m)/leni(md.nyear) for pp in md.nopp} for m in md.nmonth}
        oo['PpProd_pp_y']     = {y:{pp:sum(md.se.value*md.EeOPPROD[t,pld,pp].value for t in md.ntime for pld in md.npload if md.t_year[t]==y) for pp in md.nopp} for y in md.nyear}
        oo['PpBenefits_pp_y'] = {y:{pp:sum(md.se.value*md.EeOPPROD[t,pld,pp].value * (-1/md.se.value*md.dual[md.engy_balance[t,pld,md.op_pmarket[pp]]] - md.eOppCost[pp] - (md.Options['Fuels'] and md.eFuelCost[y,md.op_pmarket[pp],md.op_fuel[pp]]/md.eOppEff[pp])) for t in md.ntime for pld in md.npload if md.t_year[t]==y) for pp in md.nopp} for y in md.nyear}

        #Generic capacity
        oo['EnergyBalance_t_pm']['Gen power prod [GWh]'] = {t:{pm:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value for pld in md.npload for pt in md.nptech) for pm in md.npmarket} for t in md.ntime}       
        oo['EnergyBalance_t_pm']['Gen power cap [MW]']   = {t:{pm:sum(md.se.value*md.EeGENCAP[md.t_year[t],pt,pm].value for pt in md.nptech)/leni(md.nmonth) for pm in md.npmarket} for t in md.ntime}  
        for pt in md.nptech:
            oo['EnergyBalance_t_pm'][pt+' power prod [GWh]'] = {t:{pm:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value for pld in md.npload)+sum(md.se.value*md.EeOPPROD[t,pld,pp].value for pld in md.npload for pp in md.nopp if md.op_ptech[pp]==pt and md.op_pmarket[pp]==pm) for pm in md.npmarket} for t in md.ntime}       
        oo['PowerGenCap_pt_pm']  = {pt: {pm: sum(md.se.value*md.EeGENCAP[y,pt,pm].value for y in md.nyear) for pm in md.npmarket} for pt in md.nptech}
        oo['PowerGenProd_pt_pm'] = {pt: {pm: sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value for t in md.ntime for pld in md.npload)/leni(md.nyear) for pm in md.npmarket} for pt in md.nptech}
        oo['PowerGenCap_pt_y']   = {pt: {y:  sum(md.se.value*md.EeGENCAP[ky,pt,pm].value for ky in md.nyear for pm in md.npmarket if ky <= y and ky >= y-md.eLifeTime[pt,pm]) for y in md.nyear} for pt in md.nptech}
        oo['PowerGenProd_pt_y']  = {pt: {y:  sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear} for pt in md.nptech}

        #Transmission lines
        oo['TmsLineTable']  = {'Net transmission [GWh/y]' : {tl:sum(md.se.value*md.EeTRANS[t,pld,tl].value * (1-md.eTransLoss[tl]) for t in md.ntime for pld in md.npload)/leni(md.nyear) for tl in md.ntransline},
                              'Power loss [GWh/y]'        : {tl:sum(md.se.value*md.EeTRANS[t,pld,tl].value * md.eTransLoss[tl] for t in md.ntime for pld in md.npload)/leni(md.nyear) for tl in md.ntransline},
                              'Cap. Shadow [$/GWh/y]'     : {tl:sum(-md.dual[md.engy_transmissioncap[t,pld,tl]]*md.eLoadTime[pld] for t in md.ntime for pld in md.npload)/leni(md.ntime) * leni(md.nmonth) for tl in md.ntransline}}             

        #Power markets and Energy balances
        oo['EnergyBalance_t_pm']['Balance = 0 [GWh]'] = {t:{pm:+oo['EnergyBalance_t_pm']['Power demand [GWh]'][t][pm] 
                                                              - oo['EnergyBalance_t_pm']['Hydropower production [GWh]'][t][pm]
                                                              -(oo['EnergyBalance_t_pm']['Power plant prod [GWh]'][t][pm] if md.Options['Power plants'] == 1 else 0)
                                                              -(oo['EnergyBalance_t_pm']['Gen power prod [GWh]'][t][pm] if md.Options['Power technologies'] == 1 else 0)
                                                              - oo['EnergyBalance_t_pm']['Unserved demand [GWh]'][t][pm]
                                                              + oo['EnergyBalance_t_pm']['Power supply loss [GWh]'][t][pm]
                                                              + oo['EnergyBalance_t_pm']['Power transmission loss [GWh]'][t][pm]
                                                              + oo['EnergyBalance_t_pm']['Net power export [GWh]'][t][pm]
                                                              - oo['EnergyBalance_t_pm']['Net power import [GWh]'][t][pm]                                                                  
                                                              for pm in md.npmarket} for t in md.ntime}

        oo['EnergyBalance']      = {keys:sum(oo['EnergyBalance_t_pm'][keys][t][pm] for t in md.ntime for pm in md.npmarket)/leni(md.nyear) for keys in oo['EnergyBalance_t_pm'].keys()}
        oo['EnergyBalance_co']   = {keys:{co:sum(oo['EnergyBalance_t_pm'][keys][t][pm] for t in md.ntime for pm in md.npmarket if md.pmarket_country[pm]==co)/leni(md.nyear) for co in md.ncountry} for keys in oo['EnergyBalance_t_pm'].keys()}                       
        oo['EnergyBalance_m']    = {keys:{m:sum(oo['EnergyBalance_t_pm'][keys][t][pm] for t in md.ntime for pm in md.npmarket if md.t_month[t]==m)/leni(md.nyear) for m in md.nmonth} for keys in oo['EnergyBalance_t_pm'].keys()}    
        oo['EnergyBalance_y']    = {keys:{y:sum(oo['EnergyBalance_t_pm'][keys][t][pm] for t in md.ntime for pm in md.npmarket if md.t_year[t]==y) for y in md.nyear} for keys in oo['EnergyBalance_t_pm'].keys()}
        #ADD: case specific
        oo['EnergyBalance_y2']   = {keys:{y:sum(oo['EnergyBalance_t_pm'][keys][t][pm] for t in md.ntime for pm in md.npmarket if md.t_year[t]==y and pm != 'PmSouthAfrica') for y in md.nyear} for keys in oo['EnergyBalance_t_pm'].keys()}
        #Correct intensive values (energy price)
        oo['EnergyBalance']['Energy price [$/kWh]']  = oo['EnergyBalance']['Energy price [$/kWh]']/leni(md.nmonth)/leni(md.npmarket)
        oo['EnergyBalance_co']['Energy price [$/kWh]']={co: oo['EnergyBalance_co']['Energy price [$/kWh]'][co]/leni(md.nmonth)/max(1,sum(1 for pm in md.npmarket if md.pmarket_country[pm]==co)) for co in md.ncountry}
        oo['EnergyBalance_m']['Energy price [$/kWh]']= {m : oo['EnergyBalance_m']['Energy price [$/kWh]'][m]/leni(md.npmarket) for m in md.nmonth}
        oo['EnergyBalance_y']['Energy price [$/kWh]']= {y : oo['EnergyBalance_y']['Energy price [$/kWh]'][y]/leni(md.npmarket)/leni(md.nmonth) for y in md.nyear}
        oo['EnergyBalance_y2']['Energy price [$/kWh]']= {y : oo['EnergyBalance_y2']['Energy price [$/kWh]'][y]/max(1,leni(md.npmarket)-1)/leni(md.nmonth) for y in md.nyear}
        #Generic power tech
        for pt in md.nptech:
            oo['EnergyBalance_y'][pt+' production [GWh]']={y:oo['PowerGenProd_pt_y'][pt][y] for y in md.nyear}
        
        #Catchment scale energy balance
        oo['EnergyBalance_ca']={'Hydropower production [GWh]' : {c: sum(md.se.value*md.EeHPPROD[t,pld,hp].value for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_catch[hp]==c)/leni(md.nyear) for c in md.ncatch}}

        #%%PRINTS
        if PRINT ==1:
            print('*************ENERGY*************')            
            print('Hydropower production [GWh]')
            print(oo['EnergyBalance']['Hydropower production [GWh]'])
            print('Power demand [GWh]')
            print(oo['EnergyBalance']['Power demand [GWh]'])
            print('Unserved power demand [GWh]')
            print(oo['EnergyBalance']['Unserved demand [GWh]'])
            print('Net power transmition [GWh]')
            print(oo['EnergyBalance']['Net power import [GWh]'])
            print('Power transmition losses [GWh]')
            print(oo['EnergyBalance']['Power transmission loss [GWh]'])
            print('Energy balance [GWh]: PowerDem = HpProd + PpProd - PowerTransLoss - UnsPowerDem')
            print(oo['EnergyBalance']['Balance = 0 [GWh]'])
            print('PowerPlant production [GWh]')
            print(oo['EnergyBalance']['Power plant prod [GWh]'])
            print('Generic capacity investment [kWh/month]')
            print(oo['EnergyBalance']['Gen power cap [MW]'])
            print('Generic capacity production [GWh]')
            print(oo['EnergyBalance']['Gen power prod [GWh]'])
                    
#%%ECONOMICS
        #Economic balance
        oo['EconomicBalance_y_co'] = {'User benefits [M$]':      {y:{co:sum(md.wUserVal[u]*md.WwSUPPLY[t,u].value                       for t in md.ntime for u in md.nuser if md.user_country[u]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear},
                                'Water supply costs [M$]':       {y:{co:sum(md.WwSUPPLY[t,u].value/(1-md.wUserLoss[u]) * md.wSupCost[u] for t in md.ntime for u in md.nuser if md.user_country[u]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear},
                                'Crop cultivation costs [M$]':   {y:{co:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for fz in md.nfzone for cul in md.nculture if md.fzone_country[fz]==co) for co in md.ncountry} for y in md.nyear},
                                'Crop irrigation costs [M$]':    {y:{co:sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_country[fz]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear},
                                'Crop groundwater costs [M$]':   {y:{co:sum(md.AwGWSUPPLY[t,fz,cul].value * md.wGwCost[aq] for t in md.ntime for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]== md.aqui_catch[aq] and md.fzone_country[fz]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}}
        #Crop market
        oo['EconomicBalance_y_co']['Crop benefits [M$]']         = {y:{co:sum(md.kt_to_Mt*md.AcSUPPLY[y,cm,cr,cds].value * md.aStepVal[cm,cr,cds]*md.aCropVal[cm,cr] for cm in md.ncmarket for cr in md.ncrop for cds in md.ncdstep if md.cmarket_country[cm]==co) for co in md.ncountry} for y in md.nyear}
        if md.Options['Crop market'] == 0:
            oo['EconomicBalance_y_co']['Crop benefits [M$]']     = {y:{co:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * md.aFarmVal[md.fzone_country[fz],cr] for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Crop transport costs [M$]']  = {y:{co:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * md.aTransCost[ct,cr]                               for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Crop import costs [M$]']     = {y:{co:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransOut[ct]]==co) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Crop export benefits [M$]']  = {y:{co:sum(md.kt_to_Mt*md.AcTRANS[y,ct,cr].value * -1/md.kt_to_Mt*md.dual[md.agr_cropbalance[y,md.aTransIn[ct],cr]] for ct in md.nctrans for cr in md.ncrop if md.cmarket_country[md.aTransIn[ct]]==co) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Crop ext prod costs [M$]']   = {y:{co:sum(md.kt_to_Mt*md.AcEXTPROD[y,cm,cr].value * md.aCropVal[cm,cr] for cm in md.nextcmarket for cr in md.ncrop if md.cmarket_country[cm]==co) for co in md.ncountry} for y in md.nyear}
        
        #Energy markets       
        oo['EconomicBalance_y_co']['Energy benefits [M$]']       = {y:{co:sum(md.se.value*md.EeSUPPLY[t,pld,pm].value * md.eEngyVal[md.t_month[t],pm] for t in md.ntime for pld in md.npload for pm in md.npmarket if md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        if md.Options['Energy market'] == 0:
            oo['EconomicBalance_y_co']['Energy benefits [M$]']   = {y:{co: sum(md.se.value*md.EeHPPROD[t,pld,hp].value*(md.eHppVal[hp]-md.eHppCost[hp]) for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_country[hp]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy import costs [M$]']   = {y:{co:sum(md.se.value*md.EeTRANS[t,pld,tl].value * -1/md.se.value*md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy export benefits [M$]']= {y:{co:sum(md.se.value*md.EeTRANS[t,pld,tl].value * -1/md.se.value*md.dual[md.engy_balance[t,pld,md.eTransIn[tl]]] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransIn[tl]]==co  and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy trans costs [M$]']    = {y:{co:sum(md.se.value*md.EeTRANS[t,pld,tl].value * md.eTransCost[tl] for t in md.ntime for pld in md.npload for tl in md.ntransline if md.pmarket_country[md.eTransOut[tl]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy hydro O&M costs [M$]']= {y:{co:sum(md.se.value*md.EeHPPROD[t,pld,hp].value * md.eHppCost[hp] for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_country[hp]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy curtailment costs [M$]']={y:{co:sum(md.eEngyVal[md.t_month[t],pm]*(md.eEngyDem[md.t_year[t],md.t_month[t],pm]-sum(md.se.value*md.EeSUPPLY[t,pld,pm].value for pld in md.npload)) for t in md.ntime for pm in md.npmarket if md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        #Other power plants
        oo['EconomicBalance_y_co']['Energy power O&M costs [M$]']= {y:{co:sum(md.se.value*md.EeOPPROD[t,pld,pp].value*md.eOppCost[pp] for t in md.ntime for pld in md.npload for pp in md.nopp if md.pmarket_country[md.op_pmarket[pp]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy fuel costs [M$]']     = {y:{co:sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCost[y,md.op_pmarket[pp],fu] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        #Generic power technologies
        oo['EconomicBalance_y_co']['Energy Gen power O&M costs [M$]'] = {y:{co:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value*md.eVarOPEX[pt,pm] for t in md.ntime for pld in md.npload for pm in md.npmarket for pt in md.nptech if md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy Gen cap costs [M$]']  = {y:{co:sum(md.se.value*md.EeGENCAP[y,pt,pm].value*md.eCAPEX[y,pt,pm]*min(1,(md.t_year[md.Options['tfin']]-y+1)/md.eLifeTime[pt,pm]) + sum(md.se.value*md.EeGENCAP[ky,pt,pm].value for ky in md.nyear if ky <= y and ky >= y-md.eLifeTime[pt,pm])*md.eFixOPEX[pt,pm] for pm in md.npmarket for pt in md.nptech if md.pmarket_country[pm]==co) for co in md.ncountry} for y in md.nyear}
        oo['EconomicBalance_y_co']['Energy Gen fuel costs [M$]'] = {y:{co:sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCost[y,pm,fu] for t in md.ntime for pld in md.npload for pm in md.npmarket for pt in md.nptech for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co and md.t_year[t]==y) for co in md.ncountry} for y in md.nyear}                   
        #CO2 emissions
        oo['EconomicBalance_y_co']['Energy CO2 costs [M$]']      = {y:{co: sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value/md.eTechEff[pt,pm]*md.eFuelCO2[fu]*md.eCO2Val[y,pm] for t in md.ntime for pld in md.npload for pm in md.npmarket for pt in md.nptech for fu in md.nfuel if md.ptech_fuel[pt,pm]==fu and md.pmarket_country[pm]==co and md.t_year[t]==y)
                                                                         +sum(md.se.value*md.EeOPPROD[t,pld,pp].value/md.eOppEff[pp]*md.eFuelCO2[fu]*md.eCO2Val[y,md.op_pmarket[pp]] for t in md.ntime for pld in md.npload for pp in md.nopp for fu in md.nfuel if md.op_fuel[pp]==fu and md.pmarket_country[md.op_pmarket[pp]]==co and md.t_year[t]==y)
                                                                      for co in md.ncountry} for y in md.nyear}                
        #Balance
        oo['EconomicBalance_y_co'][' Balance [M$]'] = {y:{co:+ oo['EconomicBalance_y_co']['User benefits [M$]'][y][co]
                                                            - oo['EconomicBalance_y_co']['Water supply costs [M$]'][y][co]
                                                            + oo['EconomicBalance_y_co']['Crop benefits [M$]'][y][co]
                                                            - oo['EconomicBalance_y_co']['Crop cultivation costs [M$]'][y][co]
                                                            - oo['EconomicBalance_y_co']['Crop irrigation costs [M$]'][y][co]
                                                            - oo['EconomicBalance_y_co']['Crop groundwater costs [M$]'][y][co]       
                                                            - oo['EconomicBalance_y_co']['Crop transport costs [M$]'][y][co]       
                                                            - oo['EconomicBalance_y_co']['Crop import costs [M$]'][y][co]            
                                                            + oo['EconomicBalance_y_co']['Crop export benefits [M$]'][y][co]         
                                                            - oo['EconomicBalance_y_co']['Crop ext prod costs [M$]'][y][co]
                                                            + oo['EconomicBalance_y_co']['Energy benefits [M$]'][y][co]
                                                            - oo['EconomicBalance_y_co']['Energy import costs [M$]'][y][co]         
                                                            + oo['EconomicBalance_y_co']['Energy export benefits [M$]'][y][co]       
                                                            - oo['EconomicBalance_y_co']['Energy trans costs [M$]'][y][co]           
                                                            - oo['EconomicBalance_y_co']['Energy hydro O&M costs [M$]'][y][co]    
                                                            - oo['EconomicBalance_y_co']['Energy curtailment costs [M$]'][y][co]    
                                                            - oo['EconomicBalance_y_co']['Energy power O&M costs [M$]'][y][co]      
                                                            - oo['EconomicBalance_y_co']['Energy fuel costs [M$]'][y][co]            
                                                            - oo['EconomicBalance_y_co']['Energy Gen power O&M costs [M$]'][y][co]   
                                                            - oo['EconomicBalance_y_co']['Energy Gen cap costs [M$]'][y][co]       
                                                            - oo['EconomicBalance_y_co']['Energy Gen fuel costs [M$]'][y][co]      
                                                            - oo['EconomicBalance_y_co']['Energy CO2 costs [M$]'][y][co]               
                                                            for co in md.ncountry} for y in md.nyear}
        #Balance at other levels
        oo['EconomicBalance']      = {keys:     sum(oo['EconomicBalance_y_co'][keys][y][co] for y in md.nyear for co in md.ncountry)/leni(md.nyear) for keys in oo['EconomicBalance_y_co'].keys()}
        oo['EconomicBalance_co']   = {keys: {co:sum(oo['EconomicBalance_y_co'][keys][y][co] for y in md.nyear)/leni(md.nyear) for co in md.ncountry} for keys in oo['EconomicBalance_y_co'].keys()}
        oo['EconomicBalance_y']    = {keys: {y :sum(oo['EconomicBalance_y_co'][keys][y][co] for co in md.ncountry) for y in md.nyear} for keys in oo['EconomicBalance_y_co'].keys()}

        oo['EconomicBalance_ca']={'User benefits [M$]':           {c:sum(md.WwSUPPLY[t,u].value*md.wUserVal[u] for t in md.ntime for u in md.nuser if md.user_catch[u]==c)/leni(md.nyear) for c in md.ncatch},
                                 'Water supply costs [M$]':      {c:sum(md.WwSUPPLY[t,u].value/(1-md.wUserLoss[u]) * md.wSupCost[u] for t in md.ntime for u in md.nuser if md.user_catch[u]==c)/leni(md.nyear) for c in md.ncatch},                                  
                                 'Crop cultivation costs [M$]':  {c:sum(GROSSAREA[y,fz,cul] * md.aCulCost[md.fzone_type[fz],cul] for y in md.nyear for fz in md.nfzone for cul in md.nculture if md.fzone_catch[fz]==c)/leni(md.nyear)     for c in md.ncatch}, 
                                 'Crop irrigation costs [M$]':   {c:sum(md.AwSUPPLY[t,fz,cul].value * md.aIrrgCost[fz] for t in md.ntime for fz in md.nifzone for cul in md.nculture if md.fzone_catch[fz]==c)/leni(md.nyear)      for c in md.ncatch}}
        if md.Options['Crop market'] == 0:
            oo['EconomicBalance_ca']['Crop benefits [M$]']           = {c:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * md.aFarmVal[md.fzone_country[fz],cr] for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_catch[fz]==c)/leni(md.nyear) for c in md.ncatch}

        oo['EconomicBalance_ca']['Crop groundwater costs [M$]']      = {c:sum(md.AwGWSUPPLY[t,fz,cul].value * md.wGwCost[md.fzone_catch[fz]]  for t in md.ntime for fz in md.nifzone for cul in md.nculture for aq in md.naquifer if md.fzone_catch[fz]==md.aqui_catch[aq])/leni(md.nyear) for c in md.ncatch}
        if md.Options['Energy market'] == 0:
            oo['EconomicBalance_ca']['Energy benefits [M$]']         = {c:sum(md.se.value*md.EeHPPROD[t,pld,hp].value*(md.eHppVal[hp]-md.eHppCost[hp]) for t in md.ntime for pld in md.npload for hp in md.nhpp if md.hp_catch[hp]==c)/leni(md.nyear) for c in md.ncatch} 
        #%%PRINTS    
        if PRINT==1:
            print('*************ECONOMIC*************')
            print('User benefit [M$]')
            print(oo['EconomicBalance']['User benefits [M$]'])
            print('Crop benefit [M$]')
            print(oo['EconomicBalance']['Crop benefits [M$]'])
            print('Cultivation cost [M$]')
            print(oo['EconomicBalance']['Crop cultivation costs [M$]'])
            print('Irrigation cost [M$]')
            print(oo['EconomicBalance']['Crop irrigation costs [M$]'])
            print('Crop import cost [M$]')
            print(oo['EconomicBalance']['Crop import costs [M$]'])
            print('Crop export benefits [M$]')
            print(oo['EconomicBalance']['Crop export benefits [M$]'])
            print('Crop transport cost [M$]')
            print(oo['EconomicBalance']['Crop transport costs [M$]'])
            print('Energy benefit [M$]')
            print(oo['EconomicBalance']['Energy benefits [M$]'])            
            print('Energy curtailment cost [M$]')
            print(oo['EconomicBalance']['Energy curtailment costs [M$]'])
            print('Hydropower operational cost [M$]')
            print(oo['EconomicBalance']['Energy hydro O&M costs [M$]'])
            print('Other Power plants operational cost [M$]')
            print(oo['EconomicBalance']['Energy power O&M costs [M$]'])
            print('Other Power plants fuel cost [M$]')
            print(oo['EconomicBalance']['Energy fuel costs [M$]'])
            print('Generic capacity construction costs [M$]')
            print(oo['EconomicBalance']['Energy Gen cap costs [M$]'])
            print('Generic capacity operation costs [M$]')
            print(oo['EconomicBalance']['Energy Gen power O&M costs [M$]'])
            print('Generic fuel cost [M$]')
            print(oo['EconomicBalance']['Energy Gen fuel costs [M$]'])
#%% YEARLY TRENDS - exclude South Africa - #ADD CASE SPECIFIC !!!
#Capacity
        YearlyEnergy={'Power demand [GWh/y]': oo['EnergyBalance_y2']['Power demand [GWh]'],
                      'Unserved demand [GWh/y]': oo['EnergyBalance_y2']['Unserved demand [GWh]'],
                      'Power price [$/kWh]': oo['EnergyBalance_y2']['Energy price [$/kWh]'],
                      'Power trade SA [GWh/y]': {y:oo['EnergyBalance_y2']['Net power export [GWh]'][y]-oo['EnergyBalance_y2']['Net power import [GWh]'][y] for y in md.nyear},
                                           
                      }
        YearlyProd={'Hydropower [GWh/y]': oo['EnergyBalance_y2']['Hydropower production [GWh]'] }
        YearlyCapacity={'Hydropower [MW]': {y:  sum(md.eHppCap[hp] for hp in md.nhpp) 
                                              +(sum(md.IbINVEST[ip,inv].value*md.iInvCap[inv] for ip in md.ninvphase for inv in md.ninvest 
                                                   if  md.iInvType[inv]=='hydro'
                                                   and (md.invphase_t[ip]+md.iConstTime[inv] <= md.Options['tfin'] and md.t_year[md.invphase_t[ip]+md.iConstTime[inv]]<=y)
                                                   and (md.invphase_t[ip]+md.iConstTime[inv]+md.iLifeTime[inv] >= md.Options['tfin'] or md.t_year[md.invphase_t[ip]+md.iConstTime[inv]+md.iLifeTime[inv]] >= y)) 
                                              if md.Options['Investment module']==1 else 0)
                                            for y in md.nyear}
                       } 
        for pt in md.nptech:
            YearlyProd[pt+' prod [GWh/y]'] = {y: sum(md.se.value*md.EeGENPROD[t,pld,pt,pm].value for t in md.ntime for pld in md.npload for pm in md.npmarket if md.t_year[t]==y and pm != 'PmSouthAfrica')
                                                          +sum(oo['PpProd_pp_y'][y][op] for op in md.nopp if md.op_ptech[op]==pt and md.op_pmarket[op]!='PmSouthAfrica')
                                                        for y in md.nyear}
            YearlyCapacity[pt+' capacity [MW]'] = {y:  sum(md.se.value*md.EeGENCAP[ky,pt,pm].value for ky in md.nyear for pm in md.npmarket if ky <= y and ky >= y-md.eLifeTime[pt,pm] and pm != 'PmSouthAfrica')
                                                    +sum(md.eOppCap[op] for op in md.nopp if md.op_ptech[op]==pt and md.op_pmarket[op]!='PmSouthAfrica')
                                                 for y in md.nyear}
        oo['YearlyEnergy'] = YearlyEnergy
        oo['YearlyProd'] = YearlyProd
        oo['YearlyCapacity'] = YearlyCapacity
        
        oo['YearlyCrop']={'Crop demand [t/y]':{y:sum(oo['CropBalance_y_cm']['Crop demand [Mt]'][y][cm] for cm in md.ncmarket if cm != 'CmWorldMarket') for y in md.nyear},
                    'Crop ext imports [Mt/y]':oo['CropBalance_y']['Crop external prod [Mt]'],
                    'Net Cultivated land [Mha]':oo['CropBalance_y']['Net Cultivated land [Mha]'],
                    'Net Irrigated land [Mha]':oo['CropBalance_y']['Net Irrigated land [Mha]'],
                    'Crop production [Mt]':{y:oo['CropBalance_y']['Crop production [Mt]'][y]-oo['CropBalance_y']['Crop external prod [Mt]'][y] for y in md.nyear}}
        oo['YearlyWater']={'Runoff [Mm3/y]':oo['WaterBalance_y']['RunOff [Mm3]'],
                     'Agr. net cons. [Mm3/y]':{y:oo['WaterBalance_y']['Agriculture irrigation [Mm3]'][y]+oo['WaterBalance_y']['Agriculture irrig. loss [Mm3]'][y]-oo['WaterBalance_y']['Agriculture return flow [Mm3]'][y] for y in md.nyear},
                     'User net cons. [Mm3/y]':{y:oo['WaterBalance_y']['User allocation [Mm3]'][y]+oo['WaterBalance_y']['User water loss [Mm3]'][y]-oo['WaterBalance_y']['User return flow [Mm3]'][y] for y in md.nyear},
                     'Outlet flow [Mm3]':oo['WaterBalance_y']['Outlet flow [Mm3]'],
                     'Reservoir storage [Mm3]':{y:oo['WaterBalance_y']['Reservoir storage [Mm3]'][y]/12 for y in md.nyear},
                     'Reservoir ET [Mm3]':oo['WaterBalance_y']['Reservoir evaporation [Mm3]']}
#%% DEBUG VARIABLES
        #initialize
        oo['Debug']     = {'Debug is OFF':{0:'Debug is off'}}
        oo['DebugWater']= {'Debug is OFF':{0:'Debug is off'}}
        oo['DebugCrop'] = {'Debug is OFF':{0:'Debug is off'}}
        oo['DebugEflow']= {'Debug is OFF':{0:'Debug is off'}}
        #summarize debug/dummy variables
        if md.Options['Debug mode'] in [1,2]:
            oo['Debug'] = {'Water balance':sum(md.DUMYWATER[t,c].value for t in md.ntime for c in md.ncatch) if md.Options['Debug mode'] ==1 else 'NoDebugWater',
                           'Eflow relax':sum(md.DUMYEFLOW[t,ef].value for ef in md.neflow for t in md.ntime), 
                           'Crop yield':sum(md.DUMYCROP[y,fz,cr].value for cr in md.ncrop for fz in md.nfzone for y in md.nyear)}          
            oo['DebugCrop']  = {y:{fz:sum(md.DUMYCROP[y,fz,cr].value for cr in md.ncrop) for fz in md.nfzone} for y in md.nyear}
            oo['DebugEflow'] = {t:{ef:md.DUMYEFLOW[t,ef].value for ef in md.neflow} for t in md.ntime}
        if md.Options['Debug mode'] == 1:
            oo['DebugWater'] = {t:{c:md.DUMYWATER[t,c].value for c in md.ncatch} for t in md.ntime}        
#%% VALIDATION - #ADD THIS SECTION IS SPECIFIC TO ZAMBEZI STUDY CASE !!!!! (SO FAR)
        if VALIDATION == 1:
            oo['WaterAbs_fz']    = {'Sim Consumption [Mm3/y]':{fz:oo['FarmZonesTable']['Irrig net abstraction [Mm3/y]'][fz] for fz in md.nifzone},
                                  'Obs Consumption [Mm3/y]':{fz:ReadParam('ObsWaterCons')[fz] for fz in md.nifzone}}
            oo['WaterAbs_co']    = {'Sim Consumption [Mm3/y]':{co:sum(oo['FarmZonesTable']['Irrig net abstraction [Mm3/y]'][fz] for fz in md.nifzone if md.fzone_country[fz]==co) for co in md.ncountry},
                                  'Obs Consumption [Mm3/y]':ReadParam('vAgrWaterCons')}
            oo['Hydropower_hp']  = {'Hp Prod sim [GWh/y]':oo['HydropowerTable']['Power Prod [GWh/y]'],
                                  'Hp Prod obs [GWh/y]':ReadParam('vObservedProd')}
            oo['ObsCropProd_cr_co'] = {cr:{co:ReadParam('vCropProd')[co,cr]*ReadParam('vScaleFactor')[co] for co in md.ncountry} for cr in md.ncrop}
            

            oo['CropValue_co']  = {'Sim':{co:oo['EconomicBalance_co']['Crop benefits [M$]'][co]+oo['EconomicBalance_co']['Crop export benefits [M$]'][co]-oo['EconomicBalance_co']['Crop import costs [M$]'][co]-oo['EconomicBalance_co']['Crop transport costs [M$]'][co] for co in md.ncountry},
								  'Sim2':{co:sum(md.kt_to_Mt*md.AcPROD[y,fz,cr].value * (md.aCropVal[md.fzone_cmarket[fz],cr] if md.Options['Crop market']==1 else md.aFarmVal[co,cr]) for y in md.nyear for fz in md.nfzone for cr in md.ncrop if md.fzone_country[fz]==co)/leni(md.nyear) for co in md.ncountry},
                                  'Obs':{co:ReadParam('vProdValue')[co]*ReadParam('vScaleFactor')[co] for co in md.ncountry if co in ReadParam('vProdValue').keys()}}

            if md.Options['Energy market'] == 1:
                oo['EnergyImExp_co']= {'EngyExp sim [GWh/y]':oo['EnergyBalance_co']['Net power export [GWh]'],
                                      'EngyImp sim [GWh/y]':oo['EnergyBalance_co']['Net power import [GWh]'],
                                      'EngyExp obs [GWh/y]':ReadParam('PowerExp'),
                                      'EngyImp obs [GWh/y]':ReadParam('PowerImp')}
            else:
               oo['EnergyImExp_co'] = {'No Data':{'No Data':0}} 
            
#%%INVESTMENT PLANNING            
        if md.Options['Investment module'] == 1:    
            oo['InvestmentPlan_inv_ip']  = {inv:{ip:md.IbINVEST[ip,inv].value for ip in md.ninvphase} for inv in md.ninvest}
            oo['Investments_ip']         = {ip:sum(md.IbINVEST[ip,inv].value*md.iCAPEX[inv] for inv in md.ninvest) for ip in md.ninvphase}
            oo['InvestmentCap_inv_y']    = {inv: {y:sum(md.IbINVEST[ip,inv].value*md.iInvCap[inv] for ip in md.ninvphase if md.t_year[md.invphase_t[ip]]==y) for y in md.nyear} for inv in md.ninvest}
            #oo['Investments_co   
#%%RETURN
        return oo
#%% POWER BI visualization
    def export_index_mapping(self,md,expfolder):
        #mapping function - export to .csv
        def simple_map(i1,i2,param_name):
            file = os.path.join(expfolder,param_name+'.csv')
            param=md.find_component(param_name)
            if param != None and param.sparse_keys() != []:
                i1_i2=pd.DataFrame.from_dict(param,orient='index')
                i1_i2.index.names=[i1]
                i1_i2.columns=[i2]
                i1_i2.to_csv(file,sep=';',decimal=',')
            else: #generate empty file for power BI
                i1_i2 = pd.DataFrame(columns=[i1,i2])
                i1_i2.to_csv(file,sep=';',decimal=',',index=False)
        #create output folder
        if not os.path.exists(expfolder):
            os.makedirs(expfolder)
        #map indexes   
        simple_map('ncatch','ncountry','catch_country') #catch_country       
        simple_map('nres','ncatch','res_catch') #res_catch
        simple_map('nuser','ncatch','user_catch')
        simple_map('nuser','ncountry','user_country')
        simple_map('ntime','nmonth','t_month')
        simple_map('ntime','nyear','t_year')
        simple_map('neflow','ncatch','eflow_catch')
        simple_map('nfzone','ncatch','fzone_catch')
        simple_map('nfzone','ncmarket','fzone_cmarket')
        simple_map('nfzone','ncountry','fzone_country')
        simple_map('nculture','ncrop','culture_crop')
        simple_map('ncmarket','ncountry','cmarket_country')
        simple_map('nhpp','nres','hp_res')
        simple_map('nhpp','ncatch','hp_catch')
        simple_map('nhpp','npmarket','hp_pmarket')
        simple_map('nopp','npmarket','op_pmarket')
        simple_map('npmarket','ncountry','pmarket_country')

    def export_mass_balances(self,oo,expfolder,scenario='WHATIF_main'):
        #Exports masss balances in a specific format to be read by the 'WHATIF_compare_result' and exported to powerBI 
        #oo : results dictionnary
        #md : model (pyomo)
        #folder: path to result folder

        Balances={
                'EconomicBalance':oo['EconomicBalance_y_co'],
                'EnergyBalance':oo['EnergyBalance_t_pm'],
                'CropBalance':oo['CropBalance_y_cm'],
                'WaterBalance':oo['WaterBalance_t_c']
                }

        filename= scenario+'_Balances.txt'
        pickle.dump(Balances,open(os.path.join(expfolder,filename),"wb"))
#%%EXPORT RESULT FUNCTIONS    
    def exportsheet(self,writer,sheet,data,title,dataref=[],index=[],order={},total={},color={},cols=1,srowi=5,scoli=0,distance=4,nested=0):    
        scol=scoli
        srow=srowi
        rowjump=0
        for k in range(len(data)):
            #Title of the data                                    
            ptitle = pd.DataFrame([title[k]]) 
            
            #To panda dataframes
            if k in index: #data is a single dictionary and needs orient='index'
                pdata = pd.DataFrame.from_dict(data[k], orient='index')
                if dataref != []: #Substract reference data (dataref) 
                    pdataref = pd.DataFrame.from_dict(dataref[k], orient='index')
                    for obj in pdata.index:
                        if obj in pdataref.index: #only compare if object is also in reference results
                            pdata.loc[obj]=pdata.loc[obj]-pdataref.loc[obj]
            else:
                if nested == 0:
                    pdata = pd.DataFrame.from_dict(data[k])
                else: #The data is a "triple dictionnary"
                    pdata = pd.concat({dic: pd.DataFrame(v).T for dic, v in data[k].items()}, axis=0)
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
#%%
    def export_to_excel(self,oo,exppath,VALIDATION=0,NEWSHEET=0):
        ## oo = output from function read_results (processed results from model)
        ## exppath = path of the excel file to export to
        ## NEWSHEET: 1= creates a new sheet (loose formatting of existing sheet), 0=updates existing sheet (keeps formatting)
        ## VALIDATION: 1= prints out Validation outputs (ZAMBEZI specific)
        months=['oct','nov','dec','jan','feb','mar','apr','may','jun','jul','aug','sep']        
        farmingzone=['Crop Production [Mt/y]','Production Value [M$/y]','Cultivation costs [M$/y]','Irrigation costs [M$/y]',
                     'Cultivated land Net [Mha]','Cultivated land Net [%]','Cultivated land Gross [%]','Irrig net abstraction [Mm3/y]',
                     'Irrig withdrawals [Mm3/y]','Irrig losses [Mm3/y]','Irrig return flow [Mm3/y]','Crop dem. satisfaction [%]','Land shadowprice [$/ha]']

        writer = pd.ExcelWriter(exppath, engine='openpyxl') #able to save formulas (graphs die however)
        if NEWSHEET==0:
            book = load_workbook(exppath)
            writer.book = book
            writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

        #INVESTMENTS
        if oo['Options']['Investment module'] == 1:
            data=[oo['InvestmentPlan_inv_ip'],oo['Investments_ip'],oo['InvestmentCap_inv_y']]
            title=['Invest sequence','Investments [M$]','Investment Capacity [Units]']
            self.exportsheet(writer,'Investments',data,title,index=[0,1])
        
        #POWER NETWORK
        if oo['Options']['Energy market'] == 1:
            title=['Transmission lines']
            data=[oo['TmsLineTable']]
            self.exportsheet(writer,'PowerNetwork',data,title)
            
        #POWER PLANTS
        if oo['Options']['Power plants'] == 1:
            data=[oo['PowerplantTable'],oo['PpProd_pp_m'],oo['PpProd_pp_y'],oo['PpBenefits_pp_y']]
            title=['Power plant','Power production [kWh/month]','Power production [GWh/y]','Net Benefits* [M$/y]']
            self.exportsheet(writer,'PowerPlants',data,title,order={1:months})
            
        #ENERGY
        if oo['Options']['Energy market'] == 1:
            #ADD: if power tech not on some variables not defined
            data=[oo['EnergyBalance'],oo['EnergyBalance_co'],oo['EnergyBalance_ca'],oo['EnergyBalance_m'],oo['EnergyBalance_y'],oo['EnergyVal_pm_m'], oo['PowerGenCap_pt_pm'], oo['PowerGenProd_pt_pm'], oo['PowerGenCap_pt_y'], oo['PowerGenProd_pt_y']]
            title=['ENERGY BALANCE' ,'ENERGY BALANCE'     ,'ENERGY BALANCE'     ,'ENERGY BALANCE'    ,'ENERGY BALANCE'    ,'EneryShadowprice [$/kWh]', 'Generic capacity invest [MW]', 'Generic capacity prod [GWh/year]','Generic capacity invest [MW]', 'Generic capacity prod [GWh/year]']
            self.exportsheet(writer,'Energy',data,title,index=[0,1,2,3,4,5],order={3:months,5:months},total={3:1})
        
        #CROPS
        data=[oo['CropBalance_y'], oo['CropPrice_cr_y'],oo['CropPrice_cr_cm'],oo['CropProduction_cr_co'],oo['LandUse_cr_co'],oo['CropImpShare_cr_co'], oo['CropExtProd_cr_co'], oo['LandUse_cr_y'],oo['LandUse_co_ft'],oo['LandUse_cr_ft'],oo['CropTransport_cm_cm'], oo['CropLoss_cm_cm'], oo['CropTransport_ct_cr']]
        title=['Crop Balance','Crop Price [$/t]','Crop Price [$/t]','Crop Production [Mt]','Land Use [Mha]','Crop Import Share [%ofDemand]','Crop external production [Mt]','Land Use [Mha]','Land Use [Mha]','Land Use [Mha]','Crop Transport [Mt/y]', 'Crop Transport losses[Mt/y]', 'Crop transport [Mt/y]']
        self.exportsheet(writer,'Crops',data,title)
        #LAND USE
        #ADD PUT ALL LAND USE OUTPUTS HERE
        
        #FARMING ZONES
        data=[oo['FarmZonesTable'],oo['LandUse_cul_fz'],oo['LandUse_fd_ft'],oo['Yield_cul_ft'],oo['CropProduction_cr_fz'],oo['FarmBenefits_ca_ft'],oo['FarmBenefits_fz_y']]
        title=['Farming zone','Land Use [ha/y]','Land Use [ha/y]','Yield [%]','Crop production [Mt/y]','Net benefits [M$/y]','Net benefits [M$/y]']
        self.exportsheet(writer,'FarmingZones',data,title,order={0:farmingzone})        
        
        #RESERVOIRS
        data=[oo['ReservoirTable'],oo['ResStor_res_m'],oo['ResStor_res_t']]
        title=['Reservoir','Average storage [Mm3]','Reservoir storage [Mm3]']
        self.exportsheet(writer,'Reservoirs',data,title,order={1:months})
        
        #GROUNDWATER
        if oo['Options']['Groundwater']==1:
            self.exportsheet(writer,'Groundwater',[oo['GwStor_c_m'],oo['GwStor_c_t']],['Groundwater storage [Mm3]','Groundwater storage [Mm3]'],order={0:months})
        
        #HYDROPOWER
        data=[oo['HydropowerTable'],oo['HpProd_hp_m'],oo['HpStor_hp_m'],oo['HpDis_hp_m'],oo['HpProd_hp_y'],oo['HpBenefits_hp_y']]
        title=['Hydropower','Power production [GWh]','Storage* [Mm3]','Discharge [Mm3]','Power production [GWh/y]', 'Net Benefits* [M$/year]']
        self.exportsheet(writer,'Hydropower',data,title,order={1:months,2:months,3:months},total={1:2,3:2,4:2,5:2})
        
        #ENVIRONMENTAL FLOWS
        if oo['Options']['Eflows'] == 1:
            self.exportsheet(writer,'EnvFlow',[oo['EnvFlowShadow_m_ec']],['e-flow shadowprice [$/m3]'],order={0:months})
        
        #USERS
        self.exportsheet(writer,'Users',[oo['UserTable']],['Water Users'])
        
        #MASS BALANCES
        data=[oo['EconomicBalance'],oo['EconomicBalance_co'],oo['EconomicBalance_ca'],
              oo['WaterBalance'],oo['WaterBalance_co'],oo['WaterBalance_ca'],
              oo['EnergyBalance'],oo['EnergyBalance_co'],oo['EnergyBalance_ca'],
              oo['CropBalance'],oo['CropBalance_co'],oo['CropBalance_ca']]
        title=['ECONOMIC BALANCE','ECONOMIC BALANCE','ECONOMIC BALANCE',
               'WATER BALANCE','WATER BALANCE','WATER BALANCE',
               'ENERGY BALANCE','ENERGY BALANCE','ENERGY BALANCE',
               'CROP BALANCE','CROP BALANCE','CROP BALANCE']
        self.exportsheet(writer,'MassBalance',data,title,index=[0,1,2,3,4,5,6,7,8,9,10,11],cols=3)

        
        #MAIN + DEBUG
        self.exportsheet(writer,'Main',[oo['MainTable'],oo['Options'],oo['Debug'],oo['DebugWater'],oo['DebugCrop'],oo['DebugEflow']],['Configurations','Options','Debug','DebugWater','DebugCrop','DebugEflow'],index=[0,1,2,3])
        
        #VALIDATION
        if VALIDATION == 1:
            data=[oo['WaterAbs_fz'], oo['WaterAbs_co'], oo['CropProduction_cr_co'], oo['ObsCropProd_cr_co'], oo['CropValue_co'], oo['Hydropower_hp'], oo['EnergyImExp_co']]
            title=['Water Abstraction', 'Water Abstraction', 'Sim Crop production [Mt/y]', 'Obs Crop production [Mt/y]', 'Value of crop Prod [M$/y]', 'Hydropower', 'Energy']
            self.exportsheet(writer,'Validation',data,title,total={0:0,1:0,2:0,3:0,4:0,5:0,6:0})
        
        #YEARLY TRENDS
        self.exportsheet(writer,'YearlyTrends',[oo['YearlyEnergy'],oo['YearlyProd'],oo['YearlyCapacity'],oo['YearlyCrop'],oo['YearlyWater']],['Yearly Trends','Yearly production','Yearly Capacity','Yearly Crop','Yearly Water'],cols=3)
        
        #ECONOMICS
        self.exportsheet(writer,'Economics',[oo['EconomicBalance_y']],['Economics'])
        
        #Save excel files
        writer.save()
#%%
    def export_scenario_analysis(self,scenarios,ref_scen,parallelresults,expfolder):
        #Assembles results from different scenarios and exports to excel
        ##scenarios = [scen_1,...,scen_n] list of scenario names
        ##ref_scen = {scen_1:refscen_1,...,scen_n:refscen_n} reference scenario to show differential values (A-B)
        ##parallelresults = [res_1,...,res_n] where res_i are outputs of ResultAnalysis.selectedresults()
        ##expfolder = path to folder to export to
        
        #generate list of all elements present in parallelresult
        elist=[]
        for ss in range(len(scenarios)):
            for key in parallelresults[ss].keys():
                elist.append(key)        
        elist=set(elist)
        #assemble all elements from parallelresults
        oo={elem:{scenarios[ss]:parallelresults[ss][elem] for ss in range(len(scenarios)) if elem in parallelresults[ss].keys()} for elem in elist}
       
        def RelativeResults(data, ref_scen):
            relresult=[]
            for result in data:       
                relresult.append({scen+'_'+ref_scen[scen]:{elem:result[scen][elem]-result[ref_scen[scen]][elem] for elem in result[scen].keys() if elem in result[ref_scen[scen]].keys()} for scen in result.keys() if ref_scen[scen] in result.keys()})
            return relresult
        
        ScenarioOutPath = expfolder + os.sep + 'SCENARIOS_compare.xlsx'
        writer = pd.ExcelWriter(ScenarioOutPath, engine='openpyxl') #able to save formulas (graphs die however)
        
        #EconomicBalance
        data    = [oo['EconomicBalance'],oo['EconomicBalance_co'],oo['WaterEconomicBalance_co'],oo['EnergyEconomicBalance_co'],oo['AgricultureEconomicBalance_co'],oo['AgricultureEconomicBalance2_co'],
                   oo['AgricultureConsSurplus_co'], oo['AgricultureProdSurplus_co'], oo['EnergyConsSurplus_co'], oo['EnergyProdSurplus_co']]
        reldata = RelativeResults(data, ref_scen)
        title   = ['Economic Balance [$]','Economic Balance [$]','Water Economic Balance [$]',' Energy Economic Balance [$]','Agriculture Economic Balance [$]','Agriculture Economic Balance 2[$]',
                   'Agr cons surplus [$]','Agr prod surplus [$]','Egy cons surplus [$]', 'Egy prod surplus [$]']
        self.exportsheet(writer,'EconomicBalance',data,title)
        self.exportsheet(writer,'RelEconomicBalance',reldata,title)
        
        #Other indicators
        data    = [oo['OtherIndicators'], oo['HydropowerProduction_co'], oo['EnergyExportBenefit_co'], oo['EnergyExports_co'], oo['GrossCultivatedArea_co'], oo['CropExpBenefit_co'],
                   oo['NetIrrArea_co'], oo['IrrigNetCons_co'], oo['EnergyCapacityConst_co'], oo['SolarCapacityConst_co'], oo['CO2Emission_co'], oo['EnergyValue_co'], oo['CropValue_co']]
        reldata = RelativeResults(data, ref_scen)
        title   = ['Other Indicators','Hydropower Production [kWh/y]','Energy Export benefit [$/y]','Energy Exports [kWh/y]','Gross Cultivated area [ha/y]', 'Crop export benefits [$/y]',
                   'Net irrigated area [ha/y]','Agricultural water consumption [m3/y]','Additional Power Invest [kWh/month]','Solar Invest [kWh/month]','CO2 emissions [t/y]','Average energy price [$/kWh]','Average crop price [$/t]']
        self.exportsheet(writer,'OtherIndicators',data,title)
        self.exportsheet(writer,'RelOtherIndicators',reldata,title)
        
        #NILE
        data= [oo['Runoff'],oo['Benefits'],oo['Reservoir storage'], oo['Shadow value'], oo['Hydropower production'],oo['zHydropower firm'],oo['zHydropower sec'], oo['Hydropower value'], oo['Turbine capacity'], oo['Hydropower discharge'], oo['Hydropower spill'], oo['User curt.'], oo['User value']]
        title =['Runoff','Benefits','Reservoir storage','Shadow value','Hydropower production','Hydropower firm','Hydropower sec','Hydropower value','Turbine capacity','Hydropower discharge','Hydropower spill','User curt.','User value']
        self.exportsheet(writer,'TIME',data,title,index=[0,1,2,3,4,5,6,7,8,9,10,11])
        
        data= [oo['yRunoff'],oo['yBenefits'],oo['yReservoir storage'], oo['yShadow value'], oo['yHydropower production'],oo['yzHydropower firm'],oo['yzHydropower sec'], oo['yHydropower value'], oo['yTurbine capacity'], oo['yHydropower discharge'], oo['yHydropower spill'], oo['yUser curt.'], oo['yUser value']]
        title =['Runoff','Benefits','Reservoir storage','Shadow value','Hydropower production','Hydropower firm','Hydropower sec','Hydropower value','Turbine capacity','Hydropower discharge','Hydropower spill','User curt.','User value']
        self.exportsheet(writer,'YEAR',data,title,index=[0,1,2,3,4,5,6,7,8,9,10,11])
        
        self.exportsheet(writer,'SUM_NILE',[oo['SUM_NILE']],['SUM_NILE'])
        self.exportsheet(writer,'SUM_NILE_REL',RelativeResults([oo['SUM_NILE']], ref_scen),['SUM_NILE'])
        #Investments
        if sum(oo['Options'][scen]['Investment module'] for scen in scenarios)>=1:
            data    = [oo['InvestmentDecision_inv'],oo['InvestmentPhase_inv']]
            title   = ['Investment decision', 'Investment phase']
            self.exportsheet(writer,'Investments',data,title)
        
        #Options
        self.exportsheet(writer,'Options',[oo['Options']],['OPTIONS'],index=[0])
        writer.save()
        