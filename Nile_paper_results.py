# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 18:58:06 2020

@author: rapp
"""

#LIBRARIES
import sys
import os
import pickle
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
#from openpyxl import load_workbook
dirname = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(dirname, 'bin'))
#from result_analysis import ResultAnalysis

#%%OPTIONS
FOLDERNAME='investX_nx2_ph24_grid5Wed06_05_2020_16h32'
METHOD=['MPC3c_50','SDP','PF','SIM']#['PF','MPChph1','MPChph6','MPChph12','MPChph36','MPChph48','MPChph72','MPChph120']
itype='NEW_grid5_nx2_'
TRANSFER=0
power_price=50
DRANGE=(60,150,5)
RRANGE=(60,140,5)
LOADEXISTING=0#['tr_grid5_nx2_']#0#'AG_nfl24m_g5_'#'AGdev_g2_nx2_'#'TRdev_g2_n02_'#'AGdev_nf4lm_g5_nx2_'#'AGdev_50_150_g5_nx2_'#0#'AGdev_60_140_g2_n02d5' #load this aggregated result file
AGGREGATE=0#'NEWag_grid5_nx2_' #save individual results to this file
nSS=0
PHgrid=0
EXPEXCEL=0
#dirname='N://Main_syn' 

#%% Figure param
FIGSIZE=(4*7,6)
pl.rcParams['font.family'] = 'Times New Roman'
pl.rcParams['font.size'] = 26
pl.rcParams['font.style'] = 'normal'
#%% INVESTMENT TABLES
#param
path = os.path.join(dirname,'Results',FOLDERNAME)
runoff=['c'+str(cs) for cs in range(RRANGE[0],RRANGE[1]+RRANGE[2],RRANGE[2])]
demand=['d'+str(ds) for ds in range(DRANGE[0],DRANGE[1]+DRANGE[2],DRANGE[2])]

ssyn=[''] if nSS==0 else [str(kss) for kss in range(nSS)]
#runoff=['rm30','rm20','rm10','','rp10','rp20','rp30']
#demand=['dm30','dm20','dm10','d','dp10','dp20','dp30','dp40']
AGSIZE=10
def invest_table(path,opt,transfer=0,loadex=0):   
    trs=[''] if transfer==0 else ['','a']
    if EXPEXCEL==1:
        writer = pd.ExcelWriter(os.path.join(path,'SUMMARY_'+itype+opt+'.xlsx'), engine='openpyxl')    
    #load
    dload=demand if PHgrid!=1 else demand[:-int(AGSIZE/DRANGE[2])]
    scenload=[opt+r+s+d+t for r in runoff for s in ssyn for d in dload for t in trs]
    rscen={opt+r+s+d+t:opt+r+s+('d'+str(int(d[1:])-10) if t!='a' else d) for r in runoff for s in ssyn for d in demand for t in trs}
    if transfer == 0:         
        SCEN={(r,d):[opt+r+s+d for s in ssyn] for r in runoff for d in demand}
        demidx=demand[int(AGSIZE/DRANGE[2]):]
        demanda=demand[:-int(AGSIZE/DRANGE[2])]     
    else:
        SCEN={(r,d):[opt+r+s+d+'a' for s in ssyn] for d in demand for r in runoff}  
        demidx=demand[:-int(AGSIZE/DRANGE[2])]
        demanda=demand[:-int(AGSIZE/DRANGE[2])]
    if loadex == 0:
        results={scen: pickle.load(open(os.path.join(path,scen+'.txt'),"rb")) for scen in scenload} 
        if AGGREGATE != 0:
            pickle.dump(results,open(os.path.join(path,AGGREGATE+opt+'.txt'),"wb"))
    else:
        results={}
        for kpath in LOADEXISTING:
            respath=os.path.join(path,kpath+opt+'.txt')
            if os.path.exists(respath):
                temp_dic=pickle.load(open(respath,"rb"))
                results.update(temp_dic)
    #time=[t for t in results[SCEN[runoff[0],demand[0]]]['Hydropower discharge'].keys()]
    
    #functions
    def diff_result(scen,e1,e2):
        return sum(results[scenk][e1][e2]-results[rscen[scenk]][e1][e2] for scenk in scen)/len(scen)
    def data_to_xlsx(data,name):
        if EXPEXCEL==1:
            pdata = pd.DataFrame.from_dict(data,orient='index')
            pdata=pdata.reindex(runoff)
            pdata.to_excel(writer, sheet_name=name, index=runoff)#,sep=';',decimal=',')
    #Economics
    #Invest value
    if PHgrid!=1:
        Inv_value = {r:{d:diff_result(SCEN[r,d],'SUM_NILE','tBenefits') for d in demidx} for r in runoff}
        data_to_xlsx(Inv_value,'Inv_value')
    #Inv_valueXt = {r:{d:sum(diff_result(SCEN[r,d],'Benefits',t) for t in time[:-12])/(len(time)-12)*12 for d in demidx} for r in runoff}
    #data_to_xlsx(Inv_valueXt,'Inv_valueXt')
    #User val
    #User_value = {r:{d:diff_result(SCEN[r,d],'SUM_NILE','tUser value') for d in demidx} for r in runoff}
    #data_to_xlsx(User_value,'User_value')
    #Hydropwoer val
    #Dis_value = {r:{d:sum(diff_result(SCEN[r,d],'Hydropower discharge',t)*results[rscen[SCEN[r,d]]]['Head factor'][t]
    #                for t in time)/len(time)*12*power_price for d in demidx} for r in runoff}
    #data_to_xlsx(Dis_value,'Dis_value')
    #Hp_value = {r:{d:diff_result(SCEN[r,d],'SUM_NILE','tHydropower value') for d in demidx} for r in runoff}
    #data_to_xlsx(Hp_value,'Hp_value')
    #Head effet
    #Head_value = {r:{d:sum(diff_result(SCEN[r,d],'Head factor',t)*results[SCEN[r,d]]['Hydropower discharge'][t] 
    #                    for t in time)/len(time)*12*power_price for d in demidx} for r in runoff}
    #data_to_xlsx(Head_value,'Head_value')
    #Spill val
    #Spill_value = {r:{d:-diff_result(SCEN[r,d],'SUM_NILE','tHydropower spill val') for d in demidx} for r in runoff}
    #data_to_xlsx(Spill_value,'Spill_value')
    #avuservalue
    #User_AvValue = {r:{d:diff_result(SCEN[r,d],'SUM_NILE','tUser av value') for d in demidx} for r in runoff}
    #data_to_xlsx(User_AvValue,'User_AvValue')
    if PHgrid!=1:
        Shadow = {r:{d:diff_result(SCEN[r,d],'SUM_NILE','tShadow value') for d in demidx} for r in runoff}
        data_to_xlsx(Shadow,'Shadow')
    
    #Indicators
    tBenefitsm2 = {r:{d:sum(results[ks]['SUM_NILE']['tBenefitsm2'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(tBenefitsm2,'tBenefitsm2')
    FinalStorage = {r:{d:sum(results[ks]['SUM_NILE']['FinalStorage'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(FinalStorage,'FinalStorage')
    #Total benefits     
    Benefits_a = {r:{d:sum(results[ks]['SUM_NILE']['tBenefits'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(Benefits_a,'Benefits_a')
    #Hydropower
    Hpval_a = {r:{d:sum(results[ks]['SUM_NILE']['tHydropower value'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(Hpval_a,'Hpval_a')
    #Storage
    ResStor_a = {r:{d:sum(results[ks]['SUM_NILE']['tReservoir storage'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(ResStor_a,'ResStor_a')
    #Spill
    Spill_a = {r:{d:sum(results[ks]['SUM_NILE']['tHydropower spill'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(Spill_a,'Spill_a')
    #User
    UserAvValue_a = {r:{d:sum(results[ks]['SUM_NILE']['tUser av value'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(UserAvValue_a,'UserAvValue_a')
    UserCurt_a = {r:{d:sum(results[ks]['SUM_NILE']['tUser curt.'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(UserCurt_a,'UserCurt_a')
    UserVal_a = {r:{d:sum(results[ks]['SUM_NILE']['tUser value'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(UserVal_a,'UserVal_a')
    #Shadow
    Shadow_a = {r:{d:sum(results[ks]['SUM_NILE']['tShadow value'] for ks in SCEN[r,d])/len(SCEN[r,d]) for d in demanda} for r in runoff}
    data_to_xlsx(Shadow_a,'Shadow_a')
    if EXPEXCEL==1:
        writer.save()
    output={'Benefits_a':Benefits_a,'UserCurt_a':UserCurt_a, 
            'ResStor_a':ResStor_a,'Spill_a':Spill_a,'Shadow_a':Shadow_a,
            'UserAvValue_a':UserAvValue_a,'UserVal_a':UserVal_a,'Hpval_a':Hpval_a,
            'tBenefitsm2':tBenefitsm2,'FinalStorage':FinalStorage}
    if PHgrid!=1:
        output.update({'Inv_value':Inv_value})
    return output

#%% plots
def invest_plots(res,methods,refmeth=0,title=0):
    def dict_to_matrix(dic,size,thrm=0,uthr=65):
        mtrx=np.zeros(size)
        ik1=-1
        ik2=0
        for k1 in dic.keys():
            ik1+=1
            ik2=-1
            for k2 in dic[k1].keys():
                ik2+=1
                if thrm == 0 or dic[k1][k2]>uthr:
                    mtrx[ik1,ik2]=dic[k1][k2]
                else:
                    mtrx[ik1,ik2]=0
        return mtrx
    #Plot
    size=(len(runoff),len(demand)-int(AGSIZE/DRANGE[2]))
    x=range(DRANGE[0],DRANGE[1]-DRANGE[2]*(int(AGSIZE/DRANGE[2])-1),DRANGE[2])
    y=range(RRANGE[0],RRANGE[1]+RRANGE[2],RRANGE[2])
    fig,axs = pl.subplots(1, len(methods),figsize=FIGSIZE)
    if title != 0:
        fig.suptitle(title, fontsize=30)
    for k in range(len(methods)):
        m=methods[k]
        mres=dict_to_matrix(res[m],size)
        if refmeth !=0 and m!=refmeth:
            mref=dict_to_matrix(res[refmeth],size)
            mres=(mres-mref)#/mref*100
        pixim = axs[k].pcolor(x,y,mres)
        fig.colorbar(pixim,ax=axs[k])
        axs[k].set_ylabel('Runoff [%]')
        axs[k].set_xlabel('Water demand [%]')
        if methods[k] not in ['MPC3cm','MPC3c_50']:
            axs[k].set_title(methods[k])
        else:
            axs[k].set_title('MPC')
    fig.tight_layout()
    fig.subplots_adjust(left=0.05, top=0.85)
    fig.savefig(os.path.join(path,
        'absolute'+title.split(' ')[0]+'_'+'_'.join([str(k) for k in RRANGE])+'_'+'_'.join(METHOD)+'.png'),
        bbox_inches='tight')
    
#%%CREATE TABLES

refmeth=METHOD[0]

allres={m:invest_table(path,m,transfer=TRANSFER,loadex=LOADEXISTING) for m in METHOD}
benres={m:allres[m]['Benefits_a'] for m in METHOD}
if PHgrid!=1:
    invres={m:allres[m]['Inv_value'] for m in METHOD}
stores={m:allres[m]['ResStor_a'] for m in METHOD}
splres={m:allres[m]['Spill_a'] for m in METHOD}
crtres={m:allres[m]['UserCurt_a'] for m in METHOD}
shadow={m:allres[m]['Shadow_a'] for m in METHOD}
benm2={m:allres[m]['tBenefitsm2'] for m in METHOD}
fstor={m:allres[m]['FinalStorage'] for m in METHOD}
avuval={m:allres[m]['UserVal_a'] for m in METHOD}
hpval={m:allres[m]['Hpval_a'] for m in METHOD}
#%%
invest_plots(benres,METHOD,title='absTotal benefits (M$/year)')
#if PHgrid!=1:
#    invest_plots(invres,METHOD,title='absInvestment benefits [M$/year]')

invest_plots(benres,METHOD,refmeth=refmeth,title='(a) System total benefits (M$/year)')
if PHgrid!=1:
    invest_plots(invres,METHOD,refmeth=refmeth,title='(c) Transfer project benefits (M$/year)')
#
#invest_plots(stores,METHOD,refmeth=refmeth,title='Storage [BCM]')
#invest_plots(splres,METHOD,refmeth=refmeth,title='Spills [BCM/year]')
#invest_plots(crtres,METHOD,refmeth=refmeth,title='Demand curtailment [BCM/year]')
#
#invest_plots(benm2,METHOD,refmeth=refmeth,title='Benefits nolast2years [M$/year]')
#invest_plots(fstor,METHOD,refmeth=refmeth,title='Final Storage [BCM/year]')
#invest_plots(benres,METHOD,refmeth=refmeth,title='Total benefits [M$/year]')
#invest_plots(hpval,METHOD,refmeth=refmeth,title='Hydropower value [M$/year]')
#invest_plots(avuval,METHOD,refmeth=refmeth,title='Demand value [M$/year]')
    
#%%
GPLOT=1
if GPLOT==1:
    allres_t0={m:invest_table(path,m,transfer=0,loadex=LOADEXISTING) for m in METHOD}
    allres_t1={m:invest_table(path,m,transfer=1,loadex=LOADEXISTING) for m in METHOD}
    sysben={m:allres_t0[m]['Benefits_a'] for m in METHOD}
    aginvb={m:allres_t0[m]['Inv_value'] for m in METHOD}
    trinvb={m:allres_t1[m]['Inv_value'] for m in METHOD}
    invest_plots(benres,METHOD,refmeth=refmeth,title='(a) System total benefits (M$/year)')
    invest_plots(aginvb,METHOD,refmeth=refmeth,title='(b) Irrigation project benefits (M$/year)')
    invest_plots(trinvb,METHOD,refmeth=refmeth,title='(c) Transfer project benefits (M$/year)')
    