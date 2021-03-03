# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 16:35:10 2020

@author: RAPY
"""
#%% import libraries
import seaborn as sns
import matplotlib.pyplot as plt
sepin=';'
decin=','
import pandas as pd
import numpy as np
import os
dirname = os.path.abspath(os.path.dirname(__file__))

#%% Define options
FOLDERNAME='cc2C_all_loadPFx'
datapath = os.path.join(dirname,'Results',FOLDERNAME)
#%%Climate categories
#runoff data
path=os.path.join(datapath,'IFPRI_IDX','Runoff_Mm3'+'.csv')
Runoff_Mm3=pd.read_csv(path,sep=sepin,decimal=decin)
#define categories
clim=Runoff_Mm3.groupby(['scenario']).sum().drop(['nyear','ntime'],axis=1)
cat=np.percentile(clim,[0,33,67])
cat_def={0:'dry',1:'average',2:'wet'}
def define_cat(x):
    for k in range(len(cat)-1,-1,-1):
        if x>=cat[k]:
            return cat_def[k]
clim['cat']=clim.applymap(define_cat)
Runoff_Mm3['cat']=Runoff_Mm3['scenario'].apply(lambda scen: clim['cat'][scen])

#%% import data
def read_idx(name,climcat=1,XspecialcorrX=1):   
    path=os.path.join(datapath,'IFPRI_IDX',name+'.csv')
    pdata=pd.read_csv(path,sep=sepin,decimal=decin)
    if XspecialcorrX==1:#WARNING - SPECIAL CORRECTION !! WARNING !
        pdata['nyear']=pdata['nyear']-2 
    if climcat==1:
        pdata['cat']=pdata['scenario'].apply(lambda scen: clim['cat'][scen]) #define climate categories
    return pdata
#Hydrology
Runoff_Mm3=read_idx('Runoff_Mm3')
Runoff_Mm3['Runoff_Mm3']=Runoff_Mm3['Runoff_Mm3']/1000 #change unit
P_mm=read_idx('P_mm')
ET_mm=read_idx('ET_mm')
#Agriulture
Crop_area_kha=read_idx('Crop_area_kha')
Crop_price_dpt=read_idx('Crop_price_dpt')
Crop_prod_kt=read_idx('Crop_prod_kt')
Crop_val_Md=read_idx('Crop_val_Md')
Crop_yield_tpha=Crop_prod_kt.drop(['Crop_prod_kt'],axis=1)
Crop_yield_tpha['Crop_yield_tpha']=Crop_prod_kt['Crop_prod_kt']/Crop_area_kha['Crop_area_kha']
Cul_Dem_mm=read_idx('Cul_Dem_mm')
Cul_Rain_mm=read_idx('Cul_Rain_mm')
#Cul_Cons_Mm3=read_idx('Cul_Cons_Mm3')
#Energy
Hydropower_prod_GWh=read_idx('Hydropower_prod_GWh')
Hydropower_prod_GWh['Hydropower_prod_GWh']=Hydropower_prod_GWh['Hydropower_prod_GWh']/1000 #change to TWh
Hydropower_val_Md=read_idx('Hydropower_val_Md')
#Power_newprod_GWh=read_idx('Power_newprod_GWh')
#Power_newprod_GWh['Power_newprod_GWh']=Power_newprod_GWh['Power_newprod_GWh']/1000 #change to TWh
Power_price_dpkWh=read_idx('Power_price_dpkWh')
Power_price_dpkWh['Power_price_dpkWh']=Power_price_dpkWh['Power_price_dpkWh']*1000 #change to $/MWh

#%%IMPACT data
vAREAX0=read_idx('vAREAX0',climcat=0,XspecialcorrX=1)
vPCX0=read_idx('vPCX0',climcat=0,XspecialcorrX=1)
vPRODX0=read_idx('vPRODX0',climcat=0,XspecialcorrX=1)
#create yield
vYLDX0=vPRODX0.drop(['Val'],axis=1)
vYLDX0['Val']=vPRODX0['Val']/vAREAX0['Val']

#%% Plot data
def plot_data(pdat,y,vframe={},group='sum',x='nyear',
              ci='sd',hue=None,legend='auto',ax=None,
              columns=['scenario'],color=None, decade=False,
              refdat=pd.DataFrame(),reftype='rel',refy=None):
    #columns to use
    columns=columns+[x]
    if hue!=None:
        columns=columns+[hue]
    #prepare data
    def _prepare_data(data):
        #data slicer
        for var in vframe.keys(): #select specifc elements
            select=vframe[var]
            if var in data.columns:
                if type(select) is list:
                    data=data[data[var].isin(select)]
                else:
                    data=data[data[var]==select]
            else:
                print('Warning: '+var+' not found in data '+y)
        #group data   
        if group=='mean':
            data=data.groupby(by=columns).mean()
        else:
            data=data.groupby(by=columns).sum()
        #average over decades
        if decade == True:
            data=data.reset_index('nyear') #put index as column
            data['nyear']=data['nyear'].apply(lambda k:int(str(k)[0:3]+'0')) #create a column with the decade
            data.set_index('nyear', append=True, inplace=True) #put index back to indexes
            data=data.groupby(by=columns).mean()
            data=data.reset_index('nyear') #put index as column
        return data
    pdat=_prepare_data(pdat)
    if not refdat.empty:
        refdat=_prepare_data(refdat)
        if refy == None:
            refy=y
        if reftype=='rel':
            pdat[y]=(pdat[y].values-refdat[refy].values)/refdat[refy].values*100 #.sub(refdat[refy])
        elif reftype=='divide':
            pdat[y]=pdat[y].divide(refdat[refy])
        elif reftype=='multiply':
            pdat[y]=pdat[y].mul(refdat[refy])
        
    #plot average over decade   
    if decade == True:
        g1=sns.boxplot(x=x, y=y, hue=hue, data=pdat, ax=ax)
        g1.set(xlabel=None)
        g1.set(ylabel=None)
    #plot with seaborn confidence interval (std, bootstrapping or None)
    elif type(ci) is not list:
        g1=sns.lineplot(x=x, y=y, hue=hue, ci=ci, legend=legend, data=pdat, ax=ax)
        g1.set(xlabel=None)
        g1.set(ylabel=None)
        return g1
    #plot area between percentiles (e.g. ci=[0.05,0.95]) - does not work with hue !
    else:              
        lowy=pdat.groupby(level=[x]).quantile(ci[0])
        highy=pdat.groupby(level=[x]).quantile(ci[1])
        if color == None:
            #plot average
            g1=sns.lineplot(x=x, y=y, hue=hue, ci=None, legend=legend, data=pdat, ax=ax)
            g1.set(xlabel=None)
            g1.set(ylabel=None)
            #plot confidence interval based on percentiles
            g2=ax.fill_between(lowy.index,lowy[y].values,highy[y].values,alpha=0.5)
            return g1,g2
        else:
            ax.fill_between(lowy.index,lowy[y].values,highy[y].values,alpha=0.3,color=color)
#%% Runoff
sns.set_theme(style="darkgrid")
sns.lineplot(x="nyear", y="Runoff_Mm3", hue='cat',
             ci=99,
             #legend=False,#estimator=None,
             data=Runoff_Mm3.groupby(by=['nyear','scenario','cat'],axis=0).sum())
#%% Hydrology
SAVE=1
fig1, ax = plt.subplots(1,3,figsize=(12,4))
legend=False #'auto'
hue=None
ci=[0.05,0.95]
vframe={'ncountry':['Zambia','Zimbabwe','Mozambique','Malawi','Angola','Tanzania']}
vframeb=vframe.copy()
vframeb['scenario']='base'
plot_data(P_mm,'P_mm',vframe=vframe,hue=hue,legend=legend,ax=ax[0],ci=ci,group='mean')
plot_data(P_mm,'P_mm',vframe=vframeb,ax=ax[0],group='mean')
ax[0].set_title('Precipitation')
ax[0].set_ylabel('mm/month')
plot_data(ET_mm,'ET_mm',vframe=vframe,hue=hue,legend=legend,ax=ax[1],ci=ci,group='mean')
plot_data(ET_mm,'ET_mm',vframe=vframeb,ax=ax[1],group='mean')
ax[1].set_title('ET0')
ax[1].set_ylabel('mm/month')
g1,g2=plot_data(Runoff_Mm3,'Runoff_Mm3',vframe=vframe,hue=hue,legend=legend,ax=ax[2],ci=ci)
g3=plot_data(Runoff_Mm3,'Runoff_Mm3',vframe=vframeb,ax=ax[2])
ax[2].set_title('Runoff')
ax[2].set_ylabel('BCM/year')
#fig1.tight_layout()
fig1.legend([g1,g2,g3], 
            labels=['Average climate change','Historical',
                    'Climate change between 5% and 95% percentiles'], 
            loc=[0.14,0.12])
if SAVE==1:
    fig1.savefig('Hydrology_ifpri_idx',dpi=300,transparent=True)
    
#%% Hydrology by country
HYDROBYCOUNTRY=1
if HYDROBYCOUNTRY==1:
    SAVE=1
    countries=['Zambia','Zimbabwe','Mozambique','Malawi']
    fig1b, ax = plt.subplots(len(countries),3,figsize=(12,10))
    legend=False #'auto'
    hue=None
    ci=[0.05,0.95]
    for k in range(len(countries)):
        vframe={'ncountry':countries[k]}
        vframeb={'scenario':'base','ncountry':countries[k]}
        plot_data(P_mm,'P_mm',vframe=vframe,hue=hue,legend=legend,ax=ax[k,0],ci=ci)
        plot_data(P_mm,'P_mm',vframe=vframeb,ax=ax[k,0])
        
        
        plot_data(ET_mm,'ET_mm',vframe=vframe,hue=hue,legend=legend,ax=ax[k,1],ci=ci)
        plot_data(ET_mm,'ET_mm',vframe=vframeb,ax=ax[k,1])
        
        g1,g2=plot_data(Runoff_Mm3,'Runoff_Mm3',vframe=vframe,hue=hue,legend=legend,ax=ax[k,2],ci=ci)
        g3=plot_data(Runoff_Mm3,'Runoff_Mm3',vframe=vframeb,ax=ax[k,2])
        ax[k,0].set_ylabel(countries[k])
        #fig1.tight_layout()
    ax[0,0].set_title('Precipitation (mm/year)')  
    ax[0,1].set_title('ET0 (mm/year)')
    ax[0,2].set_title('Runoff (BCM/year)')
    fig1b.legend([g1,g2,g3], 
                 labels=['Average climate change','Historical',
                         'Climate change between 5% and 95% percentiles'], 
                 loc=[0.14,0.12])
    if SAVE==1:
        fig1b.savefig('Hydrology_bycountry_ifpri_idx',dpi=300,transparent=True)
#%% Agriculture
SAVE=1
crops=['sugarcane','tobacco','cotton']#['maize','cereals','cassava','roots']#['rice','pulses','oilseeds']# 
#['vegetable','fruits']#
ntype=['R','I']
ntype2=['rainfed cereals','irrigated cereals']
iterate=crops
fig2, ax = plt.subplots(5,len(iterate),figsize=(14,12))
hue=None
ci=[0.05,0.95]
legend=False

for k in range(len(iterate)):
    vframe={'ncrop':crops[k],'ncountry':['Zambia','Zimbabwe','Mozambique','Malawi']}
    #vframe={'ncrop':'cereals','ncountry':['Zambia','Zimbabwe','Mozambique','Malawi'],'ntype':ntype[k]}
    vframeb=vframe.copy()
    vframeb['scenario']='base'
    #area
    plot_data(Crop_area_kha,'Crop_area_kha',vframe=vframe,hue=hue,legend=legend,ax=ax[0,k],ci=ci)
    plot_data(Crop_area_kha,'Crop_area_kha',vframe=vframeb,legend=legend,ax=ax[0,k])
    plot_data(vAREAX0,'Val',vframe=vframe,columns=[],legend=legend,ax=ax[0,k])
    ax[0,k].set_title(iterate[k])
    #production
    plot_data(Crop_prod_kt,'Crop_prod_kt',vframe=vframe,hue=hue,legend=legend,ax=ax[1,k],ci=ci)
    plot_data(Crop_prod_kt,'Crop_prod_kt',vframe=vframeb,legend=legend,ax=ax[1,k])
    plot_data(vPRODX0,'Val',vframe=vframe,columns=[],legend=legend,ax=ax[1,k])
    #yield
    #plot_data(Crop_yield_tpha,'Crop_yield_tpha',vframe=vframe,hue=hue,group='mean',legend=legend,ax=ax[2,k],ci=ci)
    #plot_data(Crop_yield_tpha,'Crop_yield_tpha',vframe=vframeb,group='mean',legend=legend,ax=ax[2,k])
    #plot_data(vYLDX0,'Val',vframe=vframe,group='mean',columns=[],legend=legend,ax=ax[2,k])
    plot_data(Crop_prod_kt,'Crop_prod_kt',vframe=vframe,hue=hue,group='sum',legend=legend,ax=ax[2,k],ci=ci,
              refdat=Crop_area_kha,reftype='divide',refy='Crop_area_kha')
    plot_data(Crop_prod_kt,'Crop_prod_kt',vframe=vframeb,group='sum',legend=legend,ax=ax[2,k],
              refdat=Crop_area_kha,reftype='divide',refy='Crop_area_kha')
    plot_data(vPRODX0,'Val',vframe=vframe,group='mean',columns=[],legend=legend,ax=ax[2,k],
              refdat=vAREAX0,reftype='divide')
    #price
    g1,g2=plot_data(Crop_price_dpt,'Crop_price_dpt',vframe=vframe,hue=hue,group='mean',legend=legend,ax=ax[3,k],ci=ci)
    g3=plot_data(Crop_price_dpt,'Crop_price_dpt',vframe=vframeb,group='mean',legend=legend,ax=ax[3,k])
    g4=plot_data(vPCX0,'Val',vframe=vframe,group='mean',columns=[],legend=legend,ax=ax[3,k])
    #value
    plot_data(Crop_val_Md,'Crop_val_Md',vframe=vframe,hue=hue,legend=legend,ax=ax[4,k],ci=ci)
    plot_data(Crop_val_Md,'Crop_val_Md',vframe=vframeb,legend=legend,ax=ax[4,k])
    
    #set labels
    if k==0:
        ax[0,k].set_ylabel('Area (1000 ha)')
        ax[1,k].set_ylabel('Production (1000 t)')
        ax[2,k].set_ylabel('Yield (t/ha)')
        ax[3,k].set_ylabel('Price ($/t)')
        ax[4,k].set_ylabel('Value (M$)')
#fig2.tight_layout()
fig2.legend([g1,g2,g3,g4], 
            labels=['Average climate change','Historical',
                    'IMPACT model','Climate change between 5% and 95% percentiles'], 
            loc='lower center')

if SAVE==1:
    fig2.savefig('Agriculture_ifpri_idx',dpi=300)
#%% Hydropower
SAVE=1
countries=['Zambia','Zimbabwe','Mozambique','Malawi']
hue=None
decade=False
legend=False
ci=[0.05,0.95]
ci2=[0.15,0.85]
fig3, ax = plt.subplots(3,len(countries),figsize=(12,9))
for k in range(len(countries)):
    vframe={'ncountry':countries[k]}
    vframeb=vframe.copy()
    vframeb['scenario']='base'
    #Production
    plot_data(Hydropower_prod_GWh,'Hydropower_prod_GWh',vframe=vframe,hue=hue,legend=legend,ax=ax[0,k],ci=ci,decade=decade)    
    plot_data(Hydropower_prod_GWh,'Hydropower_prod_GWh',vframe=vframeb,ax=ax[0,k])
    #plot_data(Hydropower_prod_GWh,'Hydropower_prod_GWh',vframe=vframe,hue=hue,legend=legend,ax=ax[0,k],ci=ci2,color='darkblue')
    ax[0,k].set_title(countries[k])
    #Price
    plot_data(Power_price_dpkWh,'Power_price_dpkWh',vframe=vframe,hue=hue,legend=legend,group='mean',ax=ax[1,k],ci=ci,decade=decade) 
    plot_data(Power_price_dpkWh,'Power_price_dpkWh',vframe=vframeb,group='mean',ax=ax[1,k])
    #plot_data(Power_price_dpkWh,'Power_price_dpkWh',vframe=vframe,hue=hue,legend=legend,group='mean',ax=ax[1,k],ci=ci2,color='darkblue')
    #Value
    g1,g2=plot_data(Hydropower_val_Md,'Hydropower_val_Md',vframe=vframe,hue=hue,legend=legend,ax=ax[2,k],ci=ci,decade=decade)    
    g3=plot_data(Hydropower_val_Md,'Hydropower_val_Md',vframe=vframeb,ax=ax[2,k])
    #plot_data(Hydropower_val_Md,'Hydropower_val_Md',vframe=vframe,hue=hue,legend=legend,ax=ax[2,k],ci=ci2,color='darkblue')
    
    if k==0:
        ax[0,k].set_ylabel('Hydropower (Twh/year)')
        ax[1,k].set_ylabel('Price ($/MWh)')
        ax[2,k].set_ylabel('Value (M$/year)')

fig3.legend([g1,g2,g3], labels=['Average climate change','Historical','Climate change between 5% and 95% percentiles'], loc='lower center')
#fig3.tight_layout()
if SAVE==1:
    fig3.savefig('Hydropower_ifpri_idx',dpi=300)

#%% New 
NEWPRODFIGURE=0
SAVE=0
countries=['Zambia','Zimbabwe','Mozambique','Malawi','SouthAfrica']
vframe={}
hue='nptech'
legend=True
if NEWPRODFIGURE==1:
    fig4, ax = plt.subplots(1,len(countries),figsize=(10,6))
    for k in range(len(countries)):
        vframe={'ncountry':countries[k]}
        plot_data(Power_newprod_GWh,'Power_newprod_GWh',vframe=vframe,hue=hue,legend=True,ax=ax[k])
        plot_data(Power_newprod_GWh,'Power_newprod_GWh',vframe=vframe,hue=hue,legend=False,ax=ax[k])
        ax[k].set_ylabel('New Power (TWh)')
    
#%% PRIORIZATION
PRIORIZATION=0
if PRIORIZATION==1:
    #PLOT OPTIONS
    SAVE=1
    legend=False
    decade=True
    hue=None
    ci=[0.05,0.95]
    vframeag={'ntype':'I'}
    vframeeg={}
    opt=['fullag','fullpower','feb7000']
    basekw='fullag|fullpower|feb7000' #
    fig5, ax = plt.subplots(4,len(opt),figsize=(12,10))
    #DATA
    def _sort(data,keyword,b=True):
        data=data[data['scenario'].str.contains(keyword)==b]
        data['scenario']
        return data
    
    #PLOTS
    for k in range(len(opt)):
        #Hydropower value
        plot_data(_sort(Hydropower_val_Md,opt[k]),'Hydropower_val_Md',
                  vframe=vframeeg,legend=legend,ax=ax[0,k],ci=ci,decade=decade,
                  refdat=_sort(Hydropower_val_Md,basekw,b=False))
        #Irrigated value
        plot_data(_sort(Crop_val_Md,opt[k]),'Crop_val_Md',
                  vframe=vframeag,legend=legend,ax=ax[1,k],ci=ci,decade=decade,
                  refdat=_sort(Crop_val_Md,basekw,b=False))
        #Hydropower prod
        plot_data(_sort(Hydropower_prod_GWh,opt[k]),'Hydropower_prod_GWh',
                  vframe=vframeeg,legend=legend,ax=ax[2,k],ci=ci,decade=decade,
                  refdat=_sort(Hydropower_prod_GWh,basekw,b=False))
        #Irrigated consumption
        plot_data(_sort(Cul_Cons_Mm3,opt[k]),'Cul_Cons_Mm3',
                  vframe=vframeag,legend=legend,ax=ax[3,k],ci=ci,decade=decade,
                  refdat=_sort(Cul_Cons_Mm3,basekw,b=False))
    #axes
    ax[0,0].set_title('Agriculture')
    ax[0,1].set_title('Power')
    ax[0,2].set_title('Ecosystems')
    ax[0,0].set_ylabel('Hydropower value (%)')
    ax[1,0].set_ylabel('Irrigated crop value (%)')
    ax[2,0].set_ylabel('Hydropower production (%)')
    ax[3,0].set_ylabel('Irrigation consumption (%)')    
    fig5.tight_layout()
    if SAVE==1:
        fig5.savefig('Tradeoff_ifpri_idx',dpi=300)
    
#%%PLOT VERSION 2
PRIORIZATIONV2=0
if PRIORIZATIONV2==1:
    fig6, ax = plt.subplots(3,4,figsize=(12,9))
    for k in range(len(opt)):
        #Hydropower value
        plot_data(_sort(Hydropower_val_Md,opt[k]),'Hydropower_val_Md',
                  vframe=vframeeg,legend=legend,ax=ax[k,0],ci=ci,decade=decade)
        #Irrigated value
        plot_data(_sort(Crop_val_Md,opt[k]),'Crop_val_Md',
                  vframe=vframeag,legend=legend,ax=ax[k,1],ci=ci,decade=decade)
        #Hydropower prod
        plot_data(_sort(Hydropower_prod_GWh,opt[k]),'Hydropower_prod_GWh',
                  vframe=vframeeg,legend=legend,ax=ax[k,2],ci=ci,decade=decade)
        #Irrigated cons
        plot_data(_sort(Cul_Cons_Mm3,opt[k]),'Cul_Cons_Mm3',
                  vframe=vframeag,legend=legend,ax=ax[k,3],ci=ci,decade=decade)        
    #axes
    ax[0,0].set_title('Hydropower value (M$/y)')
    ax[0,1].set_title('Irrigated crop value (M$/y)')    
    ax[0,2].set_title('Hydropower production (TWh/y)')
    ax[0,3].set_title('Irrigation consumption (Mm3/y)')
    ax[0,0].set_ylabel('Agriculture')
    ax[1,0].set_ylabel('Power')
    ax[2,0].set_ylabel('Ecosystems')
    fig6.tight_layout()
#%% PRIORIZATION v3
PRIORIZATIONV3=1
if PRIORIZATIONV3==1:
    #PLOT OPTIONS
    SAVE=1
    legend=False
    decade=True
    hue=None
    ci=[0.05,0.95]
    vframeag={'ntype':'I'}
    vframeeg={}
    opt=['PF']
    fig5, ax = plt.subplots(1,4,figsize=(12,4))
    #DATA
    def _sort(pdat,opt):
        return pdat[pdat['scenario'].str.contains(opt)]
    
    #PLOTS
    for k in range(len(opt)):
        #Hydropower value
        plot_data(_sort(Hydropower_val_Md,opt[k]),'Hydropower_val_Md',
                  vframe=vframeeg,legend=legend,ax=ax[0],ci=ci,decade=decade)
        #Irrigated value
        plot_data(_sort(Crop_val_Md,opt[k]),'Crop_val_Md',
                  vframe=vframeag,legend=legend,ax=ax[1],ci=ci,decade=decade)
        #Rainfed value
        plot_data(_sort(Crop_val_Md,opt[k]),'Crop_val_Md',
                  vframe={'ntype':'R'},legend=legend,ax=ax[2],ci=ci,decade=decade)
        #Hydropower prod
        plot_data(_sort(Hydropower_prod_GWh,opt[k]),'Hydropower_prod_GWh',
                  vframe=vframeeg,legend=legend,ax=ax[3],ci=ci,decade=decade)
        #Irrigated consumption
        #plot_data(_sort(Cul_Cons_Mm3,opt[k]),'Cul_Cons_Mm3',
        #          vframe=vframeag,legend=legend,ax=ax[4],ci=ci,decade=decade)       
    #axes
    ax[0].set_ylabel('Hydropower value (M$/y)')
    ax[1].set_ylabel('Irrigated crop value (M$/y)')
    ax[2].set_ylabel('Rainfed crop value (M$/y)')
    ax[3].set_ylabel('Hydropower production (TWh/y)')
    #ax[4].set_ylabel('Irrigated consumption (Mm3/y)')
    fig5.tight_layout()
    if SAVE==1:
        fig5.savefig('Sumdecadeplot_ifpri_idx',dpi=300)   
#%% Variable hydrology
VARHYDRO=0
if VARHYDRO==1:
    SAVE=1
    crops=['maize','cereals','cassava','roots']#['sugarcane','tobacco','cotton']#['rice','pulses','oilseeds'] #
    #['vegetable','fruits']#
    iterate=crops
    fig7, ax = plt.subplots(5,len(iterate),figsize=(14,12))
    hue=None
    ci=[0.05,0.95]
    legend=False
    def _sort(pdat,keyword,b=True):
        return pdat[pdat['scenario'].str.contains(keyword)==b]
    
    for k in range(len(iterate)):
        vframe={'ncrop':crops[k],'ncountry':['Zambia','Zimbabwe','Mozambique','Malawi']}
        #area
        plot_data(_sort(Crop_area_kha,'av5',b=False),'Crop_area_kha',vframe=vframe,hue=hue,legend=legend,ax=ax[0,k],ci=ci)
        plot_data(_sort(Crop_area_kha,'av5',b=True),'Crop_area_kha',vframe=vframe,hue=hue,legend=legend,ax=ax[0,k],ci=ci)
        plot_data(vAREAX0,'Val',vframe=vframe,columns=[],legend=legend,ax=ax[0,k])
        ax[0,k].set_title(iterate[k])
        #production
        plot_data(_sort(Crop_prod_kt,'av5',b=False),'Crop_prod_kt',vframe=vframe,hue=hue,legend=legend,ax=ax[1,k],ci=ci)
        plot_data(_sort(Crop_prod_kt,'av5',b=True),'Crop_prod_kt',vframe=vframe,hue=hue,legend=legend,ax=ax[1,k],ci=ci)
        plot_data(vPRODX0,'Val',vframe=vframe,columns=[],legend=legend,ax=ax[1,k])
        #yield
        plot_data(_sort(Crop_prod_kt,'av5',b=False),'Crop_prod_kt',vframe=vframe,hue=hue,group='sum',legend=legend,ax=ax[2,k],ci=ci,
                  refdat=_sort(Crop_area_kha,'av5',b=False),reftype='divide',refy='Crop_area_kha')
        plot_data(_sort(Crop_prod_kt,'av5',b=True),'Crop_prod_kt',vframe=vframe,hue=hue,group='sum',legend=legend,ax=ax[2,k],ci=ci,
                  refdat=_sort(Crop_area_kha,'av5',b=True),reftype='divide',refy='Crop_area_kha')
        plot_data(vPRODX0,'Val',vframe=vframe,group='sum',columns=[],legend=legend,ax=ax[2,k],
                  refdat=vAREAX0,reftype='divide')
        #price
        g1,g2=plot_data(_sort(Crop_price_dpt,'av5',b=False),'Crop_price_dpt',vframe=vframe,hue=hue,group='mean',legend=legend,ax=ax[3,k],ci=ci)
        g3,g4=plot_data(_sort(Crop_price_dpt,'av5',b=True),'Crop_price_dpt',vframe=vframe,hue=hue,group='mean',legend=legend,ax=ax[3,k],ci=ci)
        g5=plot_data(vPCX0,'Val',vframe=vframe,group='mean',columns=[],legend=legend,ax=ax[3,k])
        #value
        plot_data(_sort(Crop_val_Md,'av5',b=False),'Crop_val_Md',vframe=vframe,hue=hue,legend=legend,ax=ax[4,k],ci=ci)
        plot_data(_sort(Crop_val_Md,'av5',b=True),'Crop_val_Md',vframe=vframe,hue=hue,legend=legend,ax=ax[4,k],ci=ci)
        #set labels
        if k==0:
            ax[0,k].set_ylabel('Area (1000 ha)')
            ax[1,k].set_ylabel('Production (1000 t)')
            ax[2,k].set_ylabel('Yield (t/ha)')
            ax[3,k].set_ylabel('Price ($/t)')
            ax[4,k].set_ylabel('Value (M$)')
    #fig2.tight_layout()
    fig7.legend([g1,g2,g3,g4,g5], 
                labels=['Average climate change','Average climate change with rolling average',
                        'IMPACT model','Climate change between 5% and 95% percentiles',
                        'Climate change with rolling average between 5% and 95% percentiles'], 
                loc='lower center')
    
    if SAVE==1:
        fig7.savefig('Variable_hydrology_climate_impacts',dpi=300)

#Hydrology
    SAVE=1
    fig8, ax = plt.subplots(1,3,figsize=(12,4))
    legend=False #'auto'
    hue=None
    ci=[0.05,0.95]
    plot_data(_sort(P_mm,'av5',b=False),'P_mm',vframe={},hue=hue,legend=legend,ax=ax[0],ci=ci)
    plot_data(_sort(P_mm,'av5',b=True),'P_mm',vframe={},hue=hue,legend=legend,ax=ax[0],ci=ci)
    ax[0].set_title('Precipitation')
    ax[0].set_ylabel('mm/year')
    plot_data(_sort(ET_mm,'av5',b=False),'ET_mm',vframe={},hue=hue,legend=legend,ax=ax[1],ci=ci)
    plot_data(_sort(ET_mm,'av5',b=True),'ET_mm',vframe={},hue=hue,legend=legend,ax=ax[1],ci=ci)
    ax[1].set_title('ET0')
    ax[1].set_ylabel('mm/year')
    g1,g2=plot_data(_sort(Runoff_Mm3,'av5',b=False),'Runoff_Mm3',vframe={},hue=hue,legend=legend,ax=ax[2],ci=ci)
    g3,g4=plot_data(_sort(Runoff_Mm3,'av5',b=True),'Runoff_Mm3',vframe={},hue=hue,legend=legend,ax=ax[2],ci=ci)
    ax[2].set_title('Runoff')
    ax[2].set_ylabel('BCM/year')
    #fig1.tight_layout()
    fig8.legend([g1,g2,g3,g4], 
                labels=['Average climate change','Average climate change with rolling average',
                        'Climate change between 5% and 95% percentiles',
                        'Climate change with rolling average between 5% and 95% percentiles'], 
                loc=[0.14,0.12])
    if SAVE==1:
        fig8.savefig('variability_Hydrology_ifpri_idx',dpi=300,transparent=True)