#!/usr/bin/env python

"""mech_validate.py: A python utility for benchmarking a chemical mechanism with desired experimental datasets or user-defined conditions."""

__author__      = "Bulut Tekgül"
__copyright__   = "Copyright 2019, Aalto University, Finland"

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas as pd
from pylab import figure, text, scatter, show
import matplotlib
params = {
    'text.latex.preamble': ['\\usepackage{gensymb}'],
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 1000,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'font.size': 7, # was 10
    'legend.fontsize': 7, # was 10
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'text.usetex': True,
    'figure.figsize': [3.16, 2],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)

import os
if not os.path.exists('figures'):
    os.makedirs('figures')

def exp_ndod_Mao(mechanisms,names,linestyles):
    fig,ax = plt.subplots(1)

    for mechanism,name,linestyle in zip(mechanisms,names,linestyles):
        gas = ct.Solution(mechanism)

        # Temperature range
        T_exp = np.array([725,770,809,865,911,1130,949,1001,1069,1266,1215])
        Trange = np.linspace(700,1300,30)
        p = 1500000

        Temp_arr = []
        t_arr = []
        IDT_sim = []
        for T in Trange:

            gas.TPX =  T, p, {'O2': 20.88, 'N2': 78.56, 'NC12H26': 0.56} 
            #We use contant volume reactor as the experiment was constant volume
            r = ct.IdealGasReactor(gas)
            
            #Only a reactor network can be run in time
            sim = ct.ReactorNet([r])
            
            # How long to run the reactor
            tEnd = 10e-3
            
            #Set the maximum time step applied by "sim.step" function
            sim.set_max_time_step(1e-5)

            t = 0.0
            while t<tEnd:
                t_arr = np.append(t_arr,t)
                Temp_arr = np.append(Temp_arr,r.thermo.T)
                t = sim.step()
            dTdX = np.diff(Temp_arr)/np.diff(t_arr)   
            index_max = max(range(len(dTdX)), key=dTdX.__getitem__)
            IDT_2nd = t_arr[index_max]
            IDT_sim = np.append(IDT_sim,IDT_2nd)

        ax.semilogy(1000/Trange,IDT_sim*1000,label=name,linestyle=linestyle)


    EXPndod_IDT=np.array([3239,3079,3040,3224,3479,444,2297,1572,845,132,206])/1000
    #Experimental data is already implemented
    ax.semilogy(1000/T_exp,EXPndod_IDT,'ko', mfc='none', label='Mao et al. (2019)' )
    #ax.set_xlim([5.5,7.2])
    #ax.set_ylim([0.03,None])
    #text(0.8, 0.2,'P = 1.8 atm\n$Ar$=72.7$\%$\n$O_2$=18.2$\%$\n$CH_4$=9.1$\%$', ha='center', va='center', transform=ax.transAxes)
    ax.set_xlabel('1000/T [1/K]')
    ax.set_ylabel('Ignition Delay Time (ms)')
    ax.legend(loc='best',frameon=False)
    fig.savefig('figures/ndod_mao.pdf',format='pdf', bbox_inches='tight',pad_inches=0)


def exp_methane_SEERY(mechanisms,names,linestyles):
    fig,ax = plt.subplots(1)

    for mechanism,name,linestyle in zip(mechanisms,names,linestyles):
        gas = ct.Solution(mechanism)

        # Temperature range
        Trange = np.linspace(1380,2000,15)
        p = 1.8*ct.one_atm

        Temp_arr = []
        t_arr = []
        IDT_sim = []
        for T in Trange:

            gas.TPX =  T, p, {'O2': 0.182, 'N2': 0.727, 'CH4': 0.091} 
            #We use contant volume reactor as the experiment was constant volume
            r = ct.IdealGasReactor(gas)
            
            #Only a reactor network can be run in time
            sim = ct.ReactorNet([r])
            
            # How long to run the reactor
            tEnd = 10e-3
            
            #Set the maximum time step applied by "sim.step" function
            sim.set_max_time_step(1e-5)

            t = 0.0
            while t<tEnd:
                t_arr = np.append(t_arr,t)
                Temp_arr = np.append(Temp_arr,r.thermo.T)
                t = sim.step()
            dTdX = np.diff(Temp_arr)/np.diff(t_arr)   
            index_max = max(range(len(dTdX)), key=dTdX.__getitem__)
            IDT_2nd = t_arr[index_max]
            IDT_sim = np.append(IDT_sim,IDT_2nd)

        ax.semilogy(10000/Trange,IDT_sim*1000,label=name,linestyle=linestyle)


    EXPch4_IDT=np.array([75e-3, 88e-3, 170e-3, 240e-3, 180e-3, 430e-3, 670e-3, 550e-3, 900e-3, 1800e-3])
    EXPch4_temp=np.array([5.65, 5.80, 6.02, 6.14, 6.20, 6.41, 6.53, 6.57, 6.74, 7])
    #Experimental data is already implemented
    ax.semilogy(EXPch4_temp,EXPch4_IDT,'ko', mfc='none', label='Seery & Bowman. (1970)' )
    ax.set_xlim([5.5,7.2])
    ax.set_ylim([0.03,None])
    text(0.8, 0.2,'P = 1.8 atm\n$Ar$=72.7$\%$\n$O_2$=18.2$\%$\n$CH_4$=9.1$\%$', ha='center', va='center', transform=ax.transAxes)
    ax.set_xlabel('10000/T [1/K]')
    ax.set_ylabel('Ignition Delay Time (ms)')
    ax.legend(loc=2,frameon=False)
    fig.savefig('figures/CH4_seery.pdf',format='pdf', bbox_inches='tight',pad_inches=0)


def exp_ndod_MALEWICKI(mechanisms,names,linestyles):
    fig,ax = plt.subplots(nrows=2,ncols=3)


    xls = pd.read_excel("dataset/malewicki.xls",skiprows=1,sheet_name=3)

    #Experimental data is already implemented
    ax[0,0].plot(xls['T5 /K'],xls['NC12H26'],'ko', mfc='none', label='Exp ($C_{12}H_{26}$)' )

    ax[0,1].plot(xls['T5 /K'],xls['O2'],'ko', mfc='none', label='Exp ($O_2$)' )

    ax[0,2].plot(xls['T5 /K'],xls['CO'],'ko', mfc='none', label='Exp ($CO$)' )
    ax[1,0].plot(xls['T5 /K'],xls['CO2'],'ko', mfc='none', label='Exp ($CO_2$)' )

    ax[1,1].plot(xls['T5 /K'],xls['CH4'],'ko', mfc='none', label='Exp {$CH_4$}' )

    ax[1,2].plot(xls['T5 /K'],xls['C2H4'],'ko', mfc='none', label='Exp {$C_{2}H_4$}' )


    TEXP = xls['T5 /K']
    NDOD_EXP = xls['NC12H26']
    TIME_EXP = xls['REACTION TIME /ms']
    P_EXP = xls['P5 /atm']
    #ax.set_xlim([5.5,7.2])
    #ax.set_ylim([0.03,None])
    #text(0.8, 0.2,'P = 1.8 atm\n$Ar$=72.7$\%$\n$O_2$=18.2$\%$\n$CH_4$=9.1$\%$', ha='center', va='center', transform=ax.transAxes)
    #ax.set_xlabel('10000/T [1/K]')
    #ax.set_ylabel('Ignition Delay Time (ms)')
    for mechanism,name,linestyle in zip(mechanisms,names,linestyles):
        gas = ct.Solution(mechanism)

        # Temperature range
        t_arr = []
        IDT_sim = []
        ndod_conc = []
        o2_conc = []
        co_conc = []
        co2_conc = []
        ch4_conc = []
        c2h4_conc = []
        for T_exp, t_exp, p_exp in zip(TEXP,TIME_EXP,P_EXP):
            print(T_exp)
            y_ndod = 75.85*1e-6
            y_o2 = 683.6*1e-6
            y_ar = 1-(y_ndod+y_o2)
            gas.TPX =  T_exp, p_exp*ct.one_atm, {'NC12H26': y_ndod, 'O2': y_o2 , 'AR': y_ar} 
            #We use contant volume reactor as the experiment was constant volume
            r = ct.IdealGasReactor(gas)
            
            #Only a reactor network can be run in time
            sim = ct.ReactorNet([r])
            
        
            
            #Set the maximum time step applied by "sim.step" function
            sim.set_max_time_step(1e-5)
            tEnd = t_exp*1e-3
            t = 0.0
            while t<tEnd:
                t_arr = np.append(t_arr,t)
                t = sim.step()
            ndod_conc = np.append(ndod_conc,r.thermo.X[gas.species_index('NC12H26')]*1e6)
            o2_conc = np.append(o2_conc,r.thermo.X[gas.species_index('O2')]*1e6)
            ch4_conc = np.append(ch4_conc,r.thermo.X[gas.species_index('CH4')]*1e6)
            co_conc = np.append(co_conc,r.thermo.X[gas.species_index('CO')]*1e6)
            co2_conc = np.append(co2_conc,r.thermo.X[gas.species_index('CO2')]*1e6)
            c2h4_conc = np.append(c2h4_conc,r.thermo.X[gas.species_index('C2H4')]*1e6)
        ax[0,0].plot(TEXP,ndod_conc,label=name,linestyle=linestyle)
        ax[0,1].plot(TEXP,o2_conc,linestyle=linestyle)
        ax[0,2].plot(TEXP,co_conc,linestyle=linestyle)
        ax[1,0].plot(TEXP,co2_conc,linestyle=linestyle)
        ax[1,1].plot(TEXP,ch4_conc,linestyle=linestyle)
        ax[1,2].plot(TEXP,c2h4_conc,linestyle=linestyle)

    for i, ax in enumerate(fig.axes):
        ax.legend(loc='best',frameon=False)
        ax.set_xlabel('T [K]')
        ax.set_ylabel('ppm')
    #ax[0,0].legend(loc=1,frameon=False)
    fig.savefig('figures/ndod_malewicki.pdf',format='pdf', bbox_inches='tight',pad_inches=0)

def exp_ndod_VASU(mechanisms,names,linestyles):
    fig,ax = plt.subplots(1)

    for mechanism,name,linestyle in zip(mechanisms,names,linestyles):
        gas = ct.Solution(mechanism)

        # Temperature range
        Trange = np.linspace(700,1200,30)
        p = 20*ct.one_atm

        Temp_arr = []
        t_arr = []
        IDT_sim = []
        for T in Trange:

            gas.TPX =  T, p, {'NC12H26': 0.565, 'O2': 20.89, 'N2': 78.55} 
            #We use contant volume reactor as the experiment was constant volume
            r = ct.IdealGasReactor(gas)
            
            #Only a reactor network can be run in time
            sim = ct.ReactorNet([r])
            
            # How long to run the reactor
            tEnd = 10e-3
            
            #Set the maximum time step applied by "sim.step" function
            sim.set_max_time_step(1e-5)

            t = 0.0
            while t<tEnd:
                t_arr = np.append(t_arr,t)
                Temp_arr = np.append(Temp_arr,r.thermo.T)
                t = sim.step()
            dTdX = np.diff(Temp_arr)/np.diff(t_arr)   
            index_max = max(range(len(dTdX)), key=dTdX.__getitem__)
            IDT_2nd = t_arr[index_max]
            IDT_sim = np.append(IDT_sim,IDT_2nd)

        ax.semilogy(1000/Trange,IDT_sim*1000,label=name,linestyle=linestyle)
    EXP_vasu_IDT_05=np.array([2.243,1.522,1.918,2.134,2.276,1.245,0.826,0.839,0.788,0.566,0.527,0.266,0.346,0.261,0.213,0.165])
    EXP_vasu_temp_05=np.array([747,806,873,910,943,996,1039,1049,1054,1092,1097,1117,1118,1149,1163,1177])

    EXP_vasu_temp_1 = np.array([727,773,818,822,855,869,907,942,953,957,976,978,987,991,992,1008,1015,1036])
    EXP_vasu_IDT_1 = np.array([0.809,0.556,0.881,0.805,0.875,1.040,1.081,1.116,1.141,1.064,0.8,0.912,0.699,0.645,0.739,0.570,0.508,0.432])

    #Experimental data is already implemented
    ax.semilogy(1000/EXP_vasu_temp_05,EXP_vasu_IDT_05,'ko', mfc='none', label='Vasu et al. (2009)' )
    text(0.2, 0.85,'P = 20 atm\n$\phi$=0.5', ha='center', va='center', transform=ax.transAxes)

    ax.set_xlim([0.8,1.5])
    ax.set_ylim([0.1,10])

    ax.set_xlabel('1000/T [1/K]')
    ax.set_ylabel('Ignition Delay Time (ms)')
    ax.legend(loc='best',frameon=False)
    fig.savefig('figures/ndod_vasu_1.pdf',format='pdf', bbox_inches='tight',pad_inches=0)

    fig,ax = plt.subplots(1)

    for mechanism,name,linestyle in zip(mechanisms,names,linestyles):
        gas = ct.Solution(mechanism)
        Temp_arr = []
        t_arr = []
        IDT_sim = []
        for T in Trange:

            gas.TPX =  T, p, {'NC12H26': 1.123, 'O2': 20.77, 'N2': 78.10} 
            #We use contant volume reactor as the experiment was constant volume
            r = ct.IdealGasReactor(gas)
            
            #Only a reactor network can be run in time
            sim = ct.ReactorNet([r])
            
            # How long to run the reactor
            tEnd = 10e-3
            
            #Set the maximum time step applied by "sim.step" function
            sim.set_max_time_step(1e-5)

            t = 0.0
            while t<tEnd:
                t_arr = np.append(t_arr,t)
                Temp_arr = np.append(Temp_arr,r.thermo.T)
                t = sim.step()
            dTdX = np.diff(Temp_arr)/np.diff(t_arr)   
            index_max = max(range(len(dTdX)), key=dTdX.__getitem__)
            IDT_2nd = t_arr[index_max]
            IDT_sim = np.append(IDT_sim,IDT_2nd)

        ax.semilogy(1000/Trange,IDT_sim*1000,label=name,linestyle=linestyle)
    ax.semilogy(1000/EXP_vasu_temp_1,EXP_vasu_IDT_1,'ko', mfc='none', label='Vasu et al. (2009)' )
    ax.set_xlim([0.8,1.5])
    ax.set_ylim([0.1,10])
    text(0.8, 0.2,'P = 20 atm\n$\phi$=1', ha='center', va='center', transform=ax.transAxes)

    ax.set_xlabel('1000/T [1/K]')
    ax.set_ylabel('Ignition Delay Time (ms)')

    fig.savefig('figures/ndod_vasu_2.pdf',format='pdf', bbox_inches='tight',pad_inches=0)
def main():
    mechanisms = ['mechanisms/POLIMI_TOT_1412.cti','mechanisms/Yao54.cti','mechanisms/atmadeep.cti']
    names = ['POLIMI detailed','Yao','Atma']
    linestyles = ['-','--','dotted']

    #exp_ndod_VASU(mechanisms,names,linestyles)
    #exp_methane_SEERY(mechanisms,names,linestyles)
    #exp_ndod_MALEWICKI(mechanisms,names,linestyles)
    exp_ndod_Mao(mechanisms,names,linestyles)
main()

