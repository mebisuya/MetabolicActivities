#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:53:46 2024

@author: lazaro
"""
#Import packages and set plotting parameters
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import pandas as pd




#Variables to change

#Name of the conditions
names = ["mouse", "2DG (10mM)", "Azide (1mM)", "D+A", "D+O", "human", "30C"]

#Real period, order is the same as conditions
exp_period = [152.261, 191,204, 254.515, 248.409, 324.7493333, 329.983]


#Parameters to simulate the period
Intron_delay = [9.333333333, 8.666666667, 20, 18.66666667, 11.33333333, 18, 34.66666667]
prot_deg = [20.36666667, 32.18333333, 22.42, 32.41573382, 66.5569703, 42.15489789, 40.55568795] #HalfLife
tau = [16.2489483, 21.22154159, 16.88102778	, 20.66107058, 16.56910349, 23.97703515, 36.20153438] 
mrna_deg = [0.051944147, 0.051662443, 0.053912127, 0.032581193, 0.089570337, 0.042119318, 0.010418143] 



###################


#Function to simulate periods, Based on Lewis 2003, Current Biology
#For more detailes explanations on the simulations, see Matsuda 2020, Science
def simulate(alpha,Deg_m,Deg_p,Tp,Tm, beta_0,beta_f,K,n,total_t): 
    m = [0]
    p = [0]
    for t in range(total_t):
        if (t - Tp) < 0:
            x= 0
        else:
            x= t-Tp
        if (t-Tm)  < 0:
            y = 0
        else:
            y= t-Tm        
        dp = alpha*m[x]-Deg_p*p[t]
        dm = beta_0 + beta_f/(1+p[y]**n/K**n)-Deg_m*m[t]
        if p[t] + dp<0:
            p.append(0)
        else:
            p.append(p[t] + dp)
        if m[t] + dm<0:
            m.append(0)
        else:
            m.append(m[t] + dm)
    return m,p


def sim_range(i):
    alpha= 4.5
    beta_0 = 0
    beta_f= 33
    K = 40
    n=3


    Delay_Rp = 0
    Delay_TxTl = tau[i]
    Delay_Tx = Delay_TxTl/2
    Delay_In = Intron_delay[i] 
    Deg_m = mrna_deg[i]
    Deg_p = np.log(2)/prot_deg[i]
    Delay_Tl = Delay_TxTl - Delay_Tx
    Tm= round(Delay_Rp + Delay_Tx + Delay_In)
    Tp= round(Delay_Tl)


    total_t = 4000 #time of simulation

    m,p = simulate(alpha,Deg_m,Deg_p,Tp,Tm, beta_0,beta_f,K,n,total_t)

    time = range(total_t)
    
    #We simulate the oscillations for a certain time but then crop the intiial part as we  need the oscillations to stabilize 
    #The lack of initial stability is produced by starting the expression levels at zero

    crop= 1000

    data = np.array(p[crop:])
    time = np.array(time[crop:])

    #Simulated period is estimated by measuring the distance between peaks
    max_peak = signal.argrelmax(data, order=5)
    maxpeaktime = time[max_peak]

    if len(maxpeaktime) % 2 != 0:
        maxpeaktime= np.delete(maxpeaktime, -1)
    
    periods= np.array([y - x for x,y in zip(maxpeaktime,maxpeaktime[1:])])
    period= np.average(periods)
    
    
    #Plotting
    plt.figure(figsize = (8,5))
    plt.title(names[i] + " Period: "+ str(round(exp_period[i], 1)) + " Simulated: "+ str(round(period, 1) ))
    plt.plot(data, color ="blue")
    plt.text(0.95, 0.7, "Intron delay: "+ str(round(Intron_delay[i], 1))+ "\nProt Halflife: "+ str(round(prot_deg[i], 1))+"\nTau: "+ str(round(tau[i], 1))+ "\nDelta mRNA: "+ str(round(mrna_deg[i], 4)) , fontsize=16, transform=plt.gcf().transFigure)
    plt.savefig(names[i] +".pdf", bbox_inches='tight', transparent=True)
    
    
    return period

simulated_periods = {}

for i in range(len(names)):
    simulated_periods[names[i]]= sim_range(i)
    
print(simulated_periods)
pd.DataFrame(simulated_periods, index=[0]).to_excel('Simulated_periods.xlsx', index=False)