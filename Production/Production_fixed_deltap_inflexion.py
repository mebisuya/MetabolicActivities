#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 17:46:09 2021

@author: lazaro
"""


#Import packages and set plotting parameters
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import basinhopping
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
import matplotlib as mpl

from matplotlib.lines import Line2D 

mpl.rcParams['pdf.fonttype'] = 42
plt.rc('font', family='Arial', size= 16)


##########
#Variables to change

#Names of conditions. MUST match condition number
names = ["Control", "2DG", "Azide", "Human"] 


prot_deg = [20.37, 32.19, 22.42, 42.15] #halflife, use same order as samples in file

n = 3 #repeats per experiment






#Import of data, first column is the time the rest are samples grouped by conditions
df = pd.read_csv("data_production.csv", delimiter = ",")



#################
alpha = 0.01
tau = 20
delta_m = 0.1
delta_p = 0.01

def pexp(t):
    p = alpha*(np.exp(-delta_m*(t-tau))-delta_m/delta_p *
           np.exp(-delta_p*(t-tau))) / (delta_m-delta_p) + alpha/delta_p
    p[np.where(t<tau)] = 0
    return p


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

##################

con= len(names)
names_plot = names

out = []
out_in =[]

data = []

plen = range(len(prot_deg))

k=0
m= 1

for j,i in enumerate(df.columns[1:]):
    # print(i)
    
    
    time = np.array(df["time"])
    timestep = time[1]-time[0]
    
    y = np.array(df[i])
    y = y[~np.isnan(y)]
    
    time = np.array(range(len(y)))*timestep
    
    if m <= n:       
        delta_p = np.log(2)/prot_deg[k]
        m=m+1
        
    else :
        k = k +1
        delta_p = np.log(2)/prot_deg[k]
        m=2
    
    plt.figure(figsize=(5,5))
    plt.plot(time, y, label='original data')
    plt.title("RAW data " +i)
    plt.show()
    
    y_av = moving_average(y, 5)
    time_av = time[2:-2]
    der = np.gradient(y_av)
    ip = time_av[np.argmax(der)]
    
    out_in.append(ip)
    
    plt.figure(figsize=(5,5))
    plt.plot(time_av, der, label='original data')
    plt.title("Derivative")
    plt.axvline(x=ip)
    plt.show()
    
    plt.figure(figsize=(5,5))
    plt.plot(time, y, label='original data')
    plt.title("RAW data, cropping")
    plt.axvline(x=0)
    plt.axvline(x=ip*2)
    plt.show()
    

    y = y[0:round(ip*2/timestep)]
    time = time[0:round(ip*2/timestep)]
    
    y = list(y/y[-1])
    
    data.append(np.array(y, dtype=np.float64))

    def fcost(x):
        fmod = x[0]*(np.exp(-x[2]*(time-x[1]))-x[2]/delta_p * np.exp(-delta_p*(time-x[1]))) / (x[2]-delta_p) + x[0]/delta_p
        fmod[np.where(time < x[1])] = 0
        return np.sum((fmod-y)**2)
    x0 = [0.01, 20, 0.1]
    minimizer_kwargs = {"method": "BFGS"}
    ret = basinhopping(fcost, x0, T=10, minimizer_kwargs=minimizer_kwargs, niter=500)
    alpha = ret.x[0]
    tau = ret.x[1]
    delta_m = ret.x[2]
    
    
    plt.figure(figsize=(6,6))
    plt.plot(time, pexp(time), label='fit', color='r')
    plt.plot(time, y, label='experiment', color='k')
    plt.xlabel('time (min)', fontsize=20)
    plt.ylabel('Intensity', fontsize=20)
    plt.title(i + '\nalpha, tau, delta_m, delta_p=\n' + str(ret.x) + " " + str(round(delta_p, 4)))
    plt.gca().tick_params(axis='both', labelsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig(i + "_inflexion_fixdeltaP.pdf", bbox_inches='tight', transparent=True)
   
    out.append(tau)


dict_out = {}

data2=[]
for i in range(0,len(out),3):
    av_data = np.average([np.array(out[i]),np.array(out[i+1]),np.array(out[i+2])], axis=0)
    data2.append([av_data,np.array(out[i]),np.array(out[i+1]),np.array(out[i+2]),np.array(out_in[i]),np.array(out_in[i+1]),np.array(out_in[i+2])])
    
for i in range(len(names)):
    dict_out[names[i]]= data2[i]

pd.DataFrame(dict_out,index=["Average", "tau_1", "tau_2", "tau_3", "ip_1", "ip_2", "ip_3"]).to_excel('production_tau_fixdeltaP_inflexion.xlsx', index=True)
