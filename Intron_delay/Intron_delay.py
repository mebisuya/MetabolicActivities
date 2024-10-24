#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 09:15:19 2024

@author: lazaro
"""

#Import packages and set plotting parameters

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D 

#Plot formatting
mpl.rcParams['pdf.fonttype'] = 42
plt.rc('font', family='Arial', size= 15)



##########
#Variables to change

timestep = 2 #minutes
crop_min = [25,600] #select analysis area

name = "Control"

#############

#Import of data, first FLuc then NLuc of each repeat. Only one condition at a time
#Data consists on detrended and amplitude normalised signal obtained using pyboat.


df_intron = pd.read_csv("intron_control.csv", delimiter = ",")
df_no_intron = pd.read_csv("no_intron_control.csv", delimiter = ",")



####################

def cross_corr(df, intron_presence):

    colnames= list(df.columns)

    df = df.dropna()
    # for i in df.columns:
        #     df[i] = np.array(df[i])/1000


    x = np.arange(len(df))*timestep

    fig, ax = plt.subplots(figsize=(6, 4))
    for i in df.columns:
        plt.plot(x, list(df[i]))
        plt.axvline(x=crop_min[0])
        plt.axvline(x=crop_min[1])

    crop=[]
    for i in df.columns:
        crop.append(list(df[i])[round(crop_min[0]/timestep):round(crop_min[1]/timestep)])
    df=pd.DataFrame(crop)
    df = df.transpose()
    df.columns = colnames

    x = np.arange(len(df))*timestep


    delay_list= []

    out={}
    plt.figure(figsize=(7,4))
    for i in range(0,len(df.columns),2):
        correlation = signal.correlate(df[colnames[i+1]]-np.mean(df[colnames[i+1]]), df[colnames[i]] - np.mean(df[colnames[i]]), mode="full")
        lags = signal.correlation_lags(len(df[colnames[i+1]]), len(df[colnames[i]]), mode="full")
        lag = lags[np.argmax(correlation)]
        delay = lag*timestep
        delay_list.append(delay)

    
        out[colnames[i]] = [delay]
        
    

    pd.DataFrame(out).to_csv(name + "_"+ intron_presence +'_Analysis.csv', index=False)
    return out


out_intron = cross_corr(df_intron, "intron")
out_no_intron = cross_corr(df_no_intron, "no_intron")

delay = np.concatenate(list((out_intron.values())))-np.concatenate(list((out_no_intron.values())))
out_dict = {name+"_1": delay[0],name+"_2": delay[1],name+"_3": delay[2]}
pd.DataFrame([out_dict]).to_csv(name +'_Delay_Results.csv', index= False)
    
