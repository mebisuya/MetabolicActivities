#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:00:09 2024

@author: lazaro
"""

#Import packages and set plotting parameters

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.metrics import r2_score


import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
plt.rc('font', family='Arial', size= 16)


################

#Variables to change

#Cropping of data for analysis of the most relevant part, input in minutes
crop_min = [0,450]
crop_min_plot = [0,450] #Only for plotting purposes, in some caes they have to be different, 
#especially when we want to normalise plots across conditions but only analize a certain part


################




#Import of data, first column is the time the rest are samples grouped by conditions, n=3 is assumed

df = pd.read_csv("human_degradation.csv", delimiter = ",")


#Analysis and preeliminary plotting
#Protein half-life in minutes and r-squared of fit are outputted for each sample

out = {}

data = []

for i in df.columns[1:]:
    thresh=0.05
    x = np.array(df["time"])
    timestep = x[1]-x[0]
    
    y = np.array(df[i])
    
    #Plot raw data to visualize where cropping is happening
    plt.figure(figsize=(5,5))
    plt.plot(x, y, label='original data')
    plt.axvline(x=crop_min[0])
    plt.axvline(x=crop_min[1])
    plt.show()

    y2 = y/y[0]
    data.append(np.array(y2, dtype=np.float64))
    
    y_plot = y
    x_plot = x
   
    y = y[round(crop_min[0]/timestep):round(crop_min[1]/timestep)]
    x = x[round(crop_min[0]/timestep):round(crop_min[1]/timestep)]
      
    y_plot = y_plot[round(crop_min_plot[0]/timestep):round(crop_min_plot[1]/timestep)]
    x_plot = x_plot[round(crop_min_plot[0]/timestep):round(crop_min_plot[1]/timestep)]
     
   
    y = y/y[0]
    y = np.log2(y.astype(float))
   
    y_plot = y_plot/y_plot[0]
    y_plot = np.log2(y_plot.astype(float))

    x= x.reshape(-1, 1)
    x_plot= x_plot.reshape(-1, 1)

    ransac = linear_model.RANSACRegressor(residual_threshold=thresh, random_state=5)
    ransac.fit(x, y)
    line_y_ransac = ransac.predict(x)


    inlier_mask = ransac.inlier_mask_
    positions = [i for i, x in enumerate(inlier_mask) if x]

    fig, ax = plt.subplots(figsize=(4, 4))
    plt.plot(x_plot, y_plot, ".", label='Data', color= "#cab2d6")
    plt.plot(x, line_y_ransac, 'r', label='Fitting')
    plt.ylabel("Log2(Luciferase activity)") 
    plt.xlabel("Time (min)")     
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.axvline(x=positions[0]*timestep, color='gray', ls='--', lw=2, alpha=0.7)
    plt.axvline(x=positions[-1]*timestep, color='gray', ls='--', lw=2, alpha=0.7)
    plt.text(0.72, 0.6, 'Slope = %0.4f' % ransac.estimator_.coef_.item(), horizontalalignment='center', verticalalignment='center',transform = ax.transAxes)
    plt.legend(markerscale=2.5, frameon=False)
    plt.savefig(i + ".pdf", bbox_inches='tight', transparent=True)
    plt.show()
    
    print(i)
    print(-1/ransac.estimator_.coef_)

    m, c = float(ransac.estimator_.coef_.item()), float(ransac.estimator_.intercept_.item())
    inlier_mask = ransac.inlier_mask_
    x_accpt, y_accpt = x[inlier_mask], y[inlier_mask]
    y_predict = m * x_accpt + c
    R2 = r2_score(y_accpt, y_predict)
    print(R2)
    
    out[i] = [-1/ransac.estimator_.coef_[0], R2]

    
pd.DataFrame(out).to_csv('Analysis.csv', index=False)


#Process and save the average and sd of the raw data grouped by condition (n=3)
#Bth a file with all traces and a file with the traces of each group will be generated
traces_all = {}

names= df.columns[1:]
for i in range(0,len(data),3):
    av_data = np.average([np.array(data[i]),np.array(data[i+1]),np.array(data[i+2])], axis=0)
    sd_data = np.std([np.array(data[i]),np.array(data[i+1]),np.array(data[i+2])], axis=0)
    traces = {names[i].split("_")[0]+"_Time":np.arange(len(av_data))*timestep, names[i].split("_")[0]+"_Average": av_data, names[i].split("_")[0]+"_Sd": sd_data}
    traces_all[names[i].split("_")[0]+"_Time"]= np.arange(len(av_data))*timestep
    traces_all[names[i].split("_")[0]+"_Average"]= av_data
    traces_all[names[i].split("_")[0]+"_Sd"]= sd_data
    pd.DataFrame(traces).to_csv(names[i].split("_")[0]+'_traces.csv', index=False)
    
    fig, ax = plt.subplots(figsize=(6, 4))
    plt.plot(np.arange(len(av_data))*timestep, av_data, '.')
    plt.fill_between(np.arange(len(av_data))*timestep, av_data-sd_data, av_data+sd_data,
        alpha=0.3)
    plt.ylabel("Luciferase activity") 
    plt.xlabel("Time (min)")     
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

pd.DataFrame(traces_all).to_csv('traces_all.csv', index=False)

