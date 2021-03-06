import math
from typing import ValuesView
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def import_dataset(filename):
    x=[]
    y=[]
    with open (filename) as f: #Read all the lines
        for i, line in enumerate (f): 
            #if i==0: #As you know that in the first line there is text, for the numerical exercise you bypass this line
            #    continue
            s=line.split( ) #Slit all the line and save it in a vector
            s_num=[float(val) for val in s[:]]  #takes the float value of each data until the position last-1 and assig these values to other vector 
            x.append(s_num[0])  #append the values to the general array, at each iteration adds 1 new row
            y.append(s_num[1:]) #.split('\n')[0]) does not apply since the data come from a vector which contains only floa numbers
            
    x=np.array(x,dtype=float) #generates the numpy array specifiying that there are float data
    y=np.array(y,dtype=float)#generates the numpy array specifiying that there are integral data
    #print(x[0,0])
    
    return x,y
#X: Strain Values
#Y: dependent values
X, Y= import_dataset("dataPS.csv")

import os
path='2_Fig_PS'
if not os.path.exists(path):
    os.makedirs(path)

file_name="T_Iso_Stress_strainPS.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(X[1:],Y[1:,0],linestyle='None',marker='o',markersize=5, label="s11")
ax.plot(X[1:],Y[1:,1],linestyle='None',marker='x',markersize=5, label="s22") 
ax.plot(X[1:],Y[1:,2],linestyle='None',marker='.',markersize=5, label="s33")
ax.set_xlabel("strain, e11")
ax.set_ylabel("Stress")
ax.legend()
plt.grid(True)     
fig.suptitle("Single tension")
fig.savefig(completeName)


file_name="T_Iso_plastic_strainPS.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(X[1:],Y[1:,3],linestyle='None',marker='o',markersize=5, label="epsp11")
ax.plot(X[1:],Y[1:,4],linestyle='None',marker='x',markersize=5, label="epsp22") 
ax.plot(X[1:],Y[1:,5],linestyle='None',marker='.',markersize=5, label="epsp33")
ax.set_xlabel("Strain, e11")
ax.set_ylabel("plastic strain")
ax.legend()
plt.grid(True)  
fig.suptitle("Plastic strain")
fig.savefig(completeName)


file_name="T_Iso_PorosityPS.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(X[1:],Y[1:,6],linestyle='None',marker='o',markersize=5, label="f")
ax.set_xlabel("Strain, e11")
ax.set_ylabel("porosity")
ax.legend()
plt.grid(True)  
fig.suptitle("Porosity")
fig.savefig(completeName)


file_name="T_Iso_epsp_bPS.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(X[1:],Y[1:,7],linestyle='None',marker='o',markersize=5, label="epsp_b")
ax.set_xlabel("Strain, e11")
ax.set_ylabel("microscopy equivalent plastic strain")
ax.legend()   
plt.grid(True)  
fig.suptitle("microscopy equivalent plastic strain")
fig.savefig(completeName)