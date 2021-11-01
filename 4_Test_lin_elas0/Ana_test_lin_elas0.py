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
#Von Mises Data
X1, Y1= import_dataset("dVM01.csv")
X2, Y2= import_dataset("dVM02.csv")
X3, Y3= import_dataset("dVM03.csv")
X4, Y4= import_dataset("dVM04.csv")
X5, Y5= import_dataset("dVM05.csv")

#GTN data
X1G, Y1G= import_dataset("dGN01.csv")
X2G, Y2G= import_dataset("dGN02.csv")
X3G, Y3G= import_dataset("dGN03.csv")
X4G, Y4G= import_dataset("dGN04.csv")
X5G, Y5G= import_dataset("dGN05.csv")

import os
path='Elastic'
if not os.path.exists(path):
    os.makedirs(path)

file_name="T_Iso_Stress_strain0.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(X1[1:],Y1[1:,0],linestyle='None',marker='o',markersize=5, label="s11_1VM")
ax.plot(X2[1:],Y2[1:,0],linestyle='None',marker='v',markersize=5, label="s11_2VM")
ax.plot(X3[1:],Y3[1:,0],linestyle='None',marker='v',markersize=5, label="s11_3VM") 
ax.plot(X4[1:],Y4[1:,0],linestyle='None',marker='v',markersize=5, label="s11_4VM")
ax.plot(X5[1:],Y5[1:,0],linestyle='None',marker='v',markersize=5, label="s11_5VM")

ax.plot(X1G[1:],Y1G[1:,0],linestyle='None',marker='x',markersize=5, label="s11_1GTN")
ax.plot(X2G[1:],Y2G[1:,0],linestyle='None',marker='.',markersize=5, label="s11_2GTN")
ax.plot(X3G[1:],Y3G[1:,0],linestyle='None',marker='.',markersize=5, label="s11_3GTN")
ax.plot(X4G[1:],Y4G[1:,0],linestyle='None',marker='.',markersize=5, label="s11_4GTN")
ax.plot(X5G[1:],Y5G[1:,0],linestyle='None',marker='.',markersize=5, label="s11_5GTN")
plt.grid(True)
ax.set_xlabel("strain, e11")
ax.set_ylabel("Stress")
ax.legend()     
fig.suptitle("Single tension")
fig.savefig(completeName)


