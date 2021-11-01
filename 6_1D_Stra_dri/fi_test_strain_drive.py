import math
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
    
X, Y= import_dataset("data_str_drive.csv")

import os
path='1_1D_VM'
if not os.path.exists(path):
    os.makedirs(path)

file_name="T1D_Loading_unloading.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(X,Y[:,0],linestyle='None',marker='o',markersize=5)
ax.set_ylabel("Stess, [%]")
ax.set_xlabel("time, [s]")
plt.grid(True)     
#fig.suptitle("Uniaxial_tension_(Loading/unloading)")
fig.savefig(completeName)


file_name="T1D_stress_strain_response.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(Y[:,0],Y[:,1],linestyle='None',marker='o',markersize=5)
ax.set_xlabel("e11, strain")
ax.set_ylabel("S11, stress")
plt.grid(True)
fig.suptitle("Uniaxial_tension_(stress-strain response)")
fig.savefig(completeName)


file_name="T1D_Lateral_Strain.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(Y[:,0],Y[:,2],linestyle='None',marker='o',markersize=5,label="eps22")
ax.plot(Y[:,0],Y[:,3],linestyle='None',marker='x',markersize=5,label="eps33")
ax.set_xlabel("e11, strain")
ax.set_ylabel("eps, lateral plastic strain")
fig.suptitle("Uniaxial_tension_(Lateral plastic Strain)")
plt.grid(True)
ax.legend()
fig.savefig(completeName)


file_name="T1D_Plastic_Strain.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
ax.plot(Y[:,0],Y[:,4],linestyle='None',marker='o',markersize=5,label="epsp11")
ax.plot(Y[:,0],Y[:,5],linestyle='None',marker='x',markersize=5,label="epsp22")
ax.plot(Y[:,0],Y[:,6],linestyle='None',marker='o',markersize=5,label="epsp33")
ax.set_xlabel("e11")
ax.legend()
ax.set_ylabel("epsp")     
fig.suptitle("Uniaxial_tension_(Plastic_strain)")
plt.grid(True)
fig.savefig(completeName)

#fig,ax=plt.subplots()
#ax.plot(X[1:],Y[1:,6],linestyle='None',marker='o',markersize=5, label="f")
#ax.set_ylabel("porosity")
#ax.set_xlabel("Strain, e11")
#ax.legend()     
#fig.suptitle("Porosity")
#fig.savefig("T1D_Porosity.png")
#
#fig,ax=plt.subplots()
#ax.plot(X[1:],Y[1:,7],linestyle='None',marker='o',markersize=5, label="epsp_b")
#ax.set_ylabel("microscopy equivalent plastic strain")
#ax.set_xlabel("microscopy equ pla strain, epsp_b")
#ax.legend()     
#fig.suptitle("microscopy equivalent plastic strain")
#fig.savefig("T1D_epsp_b.png")