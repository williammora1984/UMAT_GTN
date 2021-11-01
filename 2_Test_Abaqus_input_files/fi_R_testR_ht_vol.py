#Figures_test_hidrostatic test volumetric

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

def import_dataset2(filename):
    x=[]
    y=[]
    with open (filename) as f: #Read all the lines
        for i, line in enumerate (f): 
            #if i==0: #As you know that in the first line there is text, for the numerical exercise you bypass this line
            #    continue
            s=line.split( ) #Slit all the line and save it in a vector
            s_num=[float(val) for val in s[:]]  #takes the float value of each data until the position last-1 and assig these values to other vector 
            #x.append(s_num[0])  #append the values to the general array, at each iteration adds 1 new row
            y.append(s_num) #.split('\n')[0]) does not apply since the data come from a vector which contains only floa numbers
            
    #x=np.array(x,dtype=float) #generates the numpy array specifiying that there are float data
    Y1=np.array(y,dtype=float)#generates the numpy array specifiying that there are integral data
    #print(x[0,0])
    #Order the data in a tabular form
    si=Y1.shape[0]
    n_data=int(si/11)
    a=0
    for i in range (0,n_data):
        y2=Y1[a]
        for j in range(a+1,a+10):
            y2=np.append(y2,Y1[j])
        if i==0:
            y3=y2
        else:
            if i==1:
                y3=np.append([y3],[y2],0)
            else:     
                y3=np.append(y3,[y2],0)
        a=a+11

    return y3

epsv_l=np.array([0.042, 0.084, 0.126, 0.166, 0.207, 0.246, 0.285, 0.324
         ,0.362, 0.400]) 
p_l=np.array([201.4, 172.8, 148.7, 128.8, 112.1, 98.1, 86.3, 76.3,
         67.8, 60.4]) 
f_l=np.array([0.0728, 0.1130, 0.1524, 0.1911, 0.2287, 0.2645, 0.2980, 
        0.3290, 0.3575, 0.3838])

#Reference_data
X, Y= import_dataset("Reference_data/data_test_ht_vol.csv")

'''====================================================================
    Data for tension test e3 direction
======================================================================='''

#Import the data exported formm the abaqus analysis 
y3= import_dataset2("Results_abaqus/testR_te_vol.feh")
#Y1= import_dataset2("Reference_data/testR_ht_vol.feh")

#take the strain eps33 and and sig_33
strain_FE33=y3[32::8,2]
stress_FE33=y3[32::8,5]
#Compute the valumetric logaritmic strain
Lstrain_FE33=3*np.log(1+y3[32::8,2])

#Take the ISV (void volume fraction equivalent microscopic plastic strain
f_FE33=y3[32::8,22]
epsp_b33=y3[32::8,21]

'''====================================================================
    Data for tension test e2 direction
======================================================================='''
#Import the data exported formm the abaqus analysis 
y2= import_dataset2("Results_abaqus/testR_te_vol2.feh")
#take the strain eps22 and and sig_33
strain_FE22=y2[32::8,1]
stress_FE22=y2[32::8,4]
#Compute the valumetric logaritmic strain
Lstrain_FE22=3*np.log(1+y2[32::8,1])

#Take the ISV (void volume fraction equivalent microscopic plastic strain
f_FE22=y2[32::8,22]
epsp_b22=y2[32::8,21]

'''====================================================================
                    Generation of plots
======================================================================='''

import os
path='Fig_test_ht_vol'
file_name="T_ht_vol_hydsig_vs_Teps.png"
completeName = os.path.join(path, file_name)

if not os.path.exists(path):
    os.makedirs(path)

#stress a function of the logarithmic strain
fig,ax=plt.subplots()
#ax.plot(X,Y[:,0],linestyle='None',marker='o',markersize=5)
#ax.plot(epsv_l,p_l,linestyle='None',marker='s',markersize=5, label="Reference")
ax.plot(Lstrain_FE33,stress_FE33,linestyle='None',marker='o',markersize=5, label="UMAT_FEM33")
ax.plot(Lstrain_FE22,stress_FE22,linestyle='None',marker='x',markersize=2, label="UMAT_FEM22")

ax.set_xlabel("Logaritmic volumetric strain, $\epsilon_v$, [-]")   
ax.set_ylabel("-p/ $\sigma$ []")
plt.grid(True)
ax.legend()    
fig.suptitle("hydrostatic stress as a function of the logarithmic volumetric strain)")
#fig.savefig(completeName)
fig.savefig("Fig_test_ht_vol/T_ht_vol_hydsig_vs_Teps.png")


#void fraction as a function of the logarithmic strain
file_name="T_ht_vol_void_vs_epsvol.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
#ax.plot(X,Y[:,1],linestyle='None',marker='o',markersize=5)
#ax.plot(epsv_l,f_l,linestyle='None',marker='o',markersize=5, label="Reference")
ax.plot(Lstrain_FE33,f_FE33,linestyle='None',marker='o',markersize=5, label="UMAT_FEM33")
ax.plot(Lstrain_FE22,f_FE22,linestyle='None',marker='x',markersize=2, label="UMAT_FEM22")
ax.set_xlabel("Logaritmic volumetric strain, $\epsilon_v$, [-]")
ax.set_ylabel("void volume fraction, f, [%]")
plt.grid(True)
ax.legend()    
fig.suptitle("void volume fraction as a function of the logarithmic volumetric strain")
fig.savefig(completeName)

#equivalent plastic strain as a function of the logarithmic strain
file_name="T_ht_vol_epsp_b_vs_epsvol.png"
completeName = os.path.join(path, file_name)

fig,ax=plt.subplots()
#ax.plot(epsv_l,f_l,linestyle='None',marker='o',markersize=5, label="Reference")
ax.plot(Lstrain_FE33,epsp_b33,linestyle='None',marker='x',markersize=5, label="UMAT_FEM33")
ax.plot(Lstrain_FE22,epsp_b22,linestyle='None',marker='o',markersize=2, label="UMAT_FEM22")
ax.set_xlabel("Logaritmic volumetric strain, $\epsilon_v$, [-]")
ax.set_ylabel("Microscopic equivalent plastic, $bar \epsilon^p$, [-]")
plt.grid(True)
ax.legend()    
#fig.suptitle("Microscopic equivalent plastic strain as a function of the logarithmic volumetric strain")
fig.savefig(completeName)
