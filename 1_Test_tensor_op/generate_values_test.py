import numpy as np
import csv
from random import seed
from random import random
from numpy import random


#Linspace without include the end point
n_test=2
end=np.random.rand(n_test)
difference_start_end=np.random.rand(n_test)
max_n_steps=20
steps=(random.rand(n_test)*max_n_steps)//1
export=np.zeros((n_test,3+max_n_steps))
steps=np.int_(steps)
print(steps)
for i in range(0,n_test):
    export[i,0]=end[i]-difference_start_end[i]
    export[i,1]=end[i]
    export[i,2]= steps[i]
    a=np.linspace(end[i]-difference_start_end[i],end[i],steps[i],0)
    export[i,3:(3+steps[i])]=a

#print(export)

'''============================================
Write a small matrix in vector form
==============================================='''

a = np.arange(25).reshape(5,5)
#open the file in the write mode

with open('data_TT01.csv', 'w') as f:
    # create the csv writer
    writer = csv.writer(f)

    # write a row to the csv file
    for i in range (a.shape[0]):
        writer.writerow(a[i,:])

'''============================================
Diadic Product test data
==============================================='''

a = np.arange(2,6.5,0.5).reshape(3,3)
b = np.arange(2,6.5,0.5).reshape(3,3)
inv_a=np.linalg.inv(a)
inv_b=np.linalg.inv(b)
Dia2=np.einsum('ij,ij',a,b)
Trace_a=np.einsum('ii',a)

a = a.reshape(9)
b = b.reshape(9)
inv_a=inv_a.reshape(9)
inv_b=inv_b.reshape(9)
d=np.append(a,b)
export_data1=np.append(d,Trace_a)
export_data1=np.append(export_data1,Dia2)
export_data1=np.append(export_data1,inv_a)
export_data1=np.append(export_data1,inv_b)

iter=10
for i in range (iter):
    a = np.random.rand(3,3)*100
    b = np.random.rand(3,3)*200
    inv_a=np.linalg.inv(a)
    inv_b=np.linalg.inv(b)
#Diadic product using einstein summation
    Dia2=np.einsum('ij,ij',a,b)
#Trace using einstein sumation
    Trace_a=np.einsum('ii',a)
    a = a.reshape(9)
    b = b.reshape(9)
    inv_a=inv_a.reshape(9)
    inv_b=inv_b.reshape(9)
#Columns 0 to 17 , tensors A and B (3x3) in 1D (9x1).
    mat_v=np.append(a,b)
#Colum 18: Trace
    f=np.append(mat_v,Trace_a)
#Column 19: diadic product
    f=np.append(f,Dia2)
#Column 20 to 28: inverse a
    f=np.append(f,inv_a)
#Column 29 to 37: inverse b
    f=np.append(f,inv_b)
    
    if i== 0:
        export_data1=np.append([export_data1],[f],0)
    else:
        export_data1=np.append(export_data1,[f],0)

#open the file in the write mode
with open('data_TT02.csv', 'w') as f:
    # create the csv writer
    writer = csv.writer(f)

    # write a row to the csv file
    for i in range (iter+1):
        writer.writerow(export_data1[i,:])

a = np.arange(9.).reshape(3,3)

#diadic_prod_T2_T2
Dia4=np.einsum('ij,kl->ijkl',a,a)
#print (Dia4)
Dia4=Dia4.reshape(81)
#print (Dia4)

#contrac_4th_2nd
a = np.arange(81.).reshape(3,3,3,3)
b = np.arange(9.).reshape(3,3)
Con4_2=np.einsum('kl,ijkl->ij',b,a)
print(Con4_2)