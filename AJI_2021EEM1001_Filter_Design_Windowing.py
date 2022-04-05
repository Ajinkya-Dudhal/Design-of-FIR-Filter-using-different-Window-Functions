#Filter Design and Windowing code
#Ajinkya Dudhal_2021EEM1001_Assignment-4

import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import cmath
import math

hd=[int(i) for i in input("Enter samples separated by space:").split(" ")]
M=int(input("Enter window length="))

#Rectangular Window
n=[]
rw=[]
rec_h=[]
for i in range(M):                      #Implementation of Rectangular Window
    n.append(i)
    rw.append(1)
for i in range(M):                      #Windowing
    a=hd[i]*rw[i]
    rec_h.append(a)
print(rec_h)                            #FIR filter after windowing

recH=[]                                 #Calculation of DFT of filter response
for k in range(M):
    sum=0
    for i in range(M):
        sum=sum + rec_h[i]*np.exp(complex(0,(-2*np.pi*k*i)/M))
    recH.append(sum)
print(recH)
mag=[]
phase=[]
for i in range(M):
    m=abs(recH[i])
    mag.append(m)
    p=cmath.phase(recH[i])
    ph=np.degrees(p)
    phase.append(ph)
plt.subplot(2,4,1)
plt.stem(n,mag)                         #Magnitude Plot
plt.title("Rec_Mag_Spectrum")
plt.grid()
plt.subplot(2,4,5)
plt.stem(n,phase)                       #Phase Plot
plt.title("Rec_Phase_Spectrum")
plt.grid()

#Triangular Window
tw=[]
tri_h=[]
for i in range(M):
    a=1-(2*abs(i-(M-1)/2)/(M-1))
    tw.append(a)
for i in range(M):
    a=hd[i]*tw[i]
    tri_h.append(a)

triH=[]                                 #Calculation of DFT of filter response
for k in range(M):
    sum=0
    for i in range(M):
        sum=sum + tri_h[i]*np.exp(complex(0,(-2*np.pi*k*i)/M))
    triH.append(sum)
print(triH)
mag=[]
phase=[]
for i in range(M):
    m=abs(triH[i])
    mag.append(m)
    p=cmath.phase(triH[i])
    ph=np.degrees(p)
    phase.append(ph)
plt.subplot(2,4,2)
plt.stem(n,mag)                         #Magnitude Plot
plt.title("Tri_Mag_Spectrum")
plt.grid()
plt.subplot(2,4,6)
plt.stem(n,phase)                       #Phase Plot
plt.title("Tri_Phase_Spectrum")
plt.grid()

#Raised Cosine Window
Ra_co_w=[]
Ra_co_h=[]
for i in range(M):
    a=0.5*(1-np.cos((2*np.pi*i)/(M-1)))
    Ra_co_w.append(a)
for i in range(M):
    a=hd[i]*Ra_co_w[i]
    Ra_co_h.append(a)

Ra_co_H=[]                               #Calculation of DFT of filter response
for k in range(M):
    sum=0
    for i in range(M):
        sum=sum + Ra_co_h[i]*np.exp(complex(0,(-2*np.pi*k*i)/M))
    Ra_co_H.append(sum)
print(Ra_co_H)
mag=[]
phase=[]
for i in range(M):
    m=abs(Ra_co_H[i])
    mag.append(m)
    p=cmath.phase(Ra_co_H[i])
    ph=np.degrees(p)
    phase.append(ph)
plt.subplot(2,4,3)
plt.stem(n,mag)                           #Magnitude Plot
plt.title("Raised_Cosine_Mag_Spectrum")
plt.grid()
plt.subplot(2,4,7)
plt.stem(n,phase)                         #Phase Plot
plt.title("Raised_Cosine_Phase_Spectrum")
plt.grid()

#Kaiser Window
kw=[]
kaiser_h=[]
z=(M-1)/2
alpha=1.33
for i in range(M):
    x1=alpha*np.sqrt((z*z)-np.square(i-z))
    x2=alpha*np.sqrt(z)
    a=scipy.special.i0(x1)/scipy.special.i0(x2)
    kw.append(a)
for i in range(M):
    a=hd[i]*kw[i]
    kaiser_h.append(a)
print(kaiser_h)

kaiser_H=[]                                #Calculation of DFT of filter response
for k in range(M):
    sum=0
    for i in range(M):
        sum=sum + kaiser_h[i]*np.exp(complex(0,(-2*np.pi*k*i)/M))
    kaiser_H.append(sum)
print(kaiser_H)
mag=[]
phase=[]
for i in range(M):
    m=abs(kaiser_H[i])
    mag.append(m)
    p=cmath.phase(kaiser_H[i])
    ph=np.degrees(p)
    phase.append(ph)
plt.subplot(2,4,4)
plt.stem(n,mag)                            #Magnitude Plot
plt.title("Kaiser_Mag_Spectrum")
plt.grid()
plt.subplot(2,4,8)
plt.stem(n,phase)                          #Phase Plot
plt.title("Kaiser_Phase_Spectrum")
plt.grid()
plt.show()