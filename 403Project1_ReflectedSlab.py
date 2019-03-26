#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 12:03:07 2018

@author: Romy Morsy
"""

'''NUCE403 Project 1 Submission For A Reflected Slab'''

import matplotlib.pyplot as plt
import numpy as np

#Diffusion Coeff of Material 1 & 2 for both groups
DiffCoeff=[[1.5154,1.9314],[0.3438,0.2659]]
'''
DiffCoeff Example
DiffCoeff[0][0]=Group 1 Material 1
DiffCoeff[0][1]=Group 1 Material 2
DiffCoeff[1][0]=Group 2 Material 1
DiffCoeff[1][1]=Group 2 Material 2
'''

''' Out of Group Scattering XS'''
OutofGroupXS=[[0.0290,0.035],[0.0854,0.0200]]
'''Nu*FissionXS '''
nusigF = [[0.0065,0],[0.1312,0]]
SigAbs = [[0.092,0.0007],[0.0841,0.0198]]


#GroupScattering of Material 1
Mat1Scat=[[0.5425,0],[0.0198,1.3727]]
'''
MatScat Example
MatScat[0][0]=Group 1 to Group 1
MatScat[0][1]=Group 2 to Group 1
MatScat[1][0]=Group 1 to Group 2
MatScat[1][1]=Group 2 to Group 2
'''
#GroupScattering of Material 2
Mat2Scat=[[0.6060,0],[0.0343,1.9222]]
Iterations = 1001
Hfuel = 50.0
Hrefl = 25.0
Nx=100.0
dx=(Hfuel + 2*Hrefl) / (Nx)
Dis=np.arange(-(Hfuel/2. + Hrefl), (Hfuel/2. + Hrefl),dx)
NumX = len(Dis)   #Amount of nodes throughout the whole system
k = 1.0
ChiG1 = 1
ChiG2 = 0
tolerancelimit1 = 0.000001
'''Diagnal of the A-Matrix'''
a=[
   [2*DiffCoeff[0][0]/(pow(dx,2))+OutofGroupXS[0][0],2*DiffCoeff[0][1]/(pow(dx,2))+OutofGroupXS[0][1]],
   [2*DiffCoeff[1][0]/(pow(dx,2))+OutofGroupXS[1][0],2*DiffCoeff[1][1]/(pow(dx,2))+OutofGroupXS[1][1]]]
'''Right and left side of the A-Matrix'''
aside=[
       [-DiffCoeff[0][0]/(pow(dx,2)),-DiffCoeff[0][1]/(pow(dx,2))],
       [-DiffCoeff[1][0]/(pow(dx,2)),-DiffCoeff[1][1]/(pow(dx,2))]]

'''Initializing a Group 1 Source with guess values of 0'''
G1source = np.zeros(NumX,dtype=float)

'''Initializing a Group 1 Flux with guess values of 0'''
G1flux = np.zeros(NumX,dtype=float)

'''Initializing a Fission Source with guess values of 0'''
FisSource = np.zeros(NumX, dtype=float)
'''Intializing a Group 2 Flux with guess values of 0'''
G2flux = np.zeros(NumX,dtype=float)
k_graph=[] #Empty list to append to to store the divergence of delta k
Sf_graph = [] #Empty list to append to to store diveregenec of delta Sf
'''Loop Over All Iterations To Converge On k'''
for i in range(0,Iterations):
    k_old = k
    
    '''Calculating A Group 1 Flux'''
    Stop = False
    while not Stop:
        for num,i in enumerate(G1flux):
            pre = i
            if Dis[num]>=-Hrefl and Dis[num]<=Hrefl:
                if num!=0 and num!=len(G1flux) and num!=len(G1flux)-1:
                    
                    post = (1/a[0][0])*(G1source[num]-(aside[0][0]*G1flux[num-1]+aside[0][0]*G1flux[num+1]))
                    G1flux[num]=post
                    tolerance=(post-pre)/post
                    if abs(tolerance)<tolerancelimit1:
                        Stop=True
            else:
                if num!=0 and num!=len(G1flux) and num!=len(G1flux)-1:
                    post= (1/a[0][1])*(G1source[num]-(aside[0][1]*G1flux[num-1]+aside[0][1]*G1flux[num+1]))
                    G1flux[num]=post
                    tolerance=(post-pre)/post
                    if abs(tolerance)<tolerancelimit1:
                        Stop=True
    
    
    '''Calculating A Group 2 Source'''
    tempS2=[]
    for num,i in enumerate(G1flux):
        if Dis[num]>=-Hrefl and Dis[num]<=Hrefl:
            f2=i*OutofGroupXS[1][0]
        else:
            f2=i*OutofGroupXS[1][0]
        tempS2.append(f2)
    SourceG2=np.asarray(tempS2)
    
    
    '''Calculating A Group 2 Flux'''
    tolerancelimit2 = 0.000001
    Stop = False
    while not Stop:
        for num,i in enumerate(G2flux):
            pre = i
            if Dis[num]>=-Hrefl and Dis[num]<=Hrefl:
                if num!=0 and num!=len(G2flux) and num!=len(G2flux)-1:
                    
                    post= (1/a[1][0])*(SourceG2[num]-(aside[1][0]*G2flux[num-1]+aside[1][0]*G2flux[num+1]))
                    G2flux[num]=post
                    tolerance=(post-pre)/post
                    if abs(tolerance)<tolerancelimit1:
                        Stop=True
            else:
                if num!=0 and num!=len(G2flux) and num!=len(G2flux)-1:
                    post= (1/a[1][1])*(SourceG2[num]-(aside[1][1]*G2flux[num-1]+aside[1][1]*G2flux[num+1]))
                    G2flux[num]=post
                    tolerance=(post-pre)/post
                    if abs(tolerance)<tolerancelimit2:
                        Stop=True
    
    
    '''Calculating A Fission Source'''
    for val,j in enumerate(G1flux):
        pre = j
        if val >= Hrefl and val<= Hfuel + Hrefl:
            FisSource[val] = G1flux[val]*nusigF[0][0] + G2flux[val]*nusigF[1][0]
            post = FisSource[val]
            toleranceF = (post-pre)/post  
            Sf_graph.append(abs(toleranceF))
            if abs(toleranceF)<tolerancelimit1:
                Stop=True
    k = sum(FisSource)/sum(G1source)  
    
    print('k:',k)
    G1source = FisSource/k
    k_graph.append(abs(((k_old-k)/k)))

    if abs((k_old-k))/k < tolerancelimit1:
        break
'''Calculations Done For The Project'''
ScatterSourceM1 = sum(FisSource)*Mat1Scat[1][0]
ScatterSourceM2 = sum(FisSource)*Mat2Scat[1][0]
G1FissionSource = sum(G1flux)*nusigF[0][0]
G2FissionSource = sum(G2flux)*nusigF[1][0]
AbsorptionRateM1G1 = sum(G1flux)*SigAbs[0][0]*dx  
AbsorptionRateM1G2 = sum(G2flux)*SigAbs[1][0]*dx  
AbsorptionRateM2G1 = sum(G1flux)*SigAbs[0][1]*dx  
AbsorptionRateM2G2 = sum(G2flux)*SigAbs[1][1]*dx
LeakageRateG1 = (DiffCoeff[0][0]*G1flux[len(G1flux)-1]/dx)+(DiffCoeff[0][1]*G1flux[1]/dx)
LeakageRateG2 = (DiffCoeff[1][0]*G1flux[len(G1flux)-1]/dx)+(DiffCoeff[1][1]*G1flux[1]/dx)

'''Plot Of DeltaK As A Function of Iteration # Reflected'''
plt.loglog(k_graph)
plt.xlabel("Iteration Number")
plt.ylabel("DeltaK")
plt.title("Log Plot Of The Decay Of DeltaK/K Reflected")
plt.show()
'''Plot Of Scalar Fluxes Together Reflected'''
plt.plot(G1flux, color = 'g', label = "Group 1 Flux")
plt.plot(G2flux, color = "orange", label = "Group 2 Flux")
plt.xlabel("Distance across the whole element including fuel and reflector")
plt.ylabel("Group 1 and Group 2 Neutron Fluxes")
plt.title("Plot Of Scalar Fluxes Together Reflected")
plt.legend()
plt.show()
'''Plot Of Fission Source And Downscatter Source Reflected'''
plt.plot(FisSource, color = 'g', label = "Fission Source")
plt.plot(SourceG2, color = "orange", label = "Scattering Source")
plt.xlabel("Distance across the whole element including fuel and reflector")
plt.ylabel("Fission Source and Scattering Sources")
plt.title("Plot Of Fission Source And Downscatter Source Reflected")
plt.legend()
plt.show()

#Source rate = sum of the source*dx
#Group 1 scattering*flux for G2Source rate
#Removal rate is the same as absorption rate you use RemovalXS
#leakage = diffusion*flux at the edges/dx
#ALL FROM HW5 SOLUTION
