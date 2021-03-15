# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 16:36:38 2018

@author: Tiago
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import random
import os

def solucaoNumerica(U, _K1, _K2):
    #Diferenças Finitas
    for j in range (0, numpontosx):
        U[0,j] = ((j*DX)*(1-(j*DX)))
        
    U[0,numpontosx-1] = 0 
    #print(U[0,:])
        
    # iteracao
    for n in range (1, numpontost):
        U[n,0] = 0
        U[n,numpontosx-1] = 0
        for j in range (1, numpontosx -1):
            U[n,j] = S * (U[n-1,j+1] - 2*U[n-1,j] + U[n-1,j-1] 
                     + _K1*np.exp(-_K2*j*DT)) + U[n-1,j]
            #print(U[n,:])
    return U

def geraGrafico(X, T, U, nome):
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca(projection='3d')
    
    # Plot the surface.
    surf = ax.plot_surface(X, T, U, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    
    # Customize the z axis.
    #ax.set_zlim(-0.1, 1.0)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=6)
    
    
    ax.view_init(25, 45)
    ax.set_xlabel('X')
    ax.set_ylabel('T')
    ax.set_zlabel('UE')
    
    if not os.path.exists("../"+nome+"/"):
        os.makedirs("../"+nome+"/")
    
    plt.show()
    for c in range(0, 22):
        init = 0
        fig.savefig(nome+"/"+nome+str(c)+".png", dpi=100)
        ax.view_init(25, init + c*20)


k1 = 0.5
k2 = 0.5

E = 1
L = 1
TMAX = 0.2
numpontosx = 35
numpontost = 500

DX = L / numpontosx
#DT = 0.48 * np.power(DX, 2)
DT = TMAX/ numpontost



S = DT / np.power(DX, 2)
print (S)
X = np.arange(0, L, DX)
T = np.arange(0, TMAX, DT)
X, T = np.meshgrid(X, T)

UExata = X*T*0.0
URuido = X*T*0.0
UNum = X*T*0.0
UMin1 = X*T*0.0
UMin2 = X*T*0.0
UTemp = X*T*0.0
#Exata
print("Exata")
#Diferenças Finitas
for j in range (0, numpontosx):
    UExata[0,j] = ((j*DX)*(1-(j*DX)))
    
UExata[0,numpontosx-1] = 0 
#print(UExata[0,:])
    
# iteracao
for n in range (1, numpontost):
    UExata[n,0] = 0
    UExata[n,numpontosx-1] = 0
    for j in range (1, numpontosx -1):
        UExata[n,j] = S * (UExata[n-1,j+1] - 2*UExata[n-1,j] + UExata[n-1,j-1] + k1*np.exp(-k2*n*DT)) + UExata[n-1,j]
        #print(U[n,:])
        
#print(UExata)
geraGrafico(X, T, UExata, "Exata")
#Ruido
print("Com Ruido")
# iteracao
for n in range (1, numpontost):
    for j in range (1, numpontosx-1):
        URuido[n,j] = UExata[n,j] + (random.random())


geraGrafico(X, T, URuido, "Ruido")
      
#numérica
print("Numérica")

solucaoNumerica(UNum, k1, k2)
#print(UNum)
 
geraGrafico(X, T, UNum, "Numerica")

#minumos quadrados
print("Minimos A")
# iteracao
numeroTests = 100
min1 = 100000.0
k1_min = 0.0 
k2_min = 0.0
dif = 0

for c in range (1, numeroTests):
    _k1 = (1/numeroTests)*c
    _k2 = (1/numeroTests)*c
    solucaoNumerica(UTemp, _k1, _k2)
    dif = 0
    for n in range (0, numpontost):
        for j in range (0, numpontosx):
            dif = dif + np.power(UTemp[n,j] - URuido[n,j], 2)
    
    #print("k1, k2, dif", (DX)*c, (DX)*c, dif)  
    if (min1 > dif):
        min1 = dif
        k1_min = _k1 
        k2_min = _k2

    
print("Min, k1, k2", min1, k1_min, k2_min)  
solucaoNumerica(UMin1, k1_min, k2_min)
geraGrafico(X, T, UMin1, "Min1")


min1 = 100000
print("Minimos B")
# iteracao
for c in range (1, numeroTests):
    _k1 = random.random()
    _k2 = random.random()
    solucaoNumerica(UTemp, _k1, _k2)
    dif = 0
    for n in range (0, numpontost):
        for j in range (0, numpontosx):
            dif = dif + np.power(UTemp[n,j] - URuido[n,j], 2)
    
    #print("k1, k2, dif", _k1, _k2, dif)  
    if (min1 > dif):
        min1 = dif
        k1_min = _k1 
        k2_min = _k2
    
print("Min, k1, k2", min1, k1_min, k2_min)  
solucaoNumerica(UMin2, k1_min, k2_min)
geraGrafico(X, T, UMin2, "Min2")

for c in range (1, numpontosx):
    fig = plt.figure(figsize=(10,10))
    plt.plot(URuido[:, c])
    plt.plot(UNum[:, c])
    plt.plot(UMin1[:, c])
    plt.plot(UMin2[:, c])
    plt.gca().legend(("Ruido", "Numérico", "Min1", "Min Aleatório"))
    plt.show()
    fig.savefig("Grafico/Grafico"+str(c)+".png", dpi=100)

