import math
import numpy as np
import matplotlib.pyplot as plt
from riemann_scal import *

x0 = 0; y0 = 0 # On définit l'origine du repère 2D (coordonnées de la première face)
# Domaine de résolution
Lx = 1 # intervalle [x0,Lx] selon x
Ly = 1 # intervalle [y0,Ly] selon y

LL = min(Lx,Ly)
mylevel = 3 # Paramètre entier à modifier pour modifier le pas du maillage
            # Vaut 5 dans le code scilab
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de faces selon la direction x
Ny = NN * int(Ly/LL) # nombre de faces selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

neighbours = 4 # nombre de voisins par cellule

def ux(i):
    return 0
    
def uy(i):
    return 2

codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')

qini = np.zeros(codim0.shape[0])

qini[0]=1

qhist = [qini.copy()]

q0 = qini.copy()

q1 = qini.copy()

t = 0
Tfinal = 10

while t<Tfinal:
    k = 1000
    flux=np.zeros(codim0.shape[0])
    for iface in range(codim1.shape[0]):
        i=codim0to1A[iface]
        j=codim0to1B[iface]
        if i<0 or j<0:
            continue
        lambdai = ux(i)*codim0to1NX[iface]+uy(i)*codim0to1NY[iface]
        lambdaj = ux(j)*codim0to1NX[iface]+uy(j)*codim0to1NY[iface]
        statei = q0[i]
        statej = q0[j]
        [leftf,rightf,Lambda] = RIEMANN(lambdai,statei,lambdaj,statej)
        
        flux[i]+=leftf*hx
        flux[j]+=rightf*hx

        k=min(k,volume/codim0to1E[iface]/Lambda)
    
    #print(k,volume,Lambda)
    q1+=flux*k/volume
    print(q0[0::9])
    
    qhist.append(q1.copy())
    q0 = q1.copy()
    t+=k











    

