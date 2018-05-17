import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from riemann_scal import *

codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')


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
    return 1
    
def uy(i):
    return 0

## Initialisation
qini = np.zeros(codim0.shape[0])

qini[10]=1
qini[0:3]=3*[0.5]
qini[9]=0.5
qini[11]=0.5
qini[18:21]=[0.5,0.5,0.5]

qhist = [qini.copy()]

q0 = qini.copy()

q1 = qini.copy()

t = 0
Tfinal = 10

## Boucle temporelle
while t<Tfinal:
    k = 1000
    flux=np.zeros(codim0.shape[0]) # Matrice des flux
    # On parcours toutes les faces
    for iface in range(codim1.shape[0]):
        i=codim0to1A[iface] # Cellule en amont
        j=codim0to1B[iface] # Cellule en aval
        if i<0 or j<0: # Conditions de bord périodiques
            continue
        # Vitesse du fluide en i et j projeté selon la normale à la face
        lambdai = ux(i)*codim0to1NX[iface]+uy(i)*codim0to1NY[iface] 
        lambdaj = ux(j)*codim0to1NX[iface]+uy(j)*codim0to1NY[iface]
        statei = q0[i] # Concentration en i
        statej = q0[j] # Concentration en j
        # Schéma de Lax-Friedriech: calcul du flux en i et en j à travers la face
        [leftf,rightf,Lambda] = RIEMANN(lambdai,statei,lambdaj,statej) 
        # Mise à jour des flux
        flux[i]+=leftf*hx
        flux[j]+=rightf*hx
        # Calcul du pas de temps associé
        k=min(k,volume/codim0to1E[iface]/Lambda)
    # Mise à jour de la concentration
    q1+=flux*k/volume
    qhist.append(q1.copy())
    q0 = q1.copy()
    t+=k

## Plot 2D
#for q in qhist:
#    plt.matshow(q.reshape((Nx+1,Ny+1)),False,origin='lower')
#    plt.pause(0.001)



## Plot 3D
for q in qhist:
    plt.clf()
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
    Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
    x=X.reshape(Nx+1,Ny+1)
    y=Y.reshape(Nx+1,Ny+1)
    ax.plot_wireframe(x,y,q.reshape(Nx+1,Ny+1))#,False)
    plt.pause(0.1)


