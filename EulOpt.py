from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize

## Récupération des paramètres du maillage

# Domaine de resolution
Lx = 1 # intervalle [x0,Lx] selon x
Ly = 1 # intervalle [y0,Ly] selon y

LL = min(Lx,Ly)
mylevel = np.load('mylevel.npy') # Parametre entier a modifier pour modifier le pas du maillage
            # Vaut 5 dans le code scilab
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de faces selon la direction x
Ny = NN * int(Ly/LL) # nombre de faces selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

neighbours = 4 # nombre de voisins par cellule

codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')

Ur = np.load('ur_eulerien_constant.npy')
Cr = np.load('SNAPSHOTMATRIX.npy')

r = Ur.shape[1]

sol = []

## Programme d'optimisation

for t in range(50,51):
    def objective(a):
        s = 0
        for k in range(Nx):
            sk = 0
            for l in range(Ny):
                sl = 0
                for j in range(r):
                    sl += a[j]*Ur[k+l*(Nx+1),j]
                sk += (Cr[k+l*(Nx+1),t] - sl)**2
            s += sk
        return s
    
    # initial guesses
    x0 = np.zeros(r)
    
    # show initial objective
    print('Initial Objective: ' + str(objective(x0)))
    
    # optimize
    def constraint(a):
        return Ur @ a
        
    con = {'type': 'ineq', 'fun': constraint}
    
    solution = minimize(objective,x0,method='SLSQP',constraints=con)
    x = solution.x
    
    sol.append(x)
    
    # show final objective
    print('Final Objective: ' + str(objective(x)))

V = [Ur @ q for q in sol]

for q in V:
    plt.clf()
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
    Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
    x=X.reshape(Nx+1,Ny+1)
    y=Y.reshape(Nx+1,Ny+1)
    ax.plot_wireframe(x,y,q.reshape(Nx+1,Ny+1))#,False)
    plt.show()
