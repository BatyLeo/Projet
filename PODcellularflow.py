from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


r = 10 #nombre de modes pour la POD

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

codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')

v0=0#0.75
vx=0.5
vy=0.5
k=0.25
x_0,y_0=0.5,0.25


def ux(icell):
    x,y=codim0[icell]
    return  2*pi*sin(2*pi*x)*cos(2*pi*y)-2*pi*vy*v0*cos(2*pi*vx*x)*sin(2*pi*vy*y)
    
def uy(icell):
    x,y=codim0[icell]
    return -2*pi*sin(2*pi*y)*cos(2*pi*x)+2*pi*vx*v0*cos(2*pi*vy*y)*sin(2*pi*vx*x)


u=np.load('u_eulerien_cellulaire.npy')

ur=u[:,:r]
np.save('u_eulerien_cellulaire_reduit',ur)


ur=np.load('u_eulerien_cellulaire_reduit.npy')
s=np.load('s_eulerien_cellulaire.npy')





"""
# Plot des valeurs singulieres
plt.plot([i for i in range(r)],[np.log(s[i]/s[0]) for i in range(r)])
plt.show()

# Modes propres
fig=plt.figure(2)
ax=fig.gca(projection='3d')
X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
x=X.reshape(Nx+1,Ny+1)
y=Y.reshape(Nx+1,Ny+1)
ax.plot_wireframe(x,y,u[:,0].reshape(129,129))
plt.show()


"""
# matrice identite
Ir=np.zeros((r,r))
for i in range(r):
    Ir[i,i]=1

I=codim0.shape[0]

B=np.zeros((I,I))
for i in range(I):
    B[i,(i+Nx+1)%I]=1/hx/2*ux(i)
    B[i,(i+1)%I]=1/hx/2*uy(i)
    B[i,(i-1)%I]=-1/hx/2*uy(i)
    B[i,(i-Nx-1)%I]=-1/hx/2*ux(i)


Bt = ur.T @ B @ ur


P=1000
Tfinal = 0.1
dt=Tfinal/P

## Initialisation
qini = np.zeros(codim0.shape[0])

for icell in range(codim0.shape[0]):
    xi,yi=codim0[icell]
    norme=sqrt((xi-x_0)**2+(yi-y_0)**2)
    if norme<k:
        qini[icell]=exp(-1/(k**2-norme**2))
q0 = qini.copy()
q1 = qini.copy()


# Initialisation
A0 = ur.T @ q0
A = A0 

Atot = np.zeros((r, P+1));
Atot[:,0]=A0

for p in range(1,P+1):
    A=np.linalg.solve(Ir-dt*Bt,A0)
    Atot[:,p]=A
    A0=A
    
Vr = ur @ Atot # matrice des snapshots

print(Vr.shape)

## Animation 3D
for s in range(1):
    q=Vr[:,s]
    fig = plt.figure(5)
    plt.clf()
    ax = fig.gca(projection='3d')
    X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
    Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
    x=X.reshape(Nx+1,Ny+1)
    y=Y.reshape(Nx+1,Ny+1)
    ax.plot_wireframe(x,y,q.reshape(Nx+1,Ny+1))#,False)
    plt.pause(0.01)
