from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import *

from riemann_scal import *

codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')

x0 = -0.5; y0 = -0.5 # origine du repere 2D 
Lx = 2.0; Ly = 2.0 # domaine de resolution[x0,x0+Lx]*[y0,y0+Ly]

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

v0=2#0.75
vx=1/2
vy=1/2
k=0.25
x_0,y_0=0.25,0.25


"""
# Liste des conditions initiales et parametres
X0=[[0.25,0.25],[0.5,0.5],[0.5,0.25]]
V0=[0,0.5]
VX=[1,3]
CI=[]
for a in X0:
    for b in V0:
        for c in VX:
            for d in VX:
                CI.append([a[0],a[1],b,c,d])

"""

# Liste des conditions initiales et parametres
CI=[]
for i in range(64):
    X0=random()*0.5+0.25
    Y0=random()*0.5+0.25
    Vx=0.5
    Vy=0.5
    V0=random()*2.5
    CI.append([X0,Y0,V0,Vx,Vy])
    
np.save('CI',CI)

#print(CI)
CI=[[0,0,1,0.5,0.5]]

sigma=1/20

def initialisation(sigma,x1,y1):
    qini=np.zeros(codim0.shape[0])
    for i in range(len(qini)):
        x,y=codim0[i]
        qini[i]=np.exp(-((x-x1)**2+(y-y1)**2)/(2*sigma**2))
    return qini

def ux(icell,v0,vx,vy):
    x,y=codim0[icell]
    return  2*pi*sin(2*pi*x)*cos(2*pi*y)-2*pi*vy*v0*cos(2*pi*vx*x)*sin(2*pi*vy*y)
    
def uy(icell,v0,vx,vy):
    x,y=codim0[icell]
    return -2*pi*sin(2*pi*y)*cos(2*pi*x)+2*pi*vx*v0*cos(2*pi*vy*y)*sin(2*pi*vx*x)

def uux(x,y,v0,vx,vy):
    u=  2*pi*sin(2*pi*x)*cos(2*pi*y)-2*pi*vy*v0*cos(2*pi*vx*x)*sin(2*pi*vy*y)
    return u

def uuy(x,y,v0,vx,vy):
    u= -2*pi*sin(2*pi*y)*cos(2*pi*x)+2*pi*vx*v0*cos(2*pi*vy*y)*sin(2*pi*vx*x)
    return u



#qhist = [qini.copy()]

#q0 = qini.copy()

#q1 = qini.copy()

"""
fig = plt.figure(5)
ax = fig.gca(projection='3d')
X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
x=X.reshape(Nx+1,Ny+1)
y=Y.reshape(Nx+1,Ny+1)
ax.plot_wireframe(x,y,qhist[0].reshape(Nx+1,Ny+1))#,False)
plt.pause(0.1)
"""

t = 0
Tfinal = 0.5



#print(codim0.shape[0],(Nx+1)*(Ny+1))



X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
x=X.reshape(Nx+1,Ny+1)
y=Y.reshape(Nx+1,Ny+1)
plt.figure(1)
U=np.array([ux(i,v0,vx,vy) for i in range(codim0.shape[0])])
V=np.array([uy(i,v0,vx,vy) for i in range(codim0.shape[0])])
Q = plt.quiver(x,y,U.reshape(Nx+1,Ny+1),V.reshape(Nx+1,Ny+1))#, units='width')
#qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',coordinates='figure')
#plt.title('Ecoulement cellulaire')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.show()
plt.figure(2)
U=np.array([sqrt(ux(i,v0,vx,vy)**2+uy(i,v0,vx,vy)**2) for i in range(codim0.shape[0])])
CS = plt.contour(x, y, U.reshape(Nx+1,Ny+1))
plt.clabel(CS, inline=1, fontsize=10)
plt.plot([0.25,0.25,0.75,0.75,0.25],[0.25,0.75,0.75,0.25,0.25],'r')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.title('Norme de la vitesse')
plt.show()
plt.figure(3)
U=np.array([ux(i,v0,vx,vy) for i in range(codim0.shape[0])])
CS = plt.contour(x, y, U.reshape(Nx+1,Ny+1))
plt.clabel(CS, inline=1, fontsize=10)
plt.plot([0.25,0.25,0.75,0.75,0.25],[0.25,0.75,0.75,0.25,0.25],'r')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.title('ux')
plt.show()
plt.figure(4)
V=np.array([uy(i,v0,vx,vy) for i in range(codim0.shape[0])])
CS = plt.contour(x, y, V.reshape(Nx+1,Ny+1))
plt.clabel(CS, inline=1, fontsize=10)
plt.plot([0.25,0.25,0.75,0.75,0.25],[0.25,0.75,0.75,0.25,0.25],'r')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.title('uy')
plt.show()
plt.pause(0.1)



t = 0
Tfinal = 0.1
compteur=0
qhist=[]
#SNAPSHOTMATRIX=np.zeros((codim0.shape[0],1))

# boucle sur les conditions initiales
for ci in CI:
    compteur+=1
    print(compteur,ci)
    t=0
    [x_0,y_0,v0,vx,vy]=ci
    ## Initialisation
    #qini = np.zeros(codim0.shape[0])
    qini=initialisation(sigma,x_0,y_0)
    qhist.append(qini.copy())
    q0 = qini.copy()
    q1 = qini.copy()
    
    ## Boucle temporelle
    while t<Tfinal:
        k = 1000
        flux=np.zeros(codim0.shape[0]) # Matrice des flux
        # On parcours toutes les faces
        for iface in range(codim1.shape[0]):
            i=codim0to1A[iface] # Cellule en amont
            j=codim0to1B[iface] # Cellule en aval
            if i<0 or j<0: # Conditions de bord periodiques
                continue
            # Vitesse du fluide en i et j projete selon la normale a la face
            lambdai = ux(i,v0,vx,vy)*codim0to1NX[iface]+uy(i,v0,vx,vy)*codim0to1NY[iface] 
            lambdaj = ux(j,v0,vx,vy)*codim0to1NX[iface]+uy(j,v0,vx,vy)*codim0to1NY[iface]
            statei = q0[i] # Concentration en i
            statej = q0[j] # Concentration en j
            # Schema de Lax-Friedriech: calcul du flux en i et en j a travers la face
            [leftf,rightf,Lambda] = RIEMANN(lambdai,statei,lambdaj,statej) 
            # Mise a jour des flux
            flux[i]+=leftf*codim0to1E[iface]
            flux[j]+=rightf*codim0to1E[iface]
            # Calcul du pas de temps associe
            k=min(k,abs(volume/(2*(hx+hy)*Lambda)))
        # Mise a jour de la concentration
        q1+=flux*(k/volume)
        qhist.append(q1.copy())
        q0 = q1.copy()
        t+=k
        
print('fini')


SNAPSHOTMATRIX=np.zeros((codim0.shape[0],len(qhist)))
for ligne in range(codim0.shape[0]):
    for colonne in range(len(qhist)):
        SNAPSHOTMATRIX[ligne,colonne]=qhist[colonne][ligne]
        
print(SNAPSHOTMATRIX.shape[0]) # Compte le nombre de snapshots
np.save('snapshot'+str(CI[0]),SNAPSHOTMATRIX) # Sauvegarde. Attention a changer le nom pour differentes simulations


"""
M=np.load('SNAPSHOT_CELLULARFLOW.npy')

u,s,v=np.linalg.svd(M,False) # Thin SVD

np.save('u_eulerien_cellulaire',u)

np.save('s_eulerien_cellulaire',s)
"""

## Animation 3D
for q in qhist:
    fig = plt.figure(6)
    plt.clf()
    ax = fig.gca(projection='3d')
    X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
    Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
    x=X.reshape(Nx+1,Ny+1)
    y=Y.reshape(Nx+1,Ny+1)
    ax.plot_wireframe(x,y,q.reshape(Nx+1,Ny+1))#,False)
    plt.pause(0.1)
    plt.show()