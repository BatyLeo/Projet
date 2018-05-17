# -*- coding: ISO-8859-1 -*-
#---------------- Cr�ation du maillage carr� en 2 dimensions ------------------

import math
import numpy as np
import matplotlib.pyplot as plt







x0 = 0.0; y0 = 0.0 # origine du rep�re 2D 
Lx = 1.0; Ly = 1.0 # domaine de r�solution[x0,x0+Lx]*[y0,y0+Ly]

LL = min(Lx,Ly)
mylevel = 5 # Param�tre entier � modifier pour modifier le pas du maillage
            # Vaut 5 dans le code scilab
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de faces selon la direction x
Ny = NN * int(Ly/LL) # nombre de faces selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

neighbours = 4 # nombre de voisins par cellule







# Initialisation
Ncell = 0   # nombre de cellules
codim0 = np.array([[0,0]]) # chaque ligne stocke les coordonn�es (x,y) du centre de chaque cellule
Nface = 0   # nombre de fronti�res
codim1 = np.array([[0,0,0,0]]) # chaque ligne stocke les coordonn�es (x,y)
        # de deux points d�limitant une face
        # (de la droite vers la gauche pour les faces horizontales
        # du bas vers le haut pour les faces verticales)

codim0to1E = []  # liste des longueurs des fronti�res (correspond avec codim1)
codim0to1NX = [] # composantes x des normales des faces (correspond avec codim1)
codim0to1NY = [] # composantes y des normales des faces (correspond avec codim1)
Nghost = 0 # nombre de cellules fantome?? (au bord)
codim0to1A = [] # liste des num�ros des cellules en amont des fronti�res correspondant � codim1
codim0to1B = [] # liste des num�ros des cellules en aval des fronti�res correspondant � codim1


# Cr�ation du maillage
# On parcours le maillage ligne par ligne, de la droite vers la gauche, et de haut en bas
ny = 0
while ny<=Ny:
    nx = 0
    while nx<=Nx:
        #(numcell = nx + ny*Nx # num�ro de la cellule (on commence la num�rotation � 0))
        
        # -------------- On ajoute une cellule ------------------------------
        xx = x0 + nx*hx # abscisse du centre de la cellule
        yy = y0 + ny*hy # ordonn�e du centre de la cellule
        
        codim0=np.concatenate((codim0,np.array([[xx,yy]])))
        
        # --- On cr�e la face verticale � gauche de la cellule ---
        # Cas du bord x=0
        if nx==0:
            Nface+=1
            codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy-hy/2,xx-hx/2,yy+hy/2]])))
            
            ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), ou 0
            ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), ou 0
            assert(ex>=0 and ey>=0)
            E = math.sqrt(ex*ex + ey*ey) # norme de la face
            NX = ey/E; # composante x de la normale � la face (ne pas confondre avec Nx!!)
            NY = ex/E; # composante y de la normale � la face (ne pas confondre avec Ny!!)
            codim0to1E.append(E)
            codim0to1NX.append(NX)
            codim0to1NY.append(NY)
            
            # On ajoute les deux c�t�s de la face
            # Nghost+=1 # condition p�riodique
            codim0to1A.append(Ncell+Nx) # car il n'y a pas de cellule en amont
            codim0to1B.append(Ncell) # la cellule que l'on vient de cr�er est en aval
            
        
        # --- On cr�e la face verticale � droite de la cellule ---
        Nface+=1
        codim1=np.concatenate((codim1,np.array([[xx+hx/2,yy-hy/2,xx+hx/2,yy+hy/2]])))
        
        ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), or 0
        ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), or 0
        assert(ex>=0 and ey>=0)
        E = math.sqrt(ex*ex + ey*ey) # norme de la face
        NX = ey/E; # composante x de la normale � la face (ne pas confondre avec Nx!!)
        NY = ex/E; # composante y de la normale � la face (ne pas confondre avec Ny!!)
        codim0to1E.append(E)
        codim0to1NX.append(NX)
        codim0to1NY.append(NY)
        
        # Cas du bord x=Nx 
        if nx==Nx:
            codim0to1A.append(Ncell) # la cellule que l'on vient de cr�er est en amont
            #codim0to1B.append(-(Nghost+1)) # condition p�riodique
            #Nghost+=1
            codim0to1B.append(Ncell-Nx) # car il n'y a pas de cellule en aval

        else:
            codim0to1A.append(Ncell) # la cellule que l'on vient de cr�er est en amont
            codim0to1B.append(Ncell+1)

        
        # --- On cr�e la face horizontale en-dessous de la cellule ---
        # Cas o� y=0
        if ny==0:
            Nface+=1
            codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy-hy/2,xx+hx/2,yy-hy/2]])))
            
            ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), ou 0
            ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), ou 0
            assert(ex>=0 and ey>=0)
            E = math.sqrt(ex*ex + ey*ey) # norme de la face
            NX = ey/E; # composante x de la normale � la face (ne pas confondre avec Nx!!)
            NY = ex/E; # composante y de la normale � la face (ne pas confondre avec Ny!!)
            codim0to1E.append(E)
            codim0to1NX.append(NX)
            codim0to1NY.append(NY)
            #codim0to1A.append(-(Nghost+1)) # condition p�riodique           
            #Nghost+=1
            codim0to1A.append(Ncell+Ny*(Nx+1)) # il n'y a rien en amont
            codim0to1B.append(Ncell) # la cellule que l'on vient de cr�er est en aval
        
        # --- On cr�e la face horizontale au-dessus de la cellule ---
        Nface+=1
        codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy+hy/2,xx+hx/2,yy+hy/2]])))
        
        ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), ou 0
        ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), ou 0
        assert(ex>=0 and ey>=0)
        E = math.sqrt(ex*ex + ey*ey) # norme de la face
        NX = ey/E; # composante x de la normale � la face (ne pas confondre avec Nx!!)
        NY = ex/E; # composante y de la normale � la face (ne pas confondre avec Ny!!)
        codim0to1E.append(E)
        codim0to1NX.append(NX)
        codim0to1NY.append(NY)
        
        # Cas o� y=Ny-1
        if ny==Ny:           
            # Nghost+=1 # condition p�riodique
            codim0to1A.append(Ncell) # la cellule que l'on vient de cr�er est en amont
            codim0to1B.append(Ncell-Ny*(Nx+1)) # il n'y a rien en aval
        else:
            codim0to1A.append(Ncell)
            codim0to1B.append(Ncell+Nx+1)

        # On passe � la cellule suivante
        Ncell+=1 
        nx+=1
    ny+=1

codim0=codim0[1:,:]
codim1=codim1[1:,:]

'''
# Visualisation et v�rification du maillage:
    # En noir le maillage, en bleu le passage des face
for i in range(len(codim1[:,1])):
    a=[codim1[i,0],codim1[i,2]]
    b=[codim1[i,1],codim1[i,3]]
    plt.plot(a,b,"black")
    if codim0to1A[i]>=0 and codim0to1B[i]>=0:
        X=[codim0[codim0to1A[i],0],codim0[codim0to1B[i],0]]
        Y=[codim0[codim0to1A[i],1],codim0[codim0to1B[i],1]]
        plt.plot(X,Y,'b')
plt.axis('equal')
plt.show()
'''
# Sauvegarde du maillage
np.save('codim0', codim0)
np.save('codim1', codim1)
np.save('codim0to1A', codim0to1A)
np.save('codim0to1B', codim0to1B)
np.save('codim0to1E', codim0to1E)
np.save('codim0to1NX', codim0to1NX)
np.save('codim0to1NY', codim0to1NY)
