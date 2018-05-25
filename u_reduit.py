import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

u=np.load('u_eulerien_constant_combine.npy')

ur=u[:,50]

np.save('u_eulerien_constant_combine_reduit',ur)