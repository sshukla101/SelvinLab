import numpy as np
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D

ampa= np.loadtxt('AMPAR_cal.txt')
ampa=ampa[:,5:8]

nmda=np.loadtxt('NMDAR_cal.txt')
nmda=ampa[:,5:8]



fig = pylab.figure()
ax = Axes3D(fig)

ax.scatter(ampa[:,0], ampa[:,1], ampa[:,2])
ax.scatter(nmda[:,0], nmda[:,1], nmda[:,2])
plt.show()
