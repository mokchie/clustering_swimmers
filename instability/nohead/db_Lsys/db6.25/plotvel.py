from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
ts = 50000
Nsample = int(readvar('in.run','statsample'))
N = 50
fig,ax1 = plt.subplots(1,1)
nx = int(readvar('in.run','nsx'))
ny = int(readvar('in.run','nsy'))
U = 0
V = 0

U = 0
V = 0
filelist = ['data/vel-15.%d.plt'%i for i in range(ts,int(ts+Nsample*N),Nsample)]
for cn,velfile in enumerate(filelist):
    d = readstat(velfile)
    X = d['x'][0:nx]
    Y = d['y'][::ny]
    Vx = np.reshape(d['v_x'],(nx,ny),'C')
    Vy = np.reshape(d['v_y'],(nx,ny),'C')
    U+=Vx
    V+=Vy
U/=len(filelist)
V/=len(filelist)
us = 0
us = np.average(U[0])
ax1.quiver(X,Y,U-us,V)


plt.show()
