from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,10))
infile = 'in.run'
Dt = readvar(infile,'Dt')

f_ph = 'phase1-1.data'
d = readave(f_ph)

t1 = coarse(np.array(d['TimeStep'])*Dt)

ph1 = coarse(np.array(d['c_phase1min']))

omega1 = (ph1[1:]-ph1[0:-1])/(t1[1:]-t1[0:-1])
t1a = (t1[1:]+t1[0:-1])/2

fig,(ax1,ax2) = plt.subplots(2,1)
ax1.plot(t1,ph1,'b-')
ax2.plot(t1a,omega1,'b-')

plt.show()

