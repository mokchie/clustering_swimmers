from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,200))
infile = 'in.run'
Dt = readvar(infile,'Dt')

fig,(ax1,ax2) = plt.subplots(2,1)

F_ph = ['data/phase1-15.data','data/phase2-15.data']
cPh = ['c_phase1min','c_phase2min']
Ph = []
for f_ph,cph in zip(F_ph,cPh):
    d = readave(f_ph)
    t1 = coarse(np.array(d['TimeStep'])*Dt)
    ph1 = coarse(np.array(d[cph]))
    Ph.append(ph1)
    omega1 = (ph1[1:]-ph1[0:-1])/(t1[1:]-t1[0:-1])
    t1a = (t1[1:]+t1[0:-1])/2
    ax2.plot(t1a,omega1,'-')
Ph = np.array(Ph)
ax1.plot(t1,Ph[0]-Ph[1],'-')
ax1.set_xlabel('Time')
ax1.set_ylabel(r'$\Delta \phi$')
ax2.set_xlabel('Time')
ax2.set_ylabel(r'$\omega$')
plt.show()

