from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,5))
infile = 'in.run'
Dt = readvar(infile,'Dt')
tau = readvar(infile,'tau')
fig,ax1 = plt.subplots(1,1)

for cf in tau:
    for f_fvs1,c_cfv1 in zip(['data/fvs1-%s.data'%int(cf), 'data/fvs2-%s.data'%int(cf)],['c_cfvs1','c_cfvs2']):

        d1st,d1s = readrow(f_fvs1)

        t1 = coarse(np.array(d1st['TimeStep'])*Dt)

        fvs1 = np.array(d1s[c_cfv1])

        p1 = -coarse(fvs1[:,0]+fvs1[:,1])

        ax1.plot(t1[2:],p1[2:],'-')
ax1.set_xlabel('time')
ax1.set_ylabel('P')
plt.show()

