from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,5))
infile = 'in1.run'
Dt = readvar(infile,'Dt')
fn = readvar(infile,'fn')
Tau = readvar(infile,'tau')
fig,ax1 = plt.subplots(1,1)

for cf in Tau:

    f_fvs1 = 'fvs1-%d.data'%cf
    d1st,d1s = readrow(f_fvs1)

    t1 = coarse(np.array(d1st['TimeStep'])*Dt)

    fvs1 = np.array(d1s['c_cfvs1'])

    p1 = -coarse(fvs1[:,0]+fvs1[:,1])


    ax1.plot(t1[2:],p1[2:],'-')

plt.show()

