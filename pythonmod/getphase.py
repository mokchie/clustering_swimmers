from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

def getphase(infile,trjfile,typ):
    def fcos(x,b,k,theta0):
        return b*np.cos(k*x+theta0)
    kt = readvar(infile,'k')
    dt = readvar(infile,'Dt')
    Time = []
    B1 = []
    data = readxyz(trjfile,[typ,])
    Ph1 = []
    for d in data:
        timestep = d['timestep']
        X = d['x']
        Y = d['y']
        Typ = d['type']
        ID = d['id']
        X1,Y1,ID1,Typ1 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==typ))
        X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2]))
        l1 = int(len(ID1)/3)
        X1 = X1[l1:2*l1]
        Y1 = Y1[l1:2*l1]
        ID1 = ID1[l1:2*l1]
        X1,Y1 = np.transpose(sorted(zip(X1,Y1)))

        popt,pcov = curve_fit(fcos,X1,Y1-np.average(Y1),p0=[1.0,kt,0])
        b1,k1,theta01 = popt[0:3]
        if b1<0:
            b1=-b1
            theta01 -= np.pi

        B1.append(b1)
        Ph1.append(theta01)
        Time.append(timestep*dt)
    return (Time,Ph1,B1)
