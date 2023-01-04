from __future__ import division
from variable3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,root
from copy import deepcopy
from pdb import *
infile = 'in.run'
Ly = readvar(infile,'Ly')-2*readvar(infile,'dx')
N = 200
mu = readvar(infile,'mu')
rho = readvar(infile,'rho0')*readvar(infile,'m0')
print(rho)
d0 = readvar(infile,'rho0')
nu = mu/rho
Np = readvar(infile,'Np')
kBT = readvar(infile,'T')
tau = readvar(infile,'tau')
Tp = readvar(infile,'Tp')
tot = readvar(infile,'tot')*readvar(infile,'Dt')
posy = readvar(infile,'sdist')
X = np.linspace(-Ly/2,Ly/2,N+1)
X = (X[1:]+X[0:-1])/2
dx = X[1]-X[0]

Ux = np.zeros(N)
Cxx = np.ones(N)
Cxy = np.zeros(N)
Y = np.zeros(N)
y0 = []
y1 = []

#fig,(ax1,ax2) = plt.subplots(2,1)
tt = 0    
def boundary(Ux,t):
    Ux[0]  = 0.5*np.sin(2*np.pi/Tp*t)
    Ux[-1] = 0.5*np.sin(2*np.pi/Tp*t)

def f(Z,Z0,v,rho,Np,tau,kBT,dx,dt,u0,u1):
    Vx = np.append(np.append([u0,],Z[0:N-2]),[u1,])
    Cxx = Z[N-2:N-2+N]
    Cxy = Z[N-2+N:]
    Vx0 = Z0[0:N]
    Cxx0 = Z0[N:2*N]
    Cxy0 = Z0[2*N:]

    d1 = v*((Vx[2:]+Vx[0:-2]-2*Vx[1:-1])+(Vx0[2:]+Vx0[0:-2]-2*Vx0[1:-1]))/dx**2/2.0*dt + dt/2.0*Np*d0/rho*kBT*((Cxy[2:]-Cxy[0:-2])/2.0/dx + (Cxy0[2:]-Cxy0[0:-2])/2.0/dx) - (Vx[1:-1] - Vx0[1:-1])
    duxdy = np.append(np.append([0,],(Vx[2:]-Vx[0:-2])/2/dx),[0,])
    duxdy0 = np.append(np.append([0,],(Vx0[2:]-Vx0[0:-2])/2/dx),[0,])

    d2 = (2*Cxy*duxdy + (1-Cxx)/tau + 2*Cxy0*duxdy0 + (1-Cxx0)/tau)/2.0*dt - (Cxx - Cxx0)
    
    d3 = (duxdy - Cxy/tau + duxdy0 - Cxy0/tau)/2.0*dt - (Cxy - Cxy0)

    return np.append(np.append(d1,d2),d3)

Time = np.linspace(0,tot,500)
dt = Time[1]-Time[0]

for j,t in enumerate(Time):
    tt+=1
    Cxx0 = deepcopy(Cxx)
    Cxy0 = deepcopy(Cxy)
    Ux0 = deepcopy(Ux)
    print( t )
    boundary(Ux,t)
    sol = root(lambda Z: f(Z,np.append(np.append(Ux0,Cxx),Cxy),nu,rho,Np,tau,kBT,dx,dt,Ux[0],Ux[-1]),x0=np.append(np.append(Ux[1:-1],Cxx),Cxy))

    if not sol.success:
        print( sol.message )
    else:
        print( 'solution found' )
    Ux[1:-1] = sol.x[0:N-2]
    Cxx = sol.x[N-2:N-2+N]
    Cxy = sol.x[N-2+N:]

    #ax1.clear()
    #ax1.plot(X,Ux)
    #ax1.text(-0.1,0.25,'t = %.3f'%t)
    #ax1.set_ylim(-0.8,0.8)

    #ax2.clear()
    Y += (Ux+Ux0)/2*dt
    y0.append(Y[0])
    y1.append(Y[int(np.round(posy/dx))])
    #ax2.plot(Time[0:j+1],np.array(y1)-np.array(y0),'b.-')

    #fig.canvas.draw()
    #plt.pause(0.001)

with open('fd.data','w') as fp:
    for t,dy in zip(Time,np.array(y1)-np.array(y0)):
        fp.write('%f %f\n'%(t,dy))
        
