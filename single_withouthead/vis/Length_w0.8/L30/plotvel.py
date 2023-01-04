from __future__ import division
import numpy as np
from readvel import *
from scipy.optimize import curve_fit
from variable import *
import matplotlib.pyplot as plt
from copy import deepcopy

Ind = ['.']
Vy = []
Vx = []
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
def fsin(x,b,k,theta0):
    return b*np.sin(k*x+theta0)
fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2)
for ind in Ind:
    j = 0
    direct = str(ind)+'/'
    infile = direct+'in.run'
    tot = readvar(infile,'tot')
    fn = readvar(infile,'fn')
    pdiff = readvar(infile,'pdiff')
    dt = readvar(infile,'Dt')
    kt = readvar(infile,'k')
    omega = [readvar(infile,'omega'),]
    Ls = readvar(infile,'Ls')
    Ly = readvar(infile,'Ly')
    Lx = readvar(infile,'Lx')
    dx = readvar(infile,'dx')
    sdist = readvar(infile,'sdist')
    velfile = direct+'yvel-'+str(fn[j])+'.profile'
    trjfile = 'lammpstrj-%d.data'%fn[j]
    d1,d2 = readchunk(velfile,('Coord1','vx','vy'),dict_concatenate)
    data = readxyz(trjfile,[2,3])
    vTime = np.array(d1['Timestep'])*dt
    vx_vlist = d2['vx']
    vy_vlist = d2['vy']
    y_vlist = d2['Coord1']

    Time = []
    Dtheta1 = []
    Dtheta2 = []
    Dtheta = []
    Vsheet1 = []
    Vsheet2 = []
    Vsheet1off = []
    Vsheet2off = []
    count = 0
     
    for jj,d in enumerate(data):
        try:
            timestep = d['timestep']
            X = d['x']
            Y = d['y']
            Typ = d['type']
            ID = d['id']
            if jj>0:
                X1o = deepcopy(X1s)
                X2o = deepcopy(X2s)
            X1,Y1,Typ1,ID1 = np.transpose(select(zip(X,Y,Typ,ID), lambda tup: tup[2]==2))
            X2,Y2,Typ2,ID2 = np.transpose(select(zip(X,Y,Typ,ID), lambda tup: tup[2]==3))

            X1s,Y1s,ID1s = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2]))
            X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[0]))
            headx1,heady1,headid1 = sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2])[0]

            X2s,Y2s,ID2s = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2]))
            X2,Y2,ID2 = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[0]))
            headx2,heady2,headid2 = sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2])[0]

            ave1 = np.average(Y1)
            popt,pcov = curve_fit(fsin,X1,Y1-ave1,p0=[1.0,kt,0.0])
            b1,k1,theta01 = popt[0:3]
            if b1<0:
                b1=-b1
                theta01 -= np.pi
            ave2 = np.average(Y2)
            popt,pcov = curve_fit(fsin,X2,Y2-ave2,p0=[1.0,kt,0.0])
            b2,k2,theta02 = popt[0:3]
            if b2<0:
                b2=-b2
                theta02 -= np.pi
            dtheta = theta01-theta02
            dtheta -= round(dtheta/2/np.pi)*(2*np.pi)
            dtheta = dtheta/2/np.pi*360
#            if 180-dtheta<90:
#                dtheta-=360
            Dtheta.append(dtheta)
            count += 1

            ax1.clear()
            ax2.clear()
            ax4.clear()
            ax5.clear()
            ax6.clear()

            ax1.plot(X1,Y1,'b.')
            ax1.plot(X1,fsin(X1,b1,k1,theta01)+ave1,'b-')

            ax1.plot(X2,Y2,'r.')
            ax1.plot(X2,fsin(X2,b2,k2,theta02)+ave2,'r-')
            Yt1 = []
            Yt2 = []
            for i,(x1,x2,id1,id2) in enumerate(zip(X1,X2,ID1,ID2)):
               Yt1.append(-b1*np.sin(kt*(id1-np.min(ID1))*dx+omega[j]*timestep*dt)+np.average(Y1))
               Yt2.append(-b2*np.sin(kt*(id2-np.min(ID2))*dx+omega[j]*timestep*dt+pdiff)+np.average(Y2))

            popt,pcov = curve_fit(fsin,X1,Yt1-ave1,p0=[1.0,kt,0.0])
            bt1,kt1,thetat01 = popt[0:3]
            if bt1<0:
                bt1=-bt1
                thetat01 -= np.pi
            dtheta1 = thetat01-theta01
            dtheta1 -= round(dtheta1/2/np.pi)*(2*np.pi)
            dtheta1 = -dtheta1/np.pi*180
            Dtheta1.append(dtheta1)

            popt,pcov = curve_fit(fsin,X2,Yt2-ave2,p0=[1.0,kt,0.0])
            bt2,kt2,thetat02 = popt[0:3]
            if bt2<0:
                bt2=-bt2
                thetat02 -= np.pi
            dtheta2 = thetat02-theta02
            dtheta2 -= round(dtheta2/2/np.pi)*(2*np.pi)
            dtheta2 = -dtheta2/np.pi*180
            Dtheta2.append(dtheta2)
            time = timestep*dt
            Time.append(time)       

            if jj>0:
                dsheet = X1s-X1o
                dsheet -= np.round(dsheet/Lx)*Lx
                Vsheet1.append(np.average(dsheet/(Time[-1]-Time[-2])))

                dsheet = X2s-X2o
                dsheet -= np.round(dsheet/Lx)*Lx
                Vsheet2.append(np.average(dsheet/(Time[-1]-Time[-2])))
            
            ax1.plot(X1,Yt1,'b--')
            ax1.plot(X2,Yt2,'r--')
            ax1.plot([headx1,],[heady1,],'go')
            ax1.plot([headx2,],[heady2,],'go')
            if len(Dtheta)>1:
                ax2.plot(Time[1:],Dtheta[1:],'b.-')
            ax4.plot(Time,Dtheta1,'b.-')
            ax4.plot(Time,Dtheta2,'r.-')
            if jj>0:
                ax5.plot(Time[1:],Vsheet1,'b.-')
                ax5.plot(Time[1:],Vsheet2,'r.-')

            ax1.set_ylim((-sdist/2-2,sdist/2+3))
            ax1.set_xlim((-Ls/2-2,Ls/2+2))

            #ax2.set_ylim((-180,0))
            ax1.text(-Ls/2+2,sdist/2+1,'T = %.0f\n(totT = %.0f)'%(timestep*dt,tot*dt))
            ax1.set_xlabel('x')
            ax1.set_ylabel('y')
            ax5.set_xlabel('Time')
            ax5.set_ylabel('vx')

            ax4.set_xlabel('time')
            ax4.set_ylabel('phase lag')
            ax2.set_ylabel('phase diff')


            if jj>0:
                try:
                    Vsheet1off.append(Vsheet1[-1]-np.average(list(vx_v[0:5])+list(vx_v[-5:])))
                    Vsheet2off.append(Vsheet2[-1]-np.average(list(vx_v[0:5])+list(vx_v[-5:])))

                except NameError:
                    Vsheet1off.append(Vsheet1[-1])
                    Vsheet2off.append(Vsheet2[-1])
                ax6.plot(Time[1:],Vsheet1off,'b.-')
                ax6.plot(Time[1:],Vsheet2off,'r.-')
                ax6.set_xlabel('Time')
                ax6.set_ylabel('vx')
            if np.round(time) in np.round(vTime):
                ax3.clear()
                n = list(np.round(vTime)).index((np.round(time)))
                y_v = y_vlist[n]
                vx_v = vx_vlist[n]
                vy_v = vy_vlist[n]

                ax3.plot(y_v,vx_v,'.-')
                ax3.plot([-Ls,Ls],[0,0],'g--')
                ax3.text(-Ls/2+2,0.018,'T = %.0f\n(totT = %.0f)'%(time,tot*dt))
                ax3.set_xlim((-Ly/2-2,Ly/2+2))
                ax3.set_ylim((-0.01,0.025))
                ax3.set_xlabel('y')
                ax3.set_ylabel('vx')

            fig.canvas.draw()
            plt.pause(.001)
        except RuntimeError:
            pass



