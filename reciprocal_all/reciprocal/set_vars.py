from variable3 import *
from math import pi
import shutil
tau0 = 0.1
def ftau(i):
    if i==0:
        return tau0
    else:
        return round(ftau(i-1)+i*0.1,1)
nc = 6
nr = 5
N = range(nr)
infiles = ['in%s.run'%i for i in N]
for cn,(i,infile) in enumerate(zip(N,infiles)):
    tau = [ftau(j) for j in range(i*nc,(i+1)*nc)]
    shutil.copyfile('in.run',infile)
    print( 'in infile ',infile,' ...' )
    setvar(infile,'Dt',0.0005)
    setvar(infile,'t0',20)
    setvar(infile,'sample',20000)
    setvar(infile,'tot',400000)
    setlog(infile,r'log%s.lammps'%i)
    setindexvar(infile,'tau',tau)
    setvar(infile,'omega',-0.6)
    print('')


