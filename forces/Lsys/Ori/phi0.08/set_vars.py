from variable3 import *
from math import pi
from shutil import copyfile
phi = 0.08
N = range(1,10)
FN = [(0,4),(0,4),(0,2),(3,6),(0,0),(0,0),(0,0),(5,7),(8,10)]
TAU = ['0.1+${fn}*0.2','1.0+${fn}*2.0','10.0+5.0*${fn}','10.0+5.0*${fn}','2.55+${fn}','4.0+${fn}','12.0+${fn}','0.1+${fn}*0.2','0.1+${fn}*0.2']
infiles = ['in-%s.run'%i for i in N]
for cn,(i,infile) in enumerate(zip(N,infiles)):
    copyfile('in-0.run',infile)
    print( 'in infile ',infile,' ...' )
    setvar(infile,'Dt',0.002)
    setvar(infile,'sample',20000)
    setvar(infile,'fsample',500)
    setvar(infile,'tot',400000)
    setvar(infile,'start',r'${tot}+1')
    setlog(infile,r'log%s.lammps'%i)
    setloopvar(infile,'fn',FN[cn][0],FN[cn][1])
    setvar(infile,'tau',TAU[cn])
    setvar(infile,'omega',-0.8)
    setvar(infile,'Kb',3000)
    setvar(infile,'Kl',3000)
    setvar(infile,'Km',400000)
    setvar(infile,'phi1',phi)
    setvar(infile,'phi2',-phi)

    print('')


