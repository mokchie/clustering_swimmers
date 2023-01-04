from variable3 import *
from math import pi

N = range(5,6)
FN = [(9,9),]
TAU = ['1.0+${fn}',]
infiles = ['in%s.run'%i for i in N]
for cn,(i,infile) in enumerate(zip(N,infiles)):
    print( 'in infile ',infile,' ...' )
    setvar(infile,'Dt',0.0005)
    setvar(infile,'t0',20)
    setvar(infile,'sample',20000)
    setvar(infile,'tot',400000)
    setvar(infile,'start',r'${tot}+1')
    setlog(infile,r'log%s.lammps'%i)
    setloopvar(infile,'fn',FN[cn][0],FN[cn][1])
    setvar(infile,'tau',TAU[cn])
    setvar(infile,'omega',-1.0)
    print('')


