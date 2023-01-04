from variable3 import *
from math import pi

N = range(1,6)
FN = [(0,4),(5,9),(1,5),(6,10),(0,0)]
TAU = ['1.0+${fn}','1.0+${fn}','10.0+5.0*${fn}','10.0+5.0*${fn}','12.0+${fn}']
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
    print('')


