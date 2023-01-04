from variable3 import *
from math import pi
from shutil import copyfile

N = range(1,5)
FN = [(0,4),(5,9),(1,5),(6,10)]
TAU = ['1.0+${fn}','1.0+${fn}','10.0+5.0*${fn}','10.0+5.0*${fn}']
infiles = ['in%s.run'%i for i in N]
for cn,(i,infile) in enumerate(zip(N,infiles)):
    copyfile('in.run',infile)
    print( 'in infile ',infile,' ...' )
    setvar(infile,'Dt',0.002)
    setvar(infile,'t0',20)
    setvar(infile,'sample',100000)
    setvar(infile,'fsample',20000)
    setvar(infile,'tot',1000000)
    setvar(infile,'start',r'${tot}+1')
    setlog(infile,r'log%s.lammps'%i)
    setloopvar(infile,'fn',FN[cn][0],FN[cn][1])
    setvar(infile,'tau',TAU[cn])
    setvar(infile,'omega',-0.6)
    setvar(infile,'Kl',3000)
    setvar(infile,'Kth',20000)
    print('')


