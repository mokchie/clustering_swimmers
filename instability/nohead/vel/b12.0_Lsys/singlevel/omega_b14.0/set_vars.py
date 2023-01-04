from variable3 import *
from math import pi
from shutil import copyfile
N = range(1,7)
FN = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12)]
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
    setvar(infile,'tau',15.0)
    setvar(infile,'omega',r'-0.1*${fn}')
    setvar(infile,'Kb',3000)
    setvar(infile,'Kl',3000)
    setvar(infile,'b',14.0)

    print('')


