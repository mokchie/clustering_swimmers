from variable3 import *
from math import pi
import os

N = range(1,5)
FN = [(0,4),(5,9),(1,4),(5,9)]
TAU = ['0.1+${fn}*0.2','0.1+${fn}*0.2','1.0+${fn}*2.0','1.0+${fn}*2.0',]
infiles = ['in-%s.run'%i for i in N]
b0 = 12.0
IND = [11.125,]
Dirs = ['./']
for ind,direct in zip(IND,Dirs):
    for cn,(i,infile) in enumerate(zip(N,infiles)):
        infile = direct+'/'+infile
        print( 'in infile ',infile,' ...' )
        setvar(infile,'Dt',0.002)
        setvar(infile,'sample',100000)
        setvar(infile,'fsample',20000)
        setvar(infile,'tot',1000000)
        setvar(infile,'start',r'${tot}+1')
        setlog(infile,r'log%s.lammps'%i)
        setloopvar(infile,'fn',FN[cn][0],FN[cn][1])
        setvar(infile,'tau',TAU[cn])
        setvar(infile,'omega',-0.6)
        setvar(infile,'b1',b0-ind/2)
        setvar(infile,'b2',b0+ind/2)
        print('')


