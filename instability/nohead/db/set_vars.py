from variable3 import *
from math import pi
import os

N = range(1,5)
FN = [(0,9),(10,20),(2,8),(2,6)]
TAU = ['0.1+${fn}*0.2','0.1+${fn}*0.2','1.0+${fn}*2.0','10.0+5.0*${fn}']
infiles = ['in-%s.run'%i for i in N]
b0 = 12.0
IND = [8.0,8.25,8.5,8.75,9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.125,11.25,11.5,11.75,12.0]+[10.875,11.0625]
Dirs = ['db'+str(i) for i in IND]
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


