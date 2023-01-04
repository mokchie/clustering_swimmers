from create_atom import *
from math import *
from variable import *
import os
from scipy.integrate import quad
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

rho = 6.25
Lx = 20.0
Ly = 10.0
Lz = 0.1
init_phase_diff = pi
sdist = 1.0
y01 = sdist/2
y02 = -sdist/2
Num = int(rho*Lx*Ly)
xlo = -Lx/2
xhi = Lx/2
ylo = -Ly/2
yhi = Ly/2
zlo = -Lz/2
zhi = Lz/2
Ls = Lx #length of the sperm
k = 4*pi/Ls
amp = 0.2
head_pos = -Lx/2 #head position of the sperm

sperm_atoms1 = []
sperm_atoms2 = []
dx = sqrt(1.0/rho)/1.1
xcoords = [i*dx+head_pos for i in range(0,int(round(Ls/dx)))]
for x in xcoords:
    sperm_atoms1.append(make_atom(x,sdist/2-amp*sin(k*(x-head_pos)),0))
    sperm_atoms2.append(make_atom(x,-sdist/2-amp*sin(k*(x-head_pos)+init_phase_diff),0))
sperm_atoms = sperm_atoms1 + sperm_atoms2
solution = box(make_point(xlo,ylo,zlo),make_point(xhi,yhi,zhi))
box_bound = make_bound(-Lx/2, Lx/2, -Ly/2, Ly/2, -Lz/2, Lz/2)

solution_atoms = fill_region_rand(solution, box_bound, Num, dim=2, crit=sqrt(1.0/rho)/10,avoid=sperm_atoms)



fp = open('atoms.txt','w')
num = 0
n = 0
ntyp = 0

ntyp+=1
for atom in solution_atoms[0:Num]:
    num += 1
    write_atom(fp, atom, num, ntyp, 0)
    n += 1
print "%d atoms of type %d" % (n,ntyp)


print 'L=',Lx
print 'rho=',Num/(Ly*Lx)

mol_typ = 0

ntyp+=1
n = 0
sperm_atoms_index1 = []
mol_typ += 1
for atom in sperm_atoms1:
    num += 1
    write_atom(fp, atom, num, ntyp, mol_typ)
    sperm_atoms_index1.append(num)
    n += 1
print "%d atoms of type %d" % (n,ntyp)

ntyp+=1
n = 0
sperm_atoms_index2 = []
mol_typ += 1
for atom in sperm_atoms2:
    num += 1
    write_atom(fp, atom, num, ntyp, mol_typ)
    sperm_atoms_index2.append(num)
    n += 1
print "%d atoms of type %d" % (n,ntyp)


fp.close()
print "the combined density is %f" % (num/(Lx*Ly),)

bond_ntyp = 0
Bonds1 = zip(sperm_atoms_index1,shift(sperm_atoms_index1,1))
bond_ntyp += 1
Dx1 = [atom_dist(atom[0],atom[1],Lx,Ly) for atom in zip(sperm_atoms1,shift(sperm_atoms1,1))]

Bonds2 = zip(sperm_atoms_index2,shift(sperm_atoms_index2,1))
bond_ntyp += 1
Dx2 = [atom_dist(atom[0],atom[1],Lx,Ly) for atom in zip(sperm_atoms2,shift(sperm_atoms2,1))]
angle_ntyp = 0

Angles1 = zip(sperm_atoms_index1,shift(sperm_atoms_index1,1),shift(sperm_atoms_index1,2))
true_angles = [atom_angle(atom[0],atom[1],atom[2],Lx,Ly) for atom in zip(shift(sperm_atoms1,-1),shift(sperm_atoms1,0),shift(sperm_atoms1,1))]

plt.plot(xcoords,true_angles,'o')


popt,pcov = curve_fit(lambda x,bf,xf0: bf*np.sin(k*x+xf0), np.array(xcoords), np.array(true_angles))
b = abs(popt[0])
Thx1 = [(get_atom_x(atom)-head_pos)*k for atom in sperm_atoms1]
angle_ntyp += 1
plt.plot(xcoords,b*np.sin(np.array(Thx1)),'.-')
plt.show()
Angles2 = zip(sperm_atoms_index2,shift(sperm_atoms_index2,1),shift(sperm_atoms_index2,2))
Thx2 = [(get_atom_x(atom)-head_pos)*k+init_phase_diff for atom in sperm_atoms2]
angle_ntyp += 1


fb = open('bonds_and_angles.txt','w')
bond_typ = 1
nbond = 0
fb.write('\nBonds\n\n')
for i,bond in enumerate(Bonds1):
    nbond+=1
    fb.write('%d %d %d %d %f\n'%(nbond,bond_typ,bond[0],bond[1],Dx1[i]))
bond_typ += 1
for i,bond in enumerate(Bonds2):
    nbond+=1
    fb.write('%d %d %d %d %f\n'%(nbond,bond_typ,bond[0],bond[1],Dx2[i]))
angle_typ = 1
nang = 0
fb.write('\nAngles\n\n')
for i,angle in enumerate(Angles1):
    nang += 1
    thx = Thx1[i]
    fb.write('%d %d %d %d %d %f\n'%(nang,angle_typ,angle[0],angle[1],angle[2],thx))
angle_typ += 1
for i,angle in enumerate(Angles2):
    nang += 1
    thx = Thx2[i]
    fb.write('%d %d %d %d %d %f\n'%(nang,angle_typ,angle[0],angle[1],angle[2],thx))

fb.close()

mass = [1.0 for i in range(ntyp)]
moment = [1.0 for i in range(ntyp)]

number_of_atoms = num
number_of_bonds = nbond
number_of_angles = nang
number_of_dihedrals = 0
number_of_impropers = 0
number_of_atom_types = ntyp
number_of_bond_types = bond_ntyp
number_of_angle_types = angle_ntyp
number_of_dihedral_types = 0
number_of_improper_types = 0
number_of_extra_bond_peratom = 6
number_of_extra_angle_peratom = 10
number_of_extra_dihedral_peratom = 6


fp = open('atoms.txt','r')
fb = open('bonds_and_angles.txt','r')
fxyz = open('atoms.xyz','w')
fdata = open('atoms.data','w')
fdata.write('LAMMPS calculation\n\n')
fdata.write('%d atoms\n' % number_of_atoms)
fdata.write('%d bonds\n' % number_of_bonds)
fdata.write('%d angles\n' % number_of_angles)
fdata.write('%d dihedrals\n\n' % number_of_dihedrals)
#fdata.write('%d impropers\n' % number_of_impropers)
fdata.write('%d atom types\n' % number_of_atom_types)
fdata.write('%d bond types\n' % number_of_bond_types)
fdata.write('%d angle types\n' % number_of_angle_types)
fdata.write('%d dihedral types\n' % number_of_dihedral_types)
#fdata.write('%d improper types\n' % number_of_improper_types)
fdata.write('%d extra bond per atom\n' % number_of_extra_bond_peratom)
fdata.write('%d extra angle per atom\n' % number_of_extra_angle_peratom)
fdata.write('%d extra dihedral per atom\n' % number_of_extra_dihedral_peratom)

fdata.write('\n%f %f  xlo xhi\n' % (xlo,xhi))
fdata.write('%f %f  ylo yhi\n' % (ylo,yhi))
fdata.write('%f %f  zlo zhi\n' % (zlo,zhi))

fdata.write('\nMasses\n\n')
for j,m in enumerate(mass):
    fdata.write('%d %f\n' % (j+1,m))

fdata.write('\nMoments\n\n')
for j,mom in enumerate(moment):
    fdata.write('%d %f\n' % (j+1,mom))

fdata.write('\nAtoms\n\n')

fxyz.write("%d\n" % num)
fxyz.write("Atoms. Timestep: 1\n")

for line in fp:
    line_list = line.strip().split()
    n = int(line_list[0])
    mol = int(line_list[1])
    typ = int(line_list[2])
    x = float(line_list[3])
    y = float(line_list[4])
    z = float(line_list[5])
    fxyz.write("%d %f %f %f\n" % (typ, x, y, z))
    fdata.write("%d %d %d %f %f %f\n" % (n, mol, typ, x, y, z))

for line in fb:
    fdata.write(line)

fp.close()
fb.close()
fxyz.close()
fdata.close()
os.remove('atoms.txt')
os.remove('bonds_and_angles.txt')
infile = 'in.run'
setvar(infile,'rho0',num/Lx/Ly)
setvar(infile,'k',k)
setvar(infile,'Lx',Lx)
setvar(infile,'Ly',Ly)
setvar(infile,'Ls',Ls)
setvar(infile,'dx',dx)
setvar(infile,'y01',y01)
setvar(infile,'y02',y02)
setvar(infile,'b',b/pi*180)
print 'epsilon1 = bk =',amp*k
print 'epsilon2 = h/lambda =',sdist/(2*pi/k)
