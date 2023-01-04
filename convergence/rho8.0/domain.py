from create_atom import *
from math import *
from variable3 import *
import os
from scipy.integrate import quad
rho = 8.0
m0 = 6.25/rho
dx = sqrt(1.0/rho)
#rc = sqrt(1.6**2*pi*6.25/(pi*rho))

Ls = 20.0
Lx = Ls
Ns = int(round(Ls/dx))
Ls = Lx = Ns*dx
Ly = 20
Lz = 0.1
init_phase_diff = pi
sdist = Ly/2
y01 = sdist/2
y02 = -sdist/2
xlo = -Lx/2
xhi = Lx/2
ylo = -Ly/2
yhi = Ly/2
zlo = -Lz/2
zhi = Lz/2
k = 4*pi/Ls
head_pos = -Lx/2 #head position of the sperm
sperm_atoms1 = []
sperm_atoms2 = []


n_lay = 1
xcoords = [i*dx+head_pos for i in range(0,Ns)]
for x in xcoords:
    sperm_atoms1.append(make_atom(x,sdist/2,0))
    for i in range(n_lay):
        xi = x + (i%2)*dx/2
        yi0 = -i*dx
        sperm_atoms2.append(make_atom(xi,-sdist/2+(n_lay-1)/2*dx+yi0,0))

sperm_atoms = sperm_atoms1 + sperm_atoms2
solution = box(make_point(xlo,ylo,zlo),make_point(xhi,yhi,zhi))
box_bound = make_bound(-Lx/2, Lx/2, -Ly/2, Ly/2, -Lz/2, Lz/2)
sperm_region = box(make_point(xlo,-sdist/2-(n_lay-1)*dx,zlo),make_point(xhi,-sdist/2,zhi))
Num = int(rho*Lx*Ly) - len(sperm_atoms)
solution_atoms = fill_region_rand(sub_region(solution,sperm_region), box_bound, Num, dim=2, crit=sqrt(1.0/rho)/10,avoid=sperm_atoms)

fp = open('atoms.txt','w')
num = 0
n = 0
ntyp = 0

ntyp+=1
for atom in solution_atoms[0:Num]:
    num += 1
    write_atom(fp, atom, num, ntyp, 0)
    n += 1
print( "%d atoms of type %d" % (n,ntyp) )


print( 'L=',Lx )

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
print( "%d atoms of type %d" % (n,ntyp) )

ntyp+=1
n = 0
sperm_atoms_index2 = []
mol_typ += 1
for atom in sperm_atoms2:
    num += 1
    write_atom(fp, atom, num, ntyp, mol_typ)
    sperm_atoms_index2.append(num)
    n += 1
print( "%d atoms of type %d" % (n,ntyp) )


fp.close()
print( "the combined density is %f" % (num/(Lx*Ly),) )

bond_ntyp = 0
Bonds1 = zip(sperm_atoms_index1,shift(sperm_atoms_index1,1))
bond_ntyp += 1

angle_ntyp = 0

Angles1 = zip(shift(sperm_atoms_index1,-1),shift(sperm_atoms_index1,0),shift(sperm_atoms_index1,1))
Thx1 = [(get_atom_x(atom)-head_pos)*k for atom in sperm_atoms1]
angle_ntyp += 1


fb = open('bonds_and_angles.txt','w')
bond_typ = 1
nbond = 0
fb.write('\nBonds\n\n')
for i,bond in enumerate(Bonds1):
    nbond+=1
    fb.write('%d %d %d %d %f\n'%(nbond,bond_typ,bond[0],bond[1],dx))
angle_typ = 1
nang = 0
fb.write('\nAngles\n\n')
for i,angle in enumerate(Angles1):
    nang += 1
    thx = Thx1[i]
    fb.write('%d %d %d %d %d %f\n'%(nang,angle_typ,angle[0],angle[1],angle[2],thx))
fb.close()

mass = [m0 for i in range(ntyp)]
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
setvar(infile,'pdiff',init_phase_diff)
setvar(infile,'sdist',sdist)
setvar(infile,'y01',y01)
setvar(infile,'y02',y02)
setvar(infile,'m0',m0)
setindexvar(infile,'rc',[1.0,1.2,1.4,1.6,1.8,2.0])
