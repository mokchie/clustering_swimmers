from create_atom import *
from math import *
from variable import *
import os,sys
from scipy.integrate import quad

rho = 6.25
Lx = 40.0
Ly = 40.0
Lz = 0.1
rc = 1.6
if len(sys.argv)>1:
    m0 = float(sys.argv[1])
else:
    m0 = 0.125
print 'm0 =',m0
b = 1.0
xlo = -Lx/2
xhi = Lx/2
ylo = -Ly/2
yhi = Ly/2
zlo = -Lz/2
zhi = Lz/2
Ls = Lx/2 #length of the sperm
k = 2*pi/Ls
head_pos = -Ls/2 #head position of the sperm
Ns = int(round(Ls/(sqrt(1.0/rho))))
if Ns%2==1:
    Ns-=1
dx = Ls/Ns
n_lay = 3
sdist = 0
def Fw(x):
    return sqrt(1.0+b**2*k**2*cos(x)*cos(x))/2.0/pi
Q,err = quad(Fw,0,2.0*pi)

sperm_atoms1 = []
xcoords = [i*dx+head_pos for i in range(0,Ns)]
for x in xcoords:
    for i in range(n_lay):
        xi = x + (i%2)*dx/2
        yi0 = -(n_lay-1)*dx/2 + i*dx
        sperm_atoms1.append(make_atom(xi,yi0+sdist/2+b*sin(k*(xi-head_pos)),0))
sperm_atoms = sperm_atoms1
sheets_region = wave2(head_pos,head_pos+Ls,-(n_lay)*dx/2+sdist/2,(n_lay)*dx/2+sdist/2,b,k,head_pos)
box_bound = make_bound(-Lx/2, Lx/2, -Ly/2, Ly/2, -Lz/2, Lz/2)
Num = int(rho*Lx*Ly)-len(sperm_atoms)

out_n = Num
solution = sub_region(box(make_point(xlo,ylo,zlo),make_point(xhi,yhi,zhi)),sheets_region)
outer_region = solution

solution_atoms = fill_region_rand(outer_region, box_bound, Num, dim=2, crit=sqrt(1.0/rho)/10,avoid=sperm_atoms)
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

fp.close()
print "the combined density is %f" % (num/(Lx*Ly),)

mass = [m0 for i in range(ntyp)]
moment = [m0**(5.0/3.0) for i in range(ntyp)]    

number_of_atoms = num
number_of_bonds = 0
number_of_angles = 0
number_of_dihedrals = 0
number_of_impropers = 0
number_of_atom_types = ntyp
number_of_bond_types = 0
number_of_angle_types = 0
number_of_dihedral_types = 0
number_of_improper_types = 0
number_of_extra_bond_peratom = 6
number_of_extra_angle_peratom = 10
number_of_extra_dihedral_peratom = 6


fp = open('atoms.txt','r')
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

fp.close()
fxyz.close()
fdata.close()
os.remove('atoms.txt')
infile = 'in.run'
setvar(infile,'rho0',num/Lx/Ly)
setvar(infile,'k',k)
setvar(infile,'Lx',Lx)
setvar(infile,'Ly',Ly)
setvar(infile,'Ls',Ls)
setvar(infile,'dx',dx)
setvar(infile,'Q',Q)
setvar(infile,'hp',head_pos)
setvar(infile,'sdist',sdist)
setvar(infile,'b',b)
setvar(infile,'rc',rc)
