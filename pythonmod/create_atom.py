from __future__ import division
from math import pi,sin,cos,acos,ceil,floor,sqrt
from random import uniform
EPS = 1e-8
def make_atom(x,y,z):
    return (x,y,z)
def make_atom_with_r(x,y,z,r):
    return (x,y,z,r)
def get_atom_x(atom):
    return atom[0]
def get_atom_y(atom):
    return atom[1]
def get_atom_z(atom):
    return atom[2]
def get_atom_r(atom):
    return atom[3]

def make_point(x,y,z):
    return (x,y,z)
def get_point_x(point):
    return point[0]
def get_point_y(point):
    return point[1]
def get_point_z(point):
    return point[2]

def make_vect(x,y,z):
    return (x,y,z)
def get_vect_x(vect):
    return vect[0]
def get_vect_y(vect):
    return vect[1]
def get_vect_z(vect):
    return vect[2]

def add_vect(vect1, vect2):
    x1 = get_vect_x(vect1)
    x2 = get_vect_x(vect2)
    y1 = get_vect_y(vect1)
    y2 = get_vect_y(vect2)
    z1 = get_vect_z(vect1)
    z2 = get_vect_z(vect2)
    return make_vect(x1+x2, y1+y2, z1+z2)

def sub_vect(vect1, vect2):
    x1 = get_vect_x(vect1)
    x2 = get_vect_x(vect2)
    y1 = get_vect_y(vect1)
    y2 = get_vect_y(vect2)
    z1 = get_vect_z(vect1)
    z2 = get_vect_z(vect2)
    return make_vect(x1-x2, y1-y2, z1-z2)

def mul_vect(vect1, vect2):
    x1 = get_vect_x(vect1)
    x2 = get_vect_x(vect2)
    y1 = get_vect_y(vect1)
    y2 = get_vect_y(vect2)
    z1 = get_vect_z(vect1)
    z2 = get_vect_z(vect2)
    return x1*x2 + y1*y2 + z1*z2

def length_vect(vect):
    x = get_vect_x(vect)
    y = get_vect_y(vect)
    z = get_vect_z(vect)
    return sqrt(x**2 + y**2 + z**2)


def point2vect(point):
    x = get_point_x(point)
    y = get_point_y(point)
    z = get_point_z(point)
    return make_vect(x,y,z)

def rotate(axis, alpha, atom):
    x = get_atom_x(atom)
    y = get_atom_y(atom)
    z = get_atom_z(atom)
    if axis == 'x' or axis == 'X':
        yr = y*cos(alpha) - z*sin(alpha)
        zr = z*sin(alpha) + z*cos(alpha)
        return make_atom(x,yr,zr)
    elif axis == 'y' or axis == 'Y':
        xr = x*cos(alpha) + z*sin(alpha)
        zr = -x*sin(alpha) + z*cos(alpha)
        return make_atom(xr,y,zr)
    elif axis == 'z' or axis == 'Z':
        xr = x*cos(alpha) - y*sin(alpha)
        yr = x*sin(alpha) + y*cos(alpha)
        return make_atom(xr,yr,z)
    else:
        print( "Error: axis should be 'x', 'y' or 'z'\n" )
        exit(0)

def move(axis, dist, atom):
    x = get_atom_x(atom)
    y = get_atom_y(atom)
    z = get_atom_z(atom)
    if axis == 'x' or axis == 'X':
        return make_atom(x+dist,y,z)
    elif axis == 'y' or axis == 'Y':
        return make_atom(x,y+dist,z)
    elif axis == 'z' or axis == 'Z':
        return make_atom(x,y,z+dist)
    else:
        print( "Error: axis should be 'x', 'y' or 'z'\n" )
        exit(0)

def box(corner1, corner2):
    x1 = get_point_x(corner1)
    x2 = get_point_x(corner2)
    y1 = get_point_y(corner1)
    y2 = get_point_y(corner2)
    z1 = get_point_z(corner1)
    z2 = get_point_z(corner2)
    def f(atom):
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)
        return (x-x1)*(x-x2)<=0 and  (y-y1)*(y-y2)<=0 and (z-z1)*(z-z2)<=0
    return f

def cylinder(start_point, directional_vect, radius):
    def f(atom):
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)
        lateral_vect = sub_vect(make_vect(x,y,z), point2vect(start_point))
        lateral_vect_len = length_vect(lateral_vect)
        directional_vect_len = length_vect(directional_vect)
        if lateral_vect_len <= 1e-8:
            return True
        else:
            alpha = acos(mul_vect(directional_vect,lateral_vect)/directional_vect_len/lateral_vect_len)
            return alpha<=pi/2 and cos(alpha)*lateral_vect_len<=directional_vect_len and  sin(alpha)*lateral_vect_len<=radius
    return f

def sphere(center, radius):
    def f(atom):
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)

        dist = length_vect(sub_vect(make_vect(x,y,z), point2vect(center)))
        return dist<=radius
    return f

def frustum(start_point, directional_vect, start_radius, end_radius):
    def f(atom):
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)
        lateral_vect = sub_vect(make_vect(x,y,z), point2vect(start_point))
        lateral_vect_len = length_vect(lateral_vect)
        directional_vect_len = length_vect(directional_vect)
        if lateral_vect_len <= 1e-8:
            return True
        else:
            alpha = acos(mul_vect(directional_vect,lateral_vect)/directional_vect_len/lateral_vect)
            axia_length = cos(alpha)*lateral_vect_len
            radial_length = sin(alpha)*lateral_vect_len
            radial_bound = start_radius + axia_length/directional_vect_len * (end_radius - start_radius)
            return alpha<=pi/2 and axia_length <= directional_vect_length and radial_length <= radial_bound
    return f
def wave(x1,x2,y1,y2,b,k,phi1=0,phi2=None):
    def f(atom):
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)
        if x<x1 or x>x2:
            return False
        else:
            y1c = y1+b*sin(k*x+phi1)
            if phi2 is not None:
                y2c = y2+b*sin(k*x+phi2)
            else:
                y2c = y2+b*sin(k*x+phi1)
            if y<y1c or y>y2c:
                return False
            else:
                return True
    return f
def wave2(x1,x2,y1,y2,b,k,x0=0,phi1=0,phi2=None):
    def f(atom):
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)
        if x<x1 or x>x2:
            return False
        else:
            y1c = y1+b*sin(k*(x-x0)+phi1)
            if phi2 is not None:
                y2c = y2+b*sin(k*(x-x0)+phi2)
            else:
                y2c = y2+b*sin(k*(x-x0)+phi1)
            if y<y1c or y>y2c:
                return False
            else:
                return True
    return f
def empty_region():
    return lambda atom: False
def and_region(region1,region2):
    def f(atom):
        return region1(atom) and region2(atom)
    return f
def or_region(region1,region2):
    def f(atom):
        return region1(atom) or region2(atom)
    return f
def not_region(region):
    def f(atom):
        return not region(atom)
    return f
def sub_region(region1, region2):
    def f(atom):
        return region1(atom) and not region2(atom)
    return f

def make_bound(xl,xh,yl,yh,zl,zh):
    return (xl,xh,yl,yh,zl,zh)
def get_xl(bound):
    return bound[0]
def get_xh(bound):
    return bound[1]
def get_yl(bound):
    return bound[2]
def get_yh(bound):
    return bound[3]
def get_zl(bound):
    return bound[4]
def get_zh(bound):
    return bound[5]

def make_part(dx,dy,dz):
    return (dx,dy,dz)
def get_dx(part):
    return part[0]
def get_dy(part):
    return part[1]
def get_dz(part):
    return part[2]

def divide(bound, part):
    def pair_seq(start, d, n):
        result = []
        for i in range(n):
            result.append((start+i*d,start+i*d+d))
        result.reverse()
        return result
    xl = get_xl(bound)
    yl = get_yl(bound)
    zl = get_zl(bound)
    xh = get_xh(bound)
    yh = get_yh(bound)
    zh = get_zh(bound)
    dx = get_dx(part)
    dy = get_dy(part)
    dz = get_dz(part)

    x_seq = pair_seq(xl, dx, int(ceil((xh-xl)/dx)))
    y_seq = pair_seq(yl, dy, int(ceil((yh-yl)/dy)))
    z_seq = pair_seq(zl, dz, int(ceil((zh-zl)/dz)))
    result = []
    for item in z_seq:
        for jtem in y_seq:
            for ktem in x_seq:
                result.append(make_bound(ktem[0],ktem[1],jtem[0],jtem[1],item[0],item[1]))
    result.reverse()
    return result
    
def make_lattice_atom(bound,lat):
    xl = get_xl(bound)
    xh = get_xh(bound)
    yl = get_yl(bound)
    yh = get_yh(bound)
    zl = get_zl(bound)
    zh = get_zh(bound)
    dx = xh - xl
    dy = yh - yl
    dz = zh - zl
    if lat=='fcc':
        return [make_atom(xl,yl,zl), make_atom(xl,yl+dy/2,zl+dz/2), make_atom(xl+dx/2,yl,zl+dz/2), make_atom(xl+dx/2,yl+dy/2,zl)]
    elif lat=='sc':
        return [make_atom(xl,yl,zl),]

def fill_box(bound, part,lat):
    bound_list = divide(bound,part)
    result = []
    for lattice in bound_list:
        atom_list = make_lattice_atom(lattice,lat)
        for atom in atom_list:
            result.append(atom)
    return result

def fill_region(region, bound, part,lat='fcc'):
    atom_list = fill_box(bound,part,lat)
    result = []
    for atom in atom_list:
        if region(atom):
            result.append(atom)
    return result
def atom_dist(atom1, atom2,Lx=0,Ly=0,Lz=0):
    x1 = get_atom_x(atom1)
    x2 = get_atom_x(atom2)
    y1 = get_atom_y(atom1)
    y2 = get_atom_y(atom2)
    z1 = get_atom_z(atom1)
    z2 = get_atom_z(atom2)
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    if Lx>0:
        dx = dx-round(dx/Lx)*Lx
    if Ly>0:
        dy = dy-round(dy/Ly)*Ly
    if Lz>0:
        dz = dz-round(dz/Lz)*Lz
    return sqrt((dx)**2 + (dy)**2 + (dz)**2)
def atom_distsq(atom1, atom2,Lx=0,Ly=0,Lz=0):
    x1 = get_atom_x(atom1)
    x2 = get_atom_x(atom2)
    y1 = get_atom_y(atom1)
    y2 = get_atom_y(atom2)
    z1 = get_atom_z(atom1)
    z2 = get_atom_z(atom2)
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    if Lx>0:
        dx = dx-round(dx/Lx)*Lx
    if Ly>0:
        dy = dy-round(dy/Ly)*Ly
    if Lz>0:
        dz = dz-round(dz/Lz)*Lz
    return (dx)**2 + (dy)**2 + (dz)**2
def atom_angle(atom1,atom2,atom3,Lx=0,Ly=0,Lz=0):
    x1 = get_atom_x(atom1)
    y1 = get_atom_y(atom1)
    z1 = get_atom_z(atom1)
    x2 = get_atom_x(atom2)
    y2 = get_atom_y(atom2)
    z2 = get_atom_z(atom2)
    x3 = get_atom_x(atom3)
    y3 = get_atom_y(atom3)
    z3 = get_atom_z(atom3)

    dx1 = x1-x2
    dy1 = y1-y2
    dz1 = z1-z2
    dx2 = x3-x2
    dy2 = y3-y2
    dz2 = z3-z2

    if Lx>0:
        dx1 = dx1-round(dx1/Lx)*Lx
        dx2 = dx2-round(dx2/Lx)*Lx
    if Ly>0:
        dy1 = dy1-round(dy1/Ly)*Ly
        dy2 = dy2-round(dy2/Ly)*Ly
    if Lz>0:
        dz1 = dz1-round(dz1/Lz)*Lz
        dz2 = dz2-round(dz2/Lz)*Lz

    r1 = sqrt(dx1**2+dy1**2+dz1**2)
    r2 = sqrt(dx2**2+dy2**2+dz2**2)

    dot = (dx1*dx2+dy1*dy2+dz1*dz2)/r1/r2
    ang = pi-acos(dot)
    crossz = dx1*dy2-dy1*dx2
    if crossz>0:
        ang = -ang
    return ang

def fill_region_rand(region, bound, n, dim=2, crit=0.01,avoid=[],Lx=0,Ly=0,Lz=0):
    xl = get_xl(bound)
    xh = get_xh(bound)
    yl = get_yl(bound)
    yh = get_yh(bound)
    zl = get_zl(bound)
    zh = get_zh(bound)
    result = []
    i = 0
    while i<n:
        x = uniform(xl,xh)
        y = uniform(yl,yh)
        if dim==3:
            z = uniform(zl,zh)
        else:
            z = 0
        atom = make_atom(x,y,z)
        if region(atom):
            accept = True
            for atomt in result+avoid:
                if atom_dist(atomt,atom,Lx=Lx,Ly=Ly,Lz=Lz)<crit:
                    accept = False
                    break
            if accept:
                result.append(atom)
                i+=1
    return result

def fill_region_rand_bin(region, bound, n, dim=2, binnum=1, crit=0.01, crit2=None, avoid=[], func=lambda x: True,Lx=0,Ly=0,Lz=0):
    xl = get_xl(bound)
    xh = get_xh(bound)
    yl = get_yl(bound)
    yh = get_yh(bound)
    zl = get_zl(bound)
    zh = get_zh(bound)
    if crit2 is None:
        crit2 = crit
    Binc = []
    for cl,ch in zip((xl,yl,zl),(xh,yh,zh)):
        cbin = linspace(cl,ch,binnum)
        cbins = [(cb-crit*2,ce+crit*2) for cb,ce in zip(cbin[0:-1],cbin[1:])]
        Binc.append(cbins)
    Bins = []
    for i in Binc[0]:
        for j in Binc[1]:
            for k in Binc[2]:
                Bins.append([(i,j,k),[]])
    result = []
    i = 0
    while i<n:
        x = uniform(xl,xh)
        y = uniform(yl,yh)
        if dim==3:
            z = uniform(zl,zh)
        else:
            z = 0
        atom = make_atom(x,y,z)
        if region(atom) and func(atom):
            accept = True
            avoid_atoms = []
            for (xb,yb,zb),bin_atoms in Bins:
                if xb[0]<=x<=xb[1] and yb[0]<=y<=yb[1] and zb[0]<=z<=zb[1]:
                    avoid_atoms += bin_atoms
            for atomt in avoid_atoms:
                if atom_dist(atomt,atom,Lx=Lx,Ly=Ly,Lz=Lz)<crit:
                    accept = False
                    break
            if accept:
                for atomt in avoid:
                    if atom_dist(atomt,atom,Lx=Lx,Ly=Ly,Lz=Lz)<crit2:
                        accept = False
                        break
            if accept:
                if i%100==0:
                    print('number of filled atoms : ',i,'/',n)
                result.append(atom)
                for (xb,yb,zb),bin_atoms in Bins:
                    if (xb[0]<=x<=xb[1] or xb[0]<=x+Lx<=xb[1] or xb[0]<=x-Lx<=xb[1]) and (yb[0]<=y<=yb[1] or yb[0]<=y+Ly<=yb[1] or yb[0]<=y-Ly<=yb[1]) and (zb[0]<=z<=zb[1] or zb[0]<=z+Lz<=zb[1] or zb[0]<=z-Lz<=zb[1]):
                        bin_atoms.append(atom)
                i+=1
    return result

def fill_region_spheres(region, bound, n, fdist, crit, dim=2, binnum=1, Lx=0,Ly=0,Lz=0):
    xl = get_xl(bound)
    xh = get_xh(bound)
    yl = get_yl(bound)
    yh = get_yh(bound)
    zl = get_zl(bound)
    zh = get_zh(bound)
    Binc = []
    for cl,ch in zip((xl,yl,zl),(xh,yh,zh)):
        cbin = linspace(cl,ch,binnum)
        cbins = [(cb-crit*2,ce+crit*2) for cb,ce in zip(cbin[0:-1],cbin[1:])]
        Binc.append(cbins)
    Bins = []
    for i in Binc[0]:
        for j in Binc[1]:
            for k in Binc[2]:
                Bins.append([(i,j,k),[]])
    result = []
    i = 0
    accept = True
    while i<n:
        x = uniform(xl,xh)
        y = uniform(yl,yh)
        if accept:
            r = fdist()
        if dim==3:
            z = uniform(zl,zh)
        else:
            z = 0
        atom = make_atom_with_r(x,y,z,r)
        if region(atom):
            accept = True
            avoid_atoms = []
            for (xb,yb,zb),bin_atoms in Bins:
                if xb[0]<=x<=xb[1] and yb[0]<=y<=yb[1] and zb[0]<=z<=zb[1]:
                    avoid_atoms += bin_atoms
            for atomt in avoid_atoms:
                if atom_dist(atomt,atom,Lx=Lx,Ly=Ly,Lz=Lz)<=get_atom_r(atomt)+get_atom_r(atom):
                    accept = False
                    break
            if accept:
                if i%100==0:
                    print('number of filled atoms : ',i,'/',n)
                result.append(atom)
                for (xb,yb,zb),bin_atoms in Bins:
                    if (xb[0]<=x<=xb[1] or xb[0]<=x+Lx<=xb[1] or xb[0]<=x-Lx<=xb[1]) and (yb[0]<=y<=yb[1] or yb[0]<=y+Ly<=yb[1] or yb[0]<=y-Ly<=yb[1]) and (zb[0]<=z<=zb[1] or zb[0]<=z+Lz<=zb[1] or zb[0]<=z-Lz<=zb[1]):
                        bin_atoms.append(atom)
                i+=1
    return result


def fill_region_smart_avoid(region, bound, n, dim=2, binnum=1, crit=0.01, crit2=None, avoid=[], func=lambda x: True,Lx=0,Ly=0,Lz=0):
    xl = get_xl(bound)
    xh = get_xh(bound)
    yl = get_yl(bound)
    yh = get_yh(bound)
    zl = get_zl(bound)
    zh = get_zh(bound)
    if crit2 is None:
        crit2 = crit
    Binc = []
    for cl,ch in zip((xl,yl,zl),(xh,yh,zh)):
        cbin = linspace(cl,ch,binnum)
        cbins = [(cb-crit*2,ce+crit*2) for cb,ce in zip(cbin[0:-1],cbin[1:])]
        Binc.append(cbins)
    Bins = []
    for i in Binc[0]:
        for j in Binc[1]:
            for k in Binc[2]:
                Bins.append([(i,j,k),[]])
    result = []

    for atom in avoid:
        x = get_atom_x(atom)
        y = get_atom_y(atom)
        z = get_atom_z(atom)
        for (xb,yb,zb),bin_atoms in Bins:
            if (xb[0]<=x<=xb[1] or xb[0]<=x-Lx<=xb[1] or xb[0]<=x+Lx<=xb[1]) and (yb[0]<=y<=yb[1] or yb[0]<=y-Ly<=yb[1] or yb[0]<=y+Ly<=yb[1]) and (zb[0]<=z<=zb[1] or zb[0]<=z-Lz<=zb[1] or zb[0]<=z+Lz<=zb[1]):
                bin_atoms.append(atom)
                
    i = 0
    while i<n:
        x = uniform(xl,xh)
        y = uniform(yl,yh)
        if dim==3:
            z = uniform(zl,zh)
        else:
            z = 0
        atom = make_atom(x,y,z)
        if region(atom) and func(atom):
            accept = True
            avoid_atoms = []
            for (xb,yb,zb),bin_atoms in Bins:
                if xb[0]<=x<=xb[1] and yb[0]<=y<=yb[1] and zb[0]<=z<=zb[1]:
                    avoid_atoms += bin_atoms
            for atomt in avoid_atoms:
                if atom_dist(atomt,atom,Lx=Lx,Ly=Ly,Lz=Lz)<crit:
                    accept = False
                    break
            # if accept:
            #     for atomt in avoid:
            #         if atom_dist(atomt,atom,Lx=Lx,Ly=Ly,Lz=Lz)<crit2:
            #             accept = False
            #             break
            if accept:
                if i%100==0:
                    print('number of filled atoms : ',i,'/',n)
                result.append(atom)
                for (xb,yb,zb),bin_atoms in Bins:
                    if (xb[0]<=x<=xb[1] or xb[0]<=x+Lx<=xb[1] or xb[0]<=x-Lx<=xb[1]) and (yb[0]<=y<=yb[1] or yb[0]<=y+Ly<=yb[1] or yb[0]<=y-Ly<=yb[1]) and (zb[0]<=z<=zb[1] or zb[0]<=z+Lz<=zb[1] or zb[0]<=z-Lz<=zb[1]):
                        bin_atoms.append(atom)
                i+=1
    return result


def write_atom(fp,atom,num,typ,mol=-1,values=[],image=[]):
    x = get_atom_x(atom)
    y = get_atom_y(atom)
    z = get_atom_z(atom)
    if len(image)>0:
        pos = (x,y,z)+tuple(image)
    else:
        pos = (x,y,z)
    vl = len(values)
    if (mol<0):
        fp.write(("%d %d " +"%f "*vl+ "%s "*len(pos)+"\n") % ((num,typ,)+tuple(values)+pos))
    else:
        fp.write(("%d %d %d " +"%f "*vl+ "%s "*len(pos)+"\n") % ((num,mol,typ)+tuple(values)+pos))
def write_vel(fp,vel):
    n,vx,vy,vz = vel
    fp.write("%d %f %f %f\n" % (n,vx,vy,vz))
def write_bond(fp,nbond,bt,b1,b2,db=None):
    if db:
        fp.write('%d %d %d %d %f\n'%(nbond,bt,b1,b2,db))
    else:
        fp.write('%d %d %d %d\n'%(nbond,bt,b1,b2))
def write_angle(fp,nang,angle_typ,a1,a2,a3,thx=None,ang_eq=None):
    if thx:
        if ang_eq:
            fp.write('%d %d %d %d %d %f %f\n'%(nang,angle_typ,a1,a2,a3,thx,ang_eq))
        else:
            fp.write('%d %d %d %d %d %f\n'%(nang,angle_typ,a1,a2,a3,thx))
    else:
        fp.write('%d %d %d %d %d\n'%(nang,angle_typ,a1,a2,a3))

def shift(seq, n):
    if type(seq) is not list:
        seqc = list(seq)
        print( 'Warning: transforming nonlist to list' )
    else:
        seqc = seq
    n = n % len(seqc)
    return seqc[n:] + seqc[:n]

def shiftmap(x,a0,b0,a,b):
    assert x>=a0 and x<=b0
    assert abs(b0-a0)>EPS
    return (x-a0)/(b0-a0)*(b-a)+a

def get_cm(atom_list):
    n = len(atom_list)
    x = 0
    y = 0
    for atom in atom_list:
        x+=get_atom_x(atom)
        y+=get_atom_y(atom)
    x/=n
    y/=n
    return (x,y)

def linspace(x0,x1,n):
    step = (x1-x0)/n
    return [x0+i*step for i in range(n+1)]

def write_data(filedata,filexyz,Atoms,Bonds,Angles,mass,moment,xlo,xhi,ylo,yhi,zlo,zhi,image_flag=False,Vels = []):
    Lx = xhi-xlo
    Ly = yhi-ylo
    Lz = zhi-zlo
    fdata = open(filedata,'w')
    fxyz = open(filexyz,'w')
    number_of_atoms = len(Atoms)
    number_of_bonds = len(Bonds)
    number_of_angles = len(Angles)
    number_of_atom_types = len(dict.fromkeys([i[2] for i in Atoms]))
    number_of_bond_types = len(dict.fromkeys([i[1] for i in Bonds]))
    number_of_angle_types = len(dict.fromkeys([i[1] for i in Angles]))

    number_of_dihedrals = 0
    number_of_impropers = 0
    number_of_dihedral_types = 0
    number_of_improper_types = 0
    number_of_extra_bond_peratom = 6
    number_of_extra_angle_peratom = 10
    number_of_extra_dihedral_peratom = 6

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

    fxyz.write("%d\n" % number_of_atoms)
    fxyz.write("Atoms. Timestep: 1\n")

    for atom in Atoms:
        if len(atom)>4:
            atom_pos,n,typ,mol,vals = atom
        else:
            atom_pos,n,typ,mol = atom
            vals = []

        x = get_atom_x(atom_pos)
        y = get_atom_y(atom_pos)
        z = get_atom_z(atom_pos)
        fxyz.write("%d %f %f %f\n" % (typ, x, y, z))        
        if image_flag:
            image_flag_x = floor((x-xlo)/Lx)
            image_flag_y = floor((y-ylo)/Ly)
            image_flag_z = floor((z-zlo)/Lz)
            x -= image_flag_x*Lx
            y -= image_flag_y*Ly
            z -= image_flag_z*Lz
            write_atom(fdata,make_atom(x,y,z),n,typ,mol,vals,image=[image_flag_x,image_flag_y,image_flag_z])
        else:
            write_atom(fdata,atom_pos,n,typ,mol,vals)

    fxyz.close()
    if len(Vels)>0:
        fdata.write('\nVelocities\n\n')
        for vel in Vels:
            write_vel(fdata,vel)
    if (len(Bonds)>0):
        fdata.write('\nBonds\n\n')
        for bond in Bonds:
            nbond,bt,b1,b2,db = bond
            write_bond(fdata,nbond,bt,b1,b2,db)
    if (len(Angles)>0):
        fdata.write('\nAngles\n\n')
        for angle in Angles:
            nang,atyp,a1,a2,a3,thx,ang_eq = angle
            write_angle(fdata,nang,atyp,a1,a2,a3,thx,ang_eq)
    fdata.close()
    fxyz.close()
#end

def write_spheres_data(filedata,filexyz,Atoms,xlo,xhi,ylo,yhi,zlo,zhi,image_flag=False):
    Lx = xhi-xlo
    Ly = yhi-ylo
    Lz = zhi-zlo
    fdata = open(filedata,'w')
    fxyz = open(filexyz,'w')
    Bonds = []
    Angles = []
    number_of_atoms = len(Atoms)
    number_of_bonds = len(Bonds)
    number_of_angles = len(Angles)
    number_of_atom_types = len(dict.fromkeys([i[2] for i in Atoms]))
    number_of_bond_types = len(dict.fromkeys([i[1] for i in Bonds]))
    number_of_angle_types = len(dict.fromkeys([i[1] for i in Angles]))

    number_of_dihedrals = 0
    number_of_impropers = 0
    number_of_dihedral_types = 0
    number_of_improper_types = 0
    number_of_extra_bond_peratom = 6
    number_of_extra_angle_peratom = 10
    number_of_extra_dihedral_peratom = 6

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

    # fdata.write('\nMasses\n\n')
    # for j,m in enumerate(mass):
    #     fdata.write('%d %f\n' % (j+1,m))

    # fdata.write('\nMoments\n\n')
    # for j,mom in enumerate(moment):
    #     fdata.write('%d %f\n' % (j+1,mom))

    fdata.write('\nAtoms # sphere\n\n')

    fxyz.write("%d\n" % number_of_atoms)
    fxyz.write("Atoms. Timestep: 1\n")

    for atom in Atoms:
        atom_pos,n,typ,density = atom

        x = get_atom_x(atom_pos)
        y = get_atom_y(atom_pos)
        z = get_atom_z(atom_pos)
        r = get_atom_r(atom_pos)
        fxyz.write("%d %f %f %f %f\n" % (typ, x, y, z, r))        
        if image_flag:
            image_flag_x = floor((x-xlo)/Lx)
            image_flag_y = floor((y-ylo)/Ly)
            image_flag_z = floor((z-zlo)/Lz)
            x -= image_flag_x*Lx
            y -= image_flag_y*Ly
            z -= image_flag_z*Lz
            fdata.write('%d %d %f %f %f %f %f\n'%(n,typ,2*r,density,x,y,z))
        else:
            fdata.write('%d %d %f %f %f %f %f\n'%(n,typ,2*r,density,x,y,z))            
    fxyz.close()
    fdata.close()
    fxyz.close()

