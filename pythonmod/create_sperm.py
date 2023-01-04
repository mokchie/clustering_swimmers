from create_atom import *
from math import *
def make_sperm(px,py,alpha,Rh,Ls,k,dx,n_typ,mol_typ,bond_typ1,bond_typ2,angle_typ,n_shift,bn_shift,an_shift,n_layer=3,ang_eq=None):
    head_pos = 0.0 #head position of the sperm
    mil = int((n_layer-1)/2)
    Ns = int(round(Ls/dx))
    ch = head_pos - 1.2*Rh #center of the head
    Atoms_return = []
    Bonds_return = []
    Angles_return = []

    sperm_atoms1 = [[] for i in range(n_layer)]
    xcoords = [i*dx+head_pos for i in range(0,Ns)]
    for x in xcoords:
        for i in range(n_layer):
            xi = x + ((i+1)%2)*dx/2.0
            yi0 = -(n_layer-1)*dx/2 + i*dx
            sperm_atoms1[i].append(make_atom(xi,yi0,0))

    dxx = dx*sqrt(2)
    ae = ceil(2*Rh/dxx)*dxx
    head_region1 = and_region(sphere(make_point(ch,0,0),Rh*1.05),box(make_point(ch-ae,-ae,-dx/2),make_point(ch+ae,ae,dx/2)))
    head_atoms1 = fill_region(head_region1,make_bound(ch-ae,ch+ae,-ae,ae,0,dx/2),make_part(dx*sqrt(2),dx*sqrt(2),dx*sqrt(2)))

    sperm_atoms = [item for sublist in sperm_atoms1 for item in sublist]+head_atoms1
    sperm_atoms_index1 = [[] for i in range(n_layer)]


    num = n_shift
    sperm_atoms_index1 = [[] for i in range(n_layer)]
    for i in range(n_layer):
        for atom in sperm_atoms1[i]:
            num += 1
            #write_atom(fp, atom, num, ntyp, mol_typ)
            Atoms_return.append((move('y',py,move('x',px,rotate('z',alpha,atom))), num, n_typ, mol_typ))
            sperm_atoms_index1[i].append(num)


    head_atoms_index1 = []
    for atom in head_atoms1:
        num+=1
        #write_atom(fp, atom, num, ntyp, mol_typ)
        Atoms_return.append((move('y',py,move('x',px,rotate('z',alpha,atom))), num, n_typ, mol_typ))
        head_atoms_index1.append(num)

    Bonds1 = []
    for i in range(n_layer):
        Bonds1.append(list(zip(sperm_atoms_index1[i],shift(sperm_atoms_index1[i],1)))[0:-1])
    #    for atom1,atom2 in list(zip(sperm_atoms1[i],shift(sperm_atoms1[i],1)))[0:-1]:
    #        plt.plot((get_atom_x(atom1),get_atom_x(atom2)),(get_atom_y(atom1),get_atom_y(atom2)),'.-')

    for i in range(n_layer-1):
        Bonds1.append(list(zip(sperm_atoms_index1[i],shift(sperm_atoms_index1[i+1],0))))
    #    for atom1,atom2 in list(zip(sperm_atoms1[i],shift(sperm_atoms1[i+1],0))):
    #        plt.plot((get_atom_x(atom1),get_atom_x(atom2)),(get_atom_y(atom1),get_atom_y(atom2)),'.-')

        li = ((i+1)%2)*2-1
        tb = list(zip(sperm_atoms_index1[i],shift(sperm_atoms_index1[i+1],li)))
        if li==-1:
            Bonds1.append(tb[1:])
    #        for atom1,atom2 in list(zip(sperm_atoms1[i],shift(sperm_atoms1[i+1],li)))[1:]:
    #            plt.plot((get_atom_x(atom1),get_atom_x(atom2)),(get_atom_y(atom1),get_atom_y(atom2)),'.-')
        elif li==1:
            Bonds1.append(tb[0:-1])
    #        for atom1,atom2 in list(zip(sperm_atoms1[i],shift(sperm_atoms1[i+1],li)))[0:-1]:
    #            plt.plot((get_atom_x(atom1),get_atom_x(atom2)),(get_atom_y(atom1),get_atom_y(atom2)),'.-')

    head_Bonds1 = []
    head_Bonds_db1 = []
    count = 0
    d_c = sqrt(2)*dx+0.01
    zH = list(zip(head_atoms1,head_atoms_index1))
    for i,(atom1,id1) in enumerate(zH[0:-1]):
        for atom2,id2 in zH[i+1:]:
            dr = atom_dist(atom1,atom2)
            if dr<d_c:
                head_Bonds1.append([id1,id2])
                head_Bonds_db1.append(dr)
                count += 1
                #cl = (0.1,0.5,dr/0.8)
    #            plt.plot((get_atom_x(atom1),get_atom_x(atom2)),(get_atom_y(atom1),get_atom_y(atom2)),'.-')
    dr_ml = []
    q_m = make_atom(ch+Rh,0,0)
    for atom,id in zH:
        dr_m = atom_dist(atom,q_m)
        dr_ml.append(dr_m)
    dr_ml = sorted(enumerate(dr_ml),key=lambda tup: tup[1])
    U_a = []
    U_id = []
    if len(zH)>=5:
        for i in range(5):
            u_a,u_id = zH[dr_ml[i][0]]
            U_a.append(u_a)
            U_id.append(u_id)
        for i in range(n_layer):
            for u_a,u_id in zip(U_a,U_id):
                head_Bonds1.append((sperm_atoms_index1[i][0],u_id))
                head_Bonds_db1.append(atom_dist(sperm_atoms1[i][0],u_a))

                atom1 = sperm_atoms1[i][0]
                atom2 = u_a
    #            plt.plot((get_atom_x(atom1),get_atom_x(atom2)),(get_atom_y(atom1),get_atom_y(atom2)),'.-')
    Bonds1.append(head_Bonds1)

    Angles1 = list(zip(shift(sperm_atoms_index1[mil],-1),shift(sperm_atoms_index1[mil],0),shift(sperm_atoms_index1[mil],1)))[1:-1]
    Thx1 = [(get_atom_x(atom)-head_pos)*k for atom in sperm_atoms1[mil]][1:-1]


    nbond = bn_shift
    for j in range(n_layer+(n_layer-1)*2+1):
        for i,bond in enumerate(Bonds1[j]):
            nbond+=1
            if j==mil:
                bt = bond_typ2
            else:
                bt = bond_typ1
            if j>n_layer-1 and j<n_layer+(n_layer-1)*2:
                db = dx*sqrt(5.0)/2
            elif j==n_layer+(n_layer-1)*2:
                db = head_Bonds_db1[i]
            else:
                db = dx
            #fb.write('%d %d %d %d %f\n'%(nbond,bt,bond[0],bond[1],db))
            Bonds_return.append((nbond,bt,bond[0],bond[1],db))

    nang = an_shift
    for Angles,Thx in zip([Angles1,],[Thx1,]):
        for i,angle in enumerate(Angles):
            nang += 1
            thx = Thx[i]
            #fb.write('%d %d %d %d %d %f\n'%(nang,angle_typ,angle[0],angle[1],angle[2],thx))
            Angles_return.append((nang,angle_typ,angle[0],angle[1],angle[2],thx,ang_eq))

    #plt.show()
    return (Atoms_return,Bonds_return,Angles_return)
