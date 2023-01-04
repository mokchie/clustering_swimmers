from __future__ import division
import copy
import re
import pdb
from math import sqrt
def dict_replace(d,dc):
    for k,v in dc.items():
        d[k] = v
    return d
def dict_concatenate(d,dc):
    for k,v in dc.items():
        if k not in d:
            d[k] = []
        d[k].append(v)
def passline(fp,n):
    for i in range(n):
        fp.readline()
def readindex(fp):
    line = fp.readline()
    if line == '':
        return None
    line = line.strip(r" #\t")
    index_list = line.split()
    d = {}
    keys = []
    for index in index_list:
        d[index] = []
        keys.append(index)
    return (d,keys)
def readnumline(fp):
    line = fp.readline()
    if line == '':
        return None
    line = line.strip()
    num_list = line.split()
    result = []
    for num in num_list:
        if '.' in num or 'e' in num or 'nan' in num or 'inf' in num:
            result.append(float(num))
        else:
            result.append(int(num))
    return result
# def readchunk(filename,mlist=[],fd=dict_concatenate):
#     with open(filename,'r') as fp:
#         passline(fp,1)
#         d1,keys1 = readindex(fp)
#         d2,keys2 = readindex(fp)
#         if 'Number-of-chunks' not in d1:
#             print( "Error: unable to determine the number of chunks" )
#             exit(0)
#         while True:
#             line = readnumline(fp)
#             if line == None:
#                 break
#             for i,item in enumerate(line):
#                 d1[keys1[i]].append(item)
#             nor = d1["Number-of-chunks"][-1]
#             dc = {}
#             if len(mlist)!=0:
#                 for m in mlist:
#                     dc[m] = []
#             else:
#                 for m in keys2:
#                     dc[m] = []
#             for i in range(nor):
#                 line = readnumline(fp)
#                 if line == 'None':
#                     print( "Error: the chunk is not complete" )
#                     return None
#                 for k,item in enumerate(line):
#                     if keys2[k] in mlist or len(mlist)==0:
#                         dc[keys2[k]].append(item)
#             fd(d2,dc)
#     return d1,d2

def readrow(filename,mlist=[],fd=dict_concatenate):
    with open(filename,'r') as fp:
        knoc = ['Number-of-rows',"Number-of-chunks"]
        passline(fp,1)
        d1,keys1 = readindex(fp)
        #print(keys1)
        d2,keys2 = readindex(fp)
        #print(keys2)
        rkey = None
        for key in knoc:
            if key in d1:
                rkey = key
                break
        if rkey==None:
            print( "Error: unable to determine the number of rows" )
            exit(0)
        while True:
            line = readnumline(fp)
            if line == None:
                break
            for i,item in enumerate(line):
                d1[keys1[i]].append(item)
            nor = d1[rkey][-1]
            dc = {}
            if len(mlist)!=0:
                for m in mlist:
                    dc[m] = []
            else:
                for m in keys2:
                    dc[m] = []
            for i in range(nor):
                line = readnumline(fp)
                if line == 'None':
                    print( "Error: the row is not complete" )
                    return None
                for k,item in enumerate(line):
                    if keys2[k] in mlist or len(mlist)==0:
                        dc[keys2[k]].append(item)
            fd(d2,dc)
    return d1,d2
readchunk = readrow
def readave(filename,mlist=[]):
    with open(filename,'r') as fp:
        passline(fp,1)
        d,keys = readindex(fp)
        while True:
            line = readnumline(fp)
            if line==None or len(line)==0:
                break
            for i,item in enumerate(line):
                if keys[i] in mlist or len(mlist)==0:
                    d[keys[i]].append(item)
    return d
def readsimprows(filename,mlist=[]):
    with open(filename,'r') as fp:
        d,keys = readindex(fp)
        while True:
            line = readnumline(fp)
            if line==None or len(line)==0:
                break
            for i,item in enumerate(line):
                if keys[i] in mlist or len(mlist)==0:
                    d[keys[i]].append(item)
    return d
def readrdf(filename):
    with open(filename, 'r') as fp:
        passline(fp,3)
        while True:
            line = readnumline(fp)
            if line==None:
                break
            else:
                r = []
                g = []
                timestep,Nc = line
                for i in range(Nc):
                    line = readnumline(fp)
                    r.append(line[1])
                    g.append(line[2:len(line)])
    return (r,g)

def readstat(filename,mlist=[]):
    pat = re.compile(r'"([a-zA-Z0-9_ \t]+)"')
    with open(filename,'r') as fp:
        line = fp.readline()
        mat = pat.findall(line)
        d = dict()
        for i in mat:
            d[i] = []
        passline(fp,1)
        while True:
            line = readnumline(fp)
            if line:
                for i,n in enumerate(line):
                    if mat[i] in mlist or len(mlist)==0:
                        d[mat[i]].append(n)
            else:
                break
    return d
def readaverows(filename,mlist=[]):
    with open(filename,'r') as fp:
        passline(fp,1)
        d1,keys1 = readindex(fp)
        #print(keys1)
        l1 = []
        d2,keys2 = readindex(fp)
        #print(keys2)
        if 'Number-of-rows' not in d1:
            print( "Error: unable to determine the number of rows" )
            exit(0)
        while True:
            line = readnumline(fp)
            #print(line)
            if line == None:
                break
            for i,item in enumerate(line):
                d1[keys1[i]].append(item)
            nor = d1["Number-of-rows"][-1]
            dc = {}
            for m in keys2:
                dc[m] = []
            for i in range(nor):
                line = readnumline(fp)
                if line == 'None':
                    print( "Error: the row is not complete" )
                    return None
                for k,item in enumerate(line):
                    if keys2[k] in mlist or len(mlist)==0:
                        dc[keys2[k]].append(item)
            l1.append(dc)
    return d1,l1
def coarseave(lst, n):
    res = []
    r = 0
    nc = 0
    for i,v in enumerate(lst):
        r += v
        nc += 1
        if (i+1)%n== 0 or i+1==len(lst):
            res.append(r/nc)
            r = 0
            nc = 0
    return res

def readxyz(filename,typelist=[]):
    lst = []
    with open(filename, 'r') as fp:
        while True:
            line = fp.readline()
            if line == '':
                break
            if 'ITEM: TIMESTEP' in line:
                d = {}
                d['timestep'] = readnumline(fp)[0]
                lst.append(d)
            if 'ITEM: NUMBER OF ATOMS' in line:
                lst[-1]['N'] = readnumline(fp)[0]
            if 'ITEM: BOX BOUNDS' in line:
                xb = readnumline(fp)
                yb = readnumline(fp)
                zb = readnumline(fp)
                lst[-1]['bound'] = (xb,yb,zb)
            if 'ITEM: ATOMS' in line:
                ID = []
                Mol = []
                Type = []
                X = []
                Y = []
                Z = []
                n = 0
                mol = 0
                for i in range(lst[-1]['N']):
                    rnline = readnumline(fp)
                    if len(rnline)==5:
                        (aid,atype,x,y,z) = rnline
                    elif len(rnline)==6:
                        (aid,mol,atype,x,y,z) = rnline
                    else:
                        raise ValueError
                    if atype in typelist or len(typelist) == 0:
                        n+=1
                        ID.append(aid)
                        Type.append(atype)
                        X.append(x)
                        Y.append(y)
                        Z.append(z)
                        Mol.append(mol)
                lst[-1]['N']=n
                lst[-1]['id']=ID
                lst[-1]['type']=Type
                lst[-1]['x']=X
                lst[-1]['y']=Y
                lst[-1]['z']=Z
                lst[-1]['mol']=Mol
    return lst

def readxyz_custom(filename,typelist=[]):
    lst = []
    with open(filename, 'r') as fp:
        while True:
            line = fp.readline()
            if line == '':
                break
            if 'ITEM: TIMESTEP' in line:
                d = {}
                d['timestep'] = readnumline(fp)[0]
                lst.append(d)
            if 'ITEM: NUMBER OF ATOMS' in line:
                d['N'] = readnumline(fp)[0]
            if 'ITEM: BOX BOUNDS' in line:
                xb = readnumline(fp)
                yb = readnumline(fp)
                zb = readnumline(fp)
                lst[-1]['bound'] = (xb,yb,zb)

            if 'ITEM: ATOMS' in line:
                line_items = line.strip().split()[2:] 
                for item in line_items:
                    d[item] = []
                n = 0
                for i in range(d['N']):
                    rnline = readnumline(fp)
                    atype = rnline[line_items.index('type')]
                    if atype in typelist or len(typelist) == 0:
                        for j,item in enumerate(line_items):
                            d[item].append(rnline[j])
                        n += 1
                d['N'] = n
    return lst
def readdump(filename):
    lst = []
    with open(filename, 'r') as fp:
        while True:
            line = fp.readline()
            if line == '':
                break
            if 'ITEM: TIMESTEP' in line:
                d = {}
                d['timestep'] = readnumline(fp)[0]
                lst.append(d)
            if 'ITEM: NUMBER OF ATOMS' in line or 'ITEM: NUMBER OF ENTRIES' in line:
                d['N'] = readnumline(fp)[0]
            if 'ITEM: BOX BOUNDS' in line:
                xb = readnumline(fp)
                yb = readnumline(fp)
                zb = readnumline(fp)
                d['bound'] = (xb,yb,zb)

            if 'ITEM: ATOMS' in line or 'ITEM: ENTRIES' in line:
                line_items = line.strip().split()[2:] 
                for item in line_items:
                    d[item] = []
                n = 0
                for i in range(d['N']):
                    rnline = readnumline(fp)
                    for j,item in enumerate(line_items):
                        d[item].append(rnline[j])
                    n += 1
                d['N'] = n
    return lst

def select(lst,fun):
    res = []
    for i in lst:
        if fun(i):
            res.append(i)
    return res

def myround(a,n=6):
    a = round(a,n)
    if a==int(a):
        a = int(a)
    return a
def scalingunits(rc,d0):
    rho_ref = d0
    eta_ref = d0
    t_ref = rc**2*rho_ref/eta_ref
    F_ref = rho_ref*rc**3/t_ref**2
    return (rho_ref,eta_ref,t_ref,F_ref)
def scalingunits_off(rc,d0):
    rho_ref = 1
    eta_ref = 1
    t_ref = 1
    F_ref = 1
    return (rho_ref,eta_ref,t_ref,F_ref)
def savedata(filename,X,Y,header=''):
    with open(filename,'w') as fp:
        fp.write('#'+header+'\n')
        for x,y in zip(X,Y):
            fp.write('%s  %s\n'%(x,y))
def readdata(filename):
    Res = []
    with open(filename,'r') as fp:
        passline(fp,1)
        for line in fp:
            vals = [float(item) for item in line.strip().split()]
            Res.append(vals)
    return list(zip(*Res))
def fixed_or_not(filename1,filename2,eps=1e-6,typ=[]):
    dlist = readxyz_custom(filename1)
    ID0 = dlist[0]['id']
    X0 = dlist[0]['x']
    Y0 = dlist[0]['y']
    Z0 = dlist[0]['z']
    Ty0 = dlist[0]['type']
    Fix = {}
    dxyz0 = {}
    for id,ty,x,y,z in zip(ID0,Ty0,X0,Y0,Z0):
        dxyz0[id] = (x,y,z)
        if len(typ)==0 or ty in typ:
            Fix[id] = 1
        else:
            Fix[id] = 0
    with open(filename2,'w') as fxyz:
        for d in dlist:
            fxyz.write("ITEM: TIMESTEP\n%d\n"%d['timestep'])
            fxyz.write("ITEM: NUMBER OF ATOMS\n%d\n" % d['N'])
            fxyz.write("ITEM: BOX BOUNDS pp pp pp\n")
            fxyz.write("%f %f\n"%tuple(d['bound'][0]))
            fxyz.write("%f %f\n"%tuple(d['bound'][1]))
            fxyz.write("%f %f\n"%tuple(d['bound'][2]))
            vs = ['id','type','x','y','z']
            for v in d.keys():
                if v not in ['timestep','bound','N','id','type','x','y','z','fixed']:
                    vs.append(v)
            vs+=['fixed',]

            fxyz.write("ITEM: ATOMS "+' '.join(vs)+'\n')
            X = d['x']
            Y = d['y']
            Z = d['z']
            Ty = d['type']
            ID = d['id']
            for ii,(id,ty,x,y,z) in enumerate(zip(ID,Ty,X,Y,Z)):
                x0,y0,z0 = dxyz0[id]
                linelist = []
                #pdb.set_trace()
                for v in vs[0:-1]:
                    linelist.append(d[v][ii])
                formatstring = "%s "*len(vs)+"\n"
                if Fix[id] and ((x-x0)**2 + (y-y0)**2 + (z-z0)**2 < eps**2):
                    fxyz.write(formatstring  % tuple(linelist+['1',]))
                else:
                    fxyz.write(formatstring  % tuple(linelist+['0',]))
                    Fix[id] = 0


def fixed_or_not_2d(filename1,filename2,eps=1e-6,typ=[]):
    dlist = readxyz_custom(filename1)
    ID0 = dlist[0]['id']
    X0 = dlist[0]['x']
    Y0 = dlist[0]['y']
    Ty0 = dlist[0]['type']
    Fix = {}
    dxyz0 = {}
    for id,ty,x,y in zip(ID0,Ty0,X0,Y0):
        dxyz0[id] = (x,y)
        if len(typ)==0 or ty in typ:
            Fix[id] = 1
        else:
            Fix[id] = 0
    with open(filename2,'w') as fxyz:
        for d in dlist:
            fxyz.write("ITEM: TIMESTEP\n%d\n"%d['timestep'])
            fxyz.write("ITEM: NUMBER OF ATOMS\n%d\n" % d['N'])
            fxyz.write("ITEM: BOX BOUNDS pp pp pp\n")
            fxyz.write("%f %f\n"%tuple(d['bound'][0]))
            fxyz.write("%f %f\n"%tuple(d['bound'][1]))
            fxyz.write("%f %f\n"%tuple(d['bound'][2]))
            vs = ['id','type','x','y']
            for v in d.keys():
                if v not in ['timestep','bound','N','id','type','x','y','fixed']:
                    vs.append(v)
            vs+=['fixed',]

            fxyz.write("ITEM: ATOMS "+' '.join(vs)+'\n')
            X = d['x']
            Y = d['y']
            Ty = d['type']
            ID = d['id']
            for ii,(id,ty,x,y) in enumerate(zip(ID,Ty,X,Y)):
                x0,y0 = dxyz0[id]
                linelist = []
                #pdb.set_trace()
                for v in vs[0:-1]:
                    linelist.append(d[v][ii])
                formatstring = "%s "*len(vs)+"\n"
                if Fix[id] and ((x-x0)**2 + (y-y0)**2 < eps**2):
                    fxyz.write(formatstring  % tuple(linelist+['1',]))
                else:
                    fxyz.write(formatstring  % tuple(linelist+['0',]))
                    Fix[id] = 0




def append_v(filename1,filename2,eps=1e-6,typ=[]):
    dlist = readxyz_custom(filename1)
    with open(filename2,'w') as fxyz:
        for d in dlist:
            fxyz.write("ITEM: TIMESTEP\n%d\n"%d['timestep'])
            fxyz.write("ITEM: NUMBER OF ATOMS\n%d\n" % d['N'])
            fxyz.write("ITEM: BOX BOUNDS pp pp pp\n")
            fxyz.write("%f %f\n"%tuple(d['bound'][0]))
            fxyz.write("%f %f\n"%tuple(d['bound'][1]))
            fxyz.write("%f %f\n"%tuple(d['bound'][2]))
            vs = ['id','type','x','y','vx','vy']
            for v in d.keys():
                if v not in ['timestep','bound','N','id','type','x','y','vx','vy']:
                    vs.append(v)
            vs+=['v',]

            fxyz.write("ITEM: ATOMS "+' '.join(vs)+'\n')
            X = d['x']
            Y = d['y']
            Vx = d['vx']
            Vy = d['vy']
            Ty = d['type']
            ID = d['id']
            for ii,(id,ty,x,y,vx,vy) in enumerate(zip(ID,Ty,X,Y,Vx,Vy)):
                linelist = []
                for v in vs[0:-1]:
                    linelist.append(d[v][ii])
                formatstring = "%s "*len(vs)+"\n"
                fxyz.write(formatstring  % tuple(linelist+["%.8f"%(sqrt(vx**2+vy**2)),]))




                    
