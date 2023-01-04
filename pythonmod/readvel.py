from __future__ import division
import copy
import re
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
    for i in xrange(n):
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
#             print "Error: unable to determine the number of chunks"
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
#             for i in xrange(nor):
#                 line = readnumline(fp)
#                 if line == 'None':
#                     print "Error: the chunk is not complete"
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
        d2,keys2 = readindex(fp)
        rkey = None
        for key in knoc:
            if key in d1:
                rkey = key
                break
        if rkey==None:
            print "Error: unable to determine the number of rows"
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
            for i in xrange(nor):
                line = readnumline(fp)
                if line == 'None':
                    print "Error: the row is not complete"
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
    pat = re.compile(r'"(\w+)"')
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
        l1 = []
        d2,keys2 = readindex(fp)
        if 'Number-of-rows' not in d1:
            print "Error: unable to determine the number of rows"
            exit(0)
        while True:
            line = readnumline(fp)
            if line == None:
                break
            for i,item in enumerate(line):
                d1[keys1[i]].append(item)
            nor = d1["Number-of-rows"][-1]
            dc = {}
            for m in keys2:
                dc[m] = []
            for i in xrange(nor):
                line = readnumline(fp)
                if line == 'None':
                    print "Error: the row is not complete"
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
                Type = []
                X = []
                Y = []
                Z = []
                n = 0
                for i in xrange(lst[-1]['N']):
                    (aid,atype,x,y,z) = readnumline(fp)
                    if atype in typelist or len(typelist) == 0:
                        n+=1
                        ID.append(aid)
                        Type.append(atype)
                        X.append(x)
                        Y.append(y)
                        Z.append(z)
                lst[-1]['N']=n
                lst[-1]['id']=ID
                lst[-1]['type']=Type
                lst[-1]['x']=X
                lst[-1]['y']=Y
                lst[-1]['z']=Z
    return lst

def select(lst,fun):
    res = []
    for i in lst:
        if fun(i):
            res.append(i)
    return res







