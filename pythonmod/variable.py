from math import pi
import fileinput
import re
import numpy as np

PI = pi
vareq0 = re.compile(r'[ \t]*variable[ \t]+(\w+)[ \t]+equal.*') 
vareq = re.compile(r'[ \t]*variable[ \t]+(\w+)[ \t]+(equal|loop)[ \t]+([0-9\.\-]+)[ \t]*([0-9\.\-]*)[^0-9\.\+\-\*\/\{\}\(\)\$]*\n')
varexp = re.compile(r'[ \t]*variable[ \t]+(\w+)[ \t]+equal[ \t]+([0-9\.\+\-\*\/\w\{\}\(\)\$ \t]+).*')

def list2execstr(lst):
    res = 'np.array(['
    for i in lst:
        res += (str(i))+','
    res += '])'
    return res

def evaluate(exprs,varis):
    if len(exprs)==0:
        return varis
    for key,expr in exprs.items():
        try:
            v = eval(expr)
            varis[key] = v
            exprs.pop(key)
        except SyntaxError:
            pass
    for key,expr in exprs.items():
        for var, val in varis.items():
            varm = '${'+var+'}'
            if varm in expr:
                if type(val) is list or type(val) is np.ndarray:
                    #print val
                    exprs[key] = expr.replace(varm,list2execstr(val))
                else:    
                    exprs[key] = expr.replace(varm,str(val))
    return evaluate(exprs,varis)

def readallvar(filename):
    varis = dict()
    exprs = dict()
    with open(filename,'r') as fp:
        for line in fp:
            mateq = vareq.search(line)
            if mateq:
                #print line
                name = mateq.group(1)
                eq = mateq.group(2)
                val0 = mateq.group(3)
                try:
                    val1 = mateq.group(4)
                except IndexError:
                    val1 = None
                if eq == 'equal':
                    varis[name] = float(val0)
                elif eq == 'loop':
                    if val1:
                        varis[name] = np.array(range(int(val0),int(val1)+1))
                    else:
                        varis[name] = np.array(range(1, int(val0)+1))
            else:
                matexp = varexp.search(line)
                if matexp:
                    #print line
                    name = matexp.group(1)
                    expr = matexp.group(2)
                    exprs[name] = expr
    return evaluate(exprs,varis)

def readvar(filename, var):
    d = readallvar(filename)
    try:
        return d[var]
    except KeyError:
        print 'Can not find variable name ',var
        return None


def setvar(filename, var, val):
    sign = False
    for line in fileinput.input(filename,inplace=True,backup='.bak'):
        mateq = vareq0.search(line)
        if mateq:
            name = mateq.group(1)
            if name == var:
                print 'variable %s equal %s'%(var, str(val))
                sign = True
            else:
                print line.strip()
        else:
            print line.strip()
                
    if sign:
        print 'variable %s set to %s'%(var,str(val))
    else:
        print 'can not find variable %s'%var
        for line in fileinput.input(filename,inplace=True,backup='.bak'):
            mateq = vareq0.search(line)
            if mateq and not sign:
                print 'variable %s equal %s'%(var, str(val))+'\n'+line.strip()
                sign = True
            else:
                print line.strip()
        print 'variable %s has been added and set to %s'%(var,str(val))



