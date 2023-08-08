import math
from math import log
import poly
from poly import poly
import rational
from rational import rational

class Z_fun:
    def __init__(self,func):
        self.func_str = [i[0] for i in func]
        self.func_l = [i[1] for i in func]
        self.coeffs = [0 for i in func]
    
    def copy(self):
        retpoly = Z_fun([i for i in zip(self.func_str,self.func_l)])
        if type(self.coeffs[0]) is poly:
            retpoly.coeffs = [i.copy() for i in self.coeffs]
        else:
            retpoly.coeffs = [i for i in self.coeffs]
        return retpoly
    
    def __str__(self):
        retstr = "["
        for c,f in zip(self.coeffs,self.func_str):
            if c == 0 or (type(c) is rational and c.num == 0):
                continue
            if str(c) != "":
                retstr += str(c) + "*" + f + " + "
        retstr = retstr[:-3] + "]"
        return retstr
    
    def evaluate(self,x,evaltype):
        retpoly = self.copy()
        
        if evaltype != "Z":
            if type(self.coeffs[0]) is poly:
                for i in range(len(retpoly.coeffs)):
                    retpoly.coeffs[i] = retpoly.coeffs[i].evaluate(x,evaltype)
            return retpoly
        else:
            if type(self.coeffs[0]) is poly:
                print("Non-rational polynomial detected")
                exit()
            else:
                retval = self.coeffs[0].zero()
                for i,(c,f) in enumerate(zip(self.coeffs,self.func_l)):
                    #print(c,x,self.func_str[i],f(x,i),c.multiply(f(x,i)))
                    #print("\t",retval,"+",end=" ")
                    retval._add(c.scalar_multiply(f(x,i)))
                    #print(c,"*",f(x,i),"=",retval)
            return retval
    
    def integrate_unit(self,evaltype):
        retpoly = self.copy()
        if type(self.coeffs[0]) is poly:
            for i in range(len(retpoly.coeffs)):
                retpoly.coeffs[i] = retpoly.coeffs[i].integrate_unit(evaltype)
        return retpoly
    
def z_transform(np,evaltype):
    #if we don't have that evaltype, treat as constant
    if evaltype not in np.polytype:
        retZfun = Z_fun([("1/(z-1)",lambda z,c: 1 / float(z-1))])
        retZfun.coeff[0] = np.copy()
        return retZfun

    #ensure than the evaltype is leading
    if evaltype != np.polytype[0]:
        newpolytype = evaltype + re.sub(evaltype,"",retval.polytype)
        return np.convert(newpolytype).z_transform(evaltype)

    #else, evaluate immediately
    z_range = 100
    gen_func = (lambda z,c: (z.power(-(c-1))))
    retZfun = Z_fun([("1/(z-1)",lambda z,c: z.add(rational(-1,1)).recip()),
                     ("ln(z/(z-1))",lambda z,c: rational.find_close(log(float(z.multiply(z.add(rational(-1,1)).recip())))))] + \
                     [("1/z^"+str(i),gen_func) for i in range(1,z_range+1)])

    for i in range(1,len(retZfun.func_l)):
        retZfun.coeffs[i] = np.coeffs[0].zero()
    retZfun.coeffs[0] = np.coeffs[0].copy()
    if np.order > 0:
        retZfun.coeffs[1] = np.coeffs[1].copy()
    if np.order > 1:
        for o,v in enumerate(np.coeffs[2:]):
            for i in range(2,len(retZfun.func_l)):
                retZfun.coeffs[i]._add(v.scalar_multiply(rational(1,((i-1)**(o+2)))))
    return retZfun
