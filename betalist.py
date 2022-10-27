import math
from math import log
import rational
from rational import rational,gcd

import itertools
        
class betalist:
    
    ab_key = lambda p : 2*p[1]-1/(p[2]+1)
    c_key = lambda p : float(p[0])
        
    def int_solver(c,a,b):
        if type(c) is int:
            c = rational(c,1)
        if c.num == 0:
            return c
        o = rational(1,1)
        d = 1
        if a > b:
            d += a
            for b1 in range(1,b+1):
                c._multiply(rational(b1,d))
                o._multiply(rational(b1,d))
                d += 1
        else:
            d += b
            for a1 in range(1,a+1):
                c._multiply(rational(a1,d))
                o._multiply(rational(a1,d))
                d += 1
        c._multiply(rational(1,d))
        o._multiply(rational(1,d))
        return c
    
    #parts are a 3-tuple of rational,int,int (the last ints are the alpha-1 and beta-1 parameters (0s are the uniform distr.)
    def __init__(self,parts=[(rational(1,1),0,0),]):
        self.parts = []
        for p in parts:
            if type(p[0]) is rational:
                self.parts.append((p[0].copy(),p[1],p[2]))
            elif type(p[0]) is int:
                self.parts.append((rational(p[0],1),p[1],p[2]))
        
    def copy(self):
        return betalist(self.parts)
    
    def trim(self):
        self.parts = sorted(self.parts,key=betalist.ab_key)
        i = 0
        while i < len(self.parts):
            i2 = i+1
            for i2 in range(i+1,len(self.parts)):
                if (self.parts[i][1] != self.parts[i2][1]) or (self.parts[i][2] != self.parts[i2][2]):
                    break
                self.parts[i] = (self.parts[i][0].add(self.parts[i2][0]), self.parts[i][1], self.parts[i][2])
                self.parts[i2] = (rational(0,1), self.parts[i][1], self.parts[i][2])
            i = i2
        self.parts = [p for p in self.parts if not p[0].num == 0]
    
    def derivative(self,evaltype):
        if evaltype != "P":
            return self.zero()
        retlist = self.copy()
        to_add = []
        to_remove = []
        for pi,p in enumerate(retlist.parts):
            if p[1] > 0:
                retlist.parts[pi] = (p[0].scalar_multiply(p[1]+p[2]+1),p[1]-1,p[2])
            else:
                to_remove.append(p)
            if p[2] > 0:
                to_add.append((p[0].scalar_multiply(-(p[1]+p[2]+1)),p[1],p[2]-1))
        #print(retlist.parts, to_remove, to_add)
        for pr in to_remove:
            retlist.parts.remove(pr)
        retlist.parts += to_add
        retlist.trim()
        return retlist
    
    def zero(self):
        return betalist([])
    
    def __str__(self):
        if len(self.parts) == 0:
            return ""
        retstring="("
        for p in self.parts:
            retstring += str(p[0])+" B["+str(p[1])+","+str(p[2])+"] + "
        return retstring[:-3]+ ")"
    
    def str_internal(self,pre):
        if len(self.parts) == 0:
            return ""
        return str(self) + " " + pre +  " + "
    
    def evaluate(self,p,evaltype):
        if evaltype != "P":
            return self.copy()
        if type(p) is not rational:
            print("Non rational beta dist not implemented yet")
            print(p)
            exit()
            
        retval = rational(0,1)
        for pa in self.parts:
            coeff = pa[0].copy()
            alpha = pa[1]
            beta = pa[2]
            d = 1
            q = rational(1,1).subtract(p)
            for a in range(alpha):
                coeff._multiply(p)
                coeff._scalar_multiply(d)
                coeff._divide(rational(a+1,1))
                d += 1
            for b in range(beta):
                coeff._multiply(q)
                coeff._scalar_multiply(d)
                coeff._divide(rational(b+1,1))
                d += 1
            coeff._scalar_multiply(d)
            retval._add(coeff)
        return retval
    
    def integrate_unit(self,evaltype):
        if evaltype != "P":
            return self.copy()
        retval = rational(0,1)
        for p in self.parts:
            retval._add(p[0])
        return retval
    
    def power(self,k):
        if (type(k) is int) and (k > 0):
            retval = betalist()
            for i in range(k):
                retval._multiply(self)
            return retval
        print("Non-integer power")
        exit()
    
    def add(self,r):
        retval = self.copy()
        retval._add(r)
        return retval
        
    def _add(self,r):
        self.parts += [(p[0].copy(),p[1],p[2]) for p in r.parts]
        self.trim()
    
    def subtract(self,r):
        retval = self.copy()
        retval._subtract(r)
        return retval
        
    def _subtract(self,r):
        self._add(r.scalar_multiply(-1))
        
    def multiply(self,r):
        retval = self.copy()
        retval._multiply(r)
        return retval
        
    def _multiply(self,r):
        combined = []
        for p1,p2 in itertools.product(self.parts,r.parts):
            a1 = max(p1[1],p2[1])
            b1 = max(p1[2],p2[2])
            a2 = min(p1[1],p2[1])
            b2 = min(p1[2],p2[2])
            divisor = rational(1,1)
            at = a1+1
            for a2i in range(a2):
                divisor._multiply(rational(at,a2i+1))
                at += 1
            bt = b1+1
            for b2i in range(b2):
                divisor._multiply(rational(bt,b2i+1))
                bt += 1
            dt = max(p1[1]+p1[2],p2[1]+p2[2])+2
            for ci in range(min(p1[1]+p1[2],p2[1]+p2[2])+1):
                divisor._multiply(rational(ci+1,dt))
                dt += 1
            divisor._multiply(rational(dt-1,1))
            newc = divisor.multiply(p1[0]).multiply(p2[0])
            combined.append((newc,a1+a2,b1+b2))
        self.parts = combined
        self.trim()
    
    def scalar_multiply(self,k):
        retval = self.copy()
        retval._scalar_multiply(k)
        return retval
    
    def _scalar_multiply(self,k):
        if type(k) is rational:
            self.parts = [(p[0].multiply(k),p[1],p[2]) for p in self.parts]
            return
        if type(k) is int:
            self.parts = [(p[0].multiply(rational(k,1)),p[1],p[2]) for p in self.parts]
            return
        print("Shouldn't get here - nonrational scalar multiply")
        exit()

        