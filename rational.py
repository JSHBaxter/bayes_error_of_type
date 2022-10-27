import math
from math import log
        
def gcd(a,b):
    a = abs(a)
    b = abs(b)
    if b > a:
        s = a
        a = b
        b = s
    if b == 0:
        return a
    r = a % b
    while r > 0:
        a = b
        b = r
        r = a % b
    return b
    
class rational:
    
    def __init__(self,num=1,den=1):
        if type(num) is not int or type(den) is not int:
            print(num,den)
        self.num = num
        self.den = den
        
    def find_close(r,d_range=100000):
        min_diff = r+1
        min_d = 1
        min_n = 1
        for d in range(1,d_range):
            n = int(round(r * d))
            if abs(r - n/d) < min_diff:
                min_diff = abs(r - n/d)
                min_d = d
                min_n = n
        return rational(n,d)
        
    def copy(self):
        return rational(self.num,self.den)
    
    def derivative(self,evaltype):
        return self.zero()
    
    def __abs__(self):
        return rational(abs(self.num),abs(self.den))
    
    def equal(self,r):
        return ((self.num == r.num) and (self.den == r.den))
    
    def less_equal_than(self,r):
        return (self.num*r.den <= self.den*r.num)
    
    def less_than(self,r):
        return (self.num*r.den < self.den*r.num)
    
    def trim(self):
        if self.den < 0:
            self.den *= -1
            self.num *= -1
        gcdt = gcd(self.num,self.den)
        if gcdt > 1:
            self.num = self.num // gcdt
            self.den = self.den // gcdt
    
    def zero(self):
        return rational(0,1)
    
    def __str__(self):
        if self.den != 1:
            return str(self.num) + "/" + str(self.den)
        return str(self.num)
    
    def __float__(self):
        return self.num/self.den
    
    def str_internal(self,pre):
        if self.num == 0:
            return ""
        return str(self) + " " + pre +  " + "
    
    def evaluate(self,unused1=None,unused2=None):
        return float(self)
    
    def power(self,k):
        if self.num == 0:
            return rational(0,1)
        if k == 0:
            return rational(1,1)
        if k < 0:
            return rational(self.den,self.num).power(-k)
        if type(k) is int:
            return rational(self.num**k,self.den**k)
        print("Non-integer power")
        exit()
    
    def add(self,r):
        retval = self.copy()
        retval._add(r)
        return retval
        
    def _add(self,r):
        gcd_den = gcd(self.den,r.den)
        multi1 = r.den // gcd_den
        multi2 = self.den // gcd_den
        
        new_num = self.num*multi1 + r.num*multi2
        new_den = multi1 * multi2 * gcd_den
        
        self.num = new_num
        self.den = new_den
        self.trim()
    
    def subtract(self,r):
        retval = self.copy()
        retval._subtract(r)
        return retval
        
    def _subtract(self,r):
        self._add(r.scalar_multiply(-1))
        
    def divide(self,r):
        return self.multiply(r.recip())
    
    def _divide(self,r):
        self._multiply(r.recip())
        
    def multiply(self,r):
        retval = self.copy()
        retval._multiply(r)
        return retval
        
    def _multiply(self,r):
        gcd1 = gcd(self.num,r.den)
        gcd2 = gcd(self.den,r.num)
        
        new_num = (self.num // gcd1) * (r.num // gcd2)
        new_den = (self.den // gcd2) * (r.den // gcd1)
        
        self.num = new_num
        self.den = new_den
    
    def recip(self):
        return rational(self.den,self.num)
    
    def scalar_multiply(self,k):
        if type(k) is rational:
            return self.multiply(k)
        if type(k) is int:
            return self.multiply(rational(k,1))
        return float(self)*k
    
    def _scalar_multiply(self,k):
        if type(k) is rational:
            self._multiply(k)
            return
        if type(k) is int:
            self._multiply(rational(k,1))
            return
        print("Shouldn't get here - nonrational scalar multiply")
        exit()
    