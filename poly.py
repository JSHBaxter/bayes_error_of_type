import re
import rational
from rational import rational
import betalist
from betalist import betalist

class poly:
    def __init__(self, polytype, init_coeffs, base=None):
        self.polytype = polytype
        self.order = len(init_coeffs)-1
        self.base = not (type(init_coeffs[0]) is poly)
        
        init_coeffs_used = [None for i in init_coeffs]
        for i,v in enumerate(init_coeffs):
            if type(v) is int:
                init_coeffs_used[i] = rational(v,1)
            elif type(v) is float:
                init_coeffs_used[i] = rational.find_close(v)
            else:
                init_coeffs_used[i] = v
        
        if type(init_coeffs_used[0]) is poly:
            self.polytype += init_coeffs_used[0].polytype
        self.coeffs = [i.copy() for i in init_coeffs_used]
            
        self.trim()
        
    def copy(self):
        return poly(self.polytype[0], self.coeffs)
    
    def __str__(self):
        return self.str_internal("")[:-3]
    
    def str_internal(self,pre):
        retstring = ""
        for i,c in enumerate(self.coeffs):
            p = self.order-i
            if p == 0:
                pre_down = pre
            else:
                pre_down = pre + self.polytype[0] + "^" + str(p) + " "
            retstring += c.str_internal(pre_down)
        return retstring
        
    def zero(self):
        retpoly = self.copy()
        retpoly.order = 0
        retpoly.coeffs = [retpoly.coeffs[0].zero()]
        return retpoly
        
    def extend(self,k):
        if k < 1:
            return
        self.coeffs = [self.coeffs[-1].zero() for i in range(k)] + self.coeffs
        self.order += k
        
    def trim(self):
        delete_order = 0
        for e in self.coeffs:
            e.trim()
            if type(e) is rational and e.num == 0:
                delete_order += 1
            else:
                break
        if delete_order > self.order:
            self.order = 0
            self.coeffs = [self.coeffs[-1].zero()]
        elif delete_order > 0:
            self.coeffs = self.coeffs[delete_order:]
            self.order -= delete_order

    def scalar_multiply(self,k):
        if k == 0:
            return self.zero()
        retpoly = self.copy()
        retpoly._scalar_multiply(k)
        return retpoly
        
    def _scalar_multiply(self,k):
        for i in range(self.order+1):
            self.coeffs[i] = self.coeffs[i].scalar_multiply(k)
        #self.trim()
    
    def subtract(self,p):
        retpoly = self.copy()
        retpoly._subtract(p)
        return retpoly
    
    def _subtract(self,p):
        if p.polytype != self.polytype:
            print("Can only subtract polynomials of the same type for now")
            exit()
            
        if p.order <= self.order:
            for i in range(p.order+1):
                self.coeffs[-i-1]._subtract(p.coeffs[-i-1])
        else:
            self.extend(p.order-self.order)
            self._subtract(p)
        self.trim()
        
    
    def add(self,p):
        retpoly = self.copy()
        retpoly._add(p)
        return retpoly
    
    def _add(self,p):
        if p.polytype != self.polytype:
            print("Can only add polynomials of the same type for now")
            exit()
        if p.order <= self.order:
            for i in range(p.order+1):
                self.coeffs[-i-1]._add(p.coeffs[-i-1])
        else:
            self.extend(p.order-self.order)
            self._add(p)
        self.order = max([self.order,p.order])
        self.trim()
    
    def multiply(self,p):
        retpoly = self.copy()
        retpoly._multiply(p)
        return retpoly
        
        
    def _multiply(self,p):
        
        #if we have the exact same type, use a regular multiply
        if p.polytype == self.polytype:
            self._reg_multiply(p)
            return
        
        polytype_match = ""
        for c in p.polytype:
            if c in self.polytype:
                polytype_match += c
        if len(polytype_match) == 0:
            self._cyl_multiply(p)
            return
        
        print("Cannot find a functional multiplication")
        print(self)
        print(p)
        exit()

    #multiply together the polynominals noting they are of a different type
    def _cyl_multiply(self,p):
        if self.base:
            for i in range(self.order+1):
                self.coeffs[i] = p.scalar_multiply(self.coeffs[i])
        else:
            for i in range(self.order+1):
                self.coeffs[i] = self.coeffs[i].multiply(p)
        self.polytype = self.polytype[0]+self.coeffs[0].polytype
        self.base = False
            
        
    #multiply polynomials noting they have the same primary type
    def _reg_multiply(self,p):
        #create container
        torder = self.order+p.order
        if not p.base:
            new_coeffs = [p.coeffs[-1].zero() for i in range(torder+1)]
        else:
            new_coeffs = [self.coeffs[-1].zero() for i in range(torder+1)]
        
        for t in range(torder+1):
            for st in range(self.order+1):
                pt = t - st
                if pt < 0 or pt > p.order:
                    continue
                if self.base == p.base:
                    add_coeff = self.coeffs[self.order-st].multiply(p.coeffs[p.order-pt])
                    new_coeffs[torder-t] = new_coeffs[torder-t].add(add_coeff)
                elif self.base and not p.base:
                    add_coeff = p.coeffs[p.order-pt].scalar_multiply(self.coeffs[self.order-st])
                    new_coeffs[torder-t] = new_coeffs[torder-t].add(add_coeff)
                else:# not self.base and p.base:
                    add_coeff = s.coeffs[self.order-st].scalar_multiply(p.coeffs[p.order-pt])
                    new_coeffs[torder-t] = new_coeffs[torder-t].add(add_coeff)
        self.coeffs = new_coeffs
        self.order = torder
        self.base = (self.base or p.base)
        self.trim()

    def derivative(self,evaltype):
        if evaltype == self.polytype[0]:
            retval = self.copy()
            retval.coeffs = retval.coeffs[:-1]
            if retval.order == 0:
                return self.zero()
            for i,c in enumerate(retval.coeffs):
                retval.coeffs[i] = c.scalar_multiply(retval.order-i)
            retval.order -= 1
            return retval
        
        retval = self.copy();
        for c in range(retval.order+1):
            retval.coeffs[c] = retval.coeffs[c].derivative(evaltype)
        retval.trim()
        return retval
    
    def evaluate(self,x,evaltype):
        #if we are a base, evaluate immediately
        if self.base and (evaltype != self.polytype):
            retpoly = self.copy()
            if type(self.coeffs[0]) is betalist:
                retpoly.coeffs = [c.evaluate(x,evaltype) for c in self.coeffs]
            return retpoly
        
        #else, see if we can evaluate immediately and strip the outer type
        if self.base or (evaltype == self.polytype[0]):
            retval = self.coeffs[0].copy()
            if self.order == 0:
                return retval
            for v in self.coeffs[1:]:
                retval._scalar_multiply(x)
                retval._add(v)
            return retval
        
        #else, pass to a lower level and modify type name
        retval = self.copy()
        retval.polytype = re.sub(evaltype,"",retval.polytype)
        for c in range(retval.order+1):
            retval.coeffs[c] = retval.coeffs[c].evaluate(x,evaltype)
        retval.base = (len(retval.polytype) == 1)
        retval.trim()
        return retval
    
    def integrate_unit(self,evaltype):
        #if we are a base, evaluate immediately
        if self.base and (evaltype != self.polytype):
            retpoly = self.copy()
            if type(self.coeffs[0]) is betalist:
                retpoly.coeffs = [c.integrate_unit(evaltype) for c in self.coeffs]
            return retpoly
        
        #else, see if we can evaluate immediately and strip the outer type
        if self.base or (evaltype == self.polytype[0]):
            retval = self.coeffs[-1].copy()
            if self.order == 0:
                return retval
            for o,v in enumerate(reversed(self.coeffs[:-1])):
                retval._add(v.scalar_multiply(rational(1,(o+1))))
            return retval
        
        #else, pass to a lower level and modify type name
        retval = self.copy()
        retval.polytype = re.sub(evaltype,"",retval.polytype)
        for c in range(retval.order+1):
            retval.coeffs[c] = retval.coeffs[c].integrate_unit(evaltype)
        retval.base = (len(retval.polytype) == 1)
        retval.trim()
        return retval
    
    def convert(self,newtype):
        #print(self.polytype + " -> " + newtype + "\t" + str(self))
        retpoly = self.copy()
        
        #if we are basic, we cannot reorder
        if self.base:
            #print("\tImmediately return")
            return retpoly
        
        #if we match the top term, no need to reorder at this level, so just pass down
        if self.polytype == newtype:
            #print("\tImmediately return")
            return retpoly
        if self.polytype[0] == newtype[0]:
            #print("\tPass down (" + retpoly.coeffs[0].polytype + " -> " + newtype[1:] + ")")
            for i in range(len(retpoly.coeffs)):
                retpoly.coeffs[i] = retpoly.coeffs[i].convert(newtype[1:])
            retpoly.polytype = self.polytype[0]+retpoly.coeffs[0].polytype
            return retpoly
        
        #we don't match the top term, switch it down and then pass down 
        #(effectively recrusively bubbling terms, only efficient for small polys - quad  stack growth)
        
        #reorder the first two terms
        #print("\tReorder terms "+str(self)+" ",end="")
        max_order = max([p.order for p in self.coeffs])
        type0 = self.polytype[0]
        type1 = self.polytype[1]
        if self.coeffs[0].base:
            all_coeffs = [[0 for i in range(self.order+1)] for i in range(max_order+1)]
        else:
            all_coeffs = [[self.coeffs[0].coeffs[0].zero() for i in range(self.order+1)] for i in range(max_order+1)]
        for si,sv in enumerate(reversed(self.coeffs)):
            for pi,pv in enumerate(reversed(sv.coeffs)):
                all_coeffs[-pi-1][-si-1] = pv
        retpoly.coeffs = [poly(type0,c) for c in all_coeffs]
        retpoly.polytype = type1+retpoly.coeffs[0].polytype
        retpoly.order = len(retpoly.coeffs)-1
        #print(str(retpoly))
        
        for i in range(len(retpoly.coeffs)):
            retpoly.coeffs[i] = retpoly.coeffs[i].convert(newtype[1:])
        retpoly.polytype = type1+retpoly.coeffs[0].polytype
        #print(retpoly.polytype)
        
        if retpoly.polytype[1:] == newtype[1:]:
            return retpoly
        else:
            return retpoly.convert(newtype)

class N_poly(poly):
    def __init__(self, init_coeffs, base=None):
        super().__init__("N",init_coeffs,base)
        
    def copy(self):
        return N_poly(self.coeffs)
    
    def __str__(self):
        return self.str_internal("")[:-3]
    
    def trim(self):
        delete_order = 0
        for e in reversed(self.coeffs):
            e.trim()
            if type(e) is rational and e.num == 0:
                delete_order += 1
            else:
                break
        if delete_order > self.order:
            self.order = 0
            self.coeffs = [self.coeffs[-1].zero()]
        elif delete_order > 0:
            self.coeffs = self.coeffs[:-delete_order]
            self.order -= delete_order
            
    def str_internal(self,pre):
        retstring = ""
        for i,c in enumerate(self.coeffs):
            if i == 0:
                pre_down = pre
            else:
                pre_down = pre + self.polytype[0] + "^-" + str(i) + " "
            retstring += c.str_internal(pre_down)
        return retstring
    
    def evaluate(self,x,evaltype):
        #if we are a base, evaluate immediately
        if self.base and (evaltype != self.polytype):
            retpoly = self.copy()
            if type(self.coeffs[0]) is betalist:
                retpoly.coeffs = [c.evaluate(x,evaltype) for c in self.coeffs]
            return retpoly
        
        #else, see if we can evaluate immediately and strip the outer type
        if self.base or (evaltype == self.polytype[0]):
            retval = self.coeffs[0].copy()
            nden = rational(1,x)
            if self.order == 0:
                return retval
            for v in self.coeffs[1:]:
                retval._add(v.scalar_multiply(nden))
                nden._multiply(rational(1,x))
            return retval
        
        #else, pass to a lower level and modify type name
        retval = self.copy()
        retval.polytype = re.sub(evaltype,"",retval.polytype)
        for c in range(retval.order+1):
            retval.coeffs[c] = retval.coeffs[c].evaluate(x,evaltype)
        retval.base = (len(retval.polytype) == 1)
        retval.trim()
        return retval
