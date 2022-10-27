from poly import *
from betalist import *
from rational import *
import math
from math import comb
from math import exp

def list_possibility(number):
    possibilities = None
    for i in range(number):
        #i takes on the 'good' value
        old_possibilities = possibilities
        if  possibilities is None:
            possibilities = [[0,]]
        else:
            possibilities = [p+[0,] for p in possibilities]
    
        #i takes on an already used value or the next highest new one
        if old_possibilities is None:
            possibilities.append([1,])
        else:
            for p in old_possibilities:
                for n in range(1,max(p)+2):
                    possibilities.append(p+[n,])
    
    #get the possibilities for each case
    case_possibilities = {}
    for p in possibilities:
        counts = [sum([1 for pi in p if pi == item]) for item in range(max(p)+1)]
        counts = tuple(sorted([c for c in counts if c > 0],reverse=True))
        if counts in case_possibilities:
            case_possibilities[counts].append(p)
        else:
            case_possibilities[counts] = [p]

    #get the probabilities for each possibility for each case
    cases = case_possibilities.keys()
    case_probabilities = {}
    multinomials_calculated = {}
    for c in cases:
        case_probabilities[c] = N_poly([betalist([(0,0,0)])])
        for p in case_possibilities[c]:
            pcount = sum([1 for pi in p if pi == 0])
            qcount = len(p) - pcount
            coeff = betalist.int_solver(1,pcount,qcount)
            poly_coeffs = N_poly([rational(1,1)]+[rational(0,1)]*pcount)
            for i in range(max(p)):
                poly_coeffs._multiply(N_poly([rational(1,1),rational(-i,1)]))
            poly_coeffs = N_poly([betalist([(coeff,pcount,qcount)]).scalar_multiply(co) for co in poly_coeffs.coeffs])
            case_probabilities[c]._add(poly_coeffs)
    
    return case_probabilities
    
class annotators:
    def __init__(self,number):
        #self.numerator = poly("N",[poly("P",[poly("Q",[1]),]),])
        self.numerator = N_poly([betalist([(1,0,0)])])
        self.num_annotators = number
        self.case_polys = list_possibility(number)
        self.cases = tuple(self.case_polys.keys())
        self.factor = rational(number+1,1)
        for c in self.cases:
            self.case_polys[c]._scalar_multiply(self.factor)
        self.num_cases = [0 for c in self.cases]
        self.z_def = None
        self.convert()
        
    def __str__(self):
        return str(self.numerator) + " /\n" + str(self.pn_num_sum)
    
    def clear(self):
        self.numerator = N_poly([betalist([(1,0,0)])])
        self.num_cases = [0 for c in self.cases]
        self.z_def = None
        self.convert()
    
    def set_default_z(self,z):
        self.z_def = z
        self.convert()
    
    def print_case(self,number):
        if number < len(self.cases):
            print(self.cases[number])
        else:
            print("Invalid case, max is "+str(len(self.cases)))

    def add_case(self, case_number, number_add):
        self.num_cases[case_number] += number_add
        self.multiplicitous_add(self.case_polys[self.cases[case_number]],number_add)
            
    def multiplicitous_add(self, base_poly, number):
        
        mult_poly = N_poly([betalist([(1,0,0)])])
        pow_poly = base_poly.copy()
        
        power = 1
        power_appl = 0
        while number > 0:
            
            if number % 2 == 1:
                mult_poly._multiply(pow_poly)
                power_appl += power
            
            pow_poly._multiply(pow_poly)
            number = number // 2
            power *= 2
        
        if power_appl > 0:
            self.numerator = self.numerator.multiply(mult_poly)
            self.num_sum_valid = False
        
    def convert(self):
        self.num_sum_valid = True
        self.n_num_sum = self.numerator.z_transform("N")
        self.p_num_sum = self.numerator.integrate_unit("P")
        self.pn_num_sum = self.p_num_sum.z_transform("N")
        if self.z_def is not None:
            self.prob_den = self.pn_num_sum.evaluate(self.z_def,"Z")
    
    def prob(self,z=None,p=None,n=None):
        if z is None and self.z_def is None:
            print("z must be specified")
            exit()
        
        if not self.num_sum_valid:
            self.convert()
            
        prob_num = self.numerator.copy()
        if self.z_def is None or (z is not None and not z.equal(self.z_def)):
            prob_den = self.pn_num_sum.evaluate(z,"Z")
        else:
            prob_den = self.prob_den
            z = self.z_def
        
        #handle p and make sure we keep the order correct, no trimming allowed, as it messes up the z-transform
        if p is None:
            old_order = prob_num.order
            prob_num = self.p_num_sum.copy()
            new_order = prob_num.order
            if new_order < old_order:
                prob_num.extend(old_order-new_order)
                
        elif p != "all":
            old_order = prob_num.order
            if type(p) is rational:
                prob_num = prob_num.evaluate(p,"P")
            else:
                prob_num = prob_num.evaluate(rational.find_close(p),"P")
            new_order = prob_num.order
            if new_order < old_order:
                prob_num.extend(old_order-new_order)
            
        #handle n
        if n is None:
            prob_num = prob_num.z_transform("N")
            prob_num = prob_num.evaluate(z,"Z")
        else:
            multiplier = z.power(-n)
            prob_num = prob_num.evaluate(n,"N")
            prob_num = prob_num.scalar_multiply(multiplier)
        
        if type(prob_num) is not float:
            prob_num.trim()
            return prob_num.scalar_multiply(prob_den.recip())
        return prob_num / float(prob_den)
    
    def prob_mean_p(self,z=None,power=1):
        if z is None and self.z_def is None:
            print("z must be specified")
            exit()
            
        if not self.num_sum_valid:
            self.convert()
        
        p_lin = N_poly([betalist([(rational(1,power+1),power,0)])])
        prob_num = (self.numerator.multiply(p_lin)).integrate_unit("P").z_transform("N").evaluate(z,"Z")
        if self.z_def is None:
            prob_den = self.pn_num_sum.evaluate(z,"Z")
        else:
            prob_den = self.prob_den
        
        if type(prob_num) is not float:
            prob_num.trim()
            return prob_num.scalar_multiply(prob_den.recip())
        return prob_num / float(prob_den)
    
    def get_mode_p(self,z,num_modes=1):
        distri = self.prob(z,p="all",n=None)
        derivative = distri.derivative("P")
        hess = derivative.derivative("P")
        d_range_inc = 1.1
        epsilon = rational(1,20)
        
        #find candidates first via sampling at intervals of 1%
        p_max = 100
        extrema = []
        tested_values = []
        for p in range(p_max+1):
            p_test = rational(p,p_max)
            v = distri.evaluate(p_test,"P")
            tested_values.append((p_test,v))
        to_remove = []
        for pti in range(1,p_max):
            ptb = tested_values[pti-1]
            pt = tested_values[pti]
            pta = tested_values[pti+1]
            if pt[1].less_than(pta[1]) or pt[1].less_than(ptb[1]):
                to_remove.append(pt)
        if tested_values[0][1].less_than(tested_values[1][1]):
            to_remove.append(tested_values[0])
        if tested_values[-1][1].less_than(tested_values[-2][1]):
            to_remove.append(tested_values[-1])
        for t in to_remove:
            tested_values.remove(t)
            
        #sort out best candidates
        extrema = sorted(tested_values,key=(lambda e : float(e[1])))
        
        #use a few rounds of binary search to improve precision
        updated_extrema = []
        for e in extrema:
            dv = derivative.evaluate(e[0],"P")
            if e[0].num == 0:
                e = (e[0],dv,rational(0,1),rational(1,p_max))
            elif e[0].num == e[0].den:
                e = (e[0],dv,rational(p_max-1,p_max),rational(1,1))
            else:
                e = (e[0],dv,e[0].subtract(rational(1,p_max)),e[0].add(rational(1,p_max)))
            for iter in range(10):
                if e[1].num > 0:
                    new_e = rational(1,2).multiply(e[0].add(e[3]))
                    dv = derivative.evaluate(new_e,"P")
                    e = (new_e,dv,e[0],e[3])
                else:
                    new_e = rational(1,2).multiply(e[0].add(e[2]))
                    dv = derivative.evaluate(new_e,"P")
                    e = (new_e,dv,e[2],e[0])
            e = e[0]
            v = distri.evaluate(e,"P")
            h = hess.evaluate(e,"P")
            updated_extrema.append((float(e),float(v),float(h)))
        extrema = sorted(updated_extrema,key=(lambda e : float(e[1])))

        #compute heights and standard deviations
        updated_extrema = []
        for e in extrema:
            w = e[1]
            mu = e[0]
            if len(updated_extrema) > 1:
                removed_weight = sum([(g[0]/(math.sqrt(2*math.pi*g[2])))*exp(-(e[0]-g[1])**2/(2*g[2])) for g in updated_extrema])
                w -= removed_weight
            if w < 0 or e[2] >= 0.0:
                continue
            updated_extrema.append((w*math.sqrt(-2*math.pi*w/e[2]),mu,-w/e[2]))
        w_sum = sum([e[0] for e in updated_extrema])
        updated_extrema = sorted([(e[0]/w_sum,e[1],math.sqrt(e[2])) for e in updated_extrema],reverse=True)
        
        #if we don't have extrema...
        if len(updated_extrema) == 0:
            mean = self.prob_mean_p(z)
            std = math.sqrt( float(self.prob_mean_p(z,2).subtract(mean.power(2))) )
            return [(1.0,float(mean),std)]
        
        return updated_extrema
    
    
    
    
    
    
if __name__ == "__main__" :
    
    
    for k in range(3,11):
        model_base = annotators(k)
        print(model_base.cases,"\n")


