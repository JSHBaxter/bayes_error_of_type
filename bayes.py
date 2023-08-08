import poly_exp_prior
from poly import *
from betalist import *
from rational import *
import math
from math import comb
from math import exp
    
class annotators:
    
    
    def list_possibility(self,number):
        if number <= self.num_annotators:
            return
        self.list_possibility(number-1)
        factor = rational(number+1,1)
        
        self.num_annotators = number
        
        #general convenient function for updating case when appending
        def inc_list_first_val(l,v):
            ret_l = [c for c in l]
            ret_l[ret_l.index(v)] += 1
            return tuple(sorted(ret_l,reverse=True))
        
        #infrastructure for adding new cases
        new_adv_case_possibilities = {}
        def safe_add_case(k,v):
            if k not in new_adv_case_possibilities.keys():
                new_adv_case_possibilities[k] = 0
            new_adv_case_possibilities[k] += v
        
        #update cases by adding new item at the end
        for adv_case,num_configs in self.adv_case_possibilities.items():
            
            #add new G
            if adv_case[1] == 0:
                new_adv_case = (adv_case[0],1,adv_case[2],tuple(list(adv_case[3])+[1,]),adv_case[4])
            elif adv_case[0] == 0:
                new_adv_case = (adv_case[0],adv_case[1]+1,adv_case[2],inc_list_first_val(adv_case[3],adv_case[1]),adv_case[4]+1)
            else:
                new_adv_case = (adv_case[0],adv_case[1]+1,adv_case[2],inc_list_first_val(adv_case[3],adv_case[1]),adv_case[4])
            safe_add_case(new_adv_case,num_configs)
            
            #add new N
            new_adv_case = (adv_case[0],adv_case[1],adv_case[2]+1,tuple(list(adv_case[3])+[1,]),adv_case[4])
            safe_add_case(new_adv_case,num_configs)
            
            #add existing N
            if adv_case[0] == 0:
                N_c = [c for i,c in enumerate(adv_case[3]) if adv_case[1] == 0 or i != adv_case[3].index(adv_case[1])]
                for i in N_c:
                    new_adv_case = (adv_case[0],adv_case[1],adv_case[2],inc_list_first_val(adv_case[3],i),adv_case[4])
                    safe_add_case(new_adv_case,num_configs)
            else:
                N_c = [c for i,c in enumerate(adv_case[3]) if adv_case[1] == 0 or i != adv_case[3].index(adv_case[1])]
                N_c = [c for i,c in enumerate(N_c) if i != N_c.index(adv_case[4])]
                for i in N_c:
                    new_adv_case = (adv_case[0],adv_case[1],adv_case[2],inc_list_first_val(adv_case[3],i),adv_case[4])
                    safe_add_case(new_adv_case,num_configs)
                new_adv_case = (adv_case[0],adv_case[1],adv_case[2],inc_list_first_val(adv_case[3],adv_case[4]),adv_case[4]+1)
                safe_add_case(new_adv_case,num_configs)
        self.adv_case_possibilities = new_adv_case_possibilities
        
        #make new empty polys
        for adv_case,num_configs in self.adv_case_possibilities.items():
            self.case_polys[adv_case[3]] = N_poly([betalist([(0,0,0)])])
            
        #add new polys
        for adv_case,num_configs in self.adv_case_possibilities.items():
            pcount = adv_case[1]
            qcount = number-pcount
            coeff = betalist.int_solver(1,pcount,qcount)
            poly_coeffs = poly("N",[factor]+[rational(0,1)]*pcount)
            for i in range(adv_case[2]):
                poly_coeffs._multiply(poly("N",[rational(1,1),rational(-i,1)]))
            org_poly_coeffs = poly("N",[betalist([(coeff,pcount,qcount)]).scalar_multiply(co) for co in poly_coeffs.coeffs])
            poly_coeffs = N_poly([betalist([(coeff,pcount,qcount)]).scalar_multiply(co) for co in poly_coeffs.coeffs]).scalar_multiply(num_configs)
            self.case_polys[adv_case[3]]._add(poly_coeffs)
            
            if (adv_case[3],adv_case[4]) not in self.correct_point_polys_den.keys():
                self.correct_point_polys_den[(adv_case[3],adv_case[4])] = N_poly([betalist([(0,0,0)])])
            self.correct_point_polys_den[(adv_case[3],adv_case[4])]._add(poly_coeffs)
            if adv_case[0] == 0:
                if (adv_case[3],adv_case[4]) not in self.correct_point_polys_num.keys():
                    self.correct_point_polys_num[(adv_case[3],adv_case[4])] = N_poly([betalist([(0,0,0)])])
                self.correct_point_polys_num[(adv_case[3],adv_case[4])]._add(poly_coeffs)
        
        self.cases = tuple(self.case_polys.keys())

    def __init__(self):
        #create model for single annotator (can add others as necessary)
        self.n_sum_fun = lambda x : poly_exp_prior.z_transform(x,"N")
        self.num_annotators = 1
        self.adv_case_possibilities = {(0,1,0,(1,),1) : 1, (1,0,1,(1,),1) : 1}
        self.cases = (1,)
        self.case_polys = {(1,) : N_poly([betalist([(1,0,0)])])}
        self.correct_point_polys_num = {((1,),1) : N_poly([betalist([(1,1,0)])])}
        self.correct_point_polys_den = {((1,),1) : N_poly([betalist([(1,0,0)])])}
        self.clear()
        
    def increment_num_annotators(self):
        self.list_possibility(self.num_annotators+1)
        
    def __str__(self):
        return str(self.numerator) + " /\n" + str(self.pn_num_sum)
    
    def clear(self):
        self.numerator = N_poly([betalist([(1,0,0)])])
        self.z_def = None
        self.num_cases = {}
        self.convert()
    
    def set_default_z(self,z):
        self.z_def = z
        self.convert()

    def add_case(self, case, number_add):
        num_annotators = sum(case)
        while self.num_annotators < num_annotators:
            self.increment_num_annotators()
        if case in self.num_cases.keys():
            self.num_cases[case] += number_add
        else:
            self.num_cases[case] = number_add
        self.multiplicitous_add(self.case_polys[case],number_add)
            
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
        self.n_num_sum = self.n_sum_fun(self.numerator)
        self.p_num_sum = self.numerator.integrate_unit("P")
        self.pn_num_sum = self.n_sum_fun(self.p_num_sum)
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
            prob_num = self.n_sum_fun(prob_num)
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
        prob_num = self.n_sum_fun((self.numerator.multiply(p_lin)).integrate_unit("P")).evaluate(z,"Z")
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
    
    
    model_base = annotators()
    for k in range(1,7):
        #print(model_base.cases,"\n")
        
        for c in model_base.case_polys.items():
            print(c[0],c[1])
        
        for k,v1 in model_base.correct_point_polys_num.items():
            v2 = model_base.correct_point_polys_den[k]
            print(k,str(v1),"/",str(v2))
        print("\n")
        
        model_base.increment_num_annotators()


