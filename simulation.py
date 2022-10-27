import bayes
from bayes import annotators
import rational
from rational import rational

import random
import numpy as np
import math

import pickle

def single_simulation(cases,a,p,n,t):
    counts = []
    for ti in range(t):
        picked = [max(int((n/(1-p))*(random.random() - p) + 1),0) for a1 in range(a)]
        count = [sum([1 for p in picked if p == c]) for c in range(n+1)]
        count = tuple(sorted([c for c in count if c != 0],reverse=True))
        counts.append(count)
    case_counts = [sum([1 for c in counts if c == case]) for case in cases]
    return case_counts

def annotator_simulator(a,z,zn,nt,s):
    res = []
    z = rational.find_close(z)
    
    annotator = annotators(a)
    print(annotator.cases)
    for i in range(s):
    
        #pick random value of n greater than 0
        n = np.random.geometric(p=1-1/zn)
    
        #pick random value of t around
        t = np.random.poisson(nt-1)+1

        #pick random value of p between 50% and 100%
        p = 0.5 * random.random() + 0.5
    
        simulation = single_simulation(annotator.cases,a,p,n,t)
        print(a,i,t,"\t",n,"\t",p,"\t",simulation,end="")
        
        for ic,cc in enumerate(simulation):
            if cc > 0:
                annotator.add_case(ic,cc)
        
        prob_correct=annotator.prob(z,p=p,n=n)
        prob_correct_p=annotator.prob(z,p=p,n=None)
        prob_correct_n=annotator.prob(z,p=None,n=n)
        print(float(prob_correct),
              float(prob_correct_p),
              float(prob_correct_n),end="")
        modes = annotator.get_mode_p(z,num_modes=3)
        print(modes)
        res.append((a,z,zn,n,nt,t,p,simulation,prob_correct,prob_correct_p,prob_correct_n,modes))
        
        annotator.clear()
        
    return res
        
#sim_params = (k, z, z', t)
sim_type = [(3,1.25,1.25,20),
            (4,1.25,1.25,20),
            (3,1.25,2,20),
            (4,1.25,2,20),
            (3,2,1.25,20),
            (4,2,1.25,20)]

for sim in sim_types:
    res = annotator_simulator(sim[0],sim[1],sim[2],sim[3],500)

    with open('res_simulation_dump'+str(sim)+'.pkl', 'wb') as file:
        pickle.dump(res, file)
        file.close()

