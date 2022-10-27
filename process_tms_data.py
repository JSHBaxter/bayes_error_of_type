import bayes
from bayes import annotators
import rational
from rational import rational

from os.path import exists
import pickle
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

annotator = annotators(3)
z = rational(11,10)
p_max = 201
n_max = 30 

raw_data = [[32,10,2],
            [31,11,2],
            [30,12,2],
            [28,16,0],
            [30,11,3],
            [25,19,0]]
names = ["LFACEMC",
         "RFACEMC",
         "LLLIMBMC",
         "RLLIMBMC",
         "LULIMBMC",
         "RULIMBMC"]


for ind,(point_res,point_name) in enumerate(zip(raw_data,names)):
    
    pkl_filename = point_name+"_"+str(p_max)+"_"+str(n_max)+"_"+str(z.num)+"-"+str(z.den)+".pkl"
    pdf_filename = point_name+"_"+str(p_max)+"_"+str(n_max)+"_"+str(z.num)+"-"+str(z.den)+".pdf"
    
    #if we have already calculated this data, load it
    if exists(pkl_filename):
        with open(pkl_filename, 'rb') as file:
            p_evals,n_evals,res_table = pickle.load(file)
            file.close()
        
    else:
    
        #construct the model
        print("Starting " + point_name)
        annotator.clear()
        annotator.set_default_z(z)
        for ci,cn in enumerate(point_res):
            annotator.add_case(ci,cn)

        #get the results
        print("Calculating probs " + point_name)
        p_evals = np.zeros((p_max,n_max))
        n_evals = np.zeros((p_max,n_max))
        res_table = np.zeros((p_max,n_max))
        for pr in range(p_max):
            p = rational(pr,p_max-1)
            for nr in range(n_max):
                n = nr+1
                p_evals[pr,nr] = float(p)
                n_evals[pr,nr] = float(n)
                res_table[pr,nr] = float(annotator.prob(p=p,n=n))
        print("Done " + point_name)

        with open(pkl_filename, 'wb') as file:
            pickle.dump((p_evals,n_evals,res_table), file)
            file.close()
        
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    #ax.plot_surface(n_evals, p_evals, res_table, alpha=0.3)
    #ax.contourf(n_evals, p_evals, res_table, zdir='z', offset=0, levels=np.linspace(-0.01,np.max(res_table),1200), cmap=cm.coolwarm)
    ax.plot_surface(n_evals, p_evals, np.zeros_like(res_table),
                    facecolors=plt.cm.coolwarm(res_table),shade=True,
                    rstride=1, cstride=1)
    
    #plot side information
    #cmap = plt.cm.get_cmap("cividis")
    res_sum = np.zeros_like(n_evals[0,:])
    for pr in range(p_max):
        res_sum += res_table[pr,:]
    ax.plot(n_evals[0,:], n_max*res_sum/np.sum(res_sum), zdir='y',zs=1, c='k')
    res_sum = np.zeros_like(p_evals[:,0])
    for nr in range(n_max):
        res_sum += res_table[:,nr]
    ax.plot(p_evals[:,0], p_max*res_sum/np.sum(res_sum), zdir='x',zs=1, c='k')
    
    ax.set_xlabel("n")
    ax.set_xlim(1, n_max+1)
    ax.set_ylabel("p")
    ax.set_ylim(0, 1)
    ax.set_zticklabels([])
    ax.set_title("Bayesian Model Results - " + point_name)
    plt.savefig(pdf_filename)
    plt.close(fig)
    
    
    #estimate summary stats on p
    mean_p = np.sum( p_evals[int(p_max/2):,:]*res_table[int(p_max/2):,:] ) / np.sum(res_table[int(p_max/2):,:])
    mean_p2 = np.sum( (p_evals[int(p_max/2):,:]**2)*res_table[int(p_max/2):,:] ) / np.sum(res_table[int(p_max/2):,:])
    std_p = np.sqrt(mean_p2-mean_p**2)
    print(point_name,"p",mean_p,std_p)
    
    #estimate summary stats on n
    mean_n = np.sum( n_evals[int(n_max/2):,:]*res_table[int(n_max/2):,:] ) / np.sum(res_table[int(n_max/2):,:])
    mean_n2 = np.sum( (n_evals[int(n_max/2):,:]**2)*res_table[int(n_max/2):,:] ) / np.sum(res_table[int(n_max/2):,:])
    std_n = np.sqrt(mean_n2-mean_n**2)
    print(point_name,"n",mean_n,std_n)