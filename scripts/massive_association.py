#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 15:59:57 2019

@author: daniel
"""


# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:53:32 2019

@author: daniel
"""
from pathlib import Path
import os

import scipy.stats as sts
import scipy.spatial as sps
import numpy as np


from sklearn.linear_model import LinearRegression
from sklearn.linear_model import MultiTaskElasticNetCV


from parse_MFS_class import MFS_family
from parse_MFS_class import apply_environment
from Env_ball_class import Env_ball
ev =Env_ball(1000)


data_folder = os.path.join(Path(__file__).parents[2], 'data')



def expectation(data, dist=sts.gamma):
    return dist(*dist.fit(data)).mean()

def get_flux_vector(reactions, mdl):
    c=np.zeros(len(reactions))
    for i,name in enumerate(reactions):
        c[i] = mdl.reactions.get_by_id(name).flux
    return c

def one_to_all_distance(inde, matrix, metric='hamming'):
    return sps.distance.cdist(matrix[inde].reshape(1,-1),matrix[np.arange(len(matrix))!=inde], metric=metric)


def distance_to_neighbour(matrix, metric='hamming'):
    v = np.zeros(len(matrix))
    for i in range(len(matrix)):
        v[i] = min(one_to_all_distance(i, matrix, metric)[0])
    return v
        
def get_even_distances(matrix, metric='hamming'):
    '''
    1) pick a random vector and add to S
    2) define the walking distance
    3) iterate:
      a) select a random vector that has distance greater then the walking distance to all vectors in S
      b) add selected vector to set
      c) repeat until not vectors are left

    Parameters
    ----------
    matrix : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    dim = len(matrix)
    selected = [np.random.choice(np.arange(dim))]
    threshold = max(distance_to_neighbour(matrix))/2.
    a =1
    while a==1:
        k = np.array([i for i in np.arange(dim) if i not in selected])
        if len(k)==dim:
            a=0
        else:
            distance_vs = sps.distance.cdist(matrix[selected], matrix[k], metric=metric)
            m = np.ones(len(k))
            
            for i in distance_vs:
                m*=i>threshold
            
            if sum(m)==0:
                a=0
            
            else:
                selected.append(np.random.choice(k[m.astype(np.bool)]))
    print('done')
    return selected
            


def get_odds_ratio(freq_vector):
    
    fv = freq_vector + 0.5
    f = sum(fv)/len(fv)
    return np.log(f/(2-f))
    
def p_of_model(xvar, beta_0, beta_1):
    num = np.exp(beta_0 + beta_1*xvar)
    denom = 1 + num
    return num/denom
        
#get the mfs data        
fam_mfs=MFS_family('Acetobacteraceae', data_folder + '/reactomes/all_families/',data_folder + '/models/all_models' )

#get the reaction frequency data
full_freq_m = fam_mfs.freq_m.copy()
full_freq_m=full_freq_m.T[fam_mfs.include_reactome].T
av_freq_m = np.mean(full_freq_m, axis=0)


Lods_ratio = np.zeros(full_freq_m.shape)

for i in range(len(Lods_ratio)):
    Lods_ratio[i] = np.array([get_odds_ratio(i) for i in fam_mfs.mfs[str(i)]])[fam_mfs.include_reactome]

av_lods = np.mean(Lods_ratio,axis=0)
#get the model reaction frequency

model_sample = np.zeros((1000, len(av_freq_m)))

for i in range(1000):
    print(i)
    s1 = get_even_distances(fam_mfs.model_reactomes)
    mf = np.sum(fam_mfs.model_reactomes[s1], axis=0)/len(s1)
    model_sample[i] = mf[fam_mfs.include_reactome]



#get the environment
     
transporter = ev.transporters[:]
water_index = transporter.index('EX_cpd00001_e')
transporter.remove('EX_cpd00001_e')
oxygen_index  =transporter.index('EX_cpd00007_e')
transporter.remove('EX_cpd00007_e')

transporter=np.array(transporter)

envBall=[]
environment=[]

mc=fam_mfs.model.copy()
for i in range(1000):
    eb = np.delete(ev.matrix[i], water_index, axis=0)
    eb = np.delete(eb, oxygen_index, axis=0)
    d=apply_environment(mc,ev.matrix[i] ,ev.transporters)
    print(d)
    environment.append(get_flux_vector(transporter,mc))
    envBall.append(eb)

envBall = np.array(envBall)
environment = np.array(environment)
environment = -1*environment

environment = np.clip(environment, 0, 1000)


env_sum = np.sum(environment, axis=0)
env_m1 = np.mean(environment, axis=0)

transporter = transporter[env_sum>0]

environment = environment.T[env_sum>0].T




####
mfs_sets =[]


for i in range(1000):
    for z in fam_mfs.mfs[str(i)].T:
        mfs_sets.append(z[fam_mfs.include_reactome].astype(np.int))

mfs_sets = np.array(mfs_sets)



pvalues_less = np.zeros((len(mfs_sets.T), len(environment.T)))
pvalues_greater = np.zeros((len(mfs_sets.T), len(environment.T)))



for reac in range(len(mfs_sets.T)):
    
    a =mfs_sets.T[reac]==1
    
    for met in range(len(environment.T)):
        v1 = np.unique(environment.T[met][np.arange(len(a))[a==1]//1000])
        v2 =  np.unique(environment.T[met][np.arange(len(a))[a==0]//1000])
        
        pvalues_less[reac][met] = sts.mannwhitneyu(v1,v2, alternative='less')[1]
        pvalues_greater[reac][met] = sts.mannwhitneyu(v1,v2, alternative='greater')[1]
        
