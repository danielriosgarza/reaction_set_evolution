#!/usr/bin/env python2
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
    threshold = max(distance_to_neighbour(matrix))/5.
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
            
        
#get the mfs data        
fam_mfs=MFS_family('Acetobacteraceae', data_folder + '/reactomes/all_families/',data_folder + '/models/all_models' )


#get the reaction frequency data
full_freq_m = fam_mfs.freq_m.copy()
full_freq_m=full_freq_m.T[fam_mfs.include_reactome].T
av_freq_m = np.mean(full_freq_m, axis=0)

#get the residuals from the reaction frequency
residualMat = np.zeros(full_freq_m.shape)
for i in range(1000):
    reg = sts.linregress(av_freq_m, full_freq_m[i])
    t = av_freq_m*reg.slope + reg.intercept
    residualMat[i]=(full_freq_m[i] - t)


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




t1 = np.random.choice(np.arange(1000), 700, replace=0)
t2 = np.array([i for i in np.arange(1000) if i not in t1])
########Test1: Elastic Net on reaction frequencies and env ball##############



training_x = full_freq_m[t1]
training_y = envBall[t1]
test_x = full_freq_m[t2]
test_y = envBall[t2]

regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=300)
regr.fit(training_x, training_y)

predicted_envs = regr.predict(test_x)
cors_n=np.zeros(len(t2))
cors=np.zeros(len(t2))
for i in range(len(t2)):
    d=np.zeros(100)
    for z in range(100):
        d[z] = sts.pearsonr(envBall[np.random.randint(0,300)], predicted_envs[i])[0]
    cors_n[i]=np.mean(d)
    cors[i] = sts.pearsonr(envBall[t2[i]], predicted_envs[i])[0]
stats1 = np.arctanh(cors)
stats2 = np.arctanh(cors_n)
print(sts.ttest_rel(stats1,stats2))



#########test2: Elastic Net on residuals and env ball #####################


training_x = residualMat[t1]
training_y = envBall[t1]
test_x = residualMat[t2]
test_y = envBall[t2]

regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=300)
regr.fit(training_x, training_y)

predicted_envs = regr.predict(test_x)
cors_n=np.zeros(len(t2))
cors=np.zeros(len(t2))
for i in range(len(t2)):
    d=np.zeros(100)
    for z in range(100):
        d[z] = sts.pearsonr(envBall[np.random.randint(0,300)], predicted_envs[i])[0]
    cors_n[i]=np.mean(d)
    cors[i] = sts.pearsonr(envBall[t2[i]], predicted_envs[i])[0]
stats1 = np.arctanh(cors)
stats2 = np.arctanh(cors_n)
print('Elastic Net on residuals and env ball')
print(sts.ttest_rel(stats1,stats2))
plot(cors_n, color='b');plot(cors, color='r')
show()

#########test3: Elastic Net on freqMat and consumedmets #####################


training_x = full_freq_m[t1]
training_y = environment[t1]
test_x = full_freq_m[t2]
test_y = environment[t2]

regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=300)
regr.fit(training_x, training_y)

predicted_envs = regr.predict(test_x)
cors_n=np.zeros(len(t2))
cors=np.zeros(len(t2))
for i in range(len(t2)):
    d=np.zeros(100)
    for z in range(100):
        d[z] = sts.pearsonr(environment[np.random.randint(0,300)], predicted_envs[i])[0]
    cors_n[i]=np.mean(d)
    cors[i] = sts.pearsonr(environment[t2[i]], predicted_envs[i])[0]
stats1 = np.arctanh(cors)
stats2 = np.arctanh(cors_n)
print('Elastic Net on freqMat and consumedmets')
print(sts.ttest_rel(stats1,stats2))
plot(cors_n, color='b');plot(cors, color='r')
show()

#########test4: Elastic Net on residual and consumedmets #####################


training_x = residualMat[t1]
training_y = environment[t1]
test_x = residualMat[t2]
test_y = environment[t2]

regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=300)
regr.fit(training_x, training_y)

predicted_envs = regr.predict(test_x)
cors_n=np.zeros(len(t2))
cors=np.zeros(len(t2))
for i in range(len(t2)):
    d=np.zeros(100)
    for z in range(100):
        d[z] = sts.pearsonr(environment[np.random.randint(0,300)], predicted_envs[i])[0]
    cors_n[i]=np.mean(d)
    cors[i] = sts.pearsonr(environment[t2[i]], predicted_envs[i])[0]
stats1 = np.arctanh(cors)
stats2 = np.arctanh(cors_n)
print('Elastic Net on residual and consumedmets')
print(sts.ttest_rel(stats1,stats2))

plot(cors_n, color='b');plot(cors, color='r')
show()





stats = []

for i in np.arange(start=3, stop=150):
    from sklearn.decomposition import PCA
    pca=PCA(n_components=i)

    
    pc_env = pca.fit_transform(environment)
    
    
    t1 = np.random.choice(np.arange(1000), 700, replace=0)
    t2 = np.array([i for i in np.arange(1000) if i not in t1])
    
    
    training_x = residualMat[t1]
    training_y = pc_env[t1]
    
    test_x = residualMat[t2]
    test_y = pc_env[t2]
    
    regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=300)
    regr.fit(training_x, training_y)
    
    
    predicted_envs = regr.predict(test_x)
    cors_n=np.zeros(len(t2))
    cors=np.zeros(len(t2))
    for i in range(len(t2)):
        d=np.zeros(100)
        for z in range(100):
             d[z] = sts.pearsonr(pc_env[np.random.randint(0,300)], predicted_envs[i])[0]
        cors_n[i]=np.mean(d)
        cors[i] = sts.pearsonr(pc_env[t2[i]], predicted_envs[i])[0]
    t1 = np.arctanh(cors)
    t2 = np.arctanh(cors_n)
    stats.append(sts.ttest_rel(t1,t2))
    print(stats)
    

#regr.coef_
#os.system('say "your program has finished"')

pca=PCA(whiten=1, n_components=2)
pc_reac = pca.fit_transform(environment)

t1 = np.random.choice(np.arange(1000), 700, replace=0)
t2 = np.array([i for i in np.arange(1000) if i not in t1])
training_x = pc_reac[t1]
training_y = residualMat[t1]
test_x = pc_reac[t2]
test_y = residualMat[t2]

regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=300)
regr.fit(training_x, training_y)

predicted_envs = regr.predict(test_x)
cors_n=np.zeros(len(t2))
cors=np.zeros(len(t2))
for i in range(len(t2)):
    d=np.zeros(100)
    for z in range(100):
        d[z] = sts.pearsonr(residualMat[np.random.randint(0,300)], predicted_envs[i])[0]
        cors_n[i]=np.mean(d)
        cors[i] = sts.pearsonr(residualMat[t2[i]], predicted_envs[i])[0]
stats1 = np.arctanh(cors)
stats2 = np.arctanh(cors_n)
print(sts.ttest_rel(stats1,stats2))


from sklearn import linear_model
reg = linear_model.BayesianRidge()
t1 = np.random.choice(np.arange(1000), 700, replace=0)
t2 = np.array([i for i in np.arange(1000) if i not in t1])

training_x = residualMat[t1]
training_y = environment[t1]
test_x = residualMat[t2]
test_y = environment[t2]
regr.fit(training_x, training_y)