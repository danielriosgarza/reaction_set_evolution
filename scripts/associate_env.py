#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:53:32 2019

@author: daniel
"""

import scipy.stats as sts
import numpy as np

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import MultiTaskElasticNetCV
from parse_MFS_class import MFS_family
from parse_MFS_class import apply_environment
from Env_ball_class import Env_ball
ev =Env_ball(1000)






def expectation(data, dist=sts.gamma):
    return dist(*dist.fit(data)).mean()

def get_flux_vector(reactions, mdl):
    c=np.zeros(len(reactions))
    for i,name in enumerate(reactions):
        c[i] = mdl.reactions.get_by_id(name).flux
    return c



fam_mfs=MFS_family('Acetobacteraceae', '/home/daniel/generative_models/reactomes/all_families/','/home/daniel/generative_models/models/all_models' )
full_freq_m = fam_mfs.freq_m.copy()
full_freq_m=full_freq_m.T[fam_mfs.include_reactome].T
av_freq_m = np.mean(full_freq_m, axis=0)

residualMat = np.zeros(full_freq_m.shape)
for i in xrange(1000):
    reg = r=sts.linregress(av_freq_m, full_freq_m[i])
    t = av_freq_m*reg.slope + reg.intercept
    residualMat[i]=(full_freq_m[i] - t)
     



transporter = ev.transporters[:]
water_index = transporter.index('EX_cpd00001_e')
transporter.remove('EX_cpd00001_e')
oxygen_index  =transporter.index('EX_cpd00007_e')
transporter.remove('EX_cpd00007_e')

transporter=np.array(transporter)

environment=[]

mc=fam_mfs.model.copy()
for i in xrange(1000):
    apply_environment(mc,ev.matrix[i] ,ev.transporters)
    environment.append(get_flux_vector(transporter,mc))
    
environment = np.array(environment)

for i in xrange(len(environment)):
    environment[i][environment[i]>0] = 0.0
    environment[i] = np.abs(environment[i])


t1 = np.random.choice(np.arange(1000), 700, replace=0)
t2 = np.array([i for i in np.arange(1000) if i not in t1])


training_x = residualMat[t1]
training_y = environment[t1]

test_x = residualMat[t2]
test_y = environment[t2]

regr = MultiTaskElasticNetCV(cv=5, random_state=0,n_jobs=7, verbose=1, max_iter=100)
regr.fit(training_x, training_y)


predicted_envs = regr.predict(test_x)




#regr.coef_
#os.system('say "your program has finished"')
