#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:47:11 2020

@author: daniel
"""

j=env_that_sup_growth[0]
env = emet[j.astype(np.bool)]
mfs,environment = get_mfs(reactions, env, phenotype, 1000)
envs=get_average_env(emet, environment)
pop_i, pop = moran_process(10000, r, env, 1000000,.1,.1,phenotype, react_dict)
unique_env = np.unique(used_environment,axis=0)


pred = np.zeros(8)

for z,t in enumerate(reactome):
    if t in cosine_dict:
        cd = np.round(cosine_dict[t])
        pred+=(f2[z]-f3[z])*cd
