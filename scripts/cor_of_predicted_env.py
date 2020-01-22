#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:58:16 2020

@author: daniel
"""

quantile_freq=.95
quantile_met=.95
quantile_ass=.99
quantile_mod=.95
family = 'Acetobacteraceae'


test_frq = store['full_freq_m'].copy()
test_env = store['used_env'].copy()
av_test_env  = np.nanmean(used_environment,0)
inds = np.where(np.isnan(test_env))
test_env[inds] = np.take(av_test_env, inds[1])




test_dif_frq = test_frq-av_freq_m
test_dif_env = test_env - av_used_env

tfrq =test_dif_frq.T[m_diff_freq_m>.005].T
tenv = test_dif_env.T[m_diff_used_env>0.005].T


tfrq =tfrq.T[s_clss_fm!=0].T
tenv = tenv.T[s_clss_ue!=0].T
pred = np.zeros(tenv.shape)
pred.shape
    

for i in range(len(tfrq)):
    
    for z,t in enumerate(envd_reactions):
        pred[i]+=tfrq[i][z]*(cosine_dict[t])
        
pcolormesh(pred)

pc1=[]
pc2=[]


x=diff_used_env_envd.T[s_clss_ue!=0].T

for i in range(len(pred)):
    c2=[sts.pearsonr(pred[i], x[z])[0] for z in range(len(x)) if z!=i]
    c3=sts.pearsonr(pred[i], x[i])[0]
    
    pc1.append(mean(c2))
    pc2.append(c3)
    