#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:08:06 2020

@author: daniel
"""

from pathlib import Path
import os
import pickle
import gc

import scipy.stats as sts
import scipy.spatial as sps
import numpy as np


from parse_MFS_class import MFS_family
from parse_MFS_class import apply_environment
from Env_ball_class import Env_ball



data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')
    
    #obtain all data from the irreducible set
fam_mfs=MFS_family('Acetobacteraceae', data_folder + '/reactomes/all_families/',data_folder + '/models/all_models' )
    #load a pickle generated from "associate_env.py script"
store = pickle.load(open(data_folder + '/pickles/' + 'Acetobacteraceae' + '.pkl', 'rb'))



####pred data####
full_freq_mat = fam_mfs.freq_m.T[fam_mfs.include_reactome].T.copy()
full_freq_m = np.mean(full_freq_mat, axis=0)
freq_pool = full_freq_mat.flatten()
av_freq_pool = np.tile(full_freq_m, (1000,1)).flatten()
freq_resid_pool = freq_pool-av_freq_pool
scatter(av_freq_pool, freq_pool, s=5, c=np.abs(freq_resid_pool), cmap=cm.Blues, vmin=0, vmax=0.3)

used_mets = store['used_env'].copy()
used_mets_m = np.mean(used_mets, axis=0)
met_pool = used_mets.flatten()
av_met_pool = np.tile(used_mets_m, (1000,1)).flatten()
met_resid_pool = met_pool-av_met_pool
scatter(av_met_pool, met_pool, s=5, c=np.abs(met_resid_pool), cmap=cm.Reds, vmin=0, vmax=0.3)








#####convergence####

irred_sets_0 = fam_mfs.mfs['0'][fam_mfs.include_reactome].T.copy()

irred_sets_0_f = np.sum(irred_sets_0, axis=0)/len(irred_sets_0)
irred_sets_0_residual = irred_sets_0_f - full_freq_m



m_cors=[]



for it in range(2, 400, 5):
    cors=[]
    for i in range(100):
        ch = np.random.choice(np.arange(1000), size=it)
        ft = np.sum(irred_sets_0[ch], axis=0)/len(ch)
        cors.append(1-sts.pearsonr(ft, irred_sets_0_f)[0])
    m_cors.append(np.array(cors))

resid_cors = []

for it in range(2, 1000, 10):
    cors=[]
    for i in range(100):
        ch = np.random.choice(np.arange(1000), size=it)
        ft = np.sum(irred_sets_0[ch], axis=0)/len(ch)
        cors.append(1-sts.pearsonr((ft-full_freq_m), irred_sets_0_residual)[0])
    resid_cors.append(np.array(cors))
