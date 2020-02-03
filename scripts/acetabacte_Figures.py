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

def get_flux_vector(reactions, mdl):
    '''
    

    Parameters
    ----------
    reactions : TYPE list or array
        ordered list of reactions to obtain the flux from a solved model.
    mdl : cobra.model
    function assumes that the model is solved. (use model.optimize())

    Returns
    -------
    flux_vector : TYPE
        A vector containing the flux in the solved model. In the same order as the reaction 
        list.

    '''
    flux_vector=np.zeros(len(reactions))
    for i,name in enumerate(reactions):
        flux_vector[i] = mdl.reactions.get_by_id(name).flux
    return flux_vector

def get_growth_env(model_template, model, reactome, mfs_profile, transporters):
    '''
    obtain the flux of metabolites used during growth of an irreducible set
    on a specific environmte. Inputs are the ensemble model and binary vector
    containig reactions that are part of the irreducible set.
    
    Returns
    -------
    numpy array containing the flux in the same order as the 'transporters'.
    This arry is divided by its max.
    '''
    
    #model_template is added to retain information about reversibility
    
    #turn reactions on or off the reactions
    
    
    for i,name in enumerate(reactome):
        if mfs_profile[i]==0:
            model.reactions.get_by_id(name).upper_bound=0
            model.reactions.get_by_id(name).lower_bound=0
        else:
            if model_template.reactions.get_by_id(name).reversibility:
                model.reactions.get_by_id(name).upper_bound=1000
                model.reactions.get_by_id(name).lower_bound=-1000
            else:
                model.reactions.get_by_id(name).upper_bound=1000
                model.reactions.get_by_id(name).lower_bound=0
    #
    sol=model.slim_optimize()
    
    env_vec = get_flux_vector(transporters, model)
    env_vec = -1*env_vec
    env_vec = np.clip(env_vec, 0,1000)
    env_vec = env_vec/sol
    if max(env_vec)==0:
        return env_vec
    else:
        return env_vec/max(env_vec)
    
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


#########summary table#######

r_ptwy = {}

with open('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/reaction_pathway.tsv') as f:
    






#####convergence reactions####

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


####convergence metabolites ####
mc= fam_mfs.model.copy()
ge = np.zeros((1000,290))


for i in range(1000):
    print(i)
    mc2 = mc.copy()
    ge[i] = get_growth_env(mc, mc2, store['reactome'], irred_sets_0[i], store['transporter'])
    
m_cors=[]

av_ge=np.sum(ge,axis=0)

for it in range(2, 400, 5):
    cors=[]
    for i in range(100):
        ch = np.random.choice(np.arange(1000), size=it)
        ft = np.sum(ge[ch], axis=0)/len(ch)
        cors.append(1-sts.pearsonr(ft, av_ge)[0])
    m_cors.append(np.array(cors))


###comparison to model###

env_d_score1 = np.round(np.max(full_freq_mat-full_freq_m, axis=0),4)
env_d_score1 = env_d_score1/max(np.abs(env_d_score1))
env_d_score2 = np.round(np.min(full_freq_mat-full_freq_m, axis=0),4)
env_d_score2 = env_d_score2/max(np.abs(env_d_score2))
env_d_score=np.zeros(len(env_d_score1))
for i in range(len(env_d_score1)):
    if abs(env_d_score2[i])>abs(env_d_score1[i]):
        env_d_score[i] = env_d_score2[i]
    else:
        env_d_score[i] = env_d_score1[i]
    
for i in range(len(full_freq_r)):
    scatter(full_freq_m, full_freq_r[i], c='r', s=5)

scatter(full_freq_m[np.abs(env_d_score)<.9],env_d_score[np.abs(r_ma)>.9],s=10,c='r')
scatter(full_freq_m[np.abs(r_ma)>.9], r_ma[np.abs(r_ma)>.9], c='b', s=10, alpha=.5)
