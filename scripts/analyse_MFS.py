#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 08:17:28 2019

@author: daniel
"""
import os
from pathlib import Path
import pickle
from statsmodels.stats.multitest import multipletests
import networkx as nx
import numpy as np
import scipy.stats as sts
import scipy.spatial as sps

def assign_to_class(vector, sig_level=1e-150, mt=False):
    '''
    for environment_frequency - average
    -1: reaction is significantly less frequent in this environment
    0: reaction is neutral in this environment
    1: reaction is significantly more frequent in this environment
    
    '''
    SE = np.std(vector)/np.sqrt(len(vector))
    z_norm_deviate =(vector-np.mean(vector))/SE
    
    sign = z_norm_deviate/np.abs(z_norm_deviate)
    
    pvs =sts.norm.sf(np.abs(z_norm_deviate))*2
    if mt:
        pvs = multipletests(pvs, alpha=sig_level, method='Bonferroni')[1]
    
    pvs[pvs>=sig_level]=2
    pvs[(pvs<sig_level) & (sign==-1)]=3
    pvs[(pvs<sig_level) & (sign==1)]=4
    
    
    pvs[pvs==2]=0
    pvs[pvs==3]=-1
    pvs[pvs==4]=1
    
    return pvs


def cosine(a,b):
    return np.inner(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))



def build_association_network(association_d, reacs, mets):
    G=nx.MultiGraph()
    
    for i in reacs:
        G.add_node(i, type='reaction')
    for i in mets:
        G.add_node(i, type='metabolite')
    
    for i,r in enumerate(reacs):
        for z,m in enumerate(mets):
            if association_d[r][z]==-1:
                G.add_edge(r,m,type='negative')
            elif association_d[r][z]==1:
                G.add_edge(r,m,type='positive')
    return G
                
        

def main(family):
    data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')
    store = pickle.load(open(data_folder + '/pickles/' + family + '.pkl', 'rb'))
    
    used_environment = store['used_env'].copy()
    full_freq_m = store['full_freq_m'].copy()
    reactome = store['reactome'].copy()
    model_sample = store['model_sample'].copy()
    transporter = store['transporter'].copy()
    
    
    av_used_env  = np.nanmean(used_environment,0)
    inds = np.where(np.isnan(used_environment))
    used_environment[inds] = np.take(av_used_env, inds[1])


    
    
    av_freq_m = np.mean(full_freq_m, axis=0)
    
    diff_freq_m = full_freq_m-av_freq_m
    #diff_freq_m = np.arctanh(diff_freq_m)
    clss_freq_m = np.zeros(diff_freq_m.shape)
    
    for i,v in enumerate(diff_freq_m):
        clss_freq_m[i] = assign_to_class(v, mt=True)
    
    av_used_env = np.mean(used_environment, axis=0)
    
    diff_used_env = used_environment-av_used_env
    #diff_used_env = np.arctanh(diff_used_env)
    clss_used_env = np.zeros(diff_used_env.shape)
    
    for i,v in enumerate(diff_used_env):
        clss_used_env[i] = assign_to_class(v, mt=True)
        
    s_clss_fm = np.sum(np.abs(clss_freq_m), axis=0)
    s_clss_ue = np.sum(np.abs(clss_used_env), axis=0)
    
    #env_driven_reactome
    
    envd_reactions =reactome[s_clss_fm!=0]
    
    #driving_metabolites
    dm = transporter.copy()
    dm = dm[s_clss_ue!=0]
    
    #profiles
    envd_prof = clss_freq_m.T[s_clss_fm!=0].T
    dm_prof = clss_used_env.T[s_clss_ue!=0].T
    
    cosine_dict={}
    
    for i, reac in enumerate(envd_prof.T):
        cosine_dict[envd_reactions[i]] = np.array([np.arctanh(cosine(reac.flatten(), metab.flatten())) for metab in dm_prof.T])
    
    association_d={}
    for i, reac in enumerate(envd_prof.T):
        association_d[envd_reactions[i]] = assign_to_class(cosine_dict[envd_reactions[i]], mt=True)
    
    g= build_association_network(association_d, envd_reactions, dm)
    nx.write_graphml(g, os.path.join(Path(os.getcwd()).parents[0], 'files', 'networks', family) +'.graphml')

    mod_diff = model_sample - av_freq_m
    av_mod_diff = np.mean(mod_diff, axis=0)
    #av_mod_diff = np.arctanh(av_mod_diff)
    
    mod_class = assign_to_class(av_mod_diff)
    enriched_reactions=[]
    depleted_reactions = []
    p_met_drivers=[]
    n_met_drivers=[]
    
    for i,v in enumerate(mod_class[s_clss_fm!=0]):
        if v==-1:
            depleted_reactions.append(envd_reactions[i])
    for i,v in enumerate(mod_class[s_clss_fm!=0]):
        
        if v==1:
            enriched_reactions.append(envd_reactions[i])
            
            
    for i in enriched_reactions:
        k = association_d[i].copy()
        l1 = list(dm[k==1])
        l2 = list(dm[k==-1])
        p_met_drivers+=l1
        n_met_drivers+=l2
        
    for i in depleted_reactions:
        k = association_d[i].copy()
        l1 = list(dm[k==1])
        l2 = list(dm[k==-1])
        n_met_drivers+=l1
        p_met_drivers+=l2
    p_met_drivers=list(set(p_met_drivers))
    p_met_drivers.sort()
    n_met_drivers=list(set(n_met_drivers))
    n_met_drivers.sort()
    return p_met_drivers, n_met_drivers, enriched_reactions,depleted_reactions,association_d, cosine_dict,mod_class, s_clss_fm