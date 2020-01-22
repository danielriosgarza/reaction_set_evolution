#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 08:17:28 2019

@author: daniel
"""
import os
from pathlib import Path
import pickle
from pylab import *
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import networkx as nx
import numpy as np
import scipy.stats as sts
import scipy.spatial as sps


def assign_to_class(vector, sig_level=0.05, mt=False):
    '''
    for environment_frequency - average
    -1: reaction is significantly less frequent in this environment
    0: reaction is neutral in this environment
    1: reaction is significantly more frequent in this environment
    
    '''
    if sum(abs(vector))==0:
        return np.zeros(len(vector))
    SE = np.std(vector)/np.sqrt(len(vector))
    z_norm_deviate =(vector-np.mean(vector))/SE
    
    
    sign = z_norm_deviate/np.abs(z_norm_deviate)
    
    pvs =sts.norm.sf(np.abs(z_norm_deviate))*2
    if mt:
        pvs = multipletests(pvs, alpha=sig_level, method='bonferroni')[1]
    
    pvs[pvs>=sig_level]=2
    pvs[(pvs<sig_level) & (sign==-1)]=3
    pvs[(pvs<sig_level) & (sign==1)]=4
    
    
    pvs[pvs==2]=0
    pvs[pvs==3]=-1
    pvs[pvs==4]=1
    
    return pvs

def get_sig_cutoff(vector, sig_level=0.05, mt=False):
    '''
    for environment_frequency - average
    -1: reaction is significantly less frequent in this environment
    0: reaction is neutral in this environment
    1: reaction is significantly more frequent in this environment
    
    '''
    if sum(abs(vector))==0:
        return np.zeros(len(vector))
    SE = np.std(vector)/np.sqrt(len(vector))
    z_norm_deviate =(vector-np.mean(vector))/SE
    
    
    sign = z_norm_deviate/np.abs(z_norm_deviate)
    
    pvs =sts.norm.sf(np.abs(z_norm_deviate))*2
    if mt:
        pvs = multipletests(pvs, alpha=sig_level, method='bonferroni')[1]
    
    pvs[pvs>=sig_level]=2
    pvs[(pvs<sig_level) & (sign==-1)]=3
    pvs[(pvs<sig_level) & (sign==1)]=4
    
    
    pvs[pvs==2]=0
    pvs[pvs==3]=-1
    pvs[pvs==4]=1
    
    return min(vector[pvs==1]), max(vector[pvs==-1])



def get_cutoffs(cosine_dict, quantile):
    v = np.array(list(cosine_dict.values())).flatten()
    pc = np.quantile(v[v>0], quantile)
    nc = np.quantile(v[v<0], 1-quantile)
    return pc, nc
    
def assign_to_rank(vector, p_cutoff, n_cutoff):
    '''
    
    
    '''
    # SE = np.std(vector)/np.sqrt(len(vector))
    # z_norm_deviate =(vector-np.mean(vector))/SE
    ranks = np.zeros(len(vector))
    ranks[(vector>=p_cutoff)] =1
    ranks[(vector<=n_cutoff)]=-1
    return ranks

def cosine(a,b):
    return np.inner(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))


def sig_cosine(cosine_v, n_perm):
    s=int(.05*len(cosine_v))
    v= np.zeros((n_perm, s))
    sig = np.zeros(len(cosine_v))
    for i in range(n_perm):
        v[i]=np.random.choice(cosine_v,size=s, replace=True)
    v=v.flatten()
    for i in range(len(cosine_v)):
        sig[i] = min(sum(v<cosine_v[i])/len(v), sum(v>cosine_v[i])/len(v))
    return sig
    
def build_association_network(association_d, reacs, mets):
    G=nx.MultiGraph()
    
    for i in reacs:
        G.add_node(i, type='reaction')
    for i in mets:
        G.add_node(i, type='metabolite')
    
    for i,r in enumerate(reacs):
        for z,m in enumerate(mets):
            if association_d[r][z]<0:
                G.add_edge(r,m,type='negative', weight =association_d[r][z] )
            elif association_d[r][z]>0:
                G.add_edge(r,m,type='positive', weight=association_d[r][z])
    return G
                
 
def get_evolved_met_prof(evolved_env, dm, transport):
    d={dm[i]:evolved_env[i] for i in range(len(dm))}
    t=np.zeros(len(transport))
    b=transport.copy()
    b.sort()
    for i,z in enumerate(b):
        if z in dm:
            t[i] = d[z]
    return t

def main(family, quantile_ass=.99):
    data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')
    #load a pickle generated from "associate_env.py script"
    store = pickle.load(open(data_folder + '/pickles/' + family + '.pkl', 'rb'))
    
    used_environment = store['used_env'].copy()
    
    full_freq_m = store['full_freq_m'].copy()
    
    reactome = store['reactome'].copy()
    model_sample = store['model_sample'].copy()
    transporter = store['transporter'].copy()
    
    #replace nan values by the average
    av_used_env  = np.nanmean(used_environment,0)
    inds = np.where(np.isnan(used_environment))
    used_environment[inds] = np.take(av_used_env, inds[1])


    #for reaction frequency    
    av_freq_m = np.mean(full_freq_m, axis=0)
    diff_freq_m = full_freq_m-av_freq_m

    #filter out noise and find reactions that are driven by the environment
    m_diff_freq_m = np.max(np.abs(diff_freq_m), axis=0)
    
    env_driven_reactome = reactome[m_diff_freq_m>.005]
    diff_freq_m_envd = diff_freq_m.T[m_diff_freq_m>.005].T
    reaction_frequency = full_freq_m.T[m_diff_freq_m>.005].T

    clss_freq_m = np.zeros(diff_freq_m_envd.shape)
    for i,v in enumerate(diff_freq_m_envd):
        clss_freq_m[i] = v#assign_to_rank(v, fpc,fnc)
    
    #for the environment
    av_used_env = np.mean(used_environment, axis=0)
    
    diff_used_env = used_environment-av_used_env
    
    #filter out noise and find metabolites that are driven by the environment
    m_diff_used_env = np.max(np.abs(diff_used_env), axis=0)
    driving_mets = transporter[m_diff_used_env>0.005]
    diff_used_env_envd = diff_used_env.T[m_diff_used_env>0.005].T
    used_env = used_environment.T[m_diff_used_env>0.005].T
    clss_used_env = np.zeros(diff_used_env_envd.shape)
    
    for i,v in enumerate(diff_used_env_envd):
        clss_used_env[i] = v#assign_to_rank(v, epc, enc)    
     
    
    s_clss_fm = np.sum(np.abs(clss_freq_m), axis=0)
    s_clss_ue = np.sum(np.abs(clss_used_env), axis=0)
    
    #env_driven_reactome
    
    envd_reactions =env_driven_reactome[s_clss_fm!=0]
    
    #driving_metabolites
    dm = driving_mets.copy()
    dm = dm[s_clss_ue!=0]
    
    #profiles
    envd_prof = clss_freq_m.T[s_clss_fm!=0].T
    dm_prof = clss_used_env.T[s_clss_ue!=0].T
    
    #regression terms
    x=reaction_frequency.T[s_clss_fm!=0].T
    y = used_env.T[s_clss_ue!=0].T
                  
    cosine_dict={}
    
    for i, reac in enumerate(envd_prof.T):
        cosine_dict[envd_reactions[i]] = np.array([cosine(reac.flatten(), metab.flatten()) for metab in dm_prof.T])
    
    
    cosine_pool = np.array(list(cosine_dict.values())).flatten()
    pc = np.quantile(cosine_pool[cosine_pool>0], quantile_ass)
    nc = np.quantile(cosine_pool[cosine_pool<0], 1-quantile_ass)
    association_d={}
    
    for i, reac in enumerate(envd_prof.T):
        v = cosine_dict[envd_reactions[i]]
        
        association_d[envd_reactions[i]] = assign_to_rank(v, pc,nc)    
    
    g= build_association_network(association_d, envd_reactions, dm)
    nx.write_graphml(g, os.path.join(Path(os.getcwd()).parents[0], 'files', 'networks', family) +'.graphml')

    #find metabolite concentrations for models
    from sklearn.linear_model import MultiTaskElasticNetCV as EN
    enet  = EN(cv=3,verbose=1, n_jobs=7, max_iter=10000)
    print(x.shape, y.shape)
    mod=enet.fit(x, y)
    evolved_env= np.zeros((len(model_sample), len(dm)))
    
    for i,mod_prof in enumerate(model_sample):
        print(family, i)
        v = mod_prof[m_diff_freq_m>.005]
        
        p = mod.predict(v[s_clss_fm!=0].reshape(1,-1))
        p=p.flatten()
        p = p+abs(min(p))
        p=p/max(p)
        evolved_env[i] =p.copy()
    
    #av_mod_diff = np.arctanh(av_mod_diff)
    met_prof = get_evolved_met_prof(evolved_env, dm, transporter)
        
    return transporter, met_prof



# evs ={}
# path = '/home/daniel/studies/generative_models/git_rep/data/pickles'

# pick = os.listdir(path)

# for i in pick:
#     evs[i] = main(i.replace('.pkl',''))
