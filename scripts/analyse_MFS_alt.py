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


def get_residuals(matrix):
    n=len(matrix)
    av_mat = np.round(np.mean(matrix, axis=0)*n)
    mat = np.round(matrix*n)
    r1  = (mat - av_mat)
    return r1
    

def get_residual_scores(matrix):
    
    r1 = get_residuals(matrix)
    score = np.std(r1,axis=0)
    
    stand_r1 = np.zeros(r1.shape)
    
    return score/max(score)
    




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
       
    env_driv_reac_score = get_residual_scores(full_freq_m)
    
    reac_cutoff = np.std(env_driv_reac_score)
    
    env_driven_reactome = reactome[env_driv_reac_score>reac_cutoff]
    reaction_frequency = full_freq_m.T[env_driv_reac_score>reac_cutoff].T

    clss_freq_m = get_residuals(reaction_frequency)

        
    
    
    #for the environment
    
    
       
    env_driv_met_score = get_residual_scores(used_environment)
    
    met_cutoff = np.std(env_driv_met_score)
        
    
    driving_mets = transporter[env_driv_met_score>met_cutoff]
    used_env = used_environment.T[env_driv_met_score>met_cutoff].T

    clss_used_env = get_residuals(used_env)
     
    
       
    
    #regression terms
    x=reaction_frequency.copy()
    y = used_env.copy()
                  
    cosine_dict={}
    
    for i, reac in enumerate(clss_freq_m.T):
        cosine_dict[env_driven_reactome[i]] = np.array([cosine(reac.flatten(), metab.flatten()) for metab in clss_used_env.T])
    
    
    cosine_v = np.array([cosine_dict[i] for i in envd_reactions])
    
  
    #find metabolite concentrations for models
    from sklearn.linear_model import MultiTaskElasticNetCV as EN
    enet  = EN(cv=3,verbose=1, n_jobs=7, max_iter=10000)
    print(x.shape, y.shape)
    mod=enet.fit(x, y)
    evolved_env= np.zeros((len(model_sample), len(dm)))
    
    for i,mod_prof in enumerate(model_sample):
        #print(family, i)
        v = mod_prof[env_driv_reac_score>0]
        
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
