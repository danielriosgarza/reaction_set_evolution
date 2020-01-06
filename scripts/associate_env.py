#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:53:32 2019

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

#Auxiliary functions

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


def one_to_all_distance(inde, matrix, metric='hamming'):
    '''
    
    Returns
    -------
    distance of one vector in a matrix to all other vectores in the matrix.
    '''
    return sps.distance.cdist(matrix[inde].reshape(1,-1),matrix[np.arange(len(matrix))!=inde], metric=metric)


def distance_to_neighbour(matrix, metric='hamming'):
    '''
    
    Returns
    -------
    distance to closest neigbour
    '''
    v = np.zeros(len(matrix))
    for i in range(len(matrix)):
        v[i] = min(one_to_all_distance(i, matrix, metric)[0])
    return v
        
def get_even_distances(matrix, metric='hamming'):
    '''
    Buld the even_distance matrix (S)
    1) pick a random vector and add to S
    2) define the walking distance
    3) iterate:
      a) select a random vector that has distance greater then the walking distance to all vectors in S
      b) add selected vector to set
      c) repeat until not vectors are left

    Returns
    -------
    indices of vectors in the matrix that consist of a random sample that 
    don't violate the min walking distance (redundancy)

    '''
    dim = len(matrix)
    selected = [np.random.choice(np.arange(dim))]
    #redundancy threshold
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

def get_environment_sample(model, env_vector, env_transporters, reactome, mfs_matrix, interest_transporters, sample_size):
    '''
    
    Returns
    -------
    the mean flux from the import of environmental metabolites from a random sample
    of irreducible sets of reactions. 
    '''
    v = np.zeros((sample_size, len(interest_transporters)))
    mc = model.copy()
    _=apply_environment(mc, env_vector, env_transporters)
    
    for ii,i in enumerate(np.random.choice(np.arange(len(mfs_matrix)), size=sample_size)):
        v[ii] = get_growth_env(model, mc, reactome, mfs_matrix[ii], interest_transporters)
    
    return np.mean(v,axis=0)

def main(family):
    #reference from the git script
    data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')
    
    #obtain all data from the irreducible set
    fam_mfs=MFS_family(family, data_folder + '/reactomes/all_families/',data_folder + '/models/all_models' )
    
    #####get reaction frequency######
    full_freq_m = fam_mfs.freq_m.copy()
    
    #only reactions that should be included in the analysis
    full_freq_m=full_freq_m.T[fam_mfs.include_reactome].T
    av_freq_m = np.mean(full_freq_m, axis=0)
    
    #######get_model_frequency########
    
    #get the model reaction frequency
    
    model_sample = np.zeros((1000, len(av_freq_m)))
    
    for i in range(1000):
        print(i)
        s1 = get_even_distances(fam_mfs.model_reactomes)
        mf = np.sum(fam_mfs.model_reactomes[s1], axis=0)/len(s1)
        model_sample[i] = mf[fam_mfs.include_reactome]
    
    ###get_environment#####
    ev =Env_ball(1000)
         
    transporter = ev.transporters[:]
    #water_index = transporter.index('EX_cpd00001_e')
    transporter.remove('EX_cpd00001_e')
    #oxygen_index  =transporter.index('EX_cpd00007_e')
    transporter.remove('EX_cpd00007_e')
    
    #external metabolites. Water and Oxygen are excluded
    transporter=np.array(transporter)
    
       
    
    
    mc = fam_mfs.model.copy()
    
    used_environment = np.zeros((1000, 290))
    
    for i in range(1000):
        gc.collect()
        v = fam_mfs.mfs[str(i)][fam_mfs.include_reactome].T
        used_environment[i] = get_environment_sample(mc, ev.matrix[i], ev.transporters, fam_mfs.reactome[fam_mfs.include_reactome], v, transporter,200)
        print(i)  
    
    store = {'used_env':used_environment.copy(), 'model_sample':model_sample.copy(), 'full_freq_m':full_freq_m.copy(), 'reactome':fam_mfs.reactome[fam_mfs.include_reactome], 'transporter':transporter.copy()}
    pickle.dump(store, open(data_folder + '/pickles/' + family + '.pkl', 'wb'))

import sys
main(sys.argv[1])