#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:05:44 2019

@author: daniel
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import cobra
import numpy as np
from Env_ball_class import Env_ball


family=sys.argv[1]
job_state= int(sys.argv[2])


def apply_environment(mdl, env_vec,transporters):
    for i in xrange(len(transporters)):
        mdl.reactions.get_by_id(transporters[i]).lower_bound=-env_vec[i]
        mdl.reactions.get_by_id(transporters[i]).upper_bound=1000.
        
    sol_m=mdl.optimize()
    return sol_m.objective_value



def remove_reaction(orig_obj, model, reaction):
    up_orig = float(model.reactions.get_by_id(reaction).upper_bound)
    lo_orig = float(model.reactions.get_by_id(reaction).lower_bound)
    
    model.reactions.get_by_id(reaction).upper_bound=0
    model.reactions.get_by_id(reaction).lower_bound=0
    
    
    
    new_obj=np.round(float(model.slim_optimize()),decimals=12)
    
    
    if new_obj>orig_obj:
        return 0.0
    else:
        model.reactions.get_by_id(reaction).upper_bound=up_orig
        model.reactions.get_by_id(reaction).lower_bound=lo_orig
        return 1.0
    

def get_mfs(map_d, ordered_reactions, model, cuttof):
    mod= model.copy()
    
    essential_r = np.zeros(len(map_d))
    
    orig_obj= np.round(float(cuttof*model.slim_optimize()),decimals=12)
    
    for i in ordered_reactions:
        
        c=remove_reaction(orig_obj, mod,i)
        
        essential_r[map_d[i]] = c
    
    return essential_r

def get_mfs_dist(model, reactions, transporters, env, max_it=1000):
    
    rc=reactions[:]
    map_d = {reactions[i]: i for i in xrange(len(reactions))}
    
    apply_environment(model, env,transporters)
    
    all_iters=[]
    
    for i in xrange(max_it):
        np.random.shuffle(rc)
        mfs=get_mfs(map_d, rc, model, 0.01)
        all_iters.append(mfs)
    return all_iters
        


def get_results(reactions, all_iters, file_name):
    
    f =file(file_name, 'w')
    f.write('\t'.join(reactions))
    f.write('\n')
    for a in all_iters:
        f.write('\t'.join(a.astype(str)))
        f.write('\n')
    f.close()
    

def get_run(model, reactions, transporters, env, output_file, max_it=1000,seed=None):
    np.random.seed(seed)
    all_iters= get_mfs_dist(model, reactions, transporters, env, max_it)
    get_results(reactions, all_iters, output_file)


############## USE ##########

#eb=Env_ball(1000)
  
#transporters=eb.transporters[:]
#random_environments=eb.matrix.copy()

#model = cobra.io.read_sbml_model('/home/danielg/generative_models/plinc_runs/generative_models/vaginal_models/' + family + '/' + family + '.ensemble.sbml')

#model.solver='gurobi'
#reactions=[i.id for i in model.reactions if 'rxn' in i.id]

#end_env=job_state*20
#start_env=end_env-20



#[get_run(model, reactions, transporters, random_environments[i], '/home/danielg/generative_models/plinc_runs/generative_models/random_model/'+family+'/' +family+'_randiter_'+str(i)+'.tsv', 1000, job_state*3) for i in xrange(start_env, end_env,1)]

