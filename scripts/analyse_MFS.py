#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 08:17:28 2019

@author: daniel
"""
from statsmodels.stats.multitest import multipletests
import networkx as nx

def apply_environment(mdl, env_vec,transporters):
    for i in range(len(transporters)):
        mdl.reactions.get_by_id(transporters[i]).lower_bound=-env_vec[i]
        mdl.reactions.get_by_id(transporters[i]).upper_bound=1000.

    sol_m=mdl.optimize()
    return sol_m.objective_value

def assign_to_class(vector, sig_level=0.05, mt=False):
    '''
    for environment_frequency - average
    -1: reaction is significantly less frequent in this environment
    0: reaction is neutral in this environment
    1: reaction is significantly more frequent in this environment
    
    '''
    SE = np.std(vector)#/np.sqrt(len(vector))
    z_norm_deviate =(vector-np.mean(vector))/SE
    
    sign = z_norm_deviate/np.abs(z_norm_deviate)
    
    pvs =sts.norm.sf(np.abs(z_norm_deviate))
    if mt:
        pvs = multipletests(pvs, alpha=sig_level)[1]
    
    pvs[pvs>=sig_level]=2
    pvs[(pvs<sig_level) & (sign==-1)]=3
    pvs[(pvs<sig_level) & (sign==1)]=4
    
    
    pvs[pvs==2]=0
    pvs[pvs==3]=-1
    pvs[pvs==4]=1
    
    return pvs


def get_flux_vector(reactions, mdl):
    c=np.zeros(len(reactions))
    for i,name in enumerate(reactions):
        c[i] = mdl.reactions.get_by_id(name).flux
    return c

def get_growth_env(model_template, model, reactome, mfs_profile, transporters):
    
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
    return env_vec/max(env_vec)

def get_environment_sample(model, env_vector, env_transporters, reactome, mfs_matrix, interest_transporters, sample_size):
    v = np.zeros((sample_size, len(interest_transporters)))
    mc = model.copy()
    _=apply_environment(mc, env_vector, env_transporters)
    
    for ii,i in enumerate(np.random.choice(np.arange(len(mfs_matrix)), size=sample_size)):
        v[ii] = get_growth_env(model, mc, reactome, mfs_matrix[ii], interest_transporters)
    
    return np.mean(v,axis=0)


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
                
        
        


mc = fam_mfs.model.copy()

used_environment = np.zeros((1000, 290))

for i in range(1000):
    v = fam_mfs.mfs[str(i)][fam_mfs.include_reactome].T
    used_environment[i] = get_environment_sample(mc, ev.matrix[i], ev.transporters, fam_mfs.reactome[fam_mfs.include_reactome], v, transporter,200)
    print(i)    

full_freq_m = fam_mfs.freq_m.copy()
full_freq_m=full_freq_m.T[fam_mfs.include_reactome].T
av_freq_m = np.mean(full_freq_m, axis=0)


diff_freq_m = full_freq_m-av_freq_m
#transform to normal dist
diff_freq_m = np.arctanh(diff_freq_m)
clss_freq_m = np.zeros(diff_freq_m.shape)

for i,v in enumerate(diff_freq_m):
    clss_freq_m[i] = assign_to_class(v, mt=True)

av_used_env = np.mean(used_environment, axis=0)

diff_used_env = used_environment-av_used_env

clss_used_env = np.zeros(diff_used_env.shape)

for i,v in enumerate(diff_used_env):
    clss_used_env[i] = assign_to_class(v, mt=True)

#sums 
s_clss_ue = np.sum(clss_used_env, axis=0)
s_clss_fm = np.sum(clss_freq_m, axis=0)

#env_driven_reactome
envd_reactions = fam_mfs.reactome[fam_mfs.include_reactome]
envd_reactions =envd_reactions[s_clss_fm!=0]

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
nx.write_graphml(g, '/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/networks/Acetobacteraceae.graphml')


#from models

mod_diff = model_sample - av_freq_m
av_mod_diff = np.mean(mod_diff, axis=0)

av_mod_diff = np.arctanh(av_mod_diff)


mod_class = assign_to_class(av_mod_diff)

enriched_reactions=[]
depleted_reactions = []
p_met_drivers=[]
n_met_drivers=[]


for i,t in enumerate(mod_class[s_clss_fm!=0]):
    if t==-1:
        depleted_reactions.append(envd_reactions[i])
        
for i,t in enumerate(mod_class[s_clss_fm!=0]):
    if t==1:
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