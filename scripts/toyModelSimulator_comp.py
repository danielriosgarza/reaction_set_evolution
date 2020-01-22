#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 10:09:58 2019

@author: daniel
"""
from pylab import *
import pickle
import os
from pathlib import Path
import numpy as np
import networkx as nx
import scipy.stats as sts
import scipy.spatial as sps

from statsmodels.stats.multitest import multipletests
import networkx as nx

import matplotlib.pyplot as plt
plt.style.use('ggplot')

from itertools import chain, combinations
from itertools import chain
from itertools import groupby


def replication_w_error(p_in, p_del, env, model_reactions,phenotype, r_dict, pool):
    process = np.random.choice(['insert','delete', 'remain' ], p=np.array([p_in, p_del, 1-(p_in+p_del)]))
    if process =='insert':
        reaction = np.random.choice(pool)
        mod = model_reactions.copy()
        mod.add(reaction)
        reacs = [r_dict[i] for i in mod]
        m=Model(reacs)
        pc = m.is_growing(env, phenotype)
    elif process=='delete':
        reaction = np.random.choice(list(model_reactions))
        mod = model_reactions.copy()
        mod.discard(reaction)
        reacs = [r_dict[i] for i in mod]
        m=Model(reacs)
        pc = m.is_growing(env, phenotype)
    else:
        mod = model_reactions.copy()
        reacs = [r_dict[i] for i in mod]
        m=Model(reacs)
        pc = m.is_growing(env, phenotype)
    
    return mod,pc
            

def moran_process(population, mfs, env, iter_n, p_in,p_del, phenotype, r_dict):
    mfs = list(mfs.keys())
    pop = {}
    for i in range(population):
        pop[i] = set(np.random.choice(mfs))
    pop_i = pop.copy()
    for i in range(iter_n):
        
        pool = list(chain.from_iterable(list(pop.values())))
        selected = np.random.choice(list(pop.keys()), size=2)
        m,pc = replication_w_error(p_in, p_del, env, pop[selected[1]],phenotype, r_dict, pool)
        if pc:
            pop[selected[0]] = m
    
    return pop_i, pop












class Reaction:
    '''
    reaction object. Reactions have a name, reactants, and products.
    Can be reversible. Stoichiometry is one-to-one in the toy model.
    '''
    
    def __init__(self, name, reactants, products, reversibility=False):
        self.name = name
        self.reactants = set(reactants)
        self.products = set(products)
        self.reversibility = reversibility
    
    def is_active(self, metabolite_pool):
        '''
        given a metabolite pool, will the reaction be used.

        Parameters
        ----------
        metabolite_pool : TYPE set
            set of metabolites. If the reactants (in any direction)
            are in the set, the reaction is active and produces products.

        Returns
        -------
        active : TYPE bool
            

        '''
        active=True
        mp = set(metabolite_pool)
        
        if len(mp.intersection(self.reactants))<len(self.reactants):
            active=False
        if self.reversibility and active==False:
            if len(mp.intersection(self.products))==len(self.products):
                active=True
            
        return active
    def build_string(self):
        '''
        build a string representation for the reaction
        '''
        if self.reversibility:
            symbol=' <--> '
        else:
            symbol=' --> '
        
        return ' + '.join(self.reactants) + symbol + ' + '.join(self.products)
            

class Model:
    '''
    Model object. Model has intracellular and extracellular reactions and metabolites;
    Extracellular metabolites are identified by "EX" in their names.
    '''
    def __init__(self, reactions,name='model'):
        self.name = name
        self.reactions = set(reactions)
        self.ic_products, self.ic_reactants, self.ex_products, self.ex_reactants, self.reactants, self.products, self.m_to_reactions= self.parse_reactions(self.reactions)
    
    def parse_reactions(self, reactions):
        ic_prod, ic_react, ex_prod, ex_react, prod, react = set(),set(),set(),set(),set(),set()
        m_to_reactions={}
        
        for i in reactions:
            for r in i.reactants:
                if r in m_to_reactions:
                    m_to_reactions[r].append(i)
                else:
                    m_to_reactions[r]=[i]
                react.add(r)
                if 'EX' in r:
                    ex_react.add(r)
                else:
                    ic_react.add(r)
            for p in i.products:
                if p in m_to_reactions:
                    m_to_reactions[p].append(i)
                else:
                    m_to_reactions[p]=[i]
                prod.add(p)
                if 'EX' in p:
                    ex_prod.add(p)
                else:
                    ic_prod.add(p)
        return ic_prod, ic_react, ex_prod, ex_react, prod,react, m_to_reactions
    
    
        
    def is_producing_phenotype(self, metabolite_pool, phenotype):
        
        mp = set(metabolite_pool)
        
        allmets = mp.copy()
        
        allmets=allmets.union(self.products)
        
        #print allmets
        active=[]
        
        for i in self.reactions:
            if i.is_active(allmets):
                active.append(i)
        
        
        t=[]
        while set(t)!=set(active):
            
            t=active[:]
            active=[]
            
            allmets = mp.copy()
            for i in t:
                allmets=allmets.union(i.products)
            
            for i in t:
                if i.is_active(allmets):
                    active.append(i)
                
        
        allmets=mp.copy()
        for i in active:
            allmets=allmets.union(i.products)
        
        #print allmets
        #print [i.name for i in active]
        if phenotype==allmets.intersection(phenotype):
            return True
        else:
            return False
        
        
    def is_growing(self, metabolite_pool, phenotype):
        
        mp = set(metabolite_pool)
        
        allmets = mp.copy()
        
        allmets=allmets.union(self.products)
        
        #print allmets
        active=[]
        
        for i in self.reactions:
            if i.is_active(allmets):
                active.append(i)
        
        
        t=[]
        while set(t)!=set(active):
            
            t=active[:]
            active=[]
            
            allmets = mp.copy()
            for i in t:
                allmets=allmets.union(i.products)
            
            for i in t:
                if i.is_active(allmets):
                    active.append(i)
                
        
        allmets=mp.copy()
        for i in active:
            allmets=allmets.union(i.products)
        
        #print allmets
        #print [i.name for i in active]
        if phenotype==allmets.intersection(phenotype):
            #final check!
            prods=set()
            reacts=set()
            for i in active:
                for p in i.products:
                    prods.add(p)
                for r in i.reactants:
                    reacts.add(r)
            prods=prods-phenotype
            prods = prods
            b=prods.difference(reacts.intersection(prods))
            
            if b==self.ex_products.intersection(b):
                return True
            else:
                remove=[]
                
                for i in b:
                    remove+=self.m_to_reactions[i]
                
                r =self.reactions-set(remove)
                
                m=Model(r)
                if m.is_producing_phenotype(metabolite_pool, phenotype):
                    return True
                else:
                    return False
        else:
            return False
            
    
                
        
        

def autoIncrement():
 global rec
 pStart = 1 #adjust start value, if req'd 
 pInterval = 1 #adjust interval value, if req'd
 if (rec == 0): 
  rec = pStart 
 else: 
  rec = rec + pInterval 
 return str(rec).zfill(4)


def get_mfs(reactions, env, phenotype, iterations):
    
    rd = {i.name:i for i in reactions}
    
    reaction_list = list(rd.keys())
    
    rl = reaction_list[:]
    
    mfs={}
    envs={}
    for i in range(iterations):
        np.random.shuffle(rl)
        np.random.shuffle(rl)
        
        rlc=rl[:]
        
        rlt=rlc[:]
        
        for z in range(len(rl)):
            target = rlc[z]
            rlt.remove(target)
            
            r = [rd[w] for w in rlt]
            m = Model(r)
            #print rlt
            if not m.is_growing(env, phenotype):
                rlt.append(target)
        
        
        
        rlt.sort()
        r = [rd[w] for w in rlt]
        m = Model(r)
        
        
        if frozenset(rlt) in mfs:
            mfs[frozenset(rlt)]+=1
            if m.ex_reactants not in envs[frozenset(rlt)]:
                envs[frozenset(rlt)].append(m.ex_reactants)
        else:
            mfs[frozenset(rlt)]=1
            envs[frozenset(rlt)] = [m.ex_reactants]
               
    return mfs,envs
    
def get_freq(reaction_list, mfs):
    f = {i:0 for i in reaction_list}
    
    for i_set in mfs:
        
        for reaction in i_set:
            f[reaction]+=1./len(mfs)
    return np.array([f[i] for i in reaction_list])    
    
def get_average_env(env_vec, env_d):
    k = np.zeros(len(env_vec))
    denom = len(env_d)
    
    for i, t in enumerate(env_d):
        for j,m in enumerate(env_vec):
            if m in env_d[t][0]:
                k[j] += 1/denom
    return k    


def cosine(a,b):
    return np.inner(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))


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

def assign_to_rank(vector, p_cutoff, n_cutoff):
    '''
    
    
    '''
    # SE = np.std(vector)/np.sqrt(len(vector))
    # z_norm_deviate =(vector-np.mean(vector))/SE
    ranks = np.zeros(len(vector))
    ranks[(vector>0) & (vector>=p_cutoff)] =1
    ranks[(vector<0) & (vector<=n_cutoff)]=-1
    return ranks



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


rec=0
mets  =['m'+autoIncrement() for i in range(23)]

ex_mets=['EX_' + mets[i] for i in range(10)]          

rec=0
react_ids = ['r' +autoIncrement() for i in range(15)]

r1 = Reaction(react_ids[0], [ex_mets[0]], [mets[0], mets[1]])
r2 = Reaction(react_ids[1], [ex_mets[1]], [mets[0], mets[2]])
r3 = Reaction(react_ids[2], [ex_mets[2]], [mets[3]])
r4 = Reaction(react_ids[3], [ex_mets[3]], [mets[4]])
r5 = Reaction(react_ids[4], [ex_mets[4]], [mets[5]])
r6 = Reaction(react_ids[5], [ex_mets[5],ex_mets[6]], [mets[6]])
r7 = Reaction(react_ids[6], [ex_mets[7]], [mets[7]])
r8 = Reaction(react_ids[7], [ex_mets[8]], [mets[7]])
r9 = Reaction(react_ids[8], [mets[0]], [mets[8]])
r10 = Reaction(react_ids[9], [mets[1]], [ex_mets[9]])
r11 = Reaction(react_ids[10], [mets[6]], [mets[7]])
r12 = Reaction(react_ids[11], [mets[5],mets[7]], [mets[9]])
r13 = Reaction(react_ids[12], [mets[4]], [mets[10]])
r14 = Reaction(react_ids[13], [mets[2],mets[8]], [mets[11]])
r15 = Reaction(react_ids[14], [mets[3]], [mets[4], mets[2]])





reactions ={r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15}
  
react_dict = {i.name:i for i in reactions} 

k=list(react_dict.keys())
k.sort()
   
phenotype={mets[9], mets[10], mets[11]}

m=Model(reactions)

m.is_growing(ex_mets, phenotype)



r,e=get_mfs(reactions, ex_mets, phenotype, 20000)




rd={}

emet = np.array(ex_mets)

rf = np.zeros((20000, len(k)))

envs=np.zeros((20000, len(ex_mets)))

used=[]
env_that_sup_growth = []

for i in range(20000):
    j = np.round(np.random.uniform(size=10))
    if list(j) not in used:
        used.append(list(j))
        
        env = emet[j.astype(np.bool)]
        if m.is_growing(env,phenotype):
            env_that_sup_growth.append(j)
            
            mfs,environment = get_mfs(reactions, env, phenotype, 1000)
            envs[i]=get_average_env(emet, environment)
            rf[i] = get_freq(k, mfs)










#eliminate the environments with no growth




s_rf = np.sum(rf, axis=1)


quantile_freq=.90
quantile_met=.80
quantile_ass=.95


used_environment = envs[s_rf!=0].copy()

full_freq_m = rf[s_rf!=0].copy()
reactome = np.array(k)
transporter = emet.copy()  



av_freq_m = np.mean(full_freq_m, axis=0)
diff_freq_m = full_freq_m-av_freq_m



#filter out noise and find reactions that are driven by the environment
m_diff_freq_m = np.max(np.abs(diff_freq_m), axis=0)
env_driven_reactome = reactome[m_diff_freq_m>.005]
diff_freq_m_envd = diff_freq_m.T[m_diff_freq_m>.005].T

diff_freq_pool = diff_freq_m_envd.flatten()
clss_freq_m = np.zeros(diff_freq_m_envd.shape)
for i,v in enumerate(diff_freq_m_envd):
    clss_freq_m[i] = v
    


#for the environment
av_used_env = np.mean(used_environment, axis=0)
diff_used_env = used_environment-av_used_env
#filter out noise and find metabolites that are driven by the environment
m_diff_used_env = np.max(np.abs(diff_used_env), axis=0)
driving_mets = transporter[m_diff_used_env>0.005]
diff_used_env_envd = diff_used_env.T[m_diff_used_env>0.005].T
    
diff_used_env_pool = diff_used_env_envd.flatten()
clss_used_env = np.zeros(diff_used_env_envd.shape)
for i,v in enumerate(diff_used_env_envd):
    clss_used_env[i] = v

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
    
                  
cosine_dict={}

for i, reac in enumerate(envd_prof.T):
    cosine_dict[envd_reactions[i]] = np.array([cosine(reac.flatten(), metab.flatten()) for metab in dm_prof.T])

cosine_pool = np.array(list(cosine_dict.values())).flatten()
pc = np.quantile(cosine_pool[cosine_pool>0], quantile_ass)
nc = np.quantile(cosine_pool[cosine_pool<0], 1-quantile_ass)
association_d={}

for i, reac in enumerate(envd_prof.T):
    v = cosine_dict[envd_reactions[i]]
    ncp=sqrt(len(v))*(mean(v**2)-0)/std(v**2)
    alpha=ncp/sqrt(len(v))
    pa = alpha/2
    na = -pa
    association_d[envd_reactions[i]] = assign_to_rank(v, pa,na)
        
g= build_association_network(association_d, envd_reactions, dm)
nx.write_graphml(g, '/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/networks/toy_model_associations.graphml')


reacs = list(react_dict.keys())
predicted = []
average = []
predicted_environments=[]
true_used_env=[]
mfs_avg_env = []
evolved_freq = []    

# evolve conditioned to an environment
for i in range(len(env_that_sup_growth)):
    print(i)
    j=env_that_sup_growth[i]
    env = emet[j.astype(np.bool)]
    mfs,environment = get_mfs(reactions, env, phenotype, 1000)
    envs=get_average_env(emet, environment)
    mfs_avg_env.append(envs)
    pop_i, pop = moran_process(1000, r, env, 1000000,.05,.5,phenotype, react_dict)
    g1, g2 = [],[]
    
    for pp in pop_i:
        h=list(pop_i[pp].copy())
        h.sort()
        if h not in g1:
            g1.append(h)
    
    for pp in pop:
        h=list(pop[pp].copy())
        h.sort()
        if h not in g2:
            g2.append(h)
    
    
    #get the reaction_frequencies of evolved population
    
    f1,f2,f3=[],[],[]
    p1 = list(chain.from_iterable(list(g1)))
    p2 = list(chain.from_iterable(list(g2)))
    p3 = list(chain.from_iterable(list(r.keys())))
    for rr in reactome:
        f1.append(p1.count(rr)/len(g1))
    for rr in reactome:
        f2.append(p2.count(rr)/len(g2))
    for rr in reactome:
        f3.append(p3.count(rr)/9)
    f1=np.array(f1)
    f2=np.array(f2)
    f3=np.array(f3)
    
    
    evolved_freq.append(f2)
    
    #get the average metabolite usage from the evolved population
    used_mets=[]
    for mm in g2:
        reacs=[react_dict[z] for z in mm]
        m=Model(reacs)
        used_mets.append(m.ex_reactants)
    used_mets = list(chain.from_iterable(used_mets))
    mf=[]
    for mm in dm:
        mf.append(used_mets.count(mm)/len(g2))
    true_used_env.append(mf)
    
    from sklearn.linear_model import MultiTaskElasticNetCV as EN
    enet  = EN(cv=50, max_iter=100000)
    x = full_freq_m.T[m_diff_freq_m>.005].T
    y = used_environment.T[m_diff_used_env>0.005].T
    mod=enet.fit(x, y)
    p = mod.predict(f2[m_diff_freq_m>.005].reshape(1,-1))
    
    
    p=p.flatten()
    p = p+abs(min(p))
    p=p/max(p)
    
    c = [sts.pearsonr(mf,used_environment[ee][m_diff_used_env>0.005])[0] for ee in range(len(used_environment))]
    
    
    predicted.append(sts.pearsonr(p, mf)[0])
    
    average.append(mean(c))

    predicted_environments.append(p)

data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')
store = {'full_freq_m': full_freq_m, 'diff_freq_m':diff_freq_m, 'envd_prof':envd_prof, 'envd_reactions': envd_reactions,'used_environment':used_environment, 'ex_mets':ex_mets, 'dm':dm, 'diff_used_env':diff_used_env,         'dm_prof':dm_prof, 'x':x, 'y':y, 'predicted':predicted, 'average':average, 'predicted_environment':predicted_environments,'true_used_env':true_used_env, 'mfs_avg_env': mfs_avg_env, 'evolved_freq':evolved_freq, 'cosine_dict':cosine_dict, 'met_filter1':m_diff_used_env>0.005, 'met_filter2':s_clss_ue!=0, 'reac_filter1':m_diff_freq_m>.005, 'reac_filter2':s_clss_fm!=0, 'evolved_freq':evolved_freq}
pickle.dump(store, open(data_folder + '/pickles/toy_model.pkl', 'wb'))
