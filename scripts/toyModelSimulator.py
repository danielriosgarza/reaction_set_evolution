#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 10:09:58 2019

@author: daniel
"""
from pylab import *
import numpy as np
import networkx as nx
import scipy.stats as sts
import scipy.spatial as sps

from statsmodels.stats.multitest import multipletests
import networkx as nx

import matplotlib.pyplot as plt
plt.style.use('ggplot')

from itertools import chain, combinations



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
    
    for z, t in enumerate(env_d):
        l=[]
        count=0
        for s in env_d[t]:
            l+=list(s)
            count+=1
        
        for i,name in enumerate(env_vec):
            k[i] = (1./len(env_d)) * l.count(name)
    return k/max(k)    


def cosine(a,b):
    return np.inner(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))


def assign_to_class(vector, sig_level=0.05, mt=False):
    '''
    for environment_frequency - average
    -1: reaction is significantly less frequent in this environment
    0: reaction is neutral in this environment
    1: reaction is significantly more frequent in this environment
    
    '''
    SE = np.std(vector)/np.sqrt(len(vector))
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
mets  =['m'+autoIncrement() for i in range(18)]

ex_mets=['EX_' + mets[i] for i in range(10)]          

rec=0
react_ids = ['r' +autoIncrement() for i in range(12)]

r1 = Reaction(react_ids[0], [ex_mets[0]], [mets[9], mets[10]])
r2 = Reaction(react_ids[1], [ex_mets[1]], [mets[10], mets[11]])
r3 = Reaction(react_ids[2], [ex_mets[2]], [mets[11], mets[12]])
r4 = Reaction(react_ids[3], [ex_mets[3]], [mets[12]])
r5 = Reaction(react_ids[4], [ex_mets[4]], [mets[4]])
r6 = Reaction(react_ids[5], [ex_mets[5],ex_mets[6]], [mets[14]])
r7 = Reaction(react_ids[6], [ex_mets[7]], [mets[14]])
r8 = Reaction(react_ids[7], [ex_mets[8]], [mets[14]])
r9 = Reaction(react_ids[8], [mets[10]], [mets[15]])
r10 = Reaction(react_ids[9], [mets[9]], [ex_mets[9]])
r11 = Reaction(react_ids[10], [mets[11],mets[12]], [mets[16]])
r12 = Reaction(react_ids[11], [mets[4],mets[14]], [mets[17]])


reactions ={r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12}
  
react_dict = {i.name:i for i in reactions} 

k=list(react_dict.keys())
k.sort()
   
phenotype={mets[-1], mets[-2], mets[-3]}

m=Model(reactions)

m.is_growing(ex_mets, phenotype)



r,e=get_mfs(reactions, ex_mets, phenotype, 100)

rd={}

emet = np.array(ex_mets)

rf = np.zeros((10000, len(k)))

envs=np.zeros((10000, len(ex_mets)))

for i in range(10000):
    j = np.round(np.random.uniform(size=10))
    env = emet[j.astype(np.bool)]
    if m.is_growing(env,phenotype):
        
        mfs,environment = get_mfs(reactions, env, phenotype, 100)
        envs[i]=get_average_env(emet, environment)
        rf[i] = get_freq(k, mfs)

    
s_rf = np.sum(rf, axis=1)

rf = rf[s_rf!=0]
envs= envs[s_rf!=0]

av_rf = np.mean(rf, axis=0)

reactome = np.array(k)
e_driv_r = reactome[av_rf!=1]

e_rf = rf.T[av_rf!=1].T

diff_e_rf = e_rf-av_rf[av_rf!=1]

diff_e_rf  = np.arcsinh(diff_e_rf )

clss_e_rf=np.zeros(diff_e_rf.shape)
for i,t in enumerate(diff_e_rf):
    clss_e_rf[i] = assign_to_class(t, mt=True)




av_envs=np.mean(envs, axis=0)
driv_mets = emet[(av_envs!=0) & (av_envs<0.99)]
envs_driv = envs.T[(av_envs!=0) & (av_envs<0.99)].T
diff_envs = envs_driv-av_envs[(av_envs!=0) & (av_envs<0.99)]
diff_envs=np.arcsinh(diff_envs)
clss_envs = np.zeros(diff_envs.shape)
for i,t in enumerate(diff_envs):
    clss_envs[i]= assign_to_class(t, mt=True)




cosine_d={}
for i, react in enumerate(e_driv_r):
    cosine_d[react] = np.array([cosine(clss_e_rf.T[i],clss_envs.T[z]) for z in range(len(driv_mets))])
    cosine_d[react] = np.nan_to_num(cosine_d[react], nan=np.nanmean(cosine_d[react]))
    cosine_d[react] =np.arctanh(cosine_d[react])
#for i, react in enumerate(e_driv_r):
#    cosine_d[react] =np.arctanh(cosine_d[react]-.001)
    
association_d = {i: assign_to_class(cosine_d[i],mt=True) for i in cosine_d}
g=build_association_network(association_d, e_driv_r, driv_mets)
nx.write_graphml(g, '/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/networks/toy_model_associations.graphml')
