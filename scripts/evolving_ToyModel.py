#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 04:39:24 2020

@author: daniel
"""
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
    
pop_i, pop = moran_process(1000, r, em, 100000,.1,.1,phenotype, react_dict)
f1,f2,f3=[],[],[]
p1 = list(chain.from_iterable(list(pop_i.values())))
p2 = list(chain.from_iterable(list(pop.values())))
p3 = list(chain.from_iterable(list(r.keys())))
for i in reacs:
    f1.append(p1.count(i)/len(pop_i))
for i in reacs:
    f2.append(p2.count(i)/len(pop))
for i in reacs:
    f3.append(p3.count(i)/9)
f1=np.array(f1)
f2=np.array(f2)
f3=np.array(f3)