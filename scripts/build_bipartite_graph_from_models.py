#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:50:49 2019

@author: daniel
"""

import cobra
import numpy as np
import networkx as nx



def buid_bipartite_graph(model):
    G=nx.MultiDiGraph()
    
    for met in model.metabolites:
        G.add_node(met.id, name=met.name, nodeType = 'metabolite')
    
    
    for reac in model.reactions:
        if not reac.reversibility:
            G.add_node(reac.id, name=reac.name, nodeType = 'reaction')
            
            for met in reac.reactants:
                G.add_edge(met.id, reac.id, name = reac.id)
            for met in reac.products:
                G.add_edge(reac.id, met.id, name = reac.id)
        
        else:
            G.add_node(reac.id + '_r', name=reac.name, nodeType = 'reaction')
            G.add_node(reac.id + '_f', name=reac.name, nodeType = 'reaction')
            
            for met in reac.reactants:
                G.add_edge(met.id, reac.id + '_f', name = reac.id)
                G.add_edge(reac.id + '_r',met.id, name = reac.id)
            for met in reac.products:
                G.add_edge(reac.id + '_f', met.id, name = reac.id)
                G.add_edge(met.id,reac.id + '_r', name = reac.id)
    
    return G



model = cobra.io.read_sbml_model('/home/daniel/generative_models/models/all_models/Acetobacteraceae/Acetobacteraceae.ensembl.sbml')
g = buid_bipartite_graph(model)
reactions = [x for x,y in g.nodes(data=True) if y['nodeType']=='reaction']