#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:52:50 2019

@author: daniel
"""
from pylab import *
import os
import numpy as np
import cobra
import networkx as nx
from Env_ball_class import Env_ball

def apply_environment(mdl, env_vec,transporters):
    for i in range(len(transporters)):
        try:
            mdl.reactions.get_by_id(transporters[i]).lower_bound=-env_vec[i]
            mdl.reactions.get_by_id(transporters[i]).upper_bound=1000.
        except KeyError:
            pass
    sol_m=mdl.optimize()
    return sol_m.objective_value

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


class MFS_run:
    def __init__(self, file_path, file_name):
        print(file_path, file_name)
        
        self.file_name=file_name
        self.file_path = file_path
        self.file_ = os.path.join(file_path, file_name)
        
        self.environment=self.__get_environment(file_name)
        self.reactome=None
        self.binary_m=None
        self.freq_t=None
        
        self.__parse_run_file(self.file_)
    
    def __get_environment(self, file_name):
        return file_name.split(".")[0].split('_')[-1]
    
    def __parse_run_file(self, file_):
        with open(file_) as f:
        
            reactome = np.array(f.readline().strip().split('\t'))
            
            dt=[]
            sorter = np.argsort(reactome)
            
            self.reactome = reactome[sorter]
            
            for line in f:
                a=line.strip().split('\t')
                d = [*map(float, a)]
                d=np.array(d)
                dt.append(d[sorter])
            self.binary_m = np.array(dt).T
            
            self.freq_t = np.sum(self.binary_m, axis=1)/self.binary_m.shape[1]
    
    


class MFS_family:
    def __init__(self, family_name, file_path, model_path):
        self.family = family_name
        self.file_path = file_path
        self.model_path = model_path
        
        self.model = cobra.io.read_sbml_model(os.path.join(self.model_path, self.family,self.family+'.ensembl.sbml'))
        
        
        self.file_path_fam = os.path.join(self.file_path, self.family)
        self.model_path_fam = os.path.join(self.model_path, self.family, 'gapfilled')
        
        self.files = os.listdir(self.file_path_fam)
        self.model_files= os.listdir(self.model_path_fam)
        
        
        self.environments=None
        self.reactome=None
        self.freq_m=None
        
        self.model_reac_freq=None
        self.gene_counts = None
        
        self.mfs=None
        
        self.include_reactome = None
        
        self.model_reactomes = None
        
        self.__parse_run_files()
        self.__get_model_reac_freq()
        self.__get_include_reactome()
        
        
    
    def __parse_run_files(self):
        mfs1= MFS_run(self.file_path_fam, self.files[0])
        
        
        
        
        
        #sorter = np.argsort(mfs1.reactome)
        self.reactome = mfs1.reactome.copy()
        
        self.mfs={}
        self.mfs[mfs1.environment]=mfs1.binary_m.copy()
        
        
        self.environments = np.zeros(len(self.files))
        self.freq_m = np.zeros((len(self.files), len(self.reactome)))
        
        for i in range(len(self.files)):
            print (self.files[i])
            mfs= MFS_run(self.file_path_fam, self.files[i])
            
            self.environments[i] = float(mfs.environment)
            self.freq_m[i]=mfs.freq_t
            
            
            self.mfs[mfs.environment]=mfs.binary_m.copy()
        
        
        
        env_sorter = np.argsort(self.environments)
        self.environments=self.environments[env_sorter]
        self.freq_m = self.freq_m[env_sorter] 
    
    def __get_model_reac_freq(self):
        '''
        Obtain binary vectors indicating the presence of reactions in specific models.
        Summarize these vectors by their frequency.

        Returns
        -------
        None.

        '''
        models_n = len(self.model_files)
        reactions_n = len(self.reactome)
        
        v =np.zeros((models_n, reactions_n))
        freq= np.zeros(reactions_n)
        gene_counts = np.zeros(reactions_n)
        for i,name in enumerate(self.model_files):
            mod=cobra.io.read_sbml_model(os.path.join(self.model_path_fam, name))
            for r in range(reactions_n):
                if mod.reactions.has_id(self.reactome[r]):
                    freq[r]+=(1.0/models_n)
                    gene_counts[r]+=len(mod.reactions.get_by_id(self.reactome[r]).genes)
                    v[i][r] = 1.0
        self.model_reac_freq=freq
        self.gene_counts = gene_counts
        self.model_reactomes = v
    
    def __get_include_reactome(self):
        '''
        reactions that are compared. The rules are:
            1) Are associated to a gene in any of the models;
            2) Are connected to the biomass reaction in a bipartite graph;
            3) Are not lethal to the model, as defined by flux variability analysis.
        '''
        
        index = np.zeros(len(self.reactome))
        
        fva = cobra.flux_analysis.flux_variability_analysis(self.model)
        
        g = buid_bipartite_graph(self.model)
        
        
        
        for i, name in enumerate(self.reactome):
            if self.gene_counts[i]>0:#condition 1)
                
                #condition 3)    
                if fva.loc[name]['minimum']<0:
                    index[i]=1
                elif fva.loc[name]['maximum']>0:
                    index[i]=1
                
                #condition 2)
                if self.model.reactions.get_by_id(name).reversibility:
                    if nx.has_path(g, source=name+'_f', target='bio1'):
                        index[i]=1
                    elif nx.has_path(g, source=name+'_r', target='bio1'):
                        index[i]=1
                else:
                    if nx.has_path(g, source=name, target='bio1'):
                        index[i]=1
                    
                    
        self.include_reactome = index.astype(np.bool)
                    
 
ev =Env_ball(1000)
        
#families = os.listdir('/home/daniel/generative_models/reactomes/vaginal_families')

#fam_mfs={}

#for i in [families[0]]:
#    fam_mfs[i]=MFS_family(i, '/home/daniel/generative_models/reactomes/vaginal_families/','/home/daniel/generative_models/models/vaginal_models' )
                