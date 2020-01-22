#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:11:38 2019

@author: danielg
"""
import cobra
import numpy as np
import os

#based on http://compneuro.uwaterloo.ca/files/publications/voelker.2017.pd

def rv(d):
    
    u = np.random.normal(0,1,d+2)  # an array of (d+2) normally distributed random variables
    norm=np.sum(u**2) **(0.5)
    
    u = u/norm
    x = u[0:d] #take the first d coordinates
    x=x+abs(np.min(x))+1.0
    x=x/np.linalg.norm(x)
    return x



def main(simul_n):
    file_list=os.listdir('/media/danielg/4ecaf6c6-8fb4-4128-89b9-4fb4997decbf/generative_microbiome_models/models/vaginal/')
    
    transporters=[]
    
    mod1 = cobra.io.read_sbml_model('/media/danielg/4ecaf6c6-8fb4-4128-89b9-4fb4997decbf/generative_microbiome_models/models/vaginal/'+file_list[0]+'/db/all_'+file_list[0]+'_exch_reactions.sbml')
    for i in file_list:
        reactions_to_add=[]
        all_exc_reactions=cobra.io.read_sbml_model('/media/danielg/4ecaf6c6-8fb4-4128-89b9-4fb4997decbf/generative_microbiome_models/models/vaginal/'+i+'/db/all_'+i+'_exch_reactions.sbml')
        for reac in all_exc_reactions.reactions:
            if not mod1.reactions.has_id(reac.name):
                reactions_to_add.append(reac)
            
            transporters.append(reac.name)
        mod1.add_reactions(reactions_to_add)
        mod1.repair()
    
    transporters=list(set(transporters))    
    water_index=transporters.index('EX_cpd00001_e')
    
    oxygen_index=transporters.index('EX_cpd00007_e')

    
    random_environments=np.zeros((simul_n,len(transporters)))
    
    
    
    for i in xrange(simul_n):
        a=np.random.dirichlet(np.ones(len(transporters)))
        a[water_index]=1.0
        a[oxygen_index]=0.0
        a=a*1000
        random_environments[i]=a
    
    f=file('/media/danielg/4ecaf6c6-8fb4-4128-89b9-4fb4997decbf/generative_microbiome_models/files/vaginal_fams_env_'+str(simul_n)+'.tsv','w')
    
    for i in transporters:
        f.write(i+'\t')
    f.write('\n')
    
    for i in random_environments:
        for num in i:
            f.write('%.20f\t'%num)
        f.write('\n')
    f.close()
    
    cobra.io.write_sbml_model(mod1, '/media/danielg/4ecaf6c6-8fb4-4128-89b9-4fb4997decbf/generative_microbiome_models/files/vaginal_fams_exch_reacs.sbml')
main(1000)
#main(10000)
#main(100000)
