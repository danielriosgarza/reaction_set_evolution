#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 10:53:12 2020

@author: daniel
"""

from pylab import *
import os
from pathlib import Path
import numpy as np
import pickle
import seaborn as sns
import scipy.stats as sts
import matplotlib.pyplot as plt


data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')
store = pickle.load(open(data_folder + '/pickles/toy_model.pkl', 'rb'))


####labels###
reactions = np.array(['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15'])
mets = np.array(['M1_e','M2_e','M3_e','M4_e','M5_e','M6_e','M7_e','M8_e','M9_e','M10_e'])

####sorters: reaction frequency and metabolite usage ####

reaction_freq = store['full_freq_m']
s_reaction_freq = np.sum(reaction_freq, axis=1)
s0_reaction_freq = np.sum(reaction_freq, axis=0)
sorter_reac = np.argsort(s_reaction_freq)
sorter_reac = sorter_reac[::-1]
sorter_reac0 = np.argsort(s0_reaction_freq)

used_mets = store['used_environment']
s_mets = np.sum(used_mets, axis=1)
s0_mets = np.sum(used_mets, axis=0)
sorter_mets = np.argsort(s_mets)
sorter_mets = sorter_mets[::-1]
sorter_mets0 = np.argsort(s0_mets)





#####residual metabolite usage #########

residual_met_u = store['diff_used_env']
sns.heatmap(residual_met_u[sorter_mets].T[sorter_mets0].T, cmap=cm.Reds); 


xticks(np.arange(0.5, len(mets)+0.5), reactions[sorter_mets0], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(len(mets)):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/residual_met_u.svg', dpi=600)
show()


######residual reaction frequencies###

residual_reac_f = store['diff_freq_m']
sns.heatmap(residual_reac_f[sorter_mets].T[sorter_reac0].T, cmap=cm.Blues); 


xticks(np.arange(0.5, len(reactions)+0.5), reactions[sorter_reac0], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(15):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/residual_reaction_f.svg', dpi=600)
show()


#####residual metabolite usage #########
residual_met_u = store['diff_used_env']
sns.heatmap(residual_met_u[sorter_mets].T[sorter_mets0].T, cmap=cm.Reds); 


xticks(np.arange(0.5, len(mets)+0.5), mets[sorter_mets0], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(len(mets)):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/residual_met_u.svg', dpi=600)
show()





####Association Matrix #####

ass_d = store['cosine_dict']

u_reacts = reactions[store['reac_filter1']]
u_reacts = u_reacts[store['reac_filter2']]



u_mets = mets[store['met_filter1']]
u_mets = u_mets[store['met_filter2']]

sorter_u_mets=np.zeros(len(u_mets), np.int)

c=-1
for i,v in enumerate(mets[sorter_mets0]):
    if v in u_mets:
        c+=1
        sorter_u_mets[c] = np.arange(len(u_mets))[u_mets==v]

sorter_u_reac=np.zeros(len(u_reacts), np.int)

c=-1
for i,v in enumerate(reactions[sorter_reac0]):
    if v in u_reacts:
        c+=1
        sorter_u_reac[c] = np.arange(len(u_reacts))[u_reacts==v]



ass_m = np.array([ass_d[m] for m in store['envd_reactions']])


sns.heatmap(ass_m[sorter_u_reac].T[sorter_u_mets].T, cmap=cm.Purples, linewidths=0.1, linecolor='w')

xticks(np.arange(0.5, len(u_mets)+0.5), u_mets[sorter_u_mets], rotation=85, fontsize=12)
yticks(np.arange(0.5, len(u_reacts)+0.5), u_reacts[sorter_u_reac], rotation=0, fontsize=12)


savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/association_m.svg', dpi=600)
show()





####reaction frequencies #####

reac_f = store['full_freq_m'].copy()

reac_f  = reac_f.T[store['reac_filter1']].T
reac_f  = reac_f.T[store['reac_filter2']].T


sns.heatmap(reac_f[sorter_mets].T[sorter_u_reac].T, cmap=cm.Blues); 


xticks(np.arange(0.5, len(u_reacts)+0.5), u_reacts[sorter_u_reac], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(len(u_reacts)):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/reaction_f.svg', dpi=600)
show()


#####metabolite usage #########
met_usage = store['used_environment'].copy()
met_usage= met_usage.T[store['met_filter1']].T
met_usage= met_usage.T[store['met_filter2']].T

sns.heatmap(met_usage[sorter_mets].T[sorter_u_mets].T, cmap=cm.Reds)


xticks(np.arange(0.5, len(u_mets)+0.5), u_mets[sorter_u_mets], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(len(u_mets)):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/met_u.svg', dpi=600)
show()


##Evolved_freq##
evolved_freq = np.array(store['evolved_freq'])

sns.heatmap(evolved_freq[sorter_mets].T[sorter_reac0].T, cmap=cm.Blues); 


xticks(np.arange(0.5, len(reactions)+0.5), reactions[sorter_reac0], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(15):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/evolved_react_freq.svg', dpi=600)
show()

#####evolved metabolite usage #########
met_usage = store['true_used_env'].copy()

sns.heatmap(met_usage[sorter_mets].T[sorter_u_mets].T, cmap=cm.Purples)


xticks(np.arange(0.5, len(u_mets)+0.5), u_mets[sorter_u_mets], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(len(u_mets)):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/evolved_met_u.svg', dpi=600)
show()

#####predicted metabolite usage #########
met_usage = store['predicted_environment'].copy()

sns.heatmap(met_usage[sorter_mets].T[sorter_u_mets].T, cmap=cm.Purples)


xticks(np.arange(0.5, len(u_mets)+0.5), u_mets[sorter_u_mets], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(len(u_mets)):
    vlines(i,0,208,color='w', lw =.5)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/predicted_met_u.svg', dpi=600)
show()

####correlation of predicted env
plt.style.use('ggplot')
cors=[]
for i in range(208):
    cors.append(sts.pearsonr(evolved_freq[i], store['full_freq_m'][i])[0]**2)
cors=np.array(cors)
hist(cors**2, 25, density=1, color='r', alpha=0.5)
hist(cors**2, 25, density=1, histtype='step', lw=2, color='r')
xticks([0.2, 0.4, 0.6,0.8, 1.0],['0.2', '0.4', '0.6','0.8', '1.0'], fontsize=12)
savefig('/home/daniel/studies/generative_models/git_rep/reaction_set_evolution/files/Figures/Production_Figures/evolved_correlation.svg', dpi=600)
show()