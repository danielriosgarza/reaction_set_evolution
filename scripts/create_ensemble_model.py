
import os
import cobra
import numpy as np

from cobra import Model, Reaction, Metabolite
from Env_ball_class import Env_ball


class Ensemble_model:
    '''
    Create a model from a pool of reactions from related organisms.
    
    1) Take a group of models
    2) pool all of their "rxn" reactions
    3) remove all exchange reactions
    4) add all exchange reactions from the environmental ball 
    5) Even if all reactions are added, transporters are still necessary to convert met_e to met_c
    
    '''
    def __init__(self, label, path_to_models):
        self.path = path_to_models
        self.files  = os.listdir(self.path)
    
    
    def clone_model_without_transp(self, model):
        mc = Model(model.id)
        
        for reaction in model.reactions:
            
            if 'EX_' in reaction.id and '_e' in reaction.id:
                pass
            else:
                for met in reaction.metabolites:
                    if not mc.metabolites.has_id(met.id):
                        mc.add_metabolites(met.copy())
                mc.add_reaction(reaction.copy())
        
        mc.reactions.bio1.objective_coefficient=1
            
        mc.optimize()
        return mc    
            
    def make_ensemble(self):
        #open the first model
        model = cobra.io.read_sbml_model(os.path.join(self.path, self.files[0]))
        #join all non exchange reactopms from models
                
        
        mc= self.clone_model_without_transp(model)
        
        for i in self.files:
            mod= cobra.io.read_sbml_model(os.path.join(self.path, i))
            for reaction in mod.reactions:
                if 'rxn' in reaction.id:
                    if not mc.reactions.has_id(reaction.id):
                        mc.add_reaction(reaction.copy())
                        mc.repair()
        return mc
    
    def add_transporter(self, model):
        
        
        ev = Env_ball(1000)
        
        for reaction in ev.transporters:
            met_id=reaction.replace('EX_','')
            met_name = ev.metabolites_d[met_id]
            react = Reaction(reaction)
            react.name = 'export of ' + met_name
            react.lower_bound = -1000.  # This is the default
            react.upper_bound = 1000.  # This is the default
            if not model.metabolites.has_id(met_id):
                m_e = Metabolite(met_id, name=met_name,compartment='e')
                react.add_metabolites({m_e: -1.0})
                model.add_reactions([react])
            else:
                react.add_metabolites({model.metabolites.get_by_id(met_id): -1.0})
                model.add_reactions([react])
                
        model.repair()
        model.optimize()
    
    def write_model(self, model, path):
        cobra.io.write_sbml_model(model,path)    



########### usage ##############
        
#em = Ensemble_model('Actinomycetaceae', '/home/daniel/generative_models/models/vaginal_models/Actinomycetaceae/gapfilled')
#ensembl = em.make_ensemble()
#em.add_transporter(ensembl)
#em.write_model(ensembl, '/home/daniel/generative_models/models/test_ensmbl.sbml')
        
        
        
        