from ast import literal_eval
import json
from mut_framework.BTagSF_FtCM.FtCM_formula import formula_categories
from sympy import Symbol

class DataManager:
    
    def __init__(self, filename, s_name,  nEventGen, xSec): 
        with open(filename) as f:
            j = json.load(f)
        self.s_name = s_name    
        self.taggers = j["taggers"]    
        self.workPoints = j["workPoints"]    
        self.ptBins = j["ptBins"]
        self.etaBins = j["etaBins"]
        self.nEventGen = nEventGen
        self.xSec = xSec
        self.nEventPass = j["nEventPass"]
        self.tagMultiplicity = j["tagMultiplicity"]
        self.cat_counts = {literal_eval(k) : v for k,v in j["cat_counts"].items()} 

    def pretag_eff(self):     
        return float(self.nEventPass[0])/float(self.nEventGen)
        
    def tag_multiplicity(self, tagger, workPoint):
        i_tag = self.taggers.index(tagger)
        i_wp = self.workPoints[i_tag].index(workPoint)
        return self.tagMultiplicity[i_tag][i_wp]

    def cat_fractions(self):
        return {k : float(v[0])/float(self.nEventPass[0]) for k,v in self.cat_counts.items()}

    def formula(self, n_tags):
        factor = Symbol("lumi")*Symbol("xsec_"+self.s_name)*self.pretag_eff()
        return factor*formula_categories(self.cat_fractions(), n_tags)

        



