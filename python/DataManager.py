from ast import literal_eval
import json
from mut_framework.BTagSF_FtCM.FtCM_formula import formula_categories
from sympy import IndexedBase, Symbol, sympify
import pickle

class Sample:
    
    def __init__(self, filename, s_name,  nEventGen, xSec): 
        with open(filename) as f:
            j = json.load(f)
        self.s_name = s_name    
        self.taggers = j["taggers"]    
        self.workPoints = j["workPoints"]    
        self.ptBins = j["ptBins"]
        self.etaBins = j["etaBins"]
        self.nEventGen = nEventGen
        self.good_cat_jets = j["good_cat_jets"] 
        self.tag_cat_jet = j["tag_cat_jets"]
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

        
class DataManager:

    def __init__( self, samples, signal = "TTbar" ):
        self.samples = samples 
        self.signal = signal 
        self.nCat = len(samples[0].cat_counts.keys()[0])

    def full_formula(self, n_tags):    
        f = sympify(0)
        for s in self.samples:
            if s.s_name is self.signal:
                f += s.formula(n_tags)
            else:
                kappa = Symbol("kappa")
                f += kappa*s.formula(n_tags).subs(Symbol("xsec_{}".format(s.s_name)),s.xSec)
        return f

    def subs_formula(self, n_tags, extra_subs = {}):
        E = IndexedBase("E")
        subs_dict = {E[0] : Symbol("eff_x"),
                     E[1] : Symbol("eff_l"),
                     E[2] : Symbol("eff_c")}  
        b_eff_dict ={ E[i] : Symbol("eff_b{}".format(i-3)) for i in range(3,self.nCat)} 
        subs_dict = dict(subs_dict.items() + b_eff_dict.items() + extra_subs.items())
        return self.full_formula(n_tags).subs(subs_dict)

    def save_formula_list( self, filename, max_n_tags =7):
        f_list = []
        for n_tags in range(max_n_tags):
            f_list.append(self.subs_formula(n_tags))
        with open(filename,'w') as f:
            pickle.dump(f_list, f)



