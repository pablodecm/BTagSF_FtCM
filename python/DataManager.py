import json
from sympy import ccode
try:
    from symengine import Symbol, sympify
except ImportError:
    from sympy import Symbol, sympify
from mut_framework.BTagSF_FtCM.FtCM_formula import formula_categories

class Sample:
    
    def __init__(self, filename, s_name,  nEventGen, xSec, isData = False, isSignal = False): 
        with open(filename) as f:
            j = json.load(f)
        self.s_name = s_name    
        self.isData = isData
        self.isSignal = isSignal
        self.taggers = j["taggers"]    
        self.workPoints = j["workPoints"]    
        self.ptBins = j["ptBins"]
        self.etaBins = j["etaBins"]
        self.nEventGen = nEventGen
        self.good_cat_jets = j["good_cat_jets"] 
        self.tag_cat_jets = j["tag_cat_jets"]
        self.xSec = xSec
        self.nEventPass = j["nEventPass"]
        self.tagMultiplicity = j["tagMultiplicity"]
        self.cat_counts = {tuple(int(n) for n in k) : v for k,v in j["cat_counts"].items()} 

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
        if not self.isSignal: factor *= Symbol("kappa")
        return factor*formula_categories(self.cat_fractions(), n_tags)

        
class DataManager:

    def __init__( self, mc_samples, data_samples, lumi):
        self.mc_samples = mc_samples 
        self.data_samples = data_samples 
        self.lumi = lumi
        self.n_cat = len(mc_samples[0].good_cat_jets)

    def data_tag_multiplicity(self):
        tag_mul = []
        for t in range(len(self.data_samples[0].taggers)): 
            tag_mul.append([])
            for wp in range(len(self.data_samples[0].workPoints[t])):
                tag_mul_ps = [ s.tagMultiplicity[t][wp] for s in self.data_samples ] 
                tag_mul[-1].append([sum(n_jets) for n_jets in zip(*tag_mul_ps)]) 
        return tag_mul                 

    def mc_tagging_eff(self):
        g_c_j_ps = [ s.good_cat_jets for s in self.mc_samples ] 
        g_c_j = [sum(n_jets) for n_jets in zip(*g_c_j_ps)] 
        mc_eff = []
        for t in range(len(self.mc_samples[0].taggers)): 
            mc_eff.append([])
            for wp in range(len(self.mc_samples[0].workPoints[t])):
                t_c_j_ps = [ s.tag_cat_jets[t][wp] for s in self.mc_samples ] 
                t_c_j = [sum(n_jets) for n_jets in zip(*t_c_j_ps)] 
                mc_eff[-1].append([float(t_j)/g_j for t_j,g_j in zip(t_c_j,g_c_j)])
        return mc_eff                 

    def formula(self, n_tags): 
        f = sympify(0)
        for s in self.mc_samples:
            f+=s.formula(n_tags)
        return f        
    
    def formula_array(self, max_n_tags = 7, extra_subs = {}):
        f_array = []
        for n_tags in range(max_n_tags):
            f_array.append(self.formula(n_tags).subs(extra_subs))
        return f_array    
     
    def cross_sections(self):            
        xsecs = []
        for s in self.mc_samples:
            xsecs.append([])
            xsecs[-1].append("xsec_{}".format(s.s_name))
            xsecs[-1].append(False if s.isSignal else True)
            xsecs[-1].append(s.xSec)
        return xsecs

    def eff_names(self, fix = [0,1,2]):
        eff_names = []
        for i in range(self.n_cat):
            eff_names.append([])
            eff_names[-1].append("eff_{}".format(i))
            eff_names[-1].append(True if (i in fix) else False) 
        return eff_names   

    def save_to_file( self, filename, max_n_tags =7):
        save_dict = { "cross_sections": self.cross_sections(),
                      "taggers": self.mc_samples[0].taggers,
                      "workPoints": self.mc_samples[0].workPoints,
                      "data_tag_multiplicity": self.data_tag_multiplicity(),
                      "lumi": [self.lumi, True],
                      "eff_names": self.eff_names(),
                      "mc_tagging_eff" : self.mc_tagging_eff(),
                      "formula_array": [ccode(f) for f in self.formula_array()]}
        with open(filename,'w') as f:
            json.dump(save_dict, f, indent = 4)



