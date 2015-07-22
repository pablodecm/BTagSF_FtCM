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
        self.xSec = xSec
        self.nEventGen = nEventGen
        # set atributes from JSON fields
        for k,v in j.items():
            if k == "cat_counts":
                self.cat_counts = {tuple(int(n) for n in k_c) : v_c for k_c,v_c in v.items()}
            else:    
                setattr(self,k,v)

    def pretag_eff(self):     
        return float(self.nEventPass[1])/float(self.nEventGen)
        
    def tag_multiplicity(self, tagger, workPoint):
        i_tag = self.taggers.index(tagger)
        i_wp = self.workPoints[i_tag].index(workPoint)
        return self.tagMultiplicity[i_tag][i_wp]

    def cat_fractions(self):
        return {k : float(v[0])/float(self.nEventPass[1]) for k,v in self.cat_counts.items()}

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

    def set_tag_wp( self, tagger, workPoint):
        self.i_tag = self.mc_samples[0].taggers.index(tagger)
        self.i_wp = self.mc_samples[0].workPoints[self.i_tag].index(workPoint)

    def ptBins(self):
        return self.mc_samples[0].ptBins

    def good_jets_by_type(self, x_jets = False):
        gj_dict = { "data": [], "b_jets": [], "c_jets": [], "l_jets": []}
        if x_jets: gj_dict["x_jets"] = []
        for k_c in range(len(self.data_samples[0].good_x_jets)):
            gj_dict["data"].append(sum([ s.good_x_jets[k_c] for s in self.data_samples]))
        for cat in (["b_jets","c_jets","l_jets"]+["x_jets"]*x_jets):
            for k_c in range(len(self.mc_samples[0].good_b_jets)):
                j_count = [self.lumi*s.xSec*getattr(s,"good_"+cat)[k_c]/s.nEventGen \
                                                          for s in self.mc_samples]
                gj_dict[cat].append(sum(j_count))
        if not x_jets:
            for k_c in range(len(self.mc_samples[0].good_b_jets)):
                x_count = [ self.lumi*s.good_x_jets[k_c]*s.xSec/s.nEventGen for s in self.mc_samples]
                gj_dict["l_jets"][k_c] += sum(x_count)
        return gj_dict        

    def good_jets_by_sample(self):
        gj_dict = {"data" : []}
        for k_c in range(len(self.data_samples[0].good_x_jets)):
            gj_dict["data"].append(sum([ s.good_x_jets[k_c] for s in self.data_samples]))
        for s in self.mc_samples:
            gj_dict[s.s_name] = [sum(c_t) for c_t in zip(*[getattr(s, 
                                "good_{}_jets".format(t)) for t in ["b","c","l","x"]])]
            for k_c in range(len(self.mc_samples[0].good_b_jets)):
                gj_dict[s.s_name][k_c] = self.lumi*gj_dict[s.s_name][k_c]*s.xSec/s.nEventGen 
        return gj_dict     

    def tag_jets_by_type(self, x_jets = False):
        gj_dict = { "data": [], "b_jets": [], "c_jets": [], "l_jets": []}
        if x_jets: gj_dict["x_jets"] = []
        for k_c in range(len(self.data_samples[0].tag_x_jets[self.i_tag][self.i_wp])):
            gj_dict["data"].append(sum([ s.tag_x_jets[self.i_tag][self.i_wp][k_c] for s in self.data_samples]))
        for cat in (["b_jets","c_jets","l_jets"]+["x_jets"]*x_jets):
            for k_c in range(len(self.mc_samples[0].tag_b_jets)):
                j_count = [self.lumi*s.xSec*getattr(s,
                    "tag_"+cat)[self.i_tag][self.i_wp][k_c]/s.nEventGen for s in self.mc_samples]
                gj_dict[cat].append(sum(j_count))
        if not x_jets:
            for k_c in range(len(self.mc_samples[0].tag_b_jets[self.i_tag][self.i_wp] )):
                x_count = [ self.lumi*s.tag_x_jets[self.i_tag][self.i_wp][k_c] * \
                                  s.xSec/s.nEventGen for s in self.mc_samples]
                gj_dict["l_jets"][k_c] += sum(x_count)
        return gj_dict        

    def tag_jets_by_sample(self):
        gj_dict = {"data" : []}
        for k_c in range(len(self.data_samples[0].tag_x_jets[self.i_tag][self.i_wp])):
            gj_dict["data"].append(sum([ s.tag_x_jets[self.i_tag][self.i_wp][k_c] \
                                                      for s in self.data_samples]))
        for s in self.mc_samples:
            gj_dict[s.s_name] = [sum(c_t) for c_t in zip(*[getattr(s, 
                                "tag_{}_jets".format(t))[self.i_tag][self.i_wp] \
                                                  for t in ["b","c","l","x"]])]
            for k_c in range(len(self.mc_samples[0].tag_b_jets[self.i_tag][self.i_wp])):
                gj_dict[s.s_name][k_c] = self.lumi*gj_dict[s.s_name][k_c]*s.xSec/s.nEventGen 

        return gj_dict     

    def tag_mul_by_sample(self):
        tm_dict = {}
        tag_mul_ps = [ s.tagMultiplicity[self.i_tag][self.i_wp] for s in self.data_samples ] 
        tm_dict["data"] =  [sum(n_jets) for n_jets in zip(*tag_mul_ps)] 
        for s in self.mc_samples:
            tm_dict[s.s_name] = [self.lumi*b_c*s.xSec/s.nEventGen  \
                    for  b_c in s.tagMultiplicity[self.i_tag][self.i_wp]]

        return tm_dict     
 
    def jet_mul_by_sample(self):
        jm_dict = {}
        jet_mul_ps = [ s.jetMultiplicity for s in self.data_samples ] 
        jm_dict["data"] =  [sum(n_jets) for n_jets in zip(*jet_mul_ps)] 
        for s in self.mc_samples:
            jm_dict[s.s_name] = [self.lumi*b_c*s.xSec/s.nEventGen  \
                    for  b_c in s.jetMultiplicity]

        return jm_dict     
      
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



