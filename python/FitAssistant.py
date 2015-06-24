import json
from ROOT import RooRealVar, RooArgList, RooClassFactory, RooUniform, RooArgSet
from ROOT import RooCategory, RooExtendPdf, RooAddPdf, RooSimultaneous
from ROOT import RooDataHist


class FitAssistant:

    def __init__(self, filename):
        with open(filename) as f: 
            j = json.load(f)
        self.taggers = j["taggers"]    
        self.workPoints = j["workPoints"]    
        self.data_mul = j["data_tag_multiplicity"]
        self.mc_effs = j["mc_tagging_eff"] 
        self.effs = [[RooRealVar(e[0],e[0], 0.0, 1.0),
                      e[1]] for e in j["eff_names"]]    
        self.xsecs = [[RooRealVar(x[0],x[0], x[2], x[2]*0.5, x[2]*1.5),
                       x[1]] for x in j["cross_sections"]]  
        self.lumi = [RooRealVar("lumi","lumi",j["lumi"][0]), j["lumi"][1]]
        self.kappa = [RooRealVar("kappa","kappa", 1.0, 0.0, 2.0) , False]         
        self.mul_tag = RooCategory("mul_tag","Event Tag Multiplicity")
        self.arg_list = RooArgList("arg_list")
        for eff in self.effs: self.arg_list.add(eff[0])
        for xsec in self.xsecs: self.arg_list.add(xsec[0])
        self.arg_list.add(self.lumi[0])
        self.arg_list.add(self.kappa[0])
        # obtain compile expression for the FtCM formulas
        formulas = []
        for i, formula in enumerate(j["formula_array"]): 
            f_name = "formula_N{}".format(i) 
            instance = RooClassFactory.makePdfInstance( f_name,  
                                                        formula,
                                                        self.arg_list)
            formulas.append(instance) 
            self.mul_tag.defineType("{}_tag_jets".format(i))
           
        f_list = RooArgList("f_list")         
        for f in formulas: f_list.add(f)
        coefs = RooArgList("coefs")
        for f in range(len(formulas)): coefs.add(RooRealVar("coef_{}".format(f), "",
                                                            1.0/len(formulas))) 
        self.uni = RooUniform("uni","uniform",RooArgSet()) 
        self.formulas = [RooExtendPdf("ext_"+f.GetName(), "ext_"+f.GetTitle(), 
                                      self.uni , f) for f in formulas]  
        self.f_list = RooArgList("ext_f_list")
        for f in self.formulas: self.f_list.add(f)
        self.simulpdf = RooSimultaneous("simulpdf","simulpdf",
                                         self.f_list, self.mul_tag) 


    def fix_parameters(self, tagger, workPoint):
        i_tag = self.taggers.index(tagger)
        i_wp = self.workPoints[i_tag].index(workPoint)
        # set default values and fix effs
        for i,eff in enumerate(self.effs):
            eff[0].setVal(self.mc_effs[i_tag][i_wp][i])
            eff[0].setConstant(eff[1])
        for xsec in self.xsecs:    
            xsec[0].setConstant(xsec[1])
        self.lumi[0].setConstant(self.lumi[1])    
        self.kappa[0].setConstant(self.kappa[1])

    def data_hist(self, tagger, workPoint):
        i_tag = self.taggers.index(tagger)
        i_wp = self.workPoints[i_tag].index(workPoint)
        data_hist = RooDataHist("data","data", RooArgSet(self.mul_tag))
        for i, val in enumerate(self.data_mul[i_tag][i_wp]):
            self.mul_tag.setIndex(i)
            data_hist.add(RooArgSet(self.mul_tag),val)
        return data_hist

            

        



            

