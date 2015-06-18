import json
from ROOT import RooRealVar, RooArgList, RooClassFactory


class FitAssistant:

    def __init__(self, filename):
        with open(filename) as f: 
            j = json.load(f)
        self.effs = [[RooRealVar(e[0],e[0], 0.0, 1.0),
                      e[1]] for e in j["eff_names"]]    
        self.xsecs = [[RooRealVar(x[0],x[0], x[2], x[2]*0.5, x[2]*1.5),
                       e[1]] for x in j["cross_sections"]]  
        self.lumi = [RooRealVar("lumi","lumi",j["lumi"][0]), j["lumi"][1]]
        self.kappa = [RooRealVar("kappa","kappa", 1.0, 0.0, 2.0) , False]         

        self.arg_list = RooArgList("arg_list")
        for eff in self.effs: self.arg_list.add(eff[0])
        for xsec in self.xsecs: self.arg_list.add(xsec[0])
        self.arg_list.add(self.lumi[0])
        self.arg_list.add(self.kappa[0])
        # obtain compile expression for the FtCM formulas
        self.formulas = []
        for i, formula in enumerate(j["formula_array"]): 
            f_name = "formula_N{}".format(i) 
            instance = RooClassFactory.makePdfInstance( f_name,  
                                                        formula,
                                                        self.arg_list)
            self.formulas.append(instance) 



            

