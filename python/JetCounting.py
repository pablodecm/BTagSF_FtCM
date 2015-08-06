
from ROOT import TChain, vector
from ROOT import JetCounter
from os.path import split, exists
from os import makedirs


class JetCountingManager():
    

    def __init__(self, ptBins, etaBins, taggers, oldMuonSF = True, out_dir = "../output/" ):
    
        # pointer to selector
        self.jetCounter = JetCounter()

        # pt bins
        vec_ptBins = vector('double')()
        vec_ptBins += ptBins
        self.jetCounter.setPtBins(vec_ptBins)
        # eta bins
        vec_etaBins = vector('double')()
        vec_etaBins += etaBins 
        self.jetCounter.setEtaBins(vec_etaBins)
        # use old SF
        self.jetCounter.useOldMuonSF(oldMuonSF)

        for k, v in taggers.items():
            workPoints = vector('double')()
            workPoints += v
            self.jetCounter.addTagger(k, workPoints) 

        self.out_dir = out_dir    

        return None


    def process(self, samples, isData = False):
   
        if not exists(self.out_dir):
            makedirs(self.out_dir)
        out_fn = []
        for sample in samples:
            tchain = TChain("tree")
            tchain.Add(sample)
            tchain.Process(self.jetCounter, "isData"*isData)
            out_fn.append(self.out_dir+split(sample)[1].replace(".root",".json"))
            self.jetCounter.serialize(out_fn[-1])
            self.jetCounter.resetJetRegistry()

        return out_fn


            
            
            



    
