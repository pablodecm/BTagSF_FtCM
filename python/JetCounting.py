
from ROOT import TChain, vector
from ROOT import JetCounter
from os.path import split

class JetCountingManager():
    
    jetCounter = JetCounter()
    out_dir = "../output/" 

    def __init__(self, ptBins, etaBins, taggers ):
        # pt bins
        vec_ptBins = vector('double')()
        vec_ptBins += ptBins
        self.jetCounter.setPtBins(vec_ptBins)
        # eta bins
        vec_etaBins = vector('double')()
        vec_etaBins += etaBins 
        self.jetCounter.setEtaBins(vec_etaBins)

        for k, v in taggers.items():
            workPoints = vector('double')()
            workPoints += v
            self.jetCounter.addTagger(k, workPoints) 

        return None


    def process(self, samples, isData = False):
   
        out_fn = []
        for sample in samples:
            tchain = TChain("tree")
            tchain.Add(sample)
            tchain.Process(self.jetCounter, "isData"*isData)
            out_fn.append(self.out_dir+split(sample)[1].replace(".root",".json"))
            self.jetCounter.serialize(out_fn[-1])
            self.jetCounter.resetJetRegistry()

        return out_fn


            
            
            



    
